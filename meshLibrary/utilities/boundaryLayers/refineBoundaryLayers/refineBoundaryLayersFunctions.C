/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "refineBoundaryLayers.H"
#include "boundaryLayerConstraintPlanner.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGen2DEngine.H"
#include "VRWGraphList.H"
#include "meshSurfacePartitioner.H"
#include "detectBoundaryLayers.H"

#include "labelledPair.H"
#include "labelledScalar.H"

# ifdef USE_OMP
#include <omp.h>
#include <map>
#include <utility>
# endif

//#define DEBUGLayer

# ifdef DEBUGLayer
#include "OFstream.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool refineBoundaryLayers::analyseLayers()
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatch = mse.boundaryFacePatches();

    meshSurfacePartitioner mPart(mse);
    //- Partial-patch BL policy:
    //- use the existing effective requested layer count as the scope.
    //- requestedLayers > 1  => partial detected layer is allowed
    //- requestedLayers <= 1 => preserve legacy all-or-nothing behaviour
    const PtrList<boundaryPatch>& partialBoundaries = mesh_.boundaries();
    boolList allowPartialLayerPatch(partialBoundaries.size(), false);

    forAll(partialBoundaries, patchI)
    {
        label requestedLayers = globalNumLayers_;
        const word pName = partialBoundaries[patchI].patchName();

        std::map<word, label>::const_iterator it =
            numLayersForPatch_.find(pName);

        if( it != numLayersForPatch_.end() )
            requestedLayers = it->second;

        allowPartialLayerPatch[patchI] = (requestedLayers > 1);

        if( allowPartialLayerPatch[patchI] )
        {
            Info << "BLPARTIALPATCH eligible patch " << pName
                 << " requestedLayers=" << requestedLayers << endl;
        }
    }

    detectBoundaryLayers dbl
    (
        mPart,
        is2DMesh_,
        allowPartialLayerPatch
    );

    const label nGroups = dbl.nDistinctLayers();
    const labelList& faceInLayer = dbl.faceInLayer();

    //- get the hair edges
    splitEdges_ = dbl.hairEdges();

    //- Transaction diagnostic: distinguish "no layers detected" from
    //- "layers detected but hair-edge generation failed".  This runs for
    //- both the ordinary pass and a restored pass-2 mesh.
    label nDetectedLayerFaces = 0;
    forAll(faceInLayer, bfI)
        if( faceInLayer[bfI] >= 0 )
            ++nDetectedLayerFaces;

    Info << "BLANALYSE"
         << " meshPoints=" << mesh_.points().size()
         << " meshFaces=" << mesh_.faces().size()
         << " meshCells=" << mesh_.cells().size()
         << " boundaryFaces=" << bFaces.size()
         << " nGroups=" << nGroups
         << " layerFaces=" << nDetectedLayerFaces
         << " hairEdges=" << splitEdges_.size()
         << endl;

    //- Zero hair edges cannot support refineBoundaryLayers.  Previously
    //- the validation loop below could vacuously succeed when every
    //- faceInLayer entry was negative, allowing generateNewFaces() to run
    //- against empty split-edge metadata.
    //
    //- Do NOT mark refinementValid_ false here: "nothing detected to refine"
    //- and "metadata became structurally inconsistent" are different states.
    //- The caller separately checks refinementCompleted(), which remains
    //- false because done_ is set only after generateNewCells() completes.
    if( splitEdges_.size() == 0 )
    {
        Info << "BLANALYSE_NO_HAIR_EDGES"
             << " nGroups=" << nGroups
             << " layerFaces=" << nDetectedLayerFaces
             << " boundaryFaces=" << bFaces.size()
             << " -- refinement cannot execute"
             << endl;
        return false;
    }

    # ifdef DEBUGLayer
    OFstream file("hairEdges.vtk");

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << 2*splitEdges_.size() << " float\n";
    forAll(splitEdges_, seI)
    {
        const point& p = mse.mesh().points()[splitEdges_[seI].start()];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;

        const point op = mse.mesh().points()[splitEdges_[seI].end()];

        file << op.x() << ' ' << op.y() << ' ' << op.z() << nl;
    }

    //- write lines
    file << "\nLINES " << splitEdges_.size()
         << " " << 3*splitEdges_.size() << nl;
    forAll(splitEdges_, eI)
    {
        file << 2 << " " << 2*eI << " " << (2*eI+1) << nl;
    }

    file << "\n";
    # endif

    //- create point to split edges addressing
    splitEdgesAtPoint_.reverseAddressing(splitEdges_);

    //- check if the layer is valid
    bool validLayer(true);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(faceInLayer, bfI)
    {
        if( faceInLayer[bfI] < 0 )
            continue;

        const face& bf = bFaces[bfI];

        forAll(bf, pI)
            if( splitEdgesAtPoint_.sizeOfRow(bf[pI]) == 0 )
                validLayer = false;
    }

    # ifdef DEBUGLayer
    Info << "Number of independent layers in the mesh is " << nGroups << endl;
    Info << "Is valid layer " << validLayer << endl;
    # endif

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    //- create patch name to index addressing
    std::map<word, label> patchNameToIndex;
    forAll(boundaries, patchI)
        patchNameToIndex[boundaries[patchI].patchName()] = patchI;

    //- check layer labels over a patch
    layerAtPatch_.setSize(boundaries.size());
    forAll(layerAtPatch_, i)
        layerAtPatch_[i].clear();
    List<DynList<label> > groupsAtPatch(boundaries.size());
    forAll(faceInLayer, bfI)
        groupsAtPatch[facePatch[bfI]].appendIfNotIn(faceInLayer[bfI]);

    //- set the information which patches have an extruded layer
    forAll(groupsAtPatch, patchI)
    {
        const DynList<label>& layers = groupsAtPatch[patchI];

        forAll(layers, i)
        {
            if( layers[i] < 0 )
            {
                if
                (
                    patchI >= 0 &&
                    patchI < label(allowPartialLayerPatch.size()) &&
                    allowPartialLayerPatch[patchI]
                )
                {
                    //- Keep the valid detected portion of an explicitly
                    //- requested BL patch. The local unsupported faces
                    //- remain faceInLayer < 0 and are protected below.
                    continue;
                }

                //- Legacy all-or-nothing behaviour for non-BL patches.
                layerAtPatch_[patchI].clear();
                break;
            }
            else
            {
                layerAtPatch_[patchI].append(layers[i]);
            }
        }
    }

    # ifdef DEBUGLayer
    Info << "Layer at patch " << layerAtPatch_ << endl;
    # endif

    //- set the information which patches are a single boundary layer face
    patchesInLayer_.setSize(nGroups);
    forAll(layerAtPatch_, patchI)
    {
        const DynList<label>& layers = layerAtPatch_[patchI];

        forAll(layers, i)
            patchesInLayer_[layers[i]].append
            (
                boundaries[patchI].patchName()
            );
    }

    # ifdef DEBUGLayer
    Info << "Patches in layer " << patchesInLayer_ << endl;
    # endif

    //- set the number of boundary layers for each patch
    labelList nLayersAtPatch(layerAtPatch_.size(), -1);
    boolList protectedValue(layerAtPatch_.size(), false);

    forAll(patchesInLayer_, layerI)
    {
        const DynList<word>& layerPatches = patchesInLayer_[layerI];

        label maxNumLayers(1);
        bool hasLocalValue(false);

        //- find the maximum requested number of layers over the layer
        forAll(layerPatches, lpI)
        {
            const word pName = layerPatches[lpI];

            std::map<word, label>::const_iterator it =
                numLayersForPatch_.find(pName);

            if( it != numLayersForPatch_.end() )
            {
                const label ptchI = patchNameToIndex[pName];

                // CFMitch/shared-modern hard no-BL invariant:
                //
                // setNumberOfLayersForPatch() documents nLayers < 2
                // as disabling boundary layers on the patch.  The original
                // cfMesh layer-group logic allowed such a patch to inherit
                // maxNumLayers from another patch in the same detected
                // layer group, which could silently turn an explicit
                // nLayers=0 surface into a full BL-generating surface.
                //
                // Internally one layer means "retain the unsplit parent",
                // therefore normalize every explicit nLayers < 2 request
                // to one and protect it from group maximum propagation.
                //
                // This is shared by legacyEnhanced and constraintPlanner.
                // classicCfMesh uses its separate frozen implementation.
                if( it->second < 2 )
                {
                    nLayersAtPatch[ptchI] = 1;
                    protectedValue[ptchI] = true;
                    hasLocalValue = true;
                }
                //- check if the layer is interrupted at this patch
                else if
                (
                    discontinuousLayersForPatch_.find(pName) !=
                    discontinuousLayersForPatch_.end()
                )
                {
                    //- set the number of layers and lock this location
                    nLayersAtPatch[ptchI] = it->second;
                    protectedValue[ptchI] = true;
                    hasLocalValue = true;
                }
                else
                {
                    //- take the maximum number of layers
                    maxNumLayers = Foam::max(maxNumLayers, it->second);
                    hasLocalValue = true;
                }
            }
        }

        //- apply the global value if no local values exist
        if( !hasLocalValue )
            maxNumLayers = globalNumLayers_;

        //- apply the maximum number of ayer of all unprotected patches
        forAll(layerPatches, lpI)
        {
            const label ptchI = patchNameToIndex[layerPatches[lpI]];

            if( !protectedValue[ptchI] )
                nLayersAtPatch[ptchI] = maxNumLayers;
        }
    }

    if( is2DMesh_ )
    {
        polyMeshGen2DEngine mesh2DEngine(mesh_);
        const boolList& zMinPoint = mesh2DEngine.zMinPoints();
        const boolList& zMaxPoint = mesh2DEngine.zMaxPoints();

        const faceList::subList& bFaces = mse.boundaryFaces();

        boolList allZMax(mesh_.boundaries().size(), true);
        boolList allZMin(mesh_.boundaries().size(), true);

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                if( !zMinPoint[bf[pI]] )
                    allZMin[facePatch[bfI]] = false;
                if( !zMaxPoint[bf[pI]] )
                    allZMax[facePatch[bfI]] = false;
            }
        }

        //- mark empty patches as already used
        forAll(allZMin, patchI)
        {
            if( allZMin[patchI] ^ allZMax[patchI] )
            {
                nLayersAtPatch[patchI] = -1;
                layerAtPatch_[patchI].clear();
            }
        }
    }

    //- perform reduction over all processors
    reduce(nLayersAtPatch, maxOp<labelList>());

    // CFMitch/shared-modern hard no-BL audit.
    // An explicitly disabled patch must resolve to <= 1 and may never
    // inherit an active layer count from its detected layer group.
    {
        label nHardNoBLPatches = 0;
        label nHardNoBLSafe = 0;
        label nHardNoBLViolations = 0;

        forAll(nLayersAtPatch, patchI)
        {
            const word& pName = boundaries[patchI].patchName();

            const std::map<word, label>::const_iterator it =
                numLayersForPatch_.find(pName);

            if
            (
                it == numLayersForPatch_.end()
             || it->second >= 2
            )
                continue;

            ++nHardNoBLPatches;

            if( nLayersAtPatch[patchI] <= 1 )
                ++nHardNoBLSafe;
            else
                ++nHardNoBLViolations;
        }

        Info
            << "BL_HARD_NOBL_AUDIT"
            << " requestedPatches=" << nHardNoBLPatches
            << " safe=" << nHardNoBLSafe
            << " violations=" << nHardNoBLViolations
            << endl;
    }

    # ifdef DEBUGLayer
    Pout << "nLayersAtPatch " << nLayersAtPatch << endl;
    # endif

    //- set the number of boundary layers which shall be generated above
    //- each boundary face
    nLayersAtBndFace_.setSize(facePatch.size());
    nLayersAtBndFace_ = globalNumLayers_;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nLayersAtBndFace_, bfI)
    {
        const label patchI = facePatch[bfI];

        //- A locally unsupported detected-layer face must never be sent
        //- into boundary-layer refinement. This is the critical guard
        //- which lets the rest of a configured BL patch refine safely.
        if( faceInLayer[bfI] < 0 )
        {
            nLayersAtBndFace_[bfI] = 1;
            continue;
        }

        if( nLayersAtPatch[patchI] < 0 )
        {
            nLayersAtBndFace_[bfI] = 1;
        }
        else
        {
            nLayersAtBndFace_[bfI] = nLayersAtPatch[patchI];

            if( specialMode_ )
            {
                ++nLayersAtBndFace_[bfI];
            }
        }
    }

    // BL/BL junction ramp: faces touching sharp junction points get
    // reduced layer count to prevent degenerate cells at blade/hub/shroud.
    // Ring 0 (junction face)    -> 1 layer
    // Ring 1 (neighbor face)    -> max 2 layers
    // Ring 2+ resumes full nLayers
    //
    // Experimental baseline synchronization:
    // the installed test library used for the 15-layer contact-smoother
    // A/B has this historical BL/BL layer-count ramp disabled.  Keep the
    // source behavior identical while testing the contact sweep transaction.
    const bool enableBlblJunctionLayerCountRamp = false;

    if( !enableBlblJunctionLayerCountRamp )
        Info << "BL/BL junction ramp: DISABLED A/B" << endl;

    if
    (
        enableBlblJunctionLayerCountRamp
     && blblJunctionPoints_.size() > 0
    )
    {
        const meshSurfaceEngine& mseLoc = surfaceEngine();
        const VRWGraph& ptFaces = mseLoc.pointFaces();

        // Build mesh-point to boundary-point map
        labelList meshToBnd(mesh_.points().size(), -1);
        const labelList& bPts = mseLoc.boundaryPoints();
        forAll(bPts, bpI)
            meshToBnd[bPts[bpI]] = bpI;

        // Ring 0: faces directly touching junction points -> 1 layer
        // Only reduce faces that already have more than 1 layer
        boolList ring0face(nLayersAtBndFace_.size(), false);
        forAllConstIter(labelHashSet, blblJunctionPoints_, it)
        {
            const label bpI = it.key();
            if( bpI < 0 || bpI >= label(ptFaces.size()) ) continue;
            // Skip finite ramp seeds -- not hard junction caps
            if( rampSeedPoints_.size() > bpI && rampSeedPoints_[bpI] ) continue;
            forAllRow(ptFaces, bpI, pfI)
            {
                const label bfI = ptFaces(bpI, pfI);
                if( bfI < 0 || bfI >= label(nLayersAtBndFace_.size()) ) continue;
                // Virtual topology takes priority -- skip VT-handled faces
                if( bfI < label(vtFaceRing_.size()) && vtFaceRing_[bfI] >= 0 ) continue;
                if( nLayersAtBndFace_[bfI] > 1 )
                {
                    ring0face[bfI] = true;
                    nLayersAtBndFace_[bfI] = 1;
                }
            }
        }

        // Ring 1: faces sharing a point with ring0 faces -> max 2 layers
        // Only reduce faces that already have more than 2 layers
        const faceList::subList& bFacesLoc = mseLoc.boundaryFaces();
        boolList ring1face(nLayersAtBndFace_.size(), false);
        forAll(nLayersAtBndFace_, bfI)
        {
            if( !ring0face[bfI] ) continue;
            const face& f = bFacesLoc[bfI];
            forAll(f, pI)
            {
                const label meshPtI = f[pI];
                if( meshPtI < 0 || meshPtI >= label(meshToBnd.size()) ) continue;
                const label bpI = meshToBnd[meshPtI];
                if( bpI < 0 || bpI >= label(ptFaces.size()) ) continue;
                forAllRow(ptFaces, bpI, pfI)
                {
                    const label nbfI = ptFaces(bpI, pfI);
                    if( nbfI < 0 || nbfI >= label(nLayersAtBndFace_.size()) ) continue;
                    if( ring0face[nbfI] ) continue;
                    // Virtual topology takes priority -- skip VT-handled faces
                    if( nbfI < label(vtFaceRing_.size()) && vtFaceRing_[nbfI] >= 0 ) continue;
                    if( nLayersAtBndFace_[nbfI] > 2 )
                    {
                        ring1face[nbfI] = true;
                        nLayersAtBndFace_[nbfI] = 2;
                    }
                }
            }
        }

        label nRing0 = 0, nRing1 = 0;
        forAll(nLayersAtBndFace_, bfI)
        {
            if( ring0face[bfI] ) ++nRing0;
            else if( ring1face[bfI] ) ++nRing1;
        }
        Info << "BL/BL junction ramp: ring0=" << nRing0
             << " faces forced to 1 layer, ring1=" << nRing1
             << " faces capped at 2 layers" << endl;
    }

    // Acute corner face cap: controlled by acuteCornerCapLayers_ member
    // Set via setAcuteCornerCapLayers() before refineLayers()
    Info << "Acute corner face cap switch: " << acuteCornerCapLayers_
         << " acutePts=" << blblAcuteCornerPoints_.size() << endl;
    if( acuteCornerCapLayers_ && blblAcuteCornerPoints_.size() > 0 )
    {
        const meshSurfaceEngine& mseLoc = surfaceEngine();
        const VRWGraph& ptFaces = mseLoc.pointFaces();
        const faceList::subList& bFacesLoc = mseLoc.boundaryFaces();
        labelList meshToBnd2(mesh_.points().size(), -1);
        const labelList& bPts2 = mseLoc.boundaryPoints();
        forAll(bPts2, bpI)
            meshToBnd2[bPts2[bpI]] = bpI;

        boolList acRing0(nLayersAtBndFace_.size(), false);
        forAllConstIter(labelHashSet, blblJunctionPoints_, it)
        {
            const label bpI = it.key();
            if( bpI < 0 || bpI >= label(ptFaces.size()) ) continue;
            // Skip finite ramp seeds -- not hard junction caps
            if( rampSeedPoints_.size() > bpI && rampSeedPoints_[bpI] ) continue;
            forAllRow(ptFaces, bpI, pfI)
            {
                const label bfI = ptFaces(bpI, pfI);
                if( bfI < 0 || bfI >= label(nLayersAtBndFace_.size()) ) continue;
                if( nLayersAtBndFace_[bfI] > 1 )
                    nLayersAtBndFace_[bfI] = 1;
                acRing0[bfI] = true;
            }
        }

        boolList acRing1(nLayersAtBndFace_.size(), false);
        forAll(nLayersAtBndFace_, bfI)
        {
            if( !acRing0[bfI] ) continue;
            const face& f = bFacesLoc[bfI];
            forAll(f, pI)
            {
                const label meshPtI = f[pI];
                if( meshPtI < 0 || meshPtI >= label(meshToBnd2.size()) ) continue;
                const label bpI = meshToBnd2[meshPtI];
                if( bpI < 0 || bpI >= label(ptFaces.size()) ) continue;
                forAllRow(ptFaces, bpI, pfI)
                {
                    const label nbfI = ptFaces(bpI, pfI);
                    if( nbfI < 0 || nbfI >= label(nLayersAtBndFace_.size()) ) continue;
                    if( acRing0[nbfI] ) continue;
                    if( nLayersAtBndFace_[nbfI] > 1 )
                        nLayersAtBndFace_[nbfI] = 1;
                    acRing1[nbfI] = true;
                }
            }
        }

        boolList acRing2(nLayersAtBndFace_.size(), false);
        forAll(nLayersAtBndFace_, bfI)
        {
            if( !acRing1[bfI] ) continue;
            const face& f = bFacesLoc[bfI];
            forAll(f, pI)
            {
                const label meshPtI = f[pI];
                if( meshPtI < 0 || meshPtI >= label(meshToBnd2.size()) ) continue;
                const label bpI = meshToBnd2[meshPtI];
                if( bpI < 0 || bpI >= label(ptFaces.size()) ) continue;
                forAllRow(ptFaces, bpI, pfI)
                {
                    const label nbfI = ptFaces(bpI, pfI);
                    if( nbfI < 0 || nbfI >= label(nLayersAtBndFace_.size()) ) continue;
                    if( acRing0[nbfI] || acRing1[nbfI] ) continue;
                    if( nLayersAtBndFace_[nbfI] > 2 )
                        nLayersAtBndFace_[nbfI] = 2;
                    acRing2[nbfI] = true;
                }
            }
        }

        label nAC0=0, nAC1=0, nAC2=0;
        forAll(nLayersAtBndFace_, bfI)
        {
            if( acRing0[bfI] ) ++nAC0;
            else if( acRing1[bfI] ) ++nAC1;
            else if( acRing2[bfI] ) ++nAC2;
        }
        Info << "Acute corner face cap: ring0=" << nAC0
             << " ring1=" << nAC1
             << " ring2=" << nAC2
             << " faces capped" << endl;
    }

    // ============================================================
    // CFMITCH V8A NEUTRAL SEAM LAYER BACKOFF
    //
    // Preserve the macro hair and first-layer thickness.
    //
    // Only reduce the requested OUTER layer count on active BL
    // boundary faces touching a stable BL/neutral seam point.
    //
    // No coordinates are changed.
    // No h1 values are changed.
    // No faces or cells are deleted.
    // No neutral/no-BL face is reactivated.
    //
    // The resulting caps enter the existing canonical
    // forcedMaxLayersAtFace_ -> constraintPlanner path.
    // ============================================================

    if( !neutralSeamBackoffMeshPoints_.empty() )
    {
        if( neutralSeamBackoffMaxLayers_ < 1 )
        {
            FatalErrorIn
            (
                "refineBoundaryLayers::analyseLayers()"
            )
                << "CFMITCH V8A neutral seam max layers must be >= 1"
                << exit(FatalError);
        }


        const meshSurfaceEngine& seamMse =
            surfaceEngine();

        const labelList& seamBoundaryPoints =
            seamMse.boundaryPoints();

        const VRWGraph& seamPointFaces =
            seamMse.pointFaces();


        // global mesh-point label -> CURRENT boundary-point index
        labelList seamMeshToBp
        (
            mesh_.points().size(),
            -1
        );

        forAll(seamBoundaryPoints, bpI)
        {
            const label meshPtI =
                seamBoundaryPoints[bpI];

            if
            (
                meshPtI >= 0
             && meshPtI < label(seamMeshToBp.size())
            )
            {
                seamMeshToBp[meshPtI] =
                    bpI;
            }
        }


        labelHashSet seamSeedFaces;

        label nMappedPoints = 0;
        label nMissingPoints = 0;
        label nInactiveFaces = 0;
        label nOutOfRangeFaces = 0;


        forAllConstIter
        (
            labelHashSet,
            neutralSeamBackoffMeshPoints_,
            seamIt
        )
        {
            const label meshPtI =
                seamIt.key();


            if
            (
                meshPtI < 0
             || meshPtI >= label(seamMeshToBp.size())
            )
            {
                ++nMissingPoints;
                continue;
            }


            const label bpI =
                seamMeshToBp[meshPtI];


            if
            (
                bpI < 0
             || bpI >= label(seamPointFaces.size())
            )
            {
                ++nMissingPoints;
                continue;
            }


            ++nMappedPoints;


            forAllRow(seamPointFaces, bpI, pfI)
            {
                const label bfI =
                    seamPointFaces(bpI, pfI);


                if
                (
                    bfI < 0
                 || bfI >= label(nLayersAtBndFace_.size())
                )
                {
                    ++nOutOfRangeFaces;
                    continue;
                }


                // Only ACTIVE BL faces participate.
                //
                // Neutral/termination/no-BL faces are already at one
                // layer and must remain untouched.
                if( nLayersAtBndFace_[bfI] <= 1 )
                {
                    ++nInactiveFaces;
                    continue;
                }


                seamSeedFaces.insert(bfI);
            }
        }


        label nInserted = 0;
        label nLoweredExisting = 0;
        label nKeptStricter = 0;


        forAllConstIter
        (
            labelHashSet,
            seamSeedFaces,
            faceIt
        )
        {
            const label bfI =
                faceIt.key();


            if( forcedMaxLayersAtFace_.found(bfI) )
            {
                const label oldCap =
                    forcedMaxLayersAtFace_[bfI];


                if
                (
                    neutralSeamBackoffMaxLayers_
                  < oldCap
                )
                {
                    forcedMaxLayersAtFace_[bfI] =
                        neutralSeamBackoffMaxLayers_;

                    ++nLoweredExisting;
                }
                else
                {
                    // Another subsystem already requested an equal
                    // or stricter cap. Preserve the stricter policy.
                    ++nKeptStricter;
                }
            }
            else
            {
                forcedMaxLayersAtFace_.insert
                (
                    bfI,
                    neutralSeamBackoffMaxLayers_
                );

                ++nInserted;
            }
        }


        Info
            << "CFMITCH V8A NEUTRAL_SEAM_BACKOFF:"
            << " seamPoints="
            << neutralSeamBackoffMeshPoints_.size()
            << " mappedPoints="
            << nMappedPoints
            << " missingPoints="
            << nMissingPoints
            << " seedFaces="
            << seamSeedFaces.size()
            << " inactiveFaces="
            << nInactiveFaces
            << " outOfRangeFaces="
            << nOutOfRangeFaces
            << " cap="
            << neutralSeamBackoffMaxLayers_
            << " inserted="
            << nInserted
            << " loweredExisting="
            << nLoweredExisting
            << " keptStricter="
            << nKeptStricter
            << " forcedFaceTotal="
            << forcedMaxLayersAtFace_.size()
            << endl;
    }


    // Provenance/junction BL retraction: cap selected boundary faces to
    // a maximum layer count. Generalizes the old binary
    // forceSingleLayerFaces_ logic into a tapered ring cap.
    // Applied after all global/patch layer-count logic so it always wins.
    actuallyCappedFaces_.clear();
    actuallyCappedProviderPlansAtFace_.clear();

    if( forcedMaxLayersAtFace_.size() )
    {
        label nCapped = 0;
        label nAlreadyAtOrBelow = 0;
        label nOutOfRange = 0;
        label nInvalidCap = 0;
        label nCappedNoProvider = 0;

        forAllConstIter(Map<label>, forcedMaxLayersAtFace_, it)
        {
            const label bfI = it.key();
            const label maxLayers = it();

            if( bfI < 0 || bfI >= label(nLayersAtBndFace_.size()) )
            {
                ++nOutOfRange;
                continue;
            }

            if( maxLayers < 0 )
            {
                ++nInvalidCap;
                continue;
            }

            if( nLayersAtBndFace_[bfI] > maxLayers )
            {
                nLayersAtBndFace_[bfI] = maxLayers;
                ++nCapped;

                actuallyCappedFaces_.insert(bfI);

                if( forcedMaxLayerProviderPlansAtFace_.found(bfI) )
                {
                    actuallyCappedProviderPlansAtFace_.insert
                    (
                        bfI,
                        forcedMaxLayerProviderPlansAtFace_[bfI]
                    );
                }
                else
                {
                    ++nCappedNoProvider;
                }
            }
            else
            {
                ++nAlreadyAtOrBelow;
            }
        }

        Info << "refineBoundaryLayers: forced max-layer caps applied: "
             << "capped=" << nCapped
             << " alreadyAtOrBelow=" << nAlreadyAtOrBelow
             << " outOfRange=" << nOutOfRange
             << " invalidCap=" << nInvalidCap
             << " requested=" << forcedMaxLayersAtFace_.size()
             << " effectiveFaces=" << actuallyCappedFaces_.size()
             << " providerFaces="
             << actuallyCappedProviderPlansAtFace_.size()
             << " cappedNoProvider=" << nCappedNoProvider
             << endl;
    }

    # ifdef DEBUGLayer
    forAll(nLayersAtBndFace_, bfI)
    Pout << "Boundary face " << bfI << " in patch "
        << facePatch[bfI] << " num layers " << nLayersAtBndFace_[bfI] << endl;
    //::exit(1);
    # endif

    return validLayer;
}

void refineBoundaryLayers::generateNewVertices()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& pointFaces = mse.pointFaces();
    const labelList& facePatch = mse.boundaryFacePatches();
    const labelList& bp = mse.bp();

    // CFMitch v0 architecture entry point.
    // The planner is behaviour-neutral at this stage.  The
    // constraintPlanner mode intentionally retains legacyEnhanced
    // execution so exact mesh equivalence can be established before
    // responsibilities are migrated into the planner.
    boundaryLayerConstraintPlanner cfmitchPlanner
    (
        mesh_,
        boundaryLayerArchitecture_
    );

    cfmitchPlanner.report();

    //- allocate the data from storing parameters applying to a split edge
    LongList<scalar> firstLayerThickness(splitEdges_.size());
    LongList<scalar> thicknessRatio(splitEdges_.size());
    labelLongList nNodesAtEdge(splitEdges_.size());
    labelLongList nLayersAtEdge(splitEdges_.size());

    // Per-split-edge thickness scale from BLRepairPlan.
    // Min scale wins if multiple capped faces touch the same edge.
    LongList<scalar> forcedThicknessScaleAtEdge(splitEdges_.size(), scalar(1.0));

    if( forcedThicknessScaleAtFace_.size() )
    {
        const faceList::subList& bFaces = mse.boundaryFaces();

        label nScaledFaces = 0;
        label nScaledEdges = 0;
        label nOutOfRangeFaces = 0;

        forAllConstIter(Map<scalar>, forcedThicknessScaleAtFace_, it)
        {
            const label bfI = it.key();
            const scalar scale =
                Foam::max(scalar(0.0), Foam::min(scalar(1.0), it()));

            if( bfI < 0 || bfI >= label(bFaces.size()) )
            {
                ++nOutOfRangeFaces;
                continue;
            }

            ++nScaledFaces;

            const face& bf = bFaces[bfI];

            forAll(bf, fpI)
            {
                const label bpI = bf[fpI];

                if( bpI < 0 || bpI >= label(splitEdgesAtPoint_.size()) )
                    continue;

                forAllRow(splitEdgesAtPoint_, bpI, seRowI)
                {
                    const label seI = splitEdgesAtPoint_(bpI, seRowI);

                    if( seI < 0 || seI >= label(forcedThicknessScaleAtEdge.size()) )
                        continue;

                    const scalar oldScale = forcedThicknessScaleAtEdge[seI];
                    const scalar newScale = Foam::min(oldScale, scale);

                    if( newScale < oldScale - SMALL )
                        ++nScaledEdges;

                    forcedThicknessScaleAtEdge[seI] = newScale;
                }
            }
        }

        //- Do the requested faces actually own split edges? If their
        //- points have no splitEdgesAtPoint_ rows, no cap and no scale can
        //- ever reach an edge, and the face is invisible to refinement.
        //- analyseLayers() only validates faces with faceInLayer >= 0, so
        //- provenance seed faces outside a detected layer are never checked.
        {
            label nFacesNoRows = 0, nPtsNoRows = 0, nPtsWithRows = 0;
            forAllConstIter(Map<scalar>, forcedThicknessScaleAtFace_, dit)
            {
                const label dbfI = dit.key();
                if( dbfI < 0 || dbfI >= label(bFaces.size()) ) continue;
                const face& dbf = bFaces[dbfI];
                bool anyRow = false;
                forAll(dbf, dpI)
                {
                    const label gp = dbf[dpI];
                    if( gp < 0 || gp >= label(splitEdgesAtPoint_.size()) )
                        { ++nPtsNoRows; continue; }
                    if( splitEdgesAtPoint_.sizeOfRow(gp) == 0 ) ++nPtsNoRows;
                    else { ++nPtsWithRows; anyRow = true; }
                }
                if( !anyRow ) ++nFacesNoRows;
            }
            Info << "refineBoundaryLayers: seed split-edge coverage:"
                 << " facesWithNoSplitEdges=" << nFacesNoRows
                 << " pointsWithRows=" << nPtsWithRows
                 << " pointsWithoutRows=" << nPtsNoRows
                 << " splitEdgesAtPointSize=" << splitEdgesAtPoint_.size()
                 << " splitEdges=" << splitEdges_.size() << endl;
        }

        Info << "refineBoundaryLayers: forced thickness scales prepared: "
             << "scaledFaces=" << nScaledFaces
             << " scaledEdgesUpdates=" << nScaledEdges
             << " outOfRangeFaces=" << nOutOfRangeFaces
             << " requestedFaces=" << forcedThicknessScaleAtFace_.size()
             << endl;
    }

    if( cfmitchPlanner.plannerEnabled() )
    {
        const boundaryLayerPlan cfmitchPlan =
            cfmitchPlanner.solveLayerCounts
            (
                mse,
                splitEdges_,
                neutralLayerScaleAtMeshPoint_,
                nLayersAtBndFace_,
                globalThicknessRatio_,
                thicknessRatioForPatch_,
                vtFaceRing_,
                actuallyCappedFaces_,
                constraintPlannerMaxLayerStep_
            );

        nLayersAtBndFace_ =
            cfmitchPlan.faceLayers();

        Info
            << "CFMITCH: resolved layer-count plan applied"
            << " faces=" << nLayersAtBndFace_.size()
            << endl;

        // -------------------------------------------------------------
        // CFMitch V3.8 -- pre-topology edge-hex compatibility planner.
        //
        // refineEdgeHexCell can refine a parent touched by two active
        // BL boundary faces only when each generating face has one
        // well-defined opposite quadrilateral terminal face.
        //
        // A six-face polyhedron is not necessarily a topological hex.
        // In particular, a polygonal opposite face can contain an
        // edge-conforming chain of N-1 intermediate vertices while
        // remaining a single face. Assigning that whole polygon to one
        // of N terminal children leaves exactly N-1 open child edges.
        //
        // Detect that condition using only the unmodified Q0 topology.
        // Both generating boundary faces are locally terminated at one
        // layer, then the normal constraint planner is re-run so the
        // requested layer count recovers smoothly according to
        // constraintPlannerMaxLayerStep_.
        // -------------------------------------------------------------
        {
            const labelList& structuralFaceOwners =
                mse.faceOwners();

            const cellListPMG& structuralCells =
                mesh_.cells();

            const faceListPMG& structuralFaces =
                mesh_.faces();

            const PtrList<boundaryPatch>& structuralBoundaries =
                mesh_.boundaries();

            const label structuralStartBoundary =
                structuralBoundaries.size()
              ? structuralBoundaries[0].patchStart()
              : structuralFaces.size();

            auto structuralCommonEdgeCount =
            [] (const face& a, const face& b) -> label
            {
                label count = 0;

                forAll(a, aeI)
                {
                    const label a0 = a[aeI];
                    const label a1 = a[(aeI+1) % a.size()];
                    const label alo = Foam::min(a0, a1);
                    const label ahi = Foam::max(a0, a1);

                    forAll(b, beI)
                    {
                        const label b0 = b[beI];
                        const label b1 = b[(beI+1) % b.size()];

                        if
                        (
                            alo == Foam::min(b0, b1)
                         && ahi == Foam::max(b0, b1)
                        )
                        {
                            ++count;
                        }
                    }
                }

                return count;
            };

            label structuralMaxLayers = 1;

            forAll(nLayersAtBndFace_, bfI)
            {
                structuralMaxLayers =
                    Foam::max
                    (
                        structuralMaxLayers,
                        nLayersAtBndFace_[bfI]
                    );
            }

            const label structuralMaxPasses =
                structuralMaxLayers + 2;

            labelHashSet structuralAllSeeds;
            label structuralPass = 0;
            bool structuralChanged = true;

            while
            (
                structuralChanged
             && structuralPass < structuralMaxPasses
            )
            {
                ++structuralPass;
                structuralChanged = false;

                labelList activeCount
                (
                    structuralCells.size(),
                    0
                );

                labelList activeBf0
                (
                    structuralCells.size(),
                    -1
                );

                labelList activeBf1
                (
                    structuralCells.size(),
                    -1
                );

                forAll(structuralFaceOwners, bfI)
                {
                    if
                    (
                        bfI < 0
                     || bfI >= label(nLayersAtBndFace_.size())
                     || nLayersAtBndFace_[bfI] <= 1
                    )
                    {
                        continue;
                    }

                    const label cellI =
                        structuralFaceOwners[bfI];

                    if
                    (
                        cellI < 0
                     || cellI >= label(structuralCells.size())
                    )
                    {
                        continue;
                    }

                    if( activeCount[cellI] == 0 )
                    {
                        activeBf0[cellI] = bfI;
                    }
                    else if( activeCount[cellI] == 1 )
                    {
                        activeBf1[cellI] = bfI;
                    }

                    ++activeCount[cellI];
                }

                labelHashSet structuralPassSeeds;

                label nType2Parents = 0;
                label nBadParentFaceCount = 0;
                label nBadActiveFaces = 0;
                label nBadActiveCommonEdges = 0;
                label nBadOppositeCount = 0;
                label nNonQuadOpposite = 0;

                forAll(activeCount, cellI)
                {
                    if( activeCount[cellI] != 2 )
                        continue;

                    ++nType2Parents;

                    const label bf0 = activeBf0[cellI];
                    const label bf1 = activeBf1[cellI];

                    bool compatible = true;

                    const cell& parent =
                        structuralCells[cellI];

                    if( parent.size() != 6 )
                    {
                        compatible = false;
                        ++nBadParentFaceCount;
                    }

                    const label activeFaceI0 =
                        structuralStartBoundary + bf0;

                    const label activeFaceI1 =
                        structuralStartBoundary + bf1;

                    if
                    (
                        compatible
                     &&
                        (
                            activeFaceI0 < 0
                         || activeFaceI1 < 0
                         || activeFaceI0 >= label(structuralFaces.size())
                         || activeFaceI1 >= label(structuralFaces.size())
                        )
                    )
                    {
                        compatible = false;
                        ++nBadActiveFaces;
                    }

                    bool parentHasActive0 = false;
                    bool parentHasActive1 = false;

                    if( compatible )
                    {
                        forAll(parent, localFaceI)
                        {
                            parentHasActive0 =
                                parentHasActive0
                             || parent[localFaceI] == activeFaceI0;

                            parentHasActive1 =
                                parentHasActive1
                             || parent[localFaceI] == activeFaceI1;
                        }

                        if
                        (
                            !parentHasActive0
                         || !parentHasActive1
                        )
                        {
                            compatible = false;
                            ++nBadActiveFaces;
                        }
                    }

                    if
                    (
                        compatible
                     && structuralCommonEdgeCount
                        (
                            structuralFaces[activeFaceI0],
                            structuralFaces[activeFaceI1]
                        ) != 1
                    )
                    {
                        compatible = false;
                        ++nBadActiveCommonEdges;
                    }

                    label oppositeFaceI0 = -1;
                    label oppositeFaceI1 = -1;
                    label nOpposite0 = 0;
                    label nOpposite1 = 0;

                    if( compatible )
                    {
                        const face& active0 =
                            structuralFaces[activeFaceI0];

                        const face& active1 =
                            structuralFaces[activeFaceI1];

                        forAll(parent, localFaceI)
                        {
                            const label candidateFaceI =
                                parent[localFaceI];

                            if
                            (
                                candidateFaceI < 0
                             || candidateFaceI >=
                                label(structuralFaces.size())
                            )
                            {
                                compatible = false;
                                ++nBadActiveFaces;
                                break;
                            }

                            if
                            (
                                candidateFaceI != activeFaceI0
                             && structuralCommonEdgeCount
                                (
                                    active0,
                                    structuralFaces[candidateFaceI]
                                ) == 0
                            )
                            {
                                oppositeFaceI0 = candidateFaceI;
                                ++nOpposite0;
                            }

                            if
                            (
                                candidateFaceI != activeFaceI1
                             && structuralCommonEdgeCount
                                (
                                    active1,
                                    structuralFaces[candidateFaceI]
                                ) == 0
                            )
                            {
                                oppositeFaceI1 = candidateFaceI;
                                ++nOpposite1;
                            }
                        }
                    }

                    if
                    (
                        compatible
                     &&
                        (
                            nOpposite0 != 1
                         || nOpposite1 != 1
                        )
                    )
                    {
                        compatible = false;
                        ++nBadOppositeCount;
                    }

                    if
                    (
                        compatible
                     &&
                        (
                            structuralFaces[oppositeFaceI0].size() != 4
                         || structuralFaces[oppositeFaceI1].size() != 4
                        )
                    )
                    {
                        compatible = false;
                        ++nNonQuadOpposite;
                    }

                    if( !compatible )
                    {
                        structuralPassSeeds.insert(bf0);
                        structuralPassSeeds.insert(bf1);
                    }
                }

                Info
                    << "CFMITCH V3.8 EDGEHEX PLAN GUARD:"
                    << " pass=" << structuralPass
                    << " type2Parents=" << nType2Parents
                    << " seedFaces=" << structuralPassSeeds.size()
                    << " badParentFaceCount="
                    << nBadParentFaceCount
                    << " badActiveFaces="
                    << nBadActiveFaces
                    << " badActiveCommonEdges="
                    << nBadActiveCommonEdges
                    << " badOppositeCount="
                    << nBadOppositeCount
                    << " nonQuadOpposite="
                    << nNonQuadOpposite
                    << endl;

                if( structuralPassSeeds.empty() )
                    break;

                structuralChanged = true;

                forAllConstIter
                (
                    labelHashSet,
                    structuralPassSeeds,
                    seedIt
                )
                {
                    const label bfI = seedIt.key();

                    if
                    (
                        bfI < 0
                     || bfI >= label(nLayersAtBndFace_.size())
                    )
                    {
                        continue;
                    }

                    nLayersAtBndFace_[bfI] = 1;
                    structuralAllSeeds.insert(bfI);
                }

                const boundaryLayerPlan structuralPlan =
                    cfmitchPlanner.solveLayerCounts
                    (
                        mse,
                        splitEdges_,
                        neutralLayerScaleAtMeshPoint_,
                        nLayersAtBndFace_,
                        globalThicknessRatio_,
                        thicknessRatioForPatch_,
                        vtFaceRing_,
                        structuralAllSeeds,
                        constraintPlannerMaxLayerStep_
                    );

                nLayersAtBndFace_ =
                    structuralPlan.faceLayers();

                // Structural termination is a hard upper bound. The
                // planner is expected to preserve it; enforce it again
                // so a future planner implementation cannot reactivate
                // an unsupported edge-hex parent.
                forAllConstIter
                (
                    labelHashSet,
                    structuralAllSeeds,
                    seedIt
                )
                {
                    const label bfI = seedIt.key();

                    if
                    (
                        bfI >= 0
                     && bfI < label(nLayersAtBndFace_.size())
                    )
                    {
                        nLayersAtBndFace_[bfI] = 1;
                    }
                }
            }

            Info
                << "CFMITCH V3.8 EDGEHEX PLAN SUMMARY:"
                << " passes=" << structuralPass
                << " structuralSeedFaces="
                << structuralAllSeeds.size()
                << " maxPasses="
                << structuralMaxPasses
                << endl;
        }
    }
    else
    {
        // Dense stable mesh-point lookup for BL/neutral height scales.
        // Built serially; the OMP split-edge loop below is read-only.
        scalarField neutralScaleByMeshPoint
        (
            mesh_.points().size(),
            scalar(1)
        );

        label nNeutralScaleMapped = 0;

        forAllConstIter
        (
            Map<scalar>,
            neutralLayerScaleAtMeshPoint_,
            it
        )
        {
            const label meshPtI = it.key();

            if
            (
                meshPtI < 0
             || meshPtI >= label(neutralScaleByMeshPoint.size())
            )
                continue;

            neutralScaleByMeshPoint[meshPtI] =
                Foam::max
                (
                    scalar(0),
                    Foam::min(scalar(1), it())
                );

            ++nNeutralScaleMapped;
        }

        Info << "refineBoundaryLayers: neutral scale metadata:"
             << " input=" << neutralLayerScaleAtMeshPoint_.size()
             << " mapped=" << nNeutralScaleMapped
             << endl;

        // BL/neutral height-aware FACE layer-count adaptation.
        //
        // neutralScaleByMeshPoint is keyed by the ORIGINAL stable wall point.
        // After BL extrusion that point is normally the interior endpoint of a
        // split hair; the current boundary faces contain the opposite/top
        // endpoint. Transfer the scale along the hair first, then modify the
        // canonical nLayersAtBndFace_ topology.
        scalarField neutralScaleByCurrentBp
        (
            mse.boundaryPoints().size(),
            scalar(1)
        );

        boolList neutralScalePresentAtCurrentBp
        (
            mse.boundaryPoints().size(),
            false
        );

        label nNeutralScaleTransferEdges = 0;
        label nNeutralScaleTransferBpUpdates = 0;

        forAll(splitEdges_, seI)
        {
            const edge& se = splitEdges_[seI];

            scalar edgeScale = scalar(1);
            bool hasEdgeScale = false;

            const label ep0 = se.start();
            const label ep1 = se.end();

            if
            (
                ep0 >= 0
             && ep0 < label(neutralScaleByMeshPoint.size())
             && neutralScaleByMeshPoint[ep0]
                    < scalar(1) - SMALL
            )
            {
                edgeScale = neutralScaleByMeshPoint[ep0];
                hasEdgeScale = true;
            }

            if
            (
                ep1 >= 0
             && ep1 < label(neutralScaleByMeshPoint.size())
             && neutralScaleByMeshPoint[ep1]
                    < scalar(1) - SMALL
            )
            {
                edgeScale =
                    hasEdgeScale
                  ? Foam::min
                    (
                        edgeScale,
                        neutralScaleByMeshPoint[ep1]
                    )
                  : neutralScaleByMeshPoint[ep1];

                hasEdgeScale = true;
            }

            if( !hasEdgeScale )
                continue;

            ++nNeutralScaleTransferEdges;

            // Transfer to whichever split-edge endpoint is a CURRENT
            // boundary point. Usually this is e.start(), but deliberately
            // handle both endpoints so the logic is orientation-independent.
            const label endpoints[2] = {ep0, ep1};

            for( label ei=0; ei<2; ++ei )
            {
                const label meshPtI = endpoints[ei];

                if
                (
                    meshPtI < 0
                 || meshPtI >= label(bp.size())
                )
                    continue;

                const label currBpI = bp[meshPtI];

                if
                (
                    currBpI < 0
                 || currBpI >= label(neutralScaleByCurrentBp.size())
                )
                    continue;

                if( neutralScalePresentAtCurrentBp[currBpI] )
                {
                    neutralScaleByCurrentBp[currBpI] =
                        Foam::min
                        (
                            neutralScaleByCurrentBp[currBpI],
                            edgeScale
                        );
                }
                else
                {
                    neutralScaleByCurrentBp[currBpI] =
                        edgeScale;

                    neutralScalePresentAtCurrentBp[currBpI] =
                        true;
                }

                ++nNeutralScaleTransferBpUpdates;
            }
        }

        label nNeutralScaleCappedFaces = 0;
        label minNeutralFaceLayers = labelMax;
        label maxNeutralFaceLayers = 0;

        // Seeds for topology compatibility propagation.
        // Only faces DIRECTLY reduced by the BL/neutral height rule may seed
        // this ramp. Existing low-layer gap/termination/BLBL faces must not
        // initiate propagation.
        boolList neutralScaleDirectCappedFace
        (
            nLayersAtBndFace_.size(),
            false
        );

        forAll(nLayersAtBndFace_, bfI)
        {
            const label requestedNLayers =
                nLayersAtBndFace_[bfI];

            if( requestedNLayers <= 1 )
                continue;

            if( bfI < 0 || bfI >= label(bFaces.size()) )
                continue;

            const face& bf = bFaces[bfI];

            scalar faceScale = scalar(1);
            bool hasFaceScale = false;

            forAll(bf, fpI)
            {
                const label meshPtI = bf[fpI];

                if
                (
                    meshPtI < 0
                 || meshPtI >= label(bp.size())
                )
                    continue;

                const label currBpI = bp[meshPtI];

                if
                (
                    currBpI < 0
                 || currBpI >= label(neutralScaleByCurrentBp.size())
                 || !neutralScalePresentAtCurrentBp[currBpI]
                )
                    continue;

                const scalar sPt =
                    neutralScaleByCurrentBp[currBpI];

                faceScale =
                    hasFaceScale
                  ? Foam::min(faceScale, sPt)
                  : sPt;

                hasFaceScale = true;
            }

            if( !hasFaceScale )
                continue;

            scalar ratio = globalThicknessRatio_;

            if
            (
                bfI >= 0
             && bfI < label(facePatch.size())
            )
            {
                const label patchI = facePatch[bfI];

                if
                (
                    patchI >= 0
                 && patchI < label(boundaries.size())
                )
                {
                    const word& patchName =
                        boundaries[patchI].patchName();

                    const std::map<word, scalar>::const_iterator rIt =
                        thicknessRatioForPatch_.find(patchName);

                    if( rIt != thicknessRatioForPatch_.end() )
                        ratio = rIt->second;
                }
            }

            const scalar r =
                Foam::max(ratio, scalar(1e-12));

            scalar totalWeight = scalar(0);
            scalar w = scalar(1);

            for
            (
                label li=0;
                li<requestedNLayers;
                ++li
            )
            {
                totalWeight += w;
                w *= r;
            }

            const scalar availableWeight =
                faceScale * totalWeight;

            label faceCap = 1;
            scalar cumulative = scalar(0);
            w = scalar(1);

            const scalar fitTol =
                scalar(100) * SMALL
              * Foam::max(scalar(1), totalWeight);

            for
            (
                label li=1;
                li<=requestedNLayers;
                ++li
            )
            {
                cumulative += w;

                if
                (
                    cumulative
                 <= availableWeight + fitTol
                )
                {
                    faceCap = li;
                }
                else
                {
                    break;
                }

                w *= r;
            }

            faceCap =
                Foam::max
                (
                    label(1),
                    Foam::min(requestedNLayers, faceCap)
                );

            if( faceCap < nLayersAtBndFace_[bfI] )
            {
                nLayersAtBndFace_[bfI] = faceCap;
                neutralScaleDirectCappedFace[bfI] = true;

                ++nNeutralScaleCappedFaces;

                minNeutralFaceLayers =
                    Foam::min
                    (
                        minNeutralFaceLayers,
                        faceCap
                    );

                maxNeutralFaceLayers =
                    Foam::max
                    (
                        maxNeutralFaceLayers,
                        faceCap
                    );
            }
        }

        // BL/neutral FACE layer-count compatibility ramp.
        //
        // The directly height-capped faces may differ sharply from neighboring
        // full-layer faces (for example 2 -> 15). Propagate ONLY outward from
        // those direct seeds and ONLY across faces on the same boundary patch.
        //
        // Enforce:
        //
        //     neighbourLayers <= sourceLayers + 1
        //
        // one topological ring per pass. This creates a deterministic
        // 2->3->4->...->15 recovery without allowing unrelated existing
        // one-layer faces to seed a global reduction.
        boolList neutralRampFrontier
        (
            neutralScaleDirectCappedFace
        );

        boolList neutralRampAdjustedFace
        (
            nLayersAtBndFace_.size(),
            false
        );

        label nNeutralRampAdjustedFaces = 0;
        label nNeutralRampUpdates = 0;
        label nNeutralRampPasses = 0;

        // The required number of rings cannot exceed the maximum layer-count
        // difference. Add a small guard margin rather than using an arbitrary
        // fixed iteration count.
        label neutralRampMaxLayers = 1;

        forAll(nLayersAtBndFace_, bfI)
        {
            neutralRampMaxLayers =
                Foam::max
                (
                    neutralRampMaxLayers,
                    nLayersAtBndFace_[bfI]
                );
        }

        const label neutralRampMaxPasses =
            neutralRampMaxLayers + 2;

        bool neutralRampChanged = true;

        while
        (
            neutralRampChanged
         && nNeutralRampPasses < neutralRampMaxPasses
        )
        {
            neutralRampChanged = false;
            ++nNeutralRampPasses;

            boolList nextFrontier
            (
                nLayersAtBndFace_.size(),
                false
            );

            // Snapshot the layer-count field so one pass advances exactly
            // one boundary-face ring regardless of face ordering.
            const labelList layersBeforePass
            (
                nLayersAtBndFace_
            );

            forAll(neutralRampFrontier, bfI)
            {
                if( !neutralRampFrontier[bfI] )
                    continue;

                if
                (
                    bfI < 0
                 || bfI >= label(bFaces.size())
                 || bfI >= label(facePatch.size())
                )
                    continue;

                const label sourcePatch =
                    facePatch[bfI];

                const label sourceLayers =
                    layersBeforePass[bfI];

                // Permit at most a two-layer change per topological ring.
                // This keeps the neutral transition smooth while reducing the
                // spatial footprint of the compatibility ramp.
                const label neighbourCap =
                    sourceLayers + 2;

                const face& f =
                    bFaces[bfI];

                forAll(f, fpI)
                {
                    const label meshPtI =
                        f[fpI];

                    if
                    (
                        meshPtI < 0
                     || meshPtI >= label(bp.size())
                    )
                        continue;

                    const label bpI =
                        bp[meshPtI];

                    if
                    (
                        bpI < 0
                     || bpI >= label(pointFaces.size())
                    )
                        continue;

                    forAllRow(pointFaces, bpI, pfI)
                    {
                        const label nbfI =
                            pointFaces(bpI, pfI);

                        if
                        (
                            nbfI < 0
                         || nbfI >= label(nLayersAtBndFace_.size())
                         || nbfI >= label(facePatch.size())
                         || nbfI == bfI
                        )
                            continue;

                        // Never cross patch boundaries.
                        if( facePatch[nbfI] != sourcePatch )
                            continue;

                        // Existing virtual-topology treatment has priority.
                        if
                        (
                            nbfI < label(vtFaceRing_.size())
                         && vtFaceRing_[nbfI] >= 0
                        )
                            continue;

                        if
                        (
                            layersBeforePass[nbfI]
                         <= neighbourCap
                        )
                            continue;

                        if
                        (
                            nLayersAtBndFace_[nbfI]
                         > neighbourCap
                        )
                        {
                            nLayersAtBndFace_[nbfI] =
                                neighbourCap;

                            nextFrontier[nbfI] = true;
                            neutralRampChanged = true;
                            ++nNeutralRampUpdates;

                            if
                            (
                                !neutralRampAdjustedFace[nbfI]
                            )
                            {
                                neutralRampAdjustedFace[nbfI] =
                                    true;

                                ++nNeutralRampAdjustedFaces;
                            }
                        }
                    }
                }
            }

            neutralRampFrontier.transfer
            (
                nextFrontier
            );
        }

        label nNeutralRampSeeds = 0;

        forAll(neutralScaleDirectCappedFace, bfI)
        {
            if( neutralScaleDirectCappedFace[bfI] )
                ++nNeutralRampSeeds;
        }

        Info << "BL/neutral FACE layer compatibility ramp:"
             << " maxStep=2"
             << " seedFaces=" << nNeutralRampSeeds
             << " adjustedFaces=" << nNeutralRampAdjustedFaces
             << " updates=" << nNeutralRampUpdates
             << " passes=" << nNeutralRampPasses
             << " maxPasses=" << neutralRampMaxPasses
             << endl;

        Info << "BL/neutral hair-scale transfer:"
             << " sourceEdges=" << nNeutralScaleTransferEdges
             << " bpUpdates=" << nNeutralScaleTransferBpUpdates
             << endl;

        Info << "BL/neutral height-aware FACE layer caps:"
             << " cappedFaces=" << nNeutralScaleCappedFaces;

        if( nNeutralScaleCappedFaces > 0 )
        {
            Info << " minLayers=" << minNeutralFaceLayers
                 << " maxLayers=" << maxNeutralFaceLayers;
        }

        Info << endl;

    }


    // ==================================================================
    // CFMITCH V10K CROSS-PATCH TYPE1 LAYER SYNC
    //
    // Production-architecture experiment.
    //
    // V10J proved that the surviving Rotor37 internal-skew family is born
    // where independently refined refType-1 BL parents from different wall
    // patches meet through an ORIGINAL internal lateral face with unequal
    // nLayersAtBndFace_.
    //
    // The downstream topology is deterministic:
    //
    //   - nLayersAtEdge takes the maximum demand from touching BL faces;
    //   - internal lateral quads likewise use max(nSplits0,nSplits1);
    //   - generateNewCellsPrism() clamps excess derived side faces onto
    //     local child zero.
    //
    // Therefore an unequal terminal layer count at this particular
    // interface produces the complex terminal polyhedron by construction.
    //
    // Reproduce the exact V10J refType definition here, BEFORE split-edge
    // depth is derived:
    //
    //     refType[cell] =
    //         number of owned boundary faces having nLayers > 1
    //
    // An original internal face is a cross-patch type1/type1 seam iff:
    //
    //     refType[owner] == 1
    //     refType[neighbour] == 1
    //     patch(ownerBaseFace) != patch(neighbourBaseFace)
    //
    // For eligible seams enforce equality monotonically downward.  Every
    // seam participating in a mismatch also seeds the existing-style +2
    // same-patch recovery constraint.  Equality and recovery are iterated
    // to a fixed point because a propagated reduction can reach another
    // cross-patch seam.
    //
    // Existing virtual-topology faces remain authoritative and are counted
    // but excluded from modification.
    //
    // Disabled by default.
    // ==================================================================

    if
    (
        cfmitchPlanner.plannerEnabled()
     && cfmitchV10KCrossPatchLayerSync_
     && !specialMode_
    )
    {
        const labelList& v10kBoundaryOwners =
            mse.faceOwners();

        const labelList& v10kOwner =
            mesh_.owner();

        const labelList& v10kNeighbour =
            mesh_.neighbour();

        labelList v10kRefType
        (
            mesh_.cells().size(),
            0
        );

        labelList v10kCellToBfI
        (
            mesh_.cells().size(),
            -1
        );

        // --------------------------------------------------------------
        // Reproduce generateNewCells() refType/cellToBfI exactly.
        // --------------------------------------------------------------
        forAll(v10kBoundaryOwners, bfI)
        {
            if
            (
                bfI < 0
             || bfI >= label(nLayersAtBndFace_.size())
            )
                continue;

            const label cellI =
                v10kBoundaryOwners[bfI];

            if
            (
                cellI < 0
             || cellI >= label(v10kRefType.size())
            )
                continue;

            if( nLayersAtBndFace_[bfI] > 1 )
            {
                ++v10kRefType[cellI];

                if( v10kCellToBfI[cellI] < 0 )
                    v10kCellToBfI[cellI] = bfI;
            }
        }

        // Store the exact original internal interface pairs.
        labelLongList v10kPairBf0;
        labelLongList v10kPairBf1;
        labelLongList v10kPairEligible;

        label v10kCandidateInterfaces = 0;
        label v10kProtectedInterfaces = 0;
        label v10kUnequalBefore = 0;
        label v10kEligibleUnequalBefore = 0;

        label v10kDelta2 = 0;
        label v10kDelta4 = 0;
        label v10kDelta6 = 0;
        label v10kDelta8 = 0;
        label v10kDeltaOther = 0;

        label v10kMaxDeltaBefore = 0;
        label v10kMaxLayerBefore = 1;

        forAll(nLayersAtBndFace_, bfI)
        {
            v10kMaxLayerBefore =
                Foam::max
                (
                    v10kMaxLayerBefore,
                    nLayersAtBndFace_[bfI]
                );
        }

        for
        (
            label faceI=0;
            faceI<label(v10kNeighbour.size());
            ++faceI
        )
        {
            if( faceI >= label(v10kOwner.size()) )
                continue;

            const label own =
                v10kOwner[faceI];

            const label nei =
                v10kNeighbour[faceI];

            if
            (
                own < 0
             || nei < 0
             || own >= label(v10kRefType.size())
             || nei >= label(v10kRefType.size())
            )
                continue;

            if
            (
                v10kRefType[own] != 1
             || v10kRefType[nei] != 1
            )
                continue;

            const label bf0 =
                v10kCellToBfI[own];

            const label bf1 =
                v10kCellToBfI[nei];

            if
            (
                bf0 < 0
             || bf1 < 0
             || bf0 >= label(nLayersAtBndFace_.size())
             || bf1 >= label(nLayersAtBndFace_.size())
             || bf0 >= label(facePatch.size())
             || bf1 >= label(facePatch.size())
            )
                continue;

            if( facePatch[bf0] == facePatch[bf1] )
                continue;

            ++v10kCandidateInterfaces;

            const label n0 =
                nLayersAtBndFace_[bf0];

            const label n1 =
                nLayersAtBndFace_[bf1];

            const label delta =
                (n0 >= n1) ? n0-n1 : n1-n0;

            if( delta > 0 )
            {
                ++v10kUnequalBefore;

                v10kMaxDeltaBefore =
                    Foam::max
                    (
                        v10kMaxDeltaBefore,
                        delta
                    );

                if( delta == 2 )
                    ++v10kDelta2;
                else if( delta == 4 )
                    ++v10kDelta4;
                else if( delta == 6 )
                    ++v10kDelta6;
                else if( delta == 8 )
                    ++v10kDelta8;
                else
                    ++v10kDeltaOther;
            }

            bool protectedInterface = false;

            if
            (
                bf0 < label(vtFaceRing_.size())
             && vtFaceRing_[bf0] >= 0
            )
            {
                protectedInterface = true;
            }

            if
            (
                bf1 < label(vtFaceRing_.size())
             && vtFaceRing_[bf1] >= 0
            )
            {
                protectedInterface = true;
            }

            if( protectedInterface )
            {
                ++v10kProtectedInterfaces;
            }
            else if( delta > 0 )
            {
                ++v10kEligibleUnequalBefore;
            }

            v10kPairBf0.append(bf0);
            v10kPairBf1.append(bf1);
            v10kPairEligible.append
            (
                protectedInterface ? label(0) : label(1)
            );
        }

        const labelList v10kLayersBeforeSync
        (
            nLayersAtBndFace_
        );

        boolList v10kFrontier
        (
            nLayersAtBndFace_.size(),
            false
        );

        boolList v10kSeamLoweredFace
        (
            nLayersAtBndFace_.size(),
            false
        );

        boolList v10kRampAdjustedFace
        (
            nLayersAtBndFace_.size(),
            false
        );

        label v10kEqualityPairUpdates = 0;
        label v10kEqualityFaceUpdates = 0;
        label v10kUniqueSeamFacesLowered = 0;

        label v10kRampUpdates = 0;
        label v10kUniqueRampFacesAdjusted = 0;

        label v10kFixedPointPasses = 0;

        // A deliberately generous finite guard.  The solve is monotonic:
        // every legal update lowers an integer count bounded below by 1.
        const label v10kMaxPasses =
            Foam::max
            (
                label(32),
                label(v10kPairBf0.size())
              + v10kMaxLayerBefore
              + label(4)
            );

        bool v10kConverged = false;

        for
        (
            label solvePass=0;
            solvePass<v10kMaxPasses;
            ++solvePass
        )
        {
            ++v10kFixedPointPasses;

            bool changedThisPass = false;

            // ----------------------------------------------------------
            // Equality constraints on actual type1/type1 cross-patch
            // original internal interfaces.
            // ----------------------------------------------------------
            forAll(v10kPairBf0, pairI)
            {
                if( !v10kPairEligible[pairI] )
                    continue;

                const label bf0 =
                    v10kPairBf0[pairI];

                const label bf1 =
                    v10kPairBf1[pairI];

                const label n0 =
                    nLayersAtBndFace_[bf0];

                const label n1 =
                    nLayersAtBndFace_[bf1];

                if( n0 == n1 )
                    continue;

                const label target =
                    Foam::min(n0, n1);

                ++v10kEqualityPairUpdates;

                if( nLayersAtBndFace_[bf0] > target )
                {
                    nLayersAtBndFace_[bf0] =
                        target;

                    ++v10kEqualityFaceUpdates;

                    if( !v10kSeamLoweredFace[bf0] )
                    {
                        v10kSeamLoweredFace[bf0] =
                            true;

                        ++v10kUniqueSeamFacesLowered;
                    }
                }

                if( nLayersAtBndFace_[bf1] > target )
                {
                    nLayersAtBndFace_[bf1] =
                        target;

                    ++v10kEqualityFaceUpdates;

                    if( !v10kSeamLoweredFace[bf1] )
                    {
                        v10kSeamLoweredFace[bf1] =
                            true;

                        ++v10kUniqueSeamFacesLowered;
                    }
                }

                // Seed BOTH sides of the seam.  The already-short side
                // must also have a smooth same-patch recovery away from
                // the newly authoritative equal-count interface.
                v10kFrontier[bf0] = true;
                v10kFrontier[bf1] = true;

                changedThisPass = true;
            }

            // ----------------------------------------------------------
            // One same-patch +2 propagation ring.
            //
            // This deliberately follows the existing compatibility-ramp
            // adjacency: boundary faces sharing a boundary point, same
            // patch only, with VT-protected targets skipped.
            // ----------------------------------------------------------
            const labelList layersBeforePass
            (
                nLayersAtBndFace_
            );

            boolList nextFrontier
            (
                nLayersAtBndFace_.size(),
                false
            );

            forAll(v10kFrontier, bfI)
            {
                if( !v10kFrontier[bfI] )
                    continue;

                if
                (
                    bfI < 0
                 || bfI >= label(bFaces.size())
                 || bfI >= label(facePatch.size())
                )
                    continue;

                const label sourcePatch =
                    facePatch[bfI];

                const label sourceLayers =
                    layersBeforePass[bfI];

                const label neighbourCap =
                    sourceLayers + 2;

                const face& f =
                    bFaces[bfI];

                forAll(f, fpI)
                {
                    const label meshPtI =
                        f[fpI];

                    if
                    (
                        meshPtI < 0
                     || meshPtI >= label(bp.size())
                    )
                        continue;

                    const label bpI =
                        bp[meshPtI];

                    if
                    (
                        bpI < 0
                     || bpI >= label(pointFaces.size())
                    )
                        continue;

                    forAllRow(pointFaces, bpI, pfI)
                    {
                        const label nbfI =
                            pointFaces(bpI, pfI);

                        if
                        (
                            nbfI < 0
                         || nbfI >=
                            label(nLayersAtBndFace_.size())
                         || nbfI >= label(facePatch.size())
                         || nbfI == bfI
                        )
                            continue;

                        if
                        (
                            facePatch[nbfI]
                         != sourcePatch
                        )
                            continue;

                        // Preserve existing VT authority.
                        if
                        (
                            nbfI < label(vtFaceRing_.size())
                         && vtFaceRing_[nbfI] >= 0
                        )
                            continue;

                        if
                        (
                            layersBeforePass[nbfI]
                         <= neighbourCap
                        )
                            continue;

                        if
                        (
                            nLayersAtBndFace_[nbfI]
                         > neighbourCap
                        )
                        {
                            nLayersAtBndFace_[nbfI] =
                                neighbourCap;

                            nextFrontier[nbfI] =
                                true;

                            ++v10kRampUpdates;

                            if
                            (
                                !v10kRampAdjustedFace[nbfI]
                            )
                            {
                                v10kRampAdjustedFace[nbfI] =
                                    true;

                                ++v10kUniqueRampFacesAdjusted;
                            }

                            changedThisPass = true;
                        }
                    }
                }
            }

            v10kFrontier.transfer
            (
                nextFrontier
            );

            if( !changedThisPass )
            {
                v10kConverged = true;
                break;
            }
        }

        // --------------------------------------------------------------
        // Final constraint audit.
        // --------------------------------------------------------------
        label v10kUnequalAfter = 0;
        label v10kEligibleUnequalAfter = 0;
        label v10kProtectedUnequalAfter = 0;
        label v10kMaxDeltaAfter = 0;

        forAll(v10kPairBf0, pairI)
        {
            const label bf0 =
                v10kPairBf0[pairI];

            const label bf1 =
                v10kPairBf1[pairI];

            const label n0 =
                nLayersAtBndFace_[bf0];

            const label n1 =
                nLayersAtBndFace_[bf1];

            const label delta =
                (n0 >= n1) ? n0-n1 : n1-n0;

            if( delta <= 0 )
                continue;

            ++v10kUnequalAfter;

            v10kMaxDeltaAfter =
                Foam::max
                (
                    v10kMaxDeltaAfter,
                    delta
                );

            if( v10kPairEligible[pairI] )
                ++v10kEligibleUnequalAfter;
            else
                ++v10kProtectedUnequalAfter;
        }

        label v10kTotalLayerReduction = 0;
        label v10kChangedFaces = 0;
        label v10kMinFinalActiveLayers = labelMax;

        forAll(nLayersAtBndFace_, bfI)
        {
            if
            (
                nLayersAtBndFace_[bfI]
             < v10kLayersBeforeSync[bfI]
            )
            {
                ++v10kChangedFaces;

                v10kTotalLayerReduction +=
                    v10kLayersBeforeSync[bfI]
                  - nLayersAtBndFace_[bfI];
            }

            if( nLayersAtBndFace_[bfI] > 1 )
            {
                v10kMinFinalActiveLayers =
                    Foam::min
                    (
                        v10kMinFinalActiveLayers,
                        nLayersAtBndFace_[bfI]
                    );
            }
        }

        if( v10kMinFinalActiveLayers == labelMax )
            v10kMinFinalActiveLayers = 0;

        Info
            << "CFMITCH V10K CROSS_PATCH_SYNC:"
            << " candidateInterfaces="
            << v10kCandidateInterfaces
            << " protectedInterfaces="
            << v10kProtectedInterfaces
            << " unequalBefore="
            << v10kUnequalBefore
            << " eligibleUnequalBefore="
            << v10kEligibleUnequalBefore
            << " maxDeltaBefore="
            << v10kMaxDeltaBefore
            << " delta2="
            << v10kDelta2
            << " delta4="
            << v10kDelta4
            << " delta6="
            << v10kDelta6
            << " delta8="
            << v10kDelta8
            << " deltaOther="
            << v10kDeltaOther
            << " equalityPairUpdates="
            << v10kEqualityPairUpdates
            << " equalityFaceUpdates="
            << v10kEqualityFaceUpdates
            << " seamFacesLowered="
            << v10kUniqueSeamFacesLowered
            << " rampAdjustedFaces="
            << v10kUniqueRampFacesAdjusted
            << " rampUpdates="
            << v10kRampUpdates
            << " fixedPointPasses="
            << v10kFixedPointPasses
            << " maxPasses="
            << v10kMaxPasses
            << " changedFaces="
            << v10kChangedFaces
            << " totalLayerReduction="
            << v10kTotalLayerReduction
            << " minFinalActiveLayers="
            << v10kMinFinalActiveLayers
            << " unequalAfter="
            << v10kUnequalAfter
            << " eligibleUnequalAfter="
            << v10kEligibleUnequalAfter
            << " protectedUnequalAfter="
            << v10kProtectedUnequalAfter
            << " maxDeltaAfter="
            << v10kMaxDeltaAfter
            << " converged="
            << v10kConverged
            << endl;

        if
        (
            !v10kConverged
         || v10kEligibleUnequalAfter != 0
        )
        {
            refinementValid_ = false;

            FatalErrorIn
            (
                "refineBoundaryLayers::generateNewVertices()"
            )
                << "CFMITCH V10K cross-patch layer synchronization "
                << "failed to reach its fixed-point constraints."
                << " converged=" << v10kConverged
                << " eligibleUnequalAfter="
                << v10kEligibleUnequalAfter
                << " passes=" << v10kFixedPointPasses
                << " maxPasses=" << v10kMaxPasses
                << exit(FatalError);
        }
    }

    //- count the number of vertices for each split edge
    # ifdef USE_OMP
    const label nThreads = omp_get_max_threads();
    # else
    const label nThreads = 1;
    # endif

    boolList cappedEdge(splitEdges_.size(), false);
    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        //- start counting vertices at each thread
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            //- get the requested number of boundary layers
            label nLayers(1);
            scalar ratio(globalThicknessRatio_);
            scalar thickness(globalMaxThicknessFirstLayer_);
            bool overridenThickness(false);

            const label bpI = bp[e.start()];

            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                const label pos = help::positionOfEdgeInFace(e, bFaces[bfI]);
                if( pos >= 0 )
                    continue;

                const word& patchName =
                    boundaries[facePatch[bfI]].patchName();

                //- overrride the global value with the maximum number of layers
                //- at this edge
                nLayers = Foam::max(nLayers, nLayersAtBndFace_[bfI]);

                //- override with the maximum ratio
                const std::map<word, scalar>::const_iterator rIt =
                    thicknessRatioForPatch_.find(patchName);
                if( rIt != thicknessRatioForPatch_.end() )
                {
                    ratio = rIt->second;
                }

                //- override with the minimum thickness set for this edge
                const std::map<word, scalar>::const_iterator tIt =
                    maxThicknessForPatch_.find(patchName);
                if( tIt != maxThicknessForPatch_.end() )
                {
                    if( overridenThickness )
                    {
                        thickness = Foam::min(thickness, tIt->second);
                    }
                    else
                    {
                        thickness = tIt->second;
                        overridenThickness = true;
                    }
                }
            }

            // Option B split-edge cap: local gap-zone + loser-side check.
            // Checks BOTH edge endpoints (orientation-independent).
            // Uses mesh point labels (stable) + patch names (stable).
            // Only caps where edge touches gap action zone AND loser patch.
            // Minimum effective cap is 1 (ring0 already topology-suppressed).
            if( gapActionPoints_.size() > 0 && gapLoserPatchNames_.size() > 0 )
            {
                // Gap action points are always boundary surface points (e.start()).
                // Interior BL points (e.end()) are never in gapActionPoints_.
                const bool edgeInGap =
                    gapActionPoints_.found(e.start());
                if( edgeInGap )
                {
                    label edgeCap = -1;
                    forAllRow(pointFaces, bpI, pfI)
                    {
                        const label bfI = pointFaces(bpI, pfI);
                        if( bfI < 0 || bfI >= label(facePatch.size()) ) continue;
                        const word& pName =
                            boundaries[facePatch[bfI]].patchName();
                        bool isLoser = false;
                        forAll(gapLoserPatchNames_, pi)
                            if( gapLoserPatchNames_[pi] == pName )
                                { isLoser = true; break; }
                        if( !isLoser ) continue;
                        const label cap = Foam::max(label(1), gapRing1MaxLayers_);
                        edgeCap = (edgeCap < 0) ?
                            cap : Foam::min(edgeCap, cap);
                    }
                    if( edgeCap > 0 && edgeCap < nLayers )
                    {
                        nLayers = edgeCap;
                        cappedEdge[seI] = true;
                    }
                }
            }

            // BL/termination edge cap: reduce layer count at inlet/outlet junctions.
            // Checks BOTH endpoints. No patch-name check needed -- blTerminationEdgePoints_
            // already contains only BL/no-BL transition points by construction.
            // Default cap=3 means disabled (min(3,3)=3, no effect).
            // Start conservative: blTerminationRing1MaxLayers=2 reduces 3->2 only.
            if( blTerminationEdgePoints_.size() > 0
             && blTerminationRing1MaxLayers_ < 3 )
            {
                const bool edgeAtTermination =
                    blTerminationEdgePoints_.found(e.start())
                 || blTerminationEdgePoints_.found(e.end());
                if( edgeAtTermination )
                {
                    const label cap =
                        Foam::max(label(1), blTerminationRing1MaxLayers_);
                    if( cap < nLayers )
                    {
                        nLayers = cap;
                        cappedEdge[seI] = true;
                    }
                }
            }

            //- store the information
            firstLayerThickness[seI] = thickness;

            if( forcedThicknessScaleAtEdge[seI] < scalar(1.0) - SMALL )
                firstLayerThickness[seI] *= forcedThicknessScaleAtEdge[seI];

            thicknessRatio[seI] = ratio;
            nLayersAtEdge[seI] = nLayers;

            if( !specialMode_ )
            {
                nNodesAtEdge[seI] = nLayers + 1;
            }
            else
            {
                nNodesAtEdge[seI] = 3;
            }
        }
    }

    // -------------------------------------------------------------
    // CFMitch diagnostic: BL geometric-series fit audit.
    //
    // Behaviour-neutral.  This does NOT modify first-layer thickness,
    // growth ratio, layer count, split-edge topology, or point positions.
    //
    // The downstream split-edge generator uses
    //
    //   actualH1 = min(requestedH1, hairLength / seriesSum)
    //
    // where seriesSum = sum(r^i, i=0..nLayers-1).
    //
    // Quantify how often the pre-existing hair edge is too short for the
    // requested BL specification, and distinguish:
    //
    //   - requested progression already fits;
    //   - h1+nLayers fit if growth ratio were reduced toward 1;
    //   - even nLayers*h1 does not fit.
    // -------------------------------------------------------------
    if( cfmitchPlanner.plannerEnabled() && !specialMode_ )
    {
        label nAudited(0);
        label nFitsRequested(0);
        label nWouldReduceRatio(0);
        label nTooShortEvenR1(0);
        label nLegacyCompressed(0);
        label nInvalid(0);

        label nH1Ge095(0);
        label nH1075To095(0);
        label nH1050To075(0);
        label nH1025To050(0);
        label nH1Lt025(0);

        scalar minLegacyH1Ratio(GREAT);
        scalar minHairOverH1(GREAT);
        scalar maxRequestedOverHair(-GREAT);

        scalar minCandidateRatio(GREAT);
        scalar maxCandidateRatio(-GREAT);

        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            const scalar h1 = firstLayerThickness[seI];
            const scalar rRequested = thicknessRatio[seI];
            const label nLayers = nLayersAtEdge[seI];

            if
            (
                h1 <= SMALL
             || nLayers <= 0
             || e.start() < 0
             || e.end() < 0
             || e.start() >= label(points.size())
             || e.end() >= label(points.size())
            )
            {
                ++nInvalid;
                continue;
            }

            const scalar hairLength =
                mag(points[e.end()] - points[e.start()]);

            if( hairLength <= SMALL )
            {
                ++nInvalid;
                continue;
            }

            ++nAudited;

            scalar seriesSum(1.0);
            scalar term(1.0);

            if( nLayers > 1 )
            {
                for(label layerI=1; layerI<nLayers; ++layerI)
                {
                    term *= rRequested;
                    seriesSum += term;
                }
            }

            const scalar requestedTotal = h1 * seriesSum;
            const scalar minTotal = scalar(nLayers) * h1;

            const scalar legacyFitH1 =
                hairLength / Foam::max(seriesSum, VSMALL);

            const scalar legacyActualH1 =
                Foam::min(h1, legacyFitH1);

            const scalar h1Ratio =
                legacyActualH1 / Foam::max(h1, VSMALL);

            minLegacyH1Ratio =
                Foam::min(minLegacyH1Ratio, h1Ratio);

            minHairOverH1 =
                Foam::min
                (
                    minHairOverH1,
                    hairLength / h1
                );

            maxRequestedOverHair =
                Foam::max
                (
                    maxRequestedOverHair,
                    requestedTotal / hairLength
                );

            if( h1Ratio >= scalar(0.95) )
                ++nH1Ge095;
            else if( h1Ratio >= scalar(0.75) )
                ++nH1075To095;
            else if( h1Ratio >= scalar(0.50) )
                ++nH1050To075;
            else if( h1Ratio >= scalar(0.25) )
                ++nH1025To050;
            else
                ++nH1Lt025;

            const scalar fitTol =
                scalar(1e-10)
              * Foam::max(hairLength, requestedTotal);

            if( requestedTotal <= hairLength + fitTol )
            {
                ++nFitsRequested;
            }
            else
            {
                ++nLegacyCompressed;

                if( minTotal <= hairLength + fitTol )
                {
                    ++nWouldReduceRatio;

                    // Diagnostic only: determine the largest ratio in
                    // [1,rRequested] that would preserve h1+nLayers.
                    scalar candidateRatio(scalar(1.0));

                    if( rRequested > scalar(1.0) + SMALL )
                    {
                        scalar lo(scalar(1.0));
                        scalar hi(rRequested);

                        for(label iter=0; iter<50; ++iter)
                        {
                            const scalar mid =
                                scalar(0.5) * (lo + hi);

                            scalar candidateSum(1.0);
                            scalar candidateTerm(1.0);

                            for(label layerI=1; layerI<nLayers; ++layerI)
                            {
                                candidateTerm *= mid;
                                candidateSum += candidateTerm;
                            }

                            if( h1*candidateSum <= hairLength )
                                lo = mid;
                            else
                                hi = mid;
                        }

                        candidateRatio = lo;
                    }

                    minCandidateRatio =
                        Foam::min(minCandidateRatio, candidateRatio);

                    maxCandidateRatio =
                        Foam::max(maxCandidateRatio, candidateRatio);
                }
                else
                {
                    ++nTooShortEvenR1;
                }
            }
        }

        Info << "CFMITCH BL HEIGHT FIT AUDIT:"
             << " audited=" << nAudited
             << " requestedFits=" << nFitsRequested
             << " legacyCompressed=" << nLegacyCompressed
             << " wouldReduceRatio=" << nWouldReduceRatio
             << " tooShortEvenR1=" << nTooShortEvenR1
             << " invalid=" << nInvalid
             << " h1Ratio_ge0.95=" << nH1Ge095
             << " h1Ratio_0.75_0.95=" << nH1075To095
             << " h1Ratio_0.50_0.75=" << nH1050To075
             << " h1Ratio_0.25_0.50=" << nH1025To050
             << " h1Ratio_lt0.25=" << nH1Lt025;

        if( nAudited > 0 )
        {
            Info << " minLegacyH1Ratio=" << minLegacyH1Ratio
                 << " minHairOverH1=" << minHairOverH1
                 << " maxRequestedOverHair=" << maxRequestedOverHair;
        }

        if( nWouldReduceRatio > 0 )
        {
            Info << " minCandidateRatio=" << minCandidateRatio
                 << " maxCandidateRatio=" << maxCandidateRatio;
        }

        Info << endl;
    }

    // -------------------------------------------------------------
    // CFMitch V5.4 HYBRID HEIGHT POLICY PRE-SMOOTH
    //
    // Preserve the effective requested first-layer height while avoiding
    // both pathological growth-ratio flattening and unnecessary layer loss.
    //
    // Policy priority:
    //
    //   1. preserve h1
    //   2. preserve nLayers when possible
    //   3. allow growth-ratio backoff to a preferred floor
    //   4. if that floor cannot fit, reduce nLayers
    //
    // Layer-count reductions happen HERE, before the already-existing
    // split-edge compatibility smoother.  This is important: neighboring
    // columns are therefore not left with uncontrolled layer-count cliffs.
    //
    // thicknessRatio is deliberately NOT modified here.  The exact final
    // effective ratio is solved later, after layer smoothing and MPI
    // reconciliation.
    // -------------------------------------------------------------

    const scalar cfmitchHybridMinGrowthRatio(1.15);

    // -------------------------------------------------------------
    // CFMitch V5.4b -- canonical topology / height-only PRE pass.
    //
    // IMPORTANT ARCHITECTURAL CONTRACT:
    //
    //   nLayersAtBndFace_ owns boundary-layer topology.
    //
    // nLayersAtEdge is derived support geometry for that topology and
    // must not independently reduce the number of intervals after the
    // face-layer plan has been solved.
    //
    // Previous V5.4 behaviour reduced nLayersAtEdge when the requested
    // h1 + growth ratio would not fit a hair at the preferred ratio
    // floor.  That created split rows shorter than the still-canonical
    // face topology, e.g.
    //
    //     face topology = 14
    //     split row     = 8 layers / 9 points
    //
    // generateNewFaces() then legitimately failed because canonical
    // face endpoints did not exist in the shortened split row.
    //
    // V5.4b therefore performs NO topology mutation here.
    //
    // Policy:
    //
    //   1. preserve canonical layer count
    //   2. preserve requested h1 when geometrically possible
    //   3. reduce growth ratio in POST
    //   4. preferred ratio floor 1.15 is soft
    //   5. if even r=1 cannot fit n*h1, downstream legacy placement
    //      may compress h1 for that exceptional edge
    //
    // This PRE pass is classification/telemetry only.
    // -------------------------------------------------------------
    if( cfmitchPlanner.plannerEnabled() && !specialMode_ )
    {
        label nCanonicalAudited(0);
        label nCanonicalRequestedFits(0);
        label nCanonicalAboveFloorBackoff(0);
        label nCanonicalBelowFloorBackoff(0);
        label nCanonicalCannotPreserveH1(0);
        label nCanonicalInvalid(0);

        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            const scalar h1 =
                firstLayerThickness[seI];

            const scalar requestedRatio =
                thicknessRatio[seI];

            const label nLayers =
                nLayersAtEdge[seI];

            if
            (
                h1 <= SMALL
             || nLayers <= 0
             || e.start() < 0
             || e.end() < 0
             || e.start() >= label(points.size())
             || e.end() >= label(points.size())
            )
            {
                ++nCanonicalInvalid;
                continue;
            }

            const scalar hairLength =
                mag(points[e.end()] - points[e.start()]);

            if( hairLength <= SMALL )
            {
                ++nCanonicalInvalid;
                continue;
            }

            ++nCanonicalAudited;

            scalar requestedSum(1.0);
            scalar requestedTerm(1.0);

            for(label layerI=1; layerI<nLayers; ++layerI)
            {
                requestedTerm *= requestedRatio;
                requestedSum += requestedTerm;
            }

            const scalar requestedTotal =
                h1 * requestedSum;

            const scalar fitTol =
                scalar(1e-10)
              * Foam::max(hairLength, requestedTotal);

            if( requestedTotal <= hairLength + fitTol )
            {
                ++nCanonicalRequestedFits;
                continue;
            }

            // Absolute geometric lower bound for preserved h1:
            //
            // r = 1  => total thickness = nLayers * h1.
            const scalar minimumTotal =
                scalar(nLayers) * h1;

            if( minimumTotal > hairLength + fitTol )
            {
                ++nCanonicalCannotPreserveH1;
                continue;
            }

            const scalar preferredFloor =
                Foam::max
                (
                    scalar(1.0),
                    Foam::min
                    (
                        requestedRatio,
                        cfmitchHybridMinGrowthRatio
                    )
                );

            scalar floorSum(1.0);
            scalar floorTerm(1.0);

            for(label layerI=1; layerI<nLayers; ++layerI)
            {
                floorTerm *= preferredFloor;
                floorSum += floorTerm;
            }

            if( h1 * floorSum <= hairLength + fitTol )
            {
                ++nCanonicalAboveFloorBackoff;
            }
            else
            {
                // Full canonical topology and h1 still fit, but only with
                // r below the preferred 1.15 floor. POST will solve it.
                ++nCanonicalBelowFloorBackoff;
            }
        }

        Info
            << "CFMITCH V5.4b CANONICAL HEIGHT PRE:"
            << " audited=" << nCanonicalAudited
            << " requestedFits=" << nCanonicalRequestedFits
            << " aboveFloorBackoff="
            << nCanonicalAboveFloorBackoff
            << " belowFloorBackoff="
            << nCanonicalBelowFloorBackoff
            << " cannotPreserveH1="
            << nCanonicalCannotPreserveH1
            << " invalid=" << nCanonicalInvalid
            << " preferredMinRatio="
            << cfmitchHybridMinGrowthRatio
            << endl;
    }

    // Split-edge layer-count compatibility smoothing.
    // Only runs when a cap actually reduced at least one split edge.
    // Prevents 3->1 or 3->2 isolated topology cliffs.
    // Does not propagate on geometries with intentional layer differences.
    {
        label nInitiallyCapped = 0;
        forAll(cappedEdge, seI)
            if( cappedEdge[seI] )
                ++nInitiallyCapped;

        if( nInitiallyCapped > 0 )
        {
            const label nSplitEdges = splitEdges_.size();
            bool changed = true;
            label nAdjusted = 0;
            label nPasses = 0;
            const label maxPasses = 20;

            while( changed && nPasses < maxPasses )
            {
                changed = false;
                ++nPasses;

                forAll(splitEdges_, seI)
                {
                    const edge& e = splitEdges_[seI];
                    const label nL = nLayersAtEdge[seI];

                    forAllRow(splitEdgesAtPoint_, e.start(), i)
                    {
                        const label seJ = splitEdgesAtPoint_(e.start(), i);
                        if( seJ == seI || seJ < 0 || seJ >= nSplitEdges )
                            continue;
                        if( nLayersAtEdge[seJ] > nL + 1 )
                        {
                            nLayersAtEdge[seJ] = nL + 1;
                            nNodesAtEdge[seJ] =
                                specialMode_ ? 3 : nLayersAtEdge[seJ] + 1;
                            changed = true;
                            ++nAdjusted;
                        }
                    }

                    forAllRow(splitEdgesAtPoint_, e.end(), i)
                    {
                        const label seJ = splitEdgesAtPoint_(e.end(), i);
                        if( seJ == seI || seJ < 0 || seJ >= nSplitEdges )
                            continue;
                        if( nLayersAtEdge[seJ] > nL + 1 )
                        {
                            nLayersAtEdge[seJ] = nL + 1;
                            nNodesAtEdge[seJ] =
                                specialMode_ ? 3 : nLayersAtEdge[seJ] + 1;
                            changed = true;
                            ++nAdjusted;
                        }
                    }
                }
            }

            Info << "Split-edge layer-count smoothing: "
                 << "initialCapped=" << nInitiallyCapped
                 << " adjusted=" << nAdjusted
                 << " passes=" << nPasses
                 << endl;
        }
    }

    if( Pstream::parRun() )
    {
        //- transfer the information over all processor for edges
        //- at inter-processor boundaries
        const labelLongList& globalEdgeLabel =
            mesh_.addressingData().globalEdgeLabel();
        const VRWGraph& edgeAtProcs = mesh_.addressingData().edgeAtProcs();
        const Map<label>& globalToLocal =
            mesh_.addressingData().globalToLocalEdgeAddressing();
        const DynList<label>& neiProcs = mesh_.addressingData().edgeNeiProcs();
        const edgeList& edges = mesh_.addressingData().edges();
        const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();

        //- exchange point number of layers
        std::map<label, LongList<labelPair> > exchangeNumLayers;
        std::map<label, LongList<labelPair> > exchangeNumNodesAtEdge;
        std::map<label, LongList<labelledScalar> > exchangeThickness;
        std::map<label, LongList<labelledScalar> > exchangeRatio;
        forAll(neiProcs, i)
        {
            exchangeNumNodesAtEdge.insert
            (
                std::make_pair(neiProcs[i], LongList<labelPair>())
            );
            exchangeNumLayers.insert
            (
                std::make_pair(neiProcs[i], LongList<labelPair>())
            );
            exchangeThickness.insert
            (
                std::make_pair(neiProcs[i], LongList<labelledScalar>())
            );
            exchangeRatio.insert
            (
                std::make_pair(neiProcs[i], LongList<labelledScalar>())
            );
        }

        //- exchange the number of layers
        forAll(splitEdges_, seI)
        {
            const edge& se = splitEdges_[seI];

            const label s = se.start();
            label edgeI(-1);
            forAllRow(pointEdges, s, peI)
            {
                const label eI = pointEdges(s, peI);

                if( edges[eI] == se )
                {
                    edgeI = eI;
                    break;
                }
            }

            const label geI = globalEdgeLabel[edgeI];

            if( globalToLocal.found(geI) )
            {
                forAllRow(edgeAtProcs, edgeI, i)
                {
                    const label neiProc = edgeAtProcs(edgeI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeNumNodesAtEdge[neiProc].append
                    (
                        labelPair(geI, nNodesAtEdge[seI])
                    );
                    exchangeNumLayers[neiProc].append
                    (
                        labelPair(geI, nLayersAtEdge[seI])
                    );
                    exchangeThickness[neiProc].append
                    (
                        labelledScalar(geI, firstLayerThickness[seI])
                    );
                    exchangeRatio[neiProc].append
                    (
                        labelledScalar(geI, thicknessRatio[seI])
                    );
                }
            }
        }

        //- exchange number of nodes at split edge
        LongList<labelPair> receivedNumLayers;
        help::exchangeMap(exchangeNumNodesAtEdge, receivedNumLayers);

        forAll(receivedNumLayers, i)
        {
            const labelPair& lp = receivedNumLayers[i];
            if( !globalToLocal.found(lp.first()) ) continue;
            const label eI = globalToLocal[lp.first()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            nNodesAtEdge[seI] = std::max(nNodesAtEdge[seI], lp.second());
        }

        //- exchange number of layers
        receivedNumLayers.clear();
        help::exchangeMap(exchangeNumLayers, receivedNumLayers);

        forAll(receivedNumLayers, i)
        {
            const labelPair& lp = receivedNumLayers[i];
            if( !globalToLocal.found(lp.first()) ) continue;
            const label eI = globalToLocal[lp.first()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            nLayersAtEdge[seI] = std::max(nLayersAtEdge[seI], lp.second());
        }

        //- exchange thickness ratio
        LongList<labelledScalar> receivedScalar;
        help::exchangeMap(exchangeRatio, receivedScalar);

        forAll(receivedScalar, i)
        {
            const labelledScalar& ls = receivedScalar[i];
            if( !globalToLocal.found(ls.scalarLabel()) ) continue;
            const label eI = globalToLocal[ls.scalarLabel()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            thicknessRatio[seI] = std::max(thicknessRatio[seI], ls.value());
        }

        //- exchange maximum thickness of the first layer
        receivedScalar.clear();
        help::exchangeMap(exchangeThickness, receivedScalar);

        forAll(receivedScalar, i)
        {
            const labelledScalar& ls = receivedScalar[i];
            if( !globalToLocal.found(ls.scalarLabel()) ) continue;
            const label eI = globalToLocal[ls.scalarLabel()];
            const edge& e = edges[eI];
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, e.start(), i)
            {
                const label seJ = splitEdgesAtPoint_(e.start(), i);
                if( splitEdges_[seJ] == e )
                {
                    seI = seJ;
                    break;
                }
            }
            firstLayerThickness[seI] =
                std::min(firstLayerThickness[seI], ls.value());
        }
    }

    // -------------------------------------------------------------
    // CFMitch V5.4b CANONICAL HEIGHT POLICY POST-MPI
    //
    // Layer topology remains owned by the canonical face-layer plan.
    // Any pre-existing edge compatibility handling and processor-boundary
    // reconciliation have now completed.
    //
    // Solve the final largest fitting growth ratio while preserving the
    // effective first-layer height and canonical supporting layer count.
    //
    // Preferred lower bound is cfmitchHybridMinGrowthRatio.  If a shared
    // edge or other reconciliation unexpectedly requires a lower value,
    // first-layer preservation wins and the below-floor recovery is
    // reported explicitly.
    // -------------------------------------------------------------
    if( cfmitchPlanner.plannerEnabled() && !specialMode_ )
    {
        label nPostAudited(0);
        label nPostUnchanged(0);
        label nPostRatioReduced(0);
        label nPostBelowFloorRecovery(0);
        label nPostCannotPreserveH1(0);
        label nPostInvalid(0);

        label nEffLt110(0);
        label nEff110To115(0);
        label nEff115To120(0);
        label nEff120To130(0);
        label nEff130To140(0);
        label nEffGe140(0);

        scalar minEffectiveRatio(GREAT);
        scalar maxEffectiveRatio(-GREAT);

        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            const scalar h1 = firstLayerThickness[seI];
            const scalar requestedRatio = thicknessRatio[seI];
            const label nLayers = nLayersAtEdge[seI];

            if
            (
                h1 <= SMALL
             || nLayers <= 0
             || e.start() < 0
             || e.end() < 0
             || e.start() >= label(points.size())
             || e.end() >= label(points.size())
            )
            {
                ++nPostInvalid;
                continue;
            }

            const scalar hairLength =
                mag(points[e.end()] - points[e.start()]);

            if( hairLength <= SMALL )
            {
                ++nPostInvalid;
                continue;
            }

            ++nPostAudited;

            scalar requestedSum(1.0);
            scalar requestedTerm(1.0);

            for(label layerI=1; layerI<nLayers; ++layerI)
            {
                requestedTerm *= requestedRatio;
                requestedSum += requestedTerm;
            }

            const scalar requestedTotal =
                h1 * requestedSum;

            const scalar fitTol =
                scalar(1e-10)
              * Foam::max(hairLength, requestedTotal);

            if( requestedTotal <= hairLength + fitTol )
            {
                ++nPostUnchanged;
                continue;
            }

            const scalar minimumTotal =
                scalar(nLayers) * h1;

            if( minimumTotal > hairLength + fitTol )
            {
                // Final reconciled layer count cannot preserve h1 even at
                // r=1.  Keep legacy behaviour for this edge and expose it.
                ++nPostCannotPreserveH1;
                continue;
            }

            const scalar preferredFloor =
                Foam::max
                (
                    scalar(1.0),
                    Foam::min
                    (
                        requestedRatio,
                        cfmitchHybridMinGrowthRatio
                    )
                );

            // Determine whether the preferred floor still fits after all
            // smoothing/MPI reconciliation.
            scalar floorSum(1.0);
            scalar floorTerm(1.0);

            for(label layerI=1; layerI<nLayers; ++layerI)
            {
                floorTerm *= preferredFloor;
                floorSum += floorTerm;
            }

            scalar lo(preferredFloor);
            scalar hi(requestedRatio);

            if( h1 * floorSum > hairLength + fitTol )
            {
                // Unexpected but recoverable.  Preserve h1, even if doing
                // so requires dropping below the preferred growth floor.
                lo = scalar(1.0);
                hi = preferredFloor;
                ++nPostBelowFloorRecovery;
            }

            // Largest fitting ratio.
            for(label iter=0; iter<60; ++iter)
            {
                const scalar mid =
                    scalar(0.5) * (lo + hi);

                scalar candidateSum(1.0);
                scalar candidateTerm(1.0);

                for(label layerI=1; layerI<nLayers; ++layerI)
                {
                    candidateTerm *= mid;
                    candidateSum += candidateTerm;
                }

                if( h1 * candidateSum <= hairLength )
                    lo = mid;
                else
                    hi = mid;
            }

            scalar effectiveRatio = lo;

            effectiveRatio =
                Foam::max
                (
                    scalar(1.0),
                    Foam::min(requestedRatio, effectiveRatio)
                );

            thicknessRatio[seI] = effectiveRatio;

            ++nPostRatioReduced;

            minEffectiveRatio =
                Foam::min(minEffectiveRatio, effectiveRatio);

            maxEffectiveRatio =
                Foam::max(maxEffectiveRatio, effectiveRatio);

            if( effectiveRatio < scalar(1.10) )
                ++nEffLt110;
            else if( effectiveRatio < scalar(1.15) )
                ++nEff110To115;
            else if( effectiveRatio < scalar(1.20) )
                ++nEff115To120;
            else if( effectiveRatio < scalar(1.30) )
                ++nEff120To130;
            else if( effectiveRatio < scalar(1.40) )
                ++nEff130To140;
            else
                ++nEffGe140;
        }

        Info << "CFMITCH V5.4b CANONICAL HEIGHT POST:"
             << " audited=" << nPostAudited
             << " unchanged=" << nPostUnchanged
             << " ratioReduced=" << nPostRatioReduced
             << " belowFloorRecovery=" << nPostBelowFloorRecovery
             << " cannotPreserveH1=" << nPostCannotPreserveH1
             << " invalid=" << nPostInvalid;

        if( nPostRatioReduced > 0 )
        {
            Info << " minEffectiveRatio=" << minEffectiveRatio
                 << " maxEffectiveRatio=" << maxEffectiveRatio
                 << " r_lt1.10=" << nEffLt110
                 << " r_1.10_1.15=" << nEff110To115
                 << " r_1.15_1.20=" << nEff115To120
                 << " r_1.20_1.30=" << nEff120To130
                 << " r_1.30_1.40=" << nEff130To140
                 << " r_ge1.40=" << nEffGe140;
        }

        Info << endl;
    }

    //- calculate the number of additional vertices which will be generated
    //- on edges of the mesh
    DynList<label> numPointsAtThread;
    numPointsAtThread.setSize(nThreads);
    numPointsAtThread = 0;

    # ifdef USE_OMP
    # pragma omp parallel for num_threads(nThreads) schedule(static, 1)
    # endif
    forAll(nNodesAtEdge, seI)
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        numPointsAtThread[threadI] += nNodesAtEdge[seI] - 2;
    }

    //- allocate the space in a graph storing ids of points on a split edge
    newVerticesForSplitEdge_.setSizeAndRowSize(nNodesAtEdge);

    //- calculate the number of points which will be generated
    //- on split edges
    label numPoints = points.size();
    forAll(numPointsAtThread, threadI)
    {
        const label nPts = numPointsAtThread[threadI];
        numPointsAtThread[threadI] = numPoints;
        numPoints += nPts;
    }

    points.setSize(numPoints);

    # ifdef DEBUGLayer
    Info << "Generating split vertices" << endl;
    # endif

    //- generate vertices on split edges
    # ifdef USE_OMP
    # pragma omp parallel num_threads(nThreads)
    # endif
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        label& nPoints = numPointsAtThread[threadI];

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            const vector v = e.vec(points);
            const scalar magv = mag(v);

            const label nLayers = newVerticesForSplitEdge_.sizeOfRow(seI) - 1;

            scalar firstThickness = magv / nLayersAtEdge[seI];
            if( thicknessRatio[seI] > (1. + SMALL) )
            {
                firstThickness =
                    magv /
                    (
                        (1 - Foam::pow(thicknessRatio[seI], nLayersAtEdge[seI]))
                        / (1.0 - thicknessRatio[seI])
                    );

                # ifdef DEBUGLayer
                Pout << "Thread " << threadI << endl;
                Pout << "Generating vertices at split edge "
                     << " start point " << points[e.start()]
                     << " end point " << points[e.end()] << endl;
                Pout << "Edge length " << magv << endl;
                Pout << "Thickness of the first layer "
                     << firstThickness << endl;
                # endif
            }

            firstThickness =
                Foam::min
                (
                    Foam::max(firstLayerThickness[seI], SMALL),
                    firstThickness
                );

            if( specialMode_ )
            {
                scalar t = firstThickness;

                for(label i=1;i<nLayersAtEdge[seI]-1;++i)
                    t += firstThickness * Foam::pow(thicknessRatio[seI], i);

                firstThickness = t;
            }

            //- generate vertices for this edge
            newVerticesForSplitEdge_(seI, 0) = e.start();

            // C3: degenerate hair edge (zero-length, BL/BL junction point)
            // Force all intermediate vertices to surface point to produce
            // clean zero-height wedge topology instead of machine-epsilon cells
            if( magv < 1000.0*SMALL )
            {
                for(label pI=1;pI<nLayers;++pI)
                    newVerticesForSplitEdge_(seI, pI) = e.start();
                newVerticesForSplitEdge_(seI, nLayers) = e.start();
                continue;
            }

            scalar param = firstThickness;
            const vector vec = v / (magv + VSMALL);

            for(label pI=1;pI<nLayers;++pI)
            {
                //- generate the new vertex
                const point newP = points[e.start()] + param * vec;

                # ifdef DEBUGLayer
                Pout << "Split edge " << seI << " edge points " << e
                    << " start point " << points[e.start()]
                    << " end point " << points[e.end()]
                    << " param " << param
                    << " new point " << nPoints
                    << " has coordinates " << newP << endl;
                # endif

                param += firstThickness * Foam::pow(thicknessRatio[seI], pI);

                newVerticesForSplitEdge_(seI, pI) = nPoints;
                points[nPoints++] = newP;
            }

            newVerticesForSplitEdge_(seI, nLayers) = e.end();
        }
    }

    // -------------------------------------------------------------
    // CFMitch V5.4 placement audit.
    //
    // Behaviour-neutral.  Measure the ACTUAL coordinates produced by
    // split-edge subdivision after all h1/ratio/layer-count policy has
    // completed.  This distinguishes planner intent from placed geometry.
    // -------------------------------------------------------------
    if( cfmitchPlanner.plannerEnabled() && !specialMode_ )
    {
        label nPlacementAudited(0);
        label nPlacementInvalid(0);
        label nPlacementZero(0);

        label nH1Ge095(0);
        label nH1075To095(0);
        label nH1050To075(0);
        label nH1025To050(0);
        label nH1Lt025(0);

        label nL1(0);
        label nL2To3(0);
        label nL4To6(0);
        label nL7To9(0);
        label nL10To13(0);
        label nL14(0);
        label nLGt14(0);

        scalar minRequestedH1(GREAT);
        scalar maxRequestedH1(-GREAT);
        scalar minActualH1(GREAT);
        scalar maxActualH1(-GREAT);
        scalar minActualOverRequested(GREAT);
        scalar maxActualOverRequested(-GREAT);

        scalar minHair(GREAT);
        scalar maxHair(-GREAT);

        label nPrinted(0);

        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            const label rowSize =
                newVerticesForSplitEdge_.sizeOfRow(seI);

            if
            (
                rowSize < 2
             || e.start() < 0
             || e.end() < 0
             || e.start() >= label(points.size())
             || e.end() >= label(points.size())
            )
            {
                ++nPlacementInvalid;
                continue;
            }

            const label firstPoint =
                newVerticesForSplitEdge_(seI, 1);

            if
            (
                firstPoint < 0
             || firstPoint >= label(points.size())
            )
            {
                ++nPlacementInvalid;
                continue;
            }

            const scalar requestedH1 =
                firstLayerThickness[seI];

            if( requestedH1 <= SMALL )
            {
                ++nPlacementInvalid;
                continue;
            }

            const scalar actualH1 =
                mag
                (
                    points[firstPoint]
                  - points[e.start()]
                );

            const scalar hairLength =
                mag
                (
                    points[e.end()]
                  - points[e.start()]
                );

            const label placedLayers = rowSize - 1;

            const scalar h1Ratio =
                actualH1 / requestedH1;

            ++nPlacementAudited;

            minRequestedH1 =
                Foam::min(minRequestedH1, requestedH1);
            maxRequestedH1 =
                Foam::max(maxRequestedH1, requestedH1);

            minActualH1 =
                Foam::min(minActualH1, actualH1);
            maxActualH1 =
                Foam::max(maxActualH1, actualH1);

            minActualOverRequested =
                Foam::min(minActualOverRequested, h1Ratio);
            maxActualOverRequested =
                Foam::max(maxActualOverRequested, h1Ratio);

            minHair =
                Foam::min(minHair, hairLength);
            maxHair =
                Foam::max(maxHair, hairLength);

            if( actualH1 <= SMALL )
                ++nPlacementZero;

            if( h1Ratio >= scalar(0.95) )
                ++nH1Ge095;
            else if( h1Ratio >= scalar(0.75) )
                ++nH1075To095;
            else if( h1Ratio >= scalar(0.50) )
                ++nH1050To075;
            else if( h1Ratio >= scalar(0.25) )
                ++nH1025To050;
            else
                ++nH1Lt025;

            if( placedLayers <= 1 )
                ++nL1;
            else if( placedLayers <= 3 )
                ++nL2To3;
            else if( placedLayers <= 6 )
                ++nL4To6;
            else if( placedLayers <= 9 )
                ++nL7To9;
            else if( placedLayers <= 13 )
                ++nL10To13;
            else if( placedLayers == 14 )
                ++nL14;
            else
                ++nLGt14;

            // Print a small deterministic sample of genuine BL columns.
            if
            (
                nPrinted < 30
             && placedLayers >= 7
            )
            {
                ++nPrinted;

                Info << "CFMITCH V5.4 H1DIAG:"
                     << " se=" << seI
                     << " requestedH1=" << requestedH1
                     << " actualH1=" << actualH1
                     << " actualOverRequested=" << h1Ratio
                     << " ratio=" << thicknessRatio[seI]
                     << " nLayers=" << nLayersAtEdge[seI]
                     << " rowLayers=" << placedLayers
                     << " hair=" << hairLength
                     << endl;
            }
        }

        Info << "CFMITCH V5.4 PLACEMENT AUDIT:"
             << " audited=" << nPlacementAudited
             << " invalid=" << nPlacementInvalid
             << " zeroH1=" << nPlacementZero
             << " h1_ge0.95=" << nH1Ge095
             << " h1_0.75_0.95=" << nH1075To095
             << " h1_0.50_0.75=" << nH1050To075
             << " h1_0.25_0.50=" << nH1025To050
             << " h1_lt0.25=" << nH1Lt025
             << " layers1=" << nL1
             << " layers2to3=" << nL2To3
             << " layers4to6=" << nL4To6
             << " layers7to9=" << nL7To9
             << " layers10to13=" << nL10To13
             << " layers14=" << nL14
             << " layersGt14=" << nLGt14;

        if( nPlacementAudited > 0 )
        {
            Info << " minRequestedH1=" << minRequestedH1
                 << " maxRequestedH1=" << maxRequestedH1
                 << " minActualH1=" << minActualH1
                 << " maxActualH1=" << maxActualH1
                 << " minActualOverRequested="
                 << minActualOverRequested
                 << " maxActualOverRequested="
                 << maxActualOverRequested
                 << " minHair=" << minHair
                 << " maxHair=" << maxHair;
        }

        Info << endl;
    }

    // -------------------------------------------------------------
    // CFMitch face-topology authority audit.
    //
    // Behaviour-neutral.
    //
    // nLayersAtBndFace_ is the layer-count field consumed downstream by
    // generateNewFaces()/generateNewCellsPrism().  Report its distribution
    // by patch after split-edge geometry has been generated.
    // -------------------------------------------------------------
    if( cfmitchPlanner.plannerEnabled() && !specialMode_ )
    {
        const meshSurfaceEngine& topoMse =
            surfaceEngine();

        const labelList& topoFacePatch =
            topoMse.boundaryFacePatches();

        const PtrList<boundaryPatch>& topoBoundaries =
            mesh_.boundaries();

        label g0(0), g1(0), g2(0), g3to6(0);
        label g7to9(0), g10to13(0), g14(0), gGt14(0);
        label gInvalid(0);

        forAll(nLayersAtBndFace_, bfI)
        {
            const label n = nLayersAtBndFace_[bfI];

            if( n <= 0 ) ++g0;
            else if( n == 1 ) ++g1;
            else if( n == 2 ) ++g2;
            else if( n <= 6 ) ++g3to6;
            else if( n <= 9 ) ++g7to9;
            else if( n <= 13 ) ++g10to13;
            else if( n == 14 ) ++g14;
            else ++gGt14;

            if
            (
                bfI < 0
             || bfI >= label(topoFacePatch.size())
             || topoFacePatch[bfI] < 0
             || topoFacePatch[bfI] >= label(topoBoundaries.size())
            )
                ++gInvalid;
        }

        Info
            << "CFMITCH FACE TOPOLOGY AUTHORITY:"
            << " faces=" << nLayersAtBndFace_.size()
            << " n0=" << g0
            << " n1=" << g1
            << " n2=" << g2
            << " n3to6=" << g3to6
            << " n7to9=" << g7to9
            << " n10to13=" << g10to13
            << " n14=" << g14
            << " nGt14=" << gGt14
            << " invalidPatchRefs=" << gInvalid
            << endl;

        forAll(topoBoundaries, patchI)
        {
            label p0(0), p1(0), p2(0), p3to6(0);
            label p7to9(0), p10to13(0), p14(0), pGt14(0);
            label total(0);

            forAll(nLayersAtBndFace_, bfI)
            {
                if
                (
                    bfI >= label(topoFacePatch.size())
                 || topoFacePatch[bfI] != patchI
                )
                    continue;

                ++total;

                const label n = nLayersAtBndFace_[bfI];

                if( n <= 0 ) ++p0;
                else if( n == 1 ) ++p1;
                else if( n == 2 ) ++p2;
                else if( n <= 6 ) ++p3to6;
                else if( n <= 9 ) ++p7to9;
                else if( n <= 13 ) ++p10to13;
                else if( n == 14 ) ++p14;
                else ++pGt14;
            }

            if( total > 0 )
            {
                Info
                    << "CFMITCH FACE TOPOLOGY PATCH:"
                    << " patch="
                    << topoBoundaries[patchI].patchName()
                    << " faces=" << total
                    << " n0=" << p0
                    << " n1=" << p1
                    << " n2=" << p2
                    << " n3to6=" << p3to6
                    << " n7to9=" << p7to9
                    << " n10to13=" << p10to13
                    << " n14=" << p14
                    << " nGt14=" << pGt14
                    << endl;
            }
        }
    }

    if( specialMode_ )
    {
        //- set the number of layers to 2
        forAll(nLayersAtBndFace_, bfI)
            if( nLayersAtBndFace_[bfI] > 1 )
                nLayersAtBndFace_[bfI] = 2;
    }

    # ifdef DEBUGLayer
    for(label procI=0;procI<Pstream::nProcs();++procI)
    {
        if( procI == Pstream::myProcNo() )
        {
            forAll(splitEdges_, seI)
            {
                Pout << "\nSplit edge " << seI << " nodes " << splitEdges_[seI]
                    << " coordinates " << points[splitEdges_[seI][0]]
                    << " " << points[splitEdges_[seI][1]]
                    << " has new points "
                    << newVerticesForSplitEdge_[seI] << endl;

                forAllRow(newVerticesForSplitEdge_, seI, i)
                    Pout << "Point " << i << " on edge ha coordinates "
                         << points[newVerticesForSplitEdge_(seI, i)] << endl;
            }
        }

        returnReduce(1, sumOp<label>());
    }

    Info << "Finished generating vertices at split edges" << endl;
    //::exit(1);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
