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

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "OFstream.H"
#include "helperFunctions.H"
#include "helperFunctionsPar.H"
#include "demandDrivenData.H"

#include "labelledPoint.H"
#include "labelledScalar.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGLayer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::findPatchVertices
(
    const boolList& treatPatches,
    List<direction>& pVertices
) const
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    pVertices.setSize(pPatches.size());
    pVertices = NONE;

    # ifdef USE_OMP
    # pragma omp parallel for if( pPatches.size() > 1000 ) \
    schedule(dynamic, Foam::max(10, pPatches.size()/(2*omp_get_num_threads())))
    # endif
    forAll(pPatches, bpI)
    {
        bool hasTreated(false);
        bool hasNotTreated(false);

        forAllRow(pPatches, bpI, patchI)
        {
            const label patch = pPatches(bpI, patchI);
            if( treatPatches[patch] )
            {
                hasTreated = true;
            }
            else
            {
                hasNotTreated = true;
            }
        }

        if( hasTreated )
        {
            pVertices[bpI] |= PATCHNODE;

            if( hasNotTreated )
                pVertices[bpI] |= EDGENODE;
        }
    }

    if( Pstream::parRun() )
    {
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        forAll(pVertices, bpI)
            if( pVertices[bpI] && (bpAtProcs.sizeOfRow(bpI) != 0) )
                pVertices[bpI] |= PARALLELBOUNDARY;
    }
}

point boundaryLayers::createNewVertex
(
    const label bpI,
    const boolList& treatPatches,
    const List<direction>& patchVertex
) const
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const vectorField& pNormals = mse.pointNormals();
    const VRWGraph& pFaces = mse.pointFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const VRWGraph& pointPoints = mse.pointPoints();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    const pointFieldPMG& points = mesh_.points();

    # ifdef DEBUGLayer
    Info << "Creating new vertex for boundary vertex " << bpI << endl;
    Info << "Global vertex label " << bPoints[bpI] << endl;
    # endif

    vector normal(vector::zero);
    scalar dist(VGREAT);
    const point& p = points[bPoints[bpI]];

    // ============================================================
    // CFMITCH BIRTH TARGET TRACE V1B
    //
    // Diagnostic only. No geometry is modified.
    //
    // Target original Rotor37 seam points:
    //   pointI 337809
    //   pointI 337818
    // ============================================================

    const label cfmBirthPointI = bPoints[bpI];

    const bool cfmBirthTarget =
    (
        cfmBirthPointI == 337809
     || cfmBirthPointI == 337818
    );

    if( cfmBirthTarget )
    {
        Info
            << "CFMITCH BIRTHTARGET_BEGIN"
            << " bpI=" << bpI
            << " pointI=" << cfmBirthPointI
            << " p=" << p
            << " patchVertex=" << label(patchVertex[bpI])
            << " neutral="
            << (blNeutralEdgePoints_.found(bpI) ? 1 : 0)
            << " blblCorner="
            << (blblCornerPoints_.found(bpI) ? 1 : 0)
            << " blblJunction="
            << (blblJunctionPoints_.found(bpI) ? 1 : 0)
            << " pointNormal=" << pNormals[bpI]
            << endl;

        Info
            << "CFMITCH BIRTHTARGET_PATCHES"
            << " bpI=" << bpI
            << " pointI=" << cfmBirthPointI
            << " nPatches=" << pPatches.sizeOfRow(bpI);

        forAllRow(pPatches, bpI, ppI)
        {
            const label patchI = pPatches(bpI, ppI);

            const bool treated =
            (
                patchI >= 0
             && patchI < label(treatPatches.size())
             && treatPatches[patchI]
            );

            Info
                << " [patchI=" << patchI
                << " treated=" << (treated ? 1 : 0)
                << "]";
        }

        Info << endl;

        vector cfmTreatedSum(vector::zero);
        vector cfmUntreatedSum(vector::zero);

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);

            if( bfI < 0 || bfI >= label(bFaces.size()) )
                continue;

            const face& f = bFaces[bfI];
            const label patchI = boundaryFacePatches[bfI];

            vector fn(vector::zero);

            if( f.size() >= 3 )
            {
                const point& fp0 = points[f[0]];

                for(label pi=1; pi<f.size()-1; ++pi)
                {
                    fn +=
                        (points[f[pi]] - fp0)
                      ^ (points[f[pi+1]] - fp0);
                }
            }

            const scalar fnMag = mag(fn);

            vector fnUnit(vector::zero);

            if( fnMag > VSMALL )
                fnUnit = fn/fnMag;

            const bool treated =
            (
                patchI >= 0
             && patchI < label(treatPatches.size())
             && treatPatches[patchI]
            );

            if( treated )
                cfmTreatedSum += fn;
            else
                cfmUntreatedSum += fn;

            Info
                << "CFMITCH BIRTHTARGET_FACE"
                << " bpI=" << bpI
                << " pointI=" << cfmBirthPointI
                << " bfI=" << bfI
                << " patchI=" << patchI
                << " treated=" << (treated ? 1 : 0)
                << " nPts=" << f.size()
                << " areaVec=" << fn
                << " areaMag=" << fnMag
                << " unitN=" << fnUnit
                << endl;
        }

        const scalar cfmTreatedMag = mag(cfmTreatedSum);
        const scalar cfmUntreatedMag = mag(cfmUntreatedSum);

        vector cfmTreatedUnit(vector::zero);
        vector cfmUntreatedUnit(vector::zero);

        if( cfmTreatedMag > VSMALL )
            cfmTreatedUnit = cfmTreatedSum/cfmTreatedMag;

        if( cfmUntreatedMag > VSMALL )
            cfmUntreatedUnit = cfmUntreatedSum/cfmUntreatedMag;

        Info
            << "CFMITCH BIRTHTARGET_SUM"
            << " bpI=" << bpI
            << " pointI=" << cfmBirthPointI
            << " treatedSum=" << cfmTreatedSum
            << " treatedUnit=" << cfmTreatedUnit
            << " untreatedSum=" << cfmUntreatedSum
            << " untreatedUnit=" << cfmUntreatedUnit
            << " treatedDotUntreated="
            << (cfmTreatedUnit & cfmUntreatedUnit)
            << " pointNormal=" << pNormals[bpI]
            << endl;
    }

    // BIRTHDIAG trackers (function scope)
    label  bd_branch = -1;   // 0=edge size1, 2=corner size2, 3=multipatch, 9=interior
    scalar bd_edgeDotNormal = -2.0;
    // BLNOBLPATHDIAG: classify which createNewVertex path the
    // BL/no-BL termination points actually take. Diagnostic only.
    if( blNoBlEdgePoints_.found(bpI) )
    {
        DynList<label> pathDiagOtherPatches;
        forAllRow(pPatches, bpI, patchI)
        {
            const label patchLabel = pPatches(bpI, patchI);
            if( !treatPatches[patchLabel] )
                pathDiagOtherPatches.appendIfNotIn(patchLabel);
        }

        label pathDiagTermRoleCount = 0;
        forAll(pathDiagOtherPatches, opI)
        {
            const label op = pathDiagOtherPatches[opI];
            if( op >= 0
             && op < label(patchRole_.size())
             && patchRole_[op] == 1 )
            {
                ++pathDiagTermRoleCount;
            }
        }

        const scalar lsThis =
            (layerScale_.size() > bpI) ? layerScale_[bpI] : scalar(1.0);

        Info << "BLNOBLPATHDIAG bpI=" << bpI
             << " x=" << p.x() << " y=" << p.y() << " z=" << p.z()
             << " edgeNode=" << ((patchVertex[bpI] & EDGENODE) ? 1 : 0)
             << " otherPatchSize=" << pathDiagOtherPatches.size()
             << " termRoleCount=" << pathDiagTermRoleCount
             << " layerScale=" << lsThis
             << " blblCorner=" << (blblCornerPoints_.found(bpI) ? 1 : 0)
             << " blblJunction=" << (blblJunctionPoints_.found(bpI) ? 1 : 0)
             << endl;
    }

    if( patchVertex[bpI] & EDGENODE )
    {
        # ifdef DEBUGLayer
        Info << "Vertex is on the border" << endl;
        # endif

        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn
                (
                    pPatches(bpI, patchI)
                );

        if( otherPatches.size() == 1 )
        {
            //- vertex is on an edge (or BL+BL+neutral corner)
            # ifdef DEBUGLayer
            Info << "Vertex is on an edge" << endl;
            # endif

            // BL+BL+neutral corner: use per-patch conservative normal
            // Do NOT apply layerScale_ here -- common code below handles it once
            if( blblCornerPoints_.found(bpI) )
            {
                Map<vector> patchNormals;
                forAllRow(pFaces, bpI, pfI)
                {
                    const label faceI = pFaces(bpI, pfI);
                    const label patchLabel = boundaryFacePatches[faceI];
                    if( patchLabel < 0 || patchLabel >= label(treatPatches.size()) ) continue;
                    if( !treatPatches[patchLabel] ) continue;
                    const face& f = bFaces[faceI];
                    if( f.size() < 3 ) continue;
                    vector fn = vector::zero;
                    const point& p0 = points[f[0]];
                    for(label pi=1; pi<f.size()-1; ++pi)
                        fn += (points[f[pi]] - p0) ^ (points[f[pi+1]] - p0);
                    if( patchNormals.found(patchLabel) )
                        patchNormals[patchLabel] += fn;
                    else
                        patchNormals.insert(patchLabel, fn);
                }
                if( patchNormals.size() >= 2 )
                {
                    scalar bestDist = GREAT;
                    vector bestNormal = vector::zero;
                    forAllConstIter(Map<vector>, patchNormals, it)
                    {
                        vector n = it();
                        const scalar magN = mag(n);
                        if( magN < VSMALL ) continue;
                        n /= magN;
                        scalar localDist = VGREAT;
                        forAllRow(pointPoints, bpI, ppI)
                        {
                            const label bpJ = pointPoints(bpI, ppI);
                            const vector vec = points[bPoints[bpJ]] - p;
                            const scalar d = 0.5 * mag(vec & n);
                            if( d < localDist ) localDist = d;
                        }
                        if( localDist < bestDist )
                        {
                            bestDist = localDist;
                            bestNormal = n;
                        }
                    }
                    if( mag(bestNormal) > VSMALL )
                    {
                        normal = bestNormal;
                        dist = bestDist;
                        // layerScale_ applied once by common code below
                    }
                    else
                        normal = pNormals[bpI];
                }
                else
                    normal = pNormals[bpI];
            }
            else
            //- zero-dist for BL-transition edge points
            if( terminateLayersAtConcaveEdges_
             && layerScale_.size() > bpI
             && layerScale_[bpI] < 0.01 )
            {
                dist = 0.0;
            }
            else
            {
            vector v(vector::zero);

            forAllRow(pFaces, bpI, pfI)
            {
                const face& f = bFaces[pFaces(bpI, pfI)];
                const label patchLabel =
                    boundaryFacePatches[pFaces(bpI, pfI)];

                if( treatPatches[patchLabel] )
                {
                    { vector _n=vector::zero; const point& _p0=points[f[0]]; for(label _pi=1;_pi<f.size()-1;++_pi) _n+=(points[f[_pi]]-_p0)^(points[f[_pi+1]]-_p0); normal += _n; }
                }
                else
                {
                    { vector _n=vector::zero; const point& _p0=points[f[0]]; for(label _pi=1;_pi<f.size()-1;++_pi) _n+=(points[f[_pi]]-_p0)^(points[f[_pi+1]]-_p0); v += _n; }
                }
            }

            const scalar magV = mag(v) + VSMALL;
            v /= magV;

            // Scaled transition edges use the established pure wall
            // normal.  BL/neutral edges require the same direction policy
            // independently of layerScale_: neutral seam height may remain
            // unscaled while the extrusion direction must not revert to the
            // untreated-surface projection.
            const bool scaledTransition =
                terminateLayersAtConcaveEdges_
             && layerScale_.size() > bpI
             && layerScale_[bpI] < 0.99;

            const bool neutralEdgePoint =
                blNeutralEdgePoints_.found(bpI);

            if( scaledTransition || neutralEdgePoint )
            {
                const scalar magN = mag(normal) + VSMALL;
                normal /= magN;
            }
            else
            {
                normal -= (normal & v) * v;
            }

            const scalar magN = mag(normal) + VSMALL;
            normal /= magN;

            // BL/no-BL FlowBoundaryTermination ramp points still need
            // a finite local height baseline. The previous logic skipped
            // neighbor-distance clamping for all ramp points
            // (layerScale_<0.99), which could leave dist=VGREAT until the
            // final layerScale_ multiplication. On inlet/outlet contacts
            // this can create visible blowout/protrusion artifacts.
            //
            // Scoped to explicit termination patches (patchRole_==1) only
            // so periodic/neutral seam behavior is unchanged.
            const bool anyRampZone =
                terminateLayersAtConcaveEdges_
             && layerScale_.size() > bpI
             && layerScale_[bpI] < 0.99;

            bool flowRampZone = false;
            if( anyRampZone )
            {
                forAll(otherPatches, opI)
                {
                    const label op = otherPatches[opI];
                    if( op >= 0
                     && op < label(patchRole_.size())
                     && patchRole_[op] == 1 )
                    {
                        flowRampZone = true;
                        break;
                    }
                }
            }

            if( !anyRampZone )
            {
                // Original non-ramp behavior.
                const scalar distBeforeClamp = dist;

                // Diagnostic only: quantify the untreated-neighbour
                // projection distribution at BL/neutral EDGENODEs.
                // The production selector below remains unchanged.
                scalar neutralDiagMinProj = VGREAT;
                scalar neutralDiagSecondProj = VGREAT;
                scalar neutralDiagMaxProj = -VGREAT;
                scalar neutralDiagMinHalfLen = VGREAT;
                label neutralDiagMinBpJ = -1;
                label neutralDiagNCandidates = 0;

                forAllRow(pointPoints, bpI, ppI)
                {
                    const label bpJ = pointPoints(bpI, ppI);
                    if( patchVertex[bpJ] )
                        continue;

                    const vector vec = points[bPoints[bpJ]] - p;
                    const scalar prod = 0.5 * mag(vec & normal);

                    if( neutralEdgePoint )
                    {
                        const scalar halfLen = 0.5 * mag(vec);

                        ++neutralDiagNCandidates;
                        neutralDiagMaxProj =
                            Foam::max(neutralDiagMaxProj, prod);

                        if( halfLen < neutralDiagMinHalfLen )
                            neutralDiagMinHalfLen = halfLen;

                        if( prod < neutralDiagMinProj )
                        {
                            neutralDiagSecondProj = neutralDiagMinProj;
                            neutralDiagMinProj = prod;
                            neutralDiagMinBpJ = bpJ;
                        }
                        else if( prod < neutralDiagSecondProj )
                        {
                            neutralDiagSecondProj = prod;
                        }
                    }

                    if( prod < dist )
                        dist = prod;
                }

                // ----------------------------------------------------
                // CFMitch v2.8:
                // BL/neutral projected-height floor.
                //
                // Rotor37 diagnostics showed that the remaining
                // wall-adjacent bad-pyramid population is strongly
                // enriched at mixed neutral/ordinary wall faces, and
                // that the dominant discriminator is birth-height
                // collapse rather than layerScale or direction angle.
                //
                // At neutral points the established direction policy is
                // the pure treated-wall normal.  A valid neighbouring
                // boundary edge can therefore have a near-zero projection
                // onto that direction even when its physical length is
                // perfectly finite.
                //
                // Guard only that projection-derived neutral candidate.
                //
                // IMPORTANT:
                // distBeforeClamp remains a hard upper ceiling, so this
                // can never undo a stricter distance constraint computed
                // earlier in createNewVertex().
                // ----------------------------------------------------
                if
                (
                    neutralEdgePoint
                 && cfmitchNeutralProjectionFloor_ > scalar(0)
                 && neutralDiagMinProj < VGREAT
                 && neutralDiagMinHalfLen < VGREAT
                 && neutralDiagMinHalfLen > VSMALL
                )
                {
                    const scalar floorDist =
                        cfmitchNeutralProjectionFloor_
                      * neutralDiagMinHalfLen;

                    const scalar guardedNeutralDist =
                        Foam::max
                        (
                            neutralDiagMinProj,
                            floorDist
                        );

                    const scalar distBeforeV28 =
                        dist;

                    dist =
                        Foam::min
                        (
                            distBeforeClamp,
                            guardedNeutralDist
                        );

                    if( dist > distBeforeV28 )
                    {
                        static label nNeutralProjectionFloorRaised = 0;
                        ++nNeutralProjectionFloorRaised;

                        if( nNeutralProjectionFloorRaised <= 250 )
                        {
                            Info
                                << "CFMITCH V2.8 NEUTRAL PROJECTION FLOOR:"
                                << " count="
                                << nNeutralProjectionFloorRaised
                                << " bpI=" << bpI
                                << " pointI=" << bPoints[bpI]
                                << " floorFraction="
                                << cfmitchNeutralProjectionFloor_
                                << " minProj="
                                << neutralDiagMinProj
                                << " minHalfLen="
                                << neutralDiagMinHalfLen
                                << " floorDist="
                                << floorDist
                                << " distBeforeClamp="
                                << distBeforeClamp
                                << " oldDist="
                                << distBeforeV28
                                << " newDist="
                                << dist
                                << " growth="
                                << (
                                    distBeforeV28 > VSMALL
                                  ? dist / distBeforeV28
                                  : scalar(-1)
                                   )
                                << endl;
                        }

                        if( nNeutralProjectionFloorRaised == 250 )
                        {
                            Info
                                << "CFMITCH V2.8 NEUTRAL PROJECTION FLOOR:"
                                << " further raised-point diagnostics"
                                << " suppressed"
                                << endl;
                        }
                    }
                }

                if( neutralEdgePoint )
                {
                    const scalar minSecondRatio =
                        ( neutralDiagSecondProj < VGREAT
                       && neutralDiagSecondProj > VSMALL )
                      ? neutralDiagMinProj / neutralDiagSecondProj
                      : scalar(-1);

                    const scalar minLenRatio =
                        ( neutralDiagMinHalfLen < VGREAT
                       && neutralDiagMinHalfLen > VSMALL )
                      ? neutralDiagMinProj / neutralDiagMinHalfLen
                      : scalar(-1);

                    Info << "BLNEUTRALHEIGHTCAND"
                         << " bpI=" << bpI
                         << " pointI=" << bPoints[bpI]
                         << " layerScale="
                         << ((layerScale_.size() > bpI)
                             ? layerScale_[bpI] : scalar(1))
                         << " nCand=" << neutralDiagNCandidates
                         << " minBpJ=" << neutralDiagMinBpJ
                         << " minProj=" << neutralDiagMinProj
                         << " secondProj=" << neutralDiagSecondProj
                         << " maxProj=" << neutralDiagMaxProj
                         << " minHalfLen=" << neutralDiagMinHalfLen
                         << " minOverSecond=" << minSecondRatio
                         << " minOverHalfLen=" << minLenRatio
                         << " distSelected=" << dist
                         << endl;
                }
                // CLAMPDIAG
                {
                    static label nClampDiag = 0;
                    const scalar lsThis =
                        (layerScale_.size() > bpI) ? layerScale_[bpI] : scalar(1.0);
                    if( distBeforeClamp < VGREAT
                     && distBeforeClamp > VSMALL
                     && dist < scalar(0.10)*distBeforeClamp
                     && nClampDiag < 200 )
                    {
                        ++nClampDiag;
                        Info << "CLAMPDIAG bpI=" << bpI
                             << " x=" << p.x() << " y=" << p.y() << " z=" << p.z()
                             << " layerScale=" << lsThis
                             << " distBefore=" << distBeforeClamp
                             << " distAfter=" << dist
                             << " ratio=" << (dist/distBeforeClamp) << endl;
                    }
                }
            }
            else if( flowRampZone )
            {
                // FlowBoundaryTermination ramp: compute finite neighbor
                // baseline. Prefer projected spacing when meaningful;
                // otherwise fall back to local edge length.
                scalar neighProj = VGREAT;
                scalar neighLen  = VGREAT;

                forAllRow(pointPoints, bpI, ppI)
                {
                    const label bpJ = pointPoints(bpI, ppI);
                    if( bpJ < 0 || bpJ >= label(bPoints.size()) ) continue;
                    if( patchVertex[bpJ] ) continue;

                    const vector vec = points[bPoints[bpJ]] - p;
                    const scalar halfLen = 0.5 * mag(vec);
                    if( halfLen > VSMALL && halfLen < neighLen )
                        neighLen = halfLen;

                    const scalar proj = 0.5 * mag(vec & normal);
                    if( halfLen > VSMALL
                     && proj > scalar(100) * VSMALL
                     && proj > scalar(1e-4) * halfLen
                     && proj < neighProj )
                        neighProj = proj;
                }

                const scalar neighDist =
                    (neighProj < VGREAT) ? neighProj : neighLen;

                const scalar distBeforeClamp = dist;
                // Flow-termination seam guard:
                // If projected clearance is tiny relative to the local edge length,
                // a finite epsilon-height layer can create near-zero tetra/sliver
                // cells at inlet/outlet termination seams. Treat that as true local
                // termination instead.
                const bool tinyFlowTermClearance =
                    neighProj < VGREAT
                 && neighLen  < VGREAT
                 && neighProj < scalar(5e-6)
                 && neighProj < scalar(0.05) * neighLen;
                
                if( tinyFlowTermClearance )
                {
                    dist = 0.0;
                    static label nFlowTermTinyClearanceSuppress = 0;
                    if( nFlowTermTinyClearanceSuppress < 200 )
                    {
                        ++nFlowTermTinyClearanceSuppress;
                        const scalar lsThis =
                            (layerScale_.size() > bpI) ? layerScale_[bpI] : scalar(1.0);
                        Info << "FLOWTERM_TINY_CLEARANCE_SUPPRESS"
                             << " bpI=" << bpI
                             << " x=" << p.x()
                             << " y=" << p.y()
                             << " z=" << p.z()
                             << " layerScale=" << lsThis
                             << " neighProj=" << neighProj
                             << " neighLen=" << neighLen
                             << " ratio=" << (neighProj/(neighLen + VSMALL))
                             << endl;
                    }
                }
                else if( neighDist < VGREAT )
                {
                    dist = Foam::min(dist, neighDist);
                }

                // RAMPCLAMPDIAG: confirm ramp points get finite baseline.
                {
                    static label nRampClampDiag = 0;
                    const scalar lsThis =
                        (layerScale_.size() > bpI) ? layerScale_[bpI] : scalar(1.0);
                    if( nRampClampDiag < 5000 )
                    {
                        ++nRampClampDiag;
                        Info << "RAMPCLAMPDIAG bpI=" << bpI
                             << " x=" << p.x() << " y=" << p.y() << " z=" << p.z()
                             << " layerScale=" << lsThis
                             << " distBefore=" << distBeforeClamp
                             << " neighProj=" << neighProj
                             << " neighLen=" << neighLen
                             << " neighDist=" << neighDist
                             << " distAfter=" << dist
                             << endl;
                    }
                }
            }
            // else: neutral/periodic ramp zone -- behavior unchanged.
            }  // closes zero-dist else block
        }
        else if( otherPatches.size() == 2 )
        {
            # ifdef DEBUGLayer
            Info << "Vertex is a corner" << endl;
            # endif

            label otherVertex(-1);
            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ = pointPoints(bpI, ppI);

                bool found(true);
                forAll(otherPatches, opI)
                    if( !pPatches.contains(bpJ, otherPatches[opI]) )
                    {
                        found = false;
                        break;
                    }

                if( found )
                {
                    otherVertex = bpJ;
                    break;
                }
            }

            if( otherVertex == -1 )
            {
                FatalErrorIn
                (
                    "void boundaryLayers::createNewVertices"
                    "("
                        "const boolList& treatPatches,"
                        "labelList& newLabelForVertex"
                    ")"
                ) << "Cannot find moving vertex!" << exit(FatalError);
            }

            if( terminateLayersAtConcaveEdges_
             && layerScale_.size() > bpI
             && layerScale_[bpI] < 0.01 )
            {
                dist = 0.0;
            }
            else
            {
                //- normal vector is co-linear with that edge
                normal = p - points[bPoints[otherVertex]];
                dist = 0.5 * mag(normal) + VSMALL;
                normal /= 2.0 * dist;
                bd_branch = 2;
                {
                    const vector trueN = pNormals[bpI];
                    const scalar mn = mag(normal)*mag(trueN) + VSMALL;
                    bd_edgeDotNormal = (normal & trueN)/mn;
                }

                // 2-patch BL/no-BL flow-termination corner.
                // The raw edge-colinear corner distance can be too large
                // compared with the surrounding BL height field. Clamp it
                // to a finite neighboring baseline before common layerScale_
                // multiplication is applied below.
                bool cornerFlowRampZone = false;
                forAll(otherPatches, opI)
                {
                    const label op = otherPatches[opI];
                    if( op >= 0
                     && op < label(patchRole_.size())
                     && patchRole_[op] == 1 )
                    {
                        cornerFlowRampZone = true;
                        break;
                    }
                }

                if( cornerFlowRampZone )
                {
                    scalar neighProj = VGREAT;
                    scalar neighLen  = VGREAT;

                    forAllRow(pointPoints, bpI, ppI)
                    {
                        const label bpJ = pointPoints(bpI, ppI);
                        if( bpJ < 0 || bpJ >= label(bPoints.size()) ) continue;
                        if( patchVertex[bpJ] ) continue;

                        const vector vec = points[bPoints[bpJ]] - p;
                        const scalar halfLen = 0.5 * mag(vec);
                        if( halfLen > VSMALL && halfLen < neighLen )
                            neighLen = halfLen;

                        const scalar proj = 0.5 * mag(vec & normal);
                        if( halfLen > VSMALL
                         && proj > scalar(100) * VSMALL
                         && proj > scalar(1e-4) * halfLen
                         && proj < neighProj )
                        {
                            neighProj = proj;
                        }
                    }

                    const scalar neighDist =
                        (neighProj < VGREAT) ? neighProj : neighLen;

                    const scalar distBeforeCornerClamp = dist;
                    if( neighDist < VGREAT )
                        dist = Foam::min(dist, neighDist);

                    static label nCornerClampDiag = 0;
                    if( nCornerClampDiag < 200 )
                    {
                        ++nCornerClampDiag;
                        const scalar lsThis =
                            (layerScale_.size() > bpI) ? layerScale_[bpI] : scalar(1.0);

                        Info << "CORNERCLAMPDIAG bpI=" << bpI
                             << " x=" << p.x() << " y=" << p.y() << " z=" << p.z()
                             << " layerScale=" << lsThis
                             << " distBefore=" << distBeforeCornerClamp
                             << " neighProj=" << neighProj
                             << " neighLen=" << neighLen
                             << " neighDist=" << neighDist
                             << " distAfter=" << dist
                             << endl;
                    }
                }
            }
        }
        else
        {
            // Multi-patch singularity: 3+ non-treated patches meet here
            // (e.g. blade/root/periodic triple junction). A single prism
            // extrusion direction is not well-defined. Use surface normal
            // with a minimal safe distance to avoid zero-volume collapse
            // while keeping the extrusion within the local cell geometry.
            normal = pNormals[bpI];

            // Adaptive safe extrusion at multi-patch singularity.
            scalar minEdge(GREAT);
            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ = pointPoints(bpI, ppI);
                const scalar d = mag(points[bPoints[bpJ]] - p);
                if( d > VSMALL && d < minEdge )
                    minEdge = d;
            }

            scalar candidateDist =
                (minEdge < GREAT)
              ? Foam::max(scalar(0.02) * minEdge, scalar(100) * VSMALL)
              : scalar(0.0);

            bool accepted(false);

            for(label attempt=0; attempt<8 && candidateDist > VSMALL; ++attempt)
            {
                const point candidate = p - candidateDist * normal;
                bool crossesGuardPlane(false);

                forAllRow(pFaces, bpI, pfI)
                {
                    const label faceI = pFaces(bpI, pfI);
                    const label patchI = boundaryFacePatches[faceI];
                    if( patchI < 0 || patchI >= label(treatPatches.size()) )
                        continue;
                    if( treatPatches[patchI] )
                        continue;
                    const face& f = bFaces[faceI];
                    if( f.size() < 3 )
                        continue;
                    vector fn(vector::zero);
                    const point& fp0 = points[f[0]];
                    for(label pi=1; pi<f.size()-1; ++pi)
                        fn += (points[f[pi]] - fp0) ^ (points[f[pi+1]] - fp0);
                    if( mag(fn) < VSMALL )
                        continue;
                    fn /= mag(fn);
                    point fc(point::zero);
                    forAll(f, fi)
                        fc += points[f[fi]];
                    fc /= scalar(f.size());
                    const scalar s0 = (p - fc) & fn;
                    const scalar s1 = (candidate - fc) & fn;
                    if( mag(s0) > SMALL && s0 * s1 < scalar(0) )
                    {
                        crossesGuardPlane = true;
                        break;
                    }
                }

                if( !crossesGuardPlane )
                {
                    dist = candidateDist;
                    accepted = true;
                    break;
                }

                candidateDist *= scalar(0.5);
            }

            # ifdef DEBUGLayer
            if( accepted )
            {
                Info << "Multi-patch singularity at bpI=" << bpI
                     << " p=" << p
                     << ": accepted dist=" << dist
                     << " (attempt " << attempt << ")" << endl;
            }
            else
            {
                Info << "Multi-patch singularity at bpI=" << bpI
                     << " p=" << p
                     << ": adaptive bisection failed, dist=0" << endl;
            }
            # endif

            if( !accepted )
                dist = 0.0;
        }

        //- limit distances
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceLabel = pFaces(bpI, pfI);
            if( otherPatches.contains(boundaryFacePatches[faceLabel]) )
            {
                const face& f = bFaces[faceLabel];
                const label pos = f.which(bPoints[bpI]);

                if( pos != -1 )
                {
                    const point& ep1 = points[f.prevLabel(pos)];
                    const point& ep2 = points[f.nextLabel(pos)];

                    const scalar dst =
                        help::distanceOfPointFromTheEdge(ep1, ep2, p);

                    if( dst < dist )
                    {
                        if( cfmBirthTarget )
                        {
                            Info
                                << "CFMITCH BIRTHTARGET_EDGE_LIMIT"
                                << " bpI=" << bpI
                                << " pointI=" << cfmBirthPointI
                                << " faceLabel=" << faceLabel
                                << " oldDist=" << dist
                                << " dst=" << dst
                                << " newDist=" << scalar(0.9)*dst
                                << " normal=" << normal
                                << endl;
                        }

                        dist = scalar(0.9)*dst;
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                        "void boundaryLayers::createNewVertices"
                        "("
                            "const boolList& treatPatches,"
                            "labelList& newLabelForVertex"
                        ") const"
                    ) << "Face does not contains this vertex!"
                        << abort(FatalError);
                }
            }
        }
    }
    else
    {
        // BL/BL junction: per-patch normals, conservative min-dist selection
        if( blblJunctionPoints_.found(bpI) )
        {
            Map<vector> patchNormals;
            forAllRow(pFaces, bpI, pfI)
            {
                const label faceI = pFaces(bpI, pfI);
                const label patchLabel = boundaryFacePatches[faceI];
                if( patchLabel < 0 || patchLabel >= label(treatPatches.size()) ) continue;
                if( !treatPatches[patchLabel] ) continue;
                const face& f = bFaces[faceI];
                if( f.size() < 3 ) continue;
                vector fn = vector::zero;
                const point& p0 = points[f[0]];
                for(label pi=1; pi<f.size()-1; ++pi)
                    fn += (points[f[pi]] - p0) ^ (points[f[pi+1]] - p0);
                if( patchNormals.found(patchLabel) )
                    patchNormals[patchLabel] += fn;
                else
                    patchNormals.insert(patchLabel, fn);
            }
            if( patchNormals.size() >= 2 )
            {
                label junctionClass = 0; // 0=hard 1=moderate 2=mild
                Map<label>::const_iterator classIt =
                    blblJunctionClass_.find(bpI);
                if( classIt != blblJunctionClass_.end() )
                    junctionClass = classIt();

                bool usedBlendNormal = false;

                //- Moderate/mild BL/BL seams: use blended/bisector normal.
                //- This preserves boundary-to-boundary contact better than
                //- selecting one patch normal by min-distance.
                if( junctionClass > 0 )
                {
                    vector blendNormal = vector::zero;
                    label nBlend = 0;

                    forAllConstIter(Map<vector>, patchNormals, it)
                    {
                        vector n = it();
                        const scalar magN = mag(n);
                        if( magN < VSMALL ) continue;
                        blendNormal += n / magN;
                        ++nBlend;
                    }

                    if( nBlend >= 2 && mag(blendNormal) > VSMALL )
                    {
                        blendNormal /= mag(blendNormal);

                        scalar blendDist = VGREAT;
                        forAllRow(pointPoints, bpI, ppI)
                        {
                            const label bpJ = pointPoints(bpI, ppI);
                            if( bpJ < 0 || bpJ >= label(bPoints.size()) )
                                continue;
                            const vector vec = points[bPoints[bpJ]] - p;
                            const scalar d = 0.5 * mag(vec & blendNormal);
                            if( d < blendDist ) blendDist = d;
                        }

                        normal = blendNormal;
                        dist = Foam::min(dist, blendDist);
                        usedBlendNormal = true;
                    }
                }

                //- Hard corners, or failed blend: keep original conservative
                //- min-distance patch-normal selection.
                if( !usedBlendNormal )
                {
                    scalar bestDist = GREAT;
                    vector bestNormal = vector::zero;
                    forAllConstIter(Map<vector>, patchNormals, it)
                    {
                        vector n = it();
                        const scalar magN = mag(n);
                        if( magN < VSMALL ) continue;
                        n /= magN;
                        scalar localDist = VGREAT;
                        forAllRow(pointPoints, bpI, ppI)
                        {
                            const label bpJ = pointPoints(bpI, ppI);
                            if( bpJ < 0 || bpJ >= label(bPoints.size()) )
                                continue;
                            const vector vec = points[bPoints[bpJ]] - p;
                            const scalar d = 0.5 * mag(vec & n);
                            if( d < localDist ) localDist = d;
                        }
                        if( localDist < bestDist )
                        {
                            bestDist = localDist;
                            bestNormal = n;
                        }
                    }
                    if( mag(bestNormal) > VSMALL )
                    {
                        normal = bestNormal;
                        dist = Foam::min(dist, bestDist);
                    }
                    else
                        normal = pNormals[bpI];
                }
            }
            else
                normal = pNormals[bpI];
        }
        else
        {
            // Feature-aware normal: only average faces from treated patches
            vector patchNormal(vector::zero);
            forAllRow(pFaces, bpI, pfI)
            {
                const label patchLabel = boundaryFacePatches[pFaces(bpI, pfI)];
                if( treatPatches[patchLabel] )
                {
                    const face& f = bFaces[pFaces(bpI, pfI)];
                    vector _n=vector::zero;
                    const point& _p0=points[f[0]];
                    for(label _pi=1;_pi<f.size()-1;++_pi)
                        _n+=(points[f[_pi]]-_p0)^(points[f[_pi+1]]-_p0);
                    patchNormal += _n;
                }
            }
            const scalar magPN = mag(patchNormal);
            if( magPN > VSMALL )
                normal = patchNormal / magPN;
            else
                normal = pNormals[bpI];
        }

        forAllRow(pointPoints, bpI, ppI)
        {
            const scalar d =
            0.5 * mag
            (
                points[bPoints[pointPoints(bpI, ppI)]] -
                p
            );

            if( d < dist )
            {
                if( cfmBirthTarget )
                {
                    Info
                        << "CFMITCH BIRTHTARGET_GENERIC_LIMIT"
                        << " bpI=" << bpI
                        << " pointI=" << cfmBirthPointI
                        << " oldDist=" << dist
                        << " d=" << d
                        << " newDist=" << d
                        << " normal=" << normal
                        << endl;
                }

                dist = d;
            }
        }
    }

    //- create new vertex
    # ifdef DEBUGLayer
    Info << "Normal for vertex " << bpI << " is " << normal << endl;
    Info << "Distance is " << dist << endl;
    # endif

    if( cfmBirthTarget )
    {
        Info
            << "CFMITCH BIRTHTARGET_PRE_SCALE"
            << " bpI=" << bpI
            << " pointI=" << cfmBirthPointI
            << " normal=" << normal
            << " normalMag=" << mag(normal)
            << " dist=" << dist
            << " layerScale="
            << (
                   layerScale_.size() > bpI
                 ? layerScale_[bpI]
                 : scalar(1)
               )
            << endl;
    }

    // Apply layerScale_ ramp at BL/no-BL transition zones
    const scalar rawDist = dist;
    if( terminateLayersAtConcaveEdges_ && layerScale_.size() > bpI )
    {
        dist *= layerScale_[bpI];
    }
    if( dist > SMALL )
        dist = Foam::max(dist, VSMALL);
    else if( terminateLayersAtConcaveEdges_
          && layerScale_.size() > bpI
          && layerScale_[bpI] > 0.01 )
    {
        // Geometry forced a nonzero-scale point to near-zero dist.
        // Use rawDist-scaled floor to avoid collapsed layer cells.
        dist = Foam::max
        (
            scalar(1e-6) * rawDist,
            scalar(100) * VSMALL
        );
    }
    else
        dist = 0.0;

    if( cfmBirthTarget )
    {
        Info
            << "CFMITCH BIRTHTARGET_POST_SCALE"
            << " bpI=" << bpI
            << " pointI=" << cfmBirthPointI
            << " rawDist=" << rawDist
            << " finalDist=" << dist
            << " normal=" << normal
            << " predictedPMinusNewP=" << (dist*normal)
            << " predictedDispMag=" << mag(dist*normal)
            << endl;
    }

    point newP = p - dist * normal;

    if( cfmBirthTarget )
    {
        Info
            << "CFMITCH BIRTHTARGET_INITIAL_NEWP"
            << " bpI=" << bpI
            << " pointI=" << cfmBirthPointI
            << " p=" << p
            << " newP=" << newP
            << " pMinusNewP=" << (p-newP)
            << " dispMag=" << mag(p-newP)
            << endl;
    }
    if( help::isnan(newP) || help::isinf(newP) )
        return p;
    // BL/neutral crossing clamp: prevent extrusion across periodic/symmetry planes.
    // Fires for blNeutralEdgePoints_ regardless of layerScale_.
    if( !blNeutralEdgePoints_.empty() && blNeutralEdgePoints_.found(bpI) )
    {
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceI = pFaces(bpI, pfI);
            const label patchI = boundaryFacePatches[faceI];
            if( patchI < 0 || patchI >= label(patchNames_.size()) ) continue;
            // Only check true neutral guard patches (patchRole_ == 2).
            // Periodic patches are topological continuations, not guard planes.
            // Applying this crossing clamp to blade/periodic seam points
            // creates false positives and collapses valid BL extrusion.
            if( patchRole_.size() <= patchI ) continue;
            if( patchRole_[patchI] != 2 ) continue;

            bool isPeriodicPatch = false;
            if( patchI < label(patchNames_.size()) )
            {
                const word& pName = patchNames_[patchI];
                if( pName.find("periodic") != std::string::npos )
                    isPeriodicPatch = true;
            }
            if( isPeriodicPatch )
            {
                static label nSkipPeriodicGuardDiag = 0;
                if( nSkipPeriodicGuardDiag < 200 )
                {
                    ++nSkipPeriodicGuardDiag;
                    Info << "SKIPPERIODICNEUTRALGUARD bpI=" << bpI
                         << " patchI=" << patchI
                         << " x=" << p.x() << " y=" << p.y() << " z=" << p.z()
                         << endl;
                }
                continue;
            }

            const face& f = bFaces[faceI];
            vector fn = vector::zero;
            const point& fp0 = points[f[0]];
            for(label pi=1; pi<f.size()-1; ++pi)
                fn += (points[f[pi]]-fp0)^(points[f[pi+1]]-fp0);
            if( mag(fn) < VSMALL ) continue;
            fn /= mag(fn);
            point fc = point::zero;
            forAll(f, fi) fc += points[f[fi]];
            fc /= scalar(f.size());
            const scalar s0 = (p - fc) & fn;
            const scalar s1 = (newP - fc) & fn;
            if( mag(s0) > SMALL && s0*s1 < scalar(0) )
            {
                // Analytic clamp: the constraint
                // (p - d*normal - fc) & fn == 0 is linear in d.
                // Solve the exact crossing distance instead of blindly
                // halving six times. But do not accept epsilon-height
                // clamps; those become sliver cells. If the safe crossing
                // distance is too small, leave bisectAccepted=false and
                // let the existing suppress-to-surface path handle it.
                scalar localDist = dist;
                point candidate = newP;
                bool bisectAccepted = false;
                scalar dCross = -VGREAT;
                const scalar nDotFn = normal & fn;
                const scalar minCross =
                    Foam::max(scalar(0.02) * dist, scalar(100) * VSMALL);

                if( mag(nDotFn) > VSMALL )
                {
                    dCross = s0 / nDotFn;
                    if( dCross > minCross && dCross < dist )
                    {
                        localDist = scalar(0.9) * dCross;
                        candidate = p - localDist*normal;
                        const scalar sCand = (candidate - fc) & fn;
                        if( s0*sCand >= scalar(0) )
                        {
                            newP = candidate;
                            bisectAccepted = true;
                        }
                    }
                }

                static label nAnalyticClampDiag = 0;
                if( nAnalyticClampDiag < 300 )
                {
                    ++nAnalyticClampDiag;
                    Info << "ANALYTICCLAMPDIAG bpI=" << bpI
                         << " nDotFn=" << nDotFn
                         << " s0=" << s0
                         << " dist=" << dist
                         << " dCross=" << dCross
                         << " minCross=" << minCross
                         << " accepted=" << bisectAccepted
                         << " localDist=" << localDist
                         << endl;
                }
                if( !bisectAccepted )
                {
                    static label nBisectFailDiag = 0;
                    if( nBisectFailDiag < 500 )
                    {
                        ++nBisectFailDiag;
                        const scalar lsThis =
                            (layerScale_.size() > bpI) ? layerScale_[bpI] : scalar(1.0);
                        Info << "BISECTFAILSUPPRESS bpI=" << bpI
                             << " x=" << p.x() << " y=" << p.y() << " z=" << p.z()
                             << " layerScale=" << lsThis
                             << " dist=" << dist
                             << " localDist=" << localDist
                             << " finalHBeforeSuppress=" << mag(candidate - p)
                             << " -- setting newP=p"
                             << endl;
                    }
                    newP = p;
                }
                break;
            }
        }
    }
    // Robust candidate-point clamping near BL/no-BL termination patches
    // Check ALL no-BL faces -- use most restrictive constraint.
    // A point must not move further from any no-BL face than its
    // original position. This prevents BL extrusion through inlet,
    // outlet, and periodic surfaces at corner junctions.
    if( terminateLayersAtConcaveEdges_
     && layerScale_.size() > bpI
     && layerScale_[bpI] < 0.99 )
    {
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceI = pFaces(bpI, pfI);
            const label patchI = boundaryFacePatches[faceI];
            if( patchI < 0 || patchI >= label(patchRole_.size()) ) continue;
            if( patchRole_.size() <= label(patchI) || patchRole_[patchI] != 1 ) continue;  // explicit termination only
            const face& f = bFaces[faceI];
            vector fn = vector::zero;
            const point& fp0 = points[f[0]];
            for(label pi=1; pi<f.size()-1; ++pi)
                fn += (points[f[pi]]-fp0)^(points[f[pi+1]]-fp0);
            if( mag(fn) < VSMALL ) continue;
            fn /= mag(fn);
            point fc = point::zero;
            forAll(f, fi) fc += points[f[fi]];
            fc /= scalar(f.size());
            // Signed distance from no-BL face plane
            // fn points outward from no-BL patch into the domain
            // p should be on the positive side (inside domain)
            // If newP goes to negative side, it crossed through the surface
            const scalar s0 = (p - fc) & fn;
            const scalar s1 = (newP - fc) & fn;
            // Only clamp if point actually crossed the face plane
            // s0 > 0 means original point is on correct side
            // s1 < 0 means extruded point crossed to wrong side
            if( mag(s0) > SMALL && s0*s1 < scalar(0) )
            {
                // Analytic clamp: the constraint
                // (p - d*normal - fc) & fn == 0 is linear in d.
                // Solve the exact crossing distance instead of blindly
                // halving six times. But do not accept epsilon-height
                // clamps; those become sliver cells. If the safe crossing
                // distance is too small, leave bisectAccepted=false and
                // let the existing suppress-to-surface path handle it.
                scalar localDist = dist;
                point candidate = newP;
                bool bisectAccepted = false;
                scalar dCross = -VGREAT;
                const scalar nDotFn = normal & fn;
                const scalar minCross =
                    Foam::max(scalar(0.02) * dist, scalar(100) * VSMALL);

                if( mag(nDotFn) > VSMALL )
                {
                    dCross = s0 / nDotFn;
                    if( dCross > minCross && dCross < dist )
                    {
                        localDist = scalar(0.9) * dCross;
                        candidate = p - localDist*normal;
                        const scalar sCand = (candidate - fc) & fn;
                        if( s0*sCand >= scalar(0) )
                        {
                            newP = candidate;
                            bisectAccepted = true;
                        }
                    }
                }

                static label nAnalyticClampDiag = 0;
                if( nAnalyticClampDiag < 300 )
                {
                    ++nAnalyticClampDiag;
                    Info << "ANALYTICCLAMPDIAG bpI=" << bpI
                         << " nDotFn=" << nDotFn
                         << " s0=" << s0
                         << " dist=" << dist
                         << " dCross=" << dCross
                         << " minCross=" << minCross
                         << " accepted=" << bisectAccepted
                         << " localDist=" << localDist
                         << endl;
                }
                if( !bisectAccepted )
                {
                    static label nBisectFailDiag = 0;
                    if( nBisectFailDiag < 500 )
                    {
                        ++nBisectFailDiag;
                        const scalar lsThis =
                            (layerScale_.size() > bpI) ? layerScale_[bpI] : scalar(1.0);
                        Info << "BISECTFAILSUPPRESS bpI=" << bpI
                             << " x=" << p.x() << " y=" << p.y() << " z=" << p.z()
                             << " layerScale=" << lsThis
                             << " dist=" << dist
                             << " localDist=" << localDist
                             << " finalHBeforeSuppress=" << mag(candidate - p)
                             << " -- setting newP=p"
                             << endl;
                    }
                    newP = p;
                }
                break;
            }
        }
    }
    // BIRTHDIAG: log near-zero extrusion OR tangential 2-patch corner extrusion.
    {
        static label nBirthDiag = 0;
        const scalar finalH = mag(newP - p);
        const scalar lsThis =
            (layerScale_.size() > bpI) ? layerScale_[bpI] : scalar(1.0);

        const bool nearZeroBirth = finalH < scalar(1e-6);
        const bool tangentCorner =
            (bd_branch == 2 && mag(bd_edgeDotNormal) < scalar(0.10));

        if( (nearZeroBirth || tangentCorner) && nBirthDiag < 800 )
        {
            ++nBirthDiag;
            Info << "BIRTHDIAG bpI=" << bpI
                 << " branch=" << bd_branch
                 << " x=" << p.x() << " y=" << p.y() << " z=" << p.z()
                 << " finalH=" << finalH
                 << " rawDist=" << rawDist
                 << " dist=" << dist
                 << " layerScale=" << lsThis
                 << " edgeDotN=" << bd_edgeDotNormal
                 << " nearZero=" << nearZeroBirth
                 << " tangentCorner=" << tangentCorner
                 << endl;
        }
    }
    if( cfmBirthTarget )
    {
        Info
            << "CFMITCH BIRTHTARGET_FINAL"
            << " bpI=" << bpI
            << " pointI=" << cfmBirthPointI
            << " p=" << p
            << " newP=" << newP
            << " pMinusNewP=" << (p-newP)
            << " dispMag=" << mag(p-newP)
            << " normal=" << normal
            << " distVariable=" << dist
            << endl;
    }

    return newP;
}

void boundaryLayers::suppressFailedSingularityExtrusions
(
    const labelList& patchLabels
)
{
    const meshSurfaceEngine& mse = surfaceEngine();

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    boolList treatPatches(boundaries.size(), false);
    forAll(patchLabels, i)
    {
        const label patchI = patchLabels[i];
        if( patchI >= 0 && patchI < label(treatPatches.size()) )
            treatPatches[patchI] = true;
    }

    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const pointFieldPMG& points = mesh_.points();

    const VRWGraph& pFaces = mse.pointFaces();
    const VRWGraph& pointPoints = mse.pointPoints();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    if( suppressLayerAtBndFace_.size() != bFaces.size() )
        suppressLayerAtBndFace_.setSize(bFaces.size(), false);

    label nCandidates(0);
    label nFailed(0);
    label nFacesSuppressed(0);

    forAll(bPoints, bpI)
    {
        if( bpI < 0 || bpI >= label(pPatches.size()) )
            continue;

        // Only multi-patch singular points are considered here.
        // Single-patch and ordinary two-patch edge points are handled by
        // the normal BL extrusion logic.
        if( pPatches.sizeOfRow(bpI) < 3 )
            continue;

        bool hasTreated(false);
        bool hasNotTreated(false);

        forAllRow(pPatches, bpI, ppI)
        {
            const label patchI = pPatches(bpI, ppI);

            if( patchI < 0 || patchI >= label(treatPatches.size()) )
                continue;

            if( treatPatches[patchI] )
                hasTreated = true;
            else
                hasNotTreated = true;
        }

        // This mirrors the risky mixed treated / non-treated multi-patch
        // singularity case. Pure treated or pure non-treated corners are not
        // the target of this pre-pass.
        if( !hasTreated || !hasNotTreated )
            continue;

        ++nCandidates;

        const point& p = points[bPoints[bpI]];
        // Compute normal deterministically from treated face normals
        // instead of pNormals which has OMP non-determinism.
        vector normal(vector::zero);
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceI = pFaces(bpI, pfI);
            if( faceI < 0 || faceI >= label(bFaces.size()) ) continue;
            const label patchI = boundaryFacePatches[faceI];
            if( patchI < 0 || patchI >= label(treatPatches.size()) ) continue;
            if( !treatPatches[patchI] ) continue;
            const face& f = bFaces[faceI];
            if( f.size() < 3 ) continue;
            vector fn(vector::zero);
            const point& fp0 = points[f[0]];
            for(label pi=1; pi<f.size()-1; ++pi)
                fn += (points[f[pi]] - fp0) ^ (points[f[pi+1]] - fp0);
            normal += fn;
        }
        if( mag(normal) < VSMALL )
            continue;

        normal /= mag(normal);

        scalar minEdge(GREAT);

        forAllRow(pointPoints, bpI, ppI)
        {
            const label bpJ = pointPoints(bpI, ppI);
            if( bpJ < 0 || bpJ >= label(bPoints.size()) )
                continue;

            const scalar d = mag(points[bPoints[bpJ]] - p);

            if( d > VSMALL && d < minEdge )
                minEdge = d;
        }

        scalar candidateDist =
            (minEdge < GREAT)
          ? Foam::max(scalar(0.02) * minEdge, scalar(100) * VSMALL)
          : scalar(0.0);

        bool accepted(false);

        for(label attempt=0; attempt<8 && candidateDist > VSMALL; ++attempt)
        {
            const point candidate = p - candidateDist * normal;

            bool crossesGuardPlane(false);

            forAllRow(pFaces, bpI, pfI)
            {
                const label faceI = pFaces(bpI, pfI);

                if( faceI < 0 || faceI >= label(bFaces.size()) )
                    continue;

                const label patchI = boundaryFacePatches[faceI];

                if( patchI < 0 || patchI >= label(treatPatches.size()) )
                    continue;

                // Use non-treated faces as guard planes, matching the logic
                // in createNewVertex().
                if( treatPatches[patchI] )
                    continue;

                const face& f = bFaces[faceI];

                if( f.size() < 3 )
                    continue;

                vector fn(vector::zero);
                const point& fp0 = points[f[0]];

                for(label pi=1; pi<f.size()-1; ++pi)
                    fn += (points[f[pi]] - fp0) ^ (points[f[pi+1]] - fp0);

                if( mag(fn) < VSMALL )
                    continue;

                fn /= mag(fn);

                point fc(point::zero);

                forAll(f, fi)
                    fc += points[f[fi]];

                fc /= scalar(f.size());

                const scalar s0 = (p - fc) & fn;
                const scalar s1 = (candidate - fc) & fn;

                if( mag(s0) > SMALL && s0 * s1 < scalar(0) )
                {
                    crossesGuardPlane = true;
                    break;
                }
            }

            if( !crossesGuardPlane )
            {
                accepted = true;
                break;
            }

            candidateDist *= scalar(0.5);
        }

        scalar predictedDist = accepted ? candidateDist : scalar(0.0);

        // Simulate the same non-treated-face edge clamp used later in
        // createNewVertex(). If the clamp collapses the accepted bisection
        // distance to zero/near-zero, suppress the affected faces before
        // vertex/cell construction so the BL graph remains coherent.
        if( accepted )
        {
            forAllRow(pFaces, bpI, pfI)
            {
                const label faceI = pFaces(bpI, pfI);

                if( faceI < 0 || faceI >= label(bFaces.size()) )
                    continue;

                const label patchI = boundaryFacePatches[faceI];

                if( patchI < 0 || patchI >= label(treatPatches.size()) )
                    continue;

                if( treatPatches[patchI] )
                    continue;

                const face& f = bFaces[faceI];
                const label pos = f.which(bPoints[bpI]);

                if( pos == -1 )
                    continue;

                const point& ep1 = points[f.prevLabel(pos)];
                const point& ep2 = points[f.nextLabel(pos)];

                const scalar dst =
                    help::distanceOfPointFromTheEdge(ep1, ep2, p);

                if( dst < predictedDist )
                    predictedDist = scalar(0.9) * dst;
            }
        }

        const scalar minUsableDist = scalar(100) * VSMALL;

        if( predictedDist > minUsableDist )
            continue;

        ++nFailed;

        // This point would produce a zero/near-zero extrusion distance after
        // the full createNewVertex logic. Suppress touching faces BEFORE
        // vertex/cell construction.
        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);

            if( bfI < 0 || bfI >= label(suppressLayerAtBndFace_.size()) )
                continue;

            if( !suppressLayerAtBndFace_[bfI] )
            {
                suppressLayerAtBndFace_[bfI] = true;
                ++nFacesSuppressed;
            }
        }
    }

    Info << "Failed singularity extrusion pre-pass: candidates="
         << nCandidates
         << " failed=" << nFailed
         << " newlySuppressedFaces=" << nFacesSuppressed
         << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::constrainNeutralSeamCandidatesBeforeSwap()
{
    // ============================================================
    // CFMITCH CONNECTED NEUTRAL SEAM CONSTRAINT V1
    //
    // Problem being fixed:
    //
    // createNewVertex() currently chooses a BL/neutral extrusion
    // independently at each vertex from its local treated-face fan.
    //
    // Rotor37 demonstrated adjacent clean blade/periodic seam points
    // whose as-born macro hairs differed by ~45 degrees.
    //
    // V1 changes authority:
    //
    //   createNewVertex()
    //       -> local candidate / clearance provider
    //
    //   this connected pass
    //       -> seam direction authority
    //
    // For a clean wall + neutral seam:
    //
    //   tau = connected seam tangent
    //   n   = connected/smoothed neutral-manifold normal
    //   d   = n ^ tau
    //
    // Therefore d:
    //
    //   1. lies in the neutral tangent plane
    //   2. is perpendicular to the seam curve
    //   3. is computed as a connected field, not independent hairs
    //
    // V1 deliberately preserves the candidate MACRO hair length.
    //
    // refineBoundaryLayers remains authoritative for:
    //   - h1
    //   - layer count
    //   - growth ratio
    //
    // No patch-name logic.
    //
    // Only clean rank-2 BL_PATCH + NEUTRAL_PATCH seams participate.
    // BL/BL, BL/no-BL, termination and higher-order junctions stay
    // on their existing dedicated paths.
    //
    // Serial / OMP implementation first.
    // ============================================================

    if( blNeutralEdgePoints_.empty() )
        return;

    if( Pstream::parRun() )
    {
        Info
            << "CFMITCH SEAMSOLVE V1: "
            << "decomposed MPI skipped; serial/OMP only"
            << endl;

        return;
    }

    const meshSurfaceEngine& mse = surfaceEngine();

    const labelList& bPoints =
        mse.boundaryPoints();

    const faceList::subList& bFaces =
        mse.boundaryFaces();

    const labelList& bFacePatches =
        mse.boundaryFacePatches();

    const VRWGraph& pFaces =
        mse.pointFaces();

    const VRWGraph& pointPoints =
        mse.pointPoints();

    const meshSurfacePartitioner& mPart =
        surfacePartitioner();

    const VRWGraph& pPatches =
        mPart.pointPatches();

    pointFieldPMG& points =
        mesh_.points();

    const label nBP = bPoints.size();

    boolList active(nBP, false);

    labelList wallPatch(nBP, -1);
    labelList neutralPatch(nBP, -1);
    labelList topLabel(nBP, -1);

    vectorField rawNeutralNormal
    (
        nBP,
        vector::zero
    );

    vectorField neutralNormal
    (
        nBP,
        vector::zero
    );

    vectorField originalDisp
    (
        nBP,
        vector::zero
    );

    vectorField proposedDisp
    (
        nBP,
        vector::zero
    );

    scalarField originalLength
    (
        nBP,
        scalar(0)
    );

    // ------------------------------------------------------------
    // CFMITCH CONNECTED NEUTRAL SEAM CONSTRAINT V1B
    //
    // Connected vector-field state.
    // ------------------------------------------------------------

    vectorField seamTangent
    (
        nBP,
        vector::zero
    );

    vectorField preferredDir
    (
        nBP,
        vector::zero
    );

    vectorField solvedDir
    (
        nBP,
        vector::zero
    );

    boolList tangentValid
    (
        nBP,
        false
    );

    boolList directionValid
    (
        nBP,
        false
    );

    label nInitiallyEligible = 0;
    label nNoTop = 0;
    label nBadRole = 0;
    label nSpecialSkipped = 0;
    label nBadNormal = 0;


    // ============================================================
    // CLASSIFY CLEAN WALL + NEUTRAL SEAM POINTS
    // ============================================================

    forAllConstIter
    (
        labelHashSet,
        blNeutralEdgePoints_,
        iter
    )
    {
        const label bpI = iter.key();

        if
        (
            bpI < 0
         || bpI >= nBP
        )
            continue;

        if
        (
            blblCornerPoints_.found(bpI)
         || blblJunctionPoints_.found(bpI)
         || blNoBlEdgePoints_.found(bpI)
        )
        {
            ++nSpecialSkipped;
            continue;
        }

        label nWall = 0;
        label nNeutral = 0;
        label nTermination = 0;

        label wPatch = -1;
        label nPatch = -1;

        forAllRow(pPatches, bpI, ppI)
        {
            const label patchI =
                pPatches(bpI, ppI);

            if
            (
                patchI < 0
             || patchI >= label(patchRole_.size())
            )
                continue;

            const label role =
                patchRole_[patchI];

            if( role == 0 )
            {
                ++nWall;
                wPatch = patchI;
            }
            else if( role == 1 )
            {
                ++nTermination;
            }
            else if( role == 2 )
            {
                ++nNeutral;
                nPatch = patchI;
            }
        }

        if
        (
            nWall != 1
         || nNeutral != 1
         || nTermination != 0
        )
        {
            ++nBadRole;
            continue;
        }

        const label pointI =
            bPoints[bpI];

        if
        (
            pointI < 0
         || pointI >= label(newLabelForVertex_.size())
        )
        {
            ++nNoTop;
            continue;
        }

        const label topI =
            newLabelForVertex_[pointI];

        if
        (
            topI < 0
         || topI >= label(points.size())
        )
        {
            ++nNoTop;
            continue;
        }

        const vector E0 =
            points[topI]
          - points[pointI];

        const scalar L0 =
            mag(E0);

        if( L0 < VSMALL )
        {
            ++nNoTop;
            continue;
        }


        // --------------------------------------------------------
        // Neutral manifold normal from neutral-side incident faces.
        // --------------------------------------------------------

        vector nSum(vector::zero);

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI =
                pFaces(bpI, pfI);

            if
            (
                bfI < 0
             || bfI >= label(bFaces.size())
             || bfI >= label(bFacePatches.size())
            )
                continue;

            if
            (
                bFacePatches[bfI] != nPatch
            )
                continue;

            const face& f =
                bFaces[bfI];

            if( f.size() < 3 )
                continue;

            vector fn(vector::zero);

            const point& fp0 =
                points[f[0]];

            for
            (
                label pi=1;
                pi<f.size()-1;
                ++pi
            )
            {
                fn +=
                    (points[f[pi]] - fp0)
                  ^ (points[f[pi+1]] - fp0);
            }

            nSum += fn;
        }

        const scalar nMag =
            mag(nSum);

        if( nMag < VSMALL )
        {
            ++nBadNormal;
            continue;
        }

        active[bpI] = true;

        wallPatch[bpI] =
            wPatch;

        neutralPatch[bpI] =
            nPatch;

        topLabel[bpI] =
            topI;

        rawNeutralNormal[bpI] =
            nSum/nMag;

        neutralNormal[bpI] =
            rawNeutralNormal[bpI];

        originalDisp[bpI] =
            E0;

        originalLength[bpI] =
            L0;

        ++nInitiallyEligible;
    }


    // ============================================================
    // GRAPH VALIDATION
    //
    // Clean seam topology:
    //
    // degree 1 -> open endpoint
    // degree 2 -> ordinary interior / closed loop
    //
    // degree >2 is a junction and belongs elsewhere.
    // ============================================================

    label nGraphPruned = 0;

    for
    (
        label prunePass=0;
        prunePass<4;
        ++prunePass
    )
    {
        boolList keep(active);

        bool changed = false;

        forAll(active, bpI)
        {
            if( !active[bpI] )
                continue;

            label degree = 0;

            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ =
                    pointPoints(bpI, ppI);

                if
                (
                    bpJ < 0
                 || bpJ >= nBP
                 || !active[bpJ]
                )
                    continue;

                if
                (
                    wallPatch[bpJ]
                        != wallPatch[bpI]
                 || neutralPatch[bpJ]
                        != neutralPatch[bpI]
                )
                    continue;

                ++degree;
            }

            if
            (
                degree < 1
             || degree > 2
            )
            {
                keep[bpI] = false;
                changed = true;
                ++nGraphPruned;
            }
        }

        active = keep;

        if( !changed )
            break;
    }


    // ============================================================
    // CONNECTED NEUTRAL NORMAL FIELD
    //
    // Local face normals are noisy at clipped seams.
    //
    // Smooth ONLY along a seam having the same wall/neutral pair.
    // Sign-align neighbors before averaging.
    // ============================================================

    const label normalSmoothPasses = 4;

    for
    (
        label pass=0;
        pass<normalSmoothPasses;
        ++pass
    )
    {
        vectorField nextNormal
        (
            neutralNormal
        );

        forAll(active, bpI)
        {
            if( !active[bpI] )
                continue;

            vector acc =
                scalar(2.0)
               *neutralNormal[bpI];

            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ =
                    pointPoints(bpI, ppI);

                if
                (
                    bpJ < 0
                 || bpJ >= nBP
                 || !active[bpJ]
                )
                    continue;

                if
                (
                    wallPatch[bpJ]
                        != wallPatch[bpI]
                 || neutralPatch[bpJ]
                        != neutralPatch[bpI]
                )
                    continue;

                vector nJ =
                    neutralNormal[bpJ];

                if
                (
                    (neutralNormal[bpI] & nJ)
                    < scalar(0)
                )
                {
                    nJ = -nJ;
                }

                acc += nJ;
            }

            const scalar aMag =
                mag(acc);

            if( aMag > VSMALL )
            {
                nextNormal[bpI] =
                    acc/aMag;
            }
        }

        neutralNormal = nextNormal;
    }


    // ============================================================
    // CFMITCH NEUTRAL MANIFOLD STENCIL V6
    //
    // Build the neutral constraint plane from a small interior
    // stencil on the same neutral patch rather than only from the
    // clipped contact-line face fan.
    //
    // depth 0 = seam point
    // depth 1 = immediate neutral-surface neighbours
    // depth 2 = one additional neutral-surface ring
    //
    // The resulting area-weighted surface normal becomes the
    // authoritative neutral constraint normal used by V5's
    // projected-native displacement architecture.
    // ============================================================

    const label neutralStencilDepth = 2;

    vectorField neutralBeforeStencil
    (
        neutralNormal
    );

    label nStencilSolved = 0;
    label nStencilFallback = 0;
    label nStencilFacesTotal = 0;

    label minStencilFaces = 1000000000;
    label maxStencilFaces = 0;


    forAll(active, bpI)
    {
        if( !active[bpI] )
            continue;

        const label nPatch =
            neutralPatch[bpI];

        if( nPatch < 0 )
        {
            ++nStencilFallback;
            continue;
        }


        // --------------------------------------------------------
        // Build point stencil constrained to the same neutral patch.
        // --------------------------------------------------------

        labelHashSet visitedPoints;
        labelLongList frontier;

        visitedPoints.insert(bpI);
        frontier.append(bpI);


        for
        (
            label depth=0;
            depth<neutralStencilDepth;
            ++depth
        )
        {
            labelLongList nextFrontier;

            forAll(frontier, fI)
            {
                const label bpK =
                    frontier[fI];

                forAllRow(pointPoints, bpK, ppI)
                {
                    const label bpJ =
                        pointPoints(bpK, ppI);

                    if
                    (
                        bpJ < 0
                     || bpJ >= nBP
                     || visitedPoints.found(bpJ)
                    )
                        continue;


                    bool onNeutralPatch = false;

                    forAllRow(pFaces, bpJ, pfI)
                    {
                        const label bfI =
                            pFaces(bpJ, pfI);

                        if
                        (
                            bfI < 0
                         || bfI >= label(bFacePatches.size())
                        )
                            continue;

                        if
                        (
                            bFacePatches[bfI] == nPatch
                        )
                        {
                            onNeutralPatch = true;
                            break;
                        }
                    }


                    if( !onNeutralPatch )
                        continue;

                    visitedPoints.insert(bpJ);
                    nextFrontier.append(bpJ);
                }
            }


            frontier = nextFrontier;

            if( frontier.size() == 0 )
                break;
        }


        // --------------------------------------------------------
        // Collect unique neutral-patch faces touched by stencil.
        // --------------------------------------------------------

        labelHashSet stencilFaces;

        forAllConstIter
        (
            labelHashSet,
            visitedPoints,
            pIter
        )
        {
            const label bpK =
                pIter.key();

            forAllRow(pFaces, bpK, pfI)
            {
                const label bfI =
                    pFaces(bpK, pfI);

                if
                (
                    bfI < 0
                 || bfI >= label(bFaces.size())
                 || bfI >= label(bFacePatches.size())
                )
                    continue;

                if
                (
                    bFacePatches[bfI] != nPatch
                )
                    continue;

                stencilFaces.insert(bfI);
            }
        }


        // --------------------------------------------------------
        // Area-weighted normal over the stencil.
        // --------------------------------------------------------

        vector nSum(vector::zero);

        forAllConstIter
        (
            labelHashSet,
            stencilFaces,
            fIter
        )
        {
            const label bfI =
                fIter.key();

            const face& f =
                bFaces[bfI];

            if( f.size() < 3 )
                continue;

            vector fn(vector::zero);

            const point& fp0 =
                points[f[0]];

            for
            (
                label pi=1;
                pi<f.size()-1;
                ++pi
            )
            {
                fn +=
                    (points[f[pi]] - fp0)
                  ^ (points[f[pi+1]] - fp0);
            }

            if( mag(fn) < VSMALL )
                continue;


            // Sign-align against original neutral-side normal.
            if
            (
                (fn & rawNeutralNormal[bpI])
                < scalar(0)
            )
            {
                fn = -fn;
            }

            nSum += fn;
        }


        const label nFaces =
            stencilFaces.size();

        nStencilFacesTotal +=
            nFaces;

        minStencilFaces =
            Foam::min
            (
                minStencilFaces,
                nFaces
            );

        maxStencilFaces =
            Foam::max
            (
                maxStencilFaces,
                nFaces
            );


        const scalar nMag =
            mag(nSum);

        if( nMag < VSMALL )
        {
            ++nStencilFallback;
            continue;
        }


        neutralNormal[bpI] =
            nSum/nMag;

        ++nStencilSolved;
    }


    // ------------------------------------------------------------
    // Light connected regularisation ALONG each same-patch seam.
    // Strong self-weight preserves real neutral-surface curvature.
    // ------------------------------------------------------------

    const label manifoldSmoothPasses = 6;

    for
    (
        label pass=0;
        pass<manifoldSmoothPasses;
        ++pass
    )
    {
        vectorField nextNormal
        (
            neutralNormal
        );


        forAll(active, bpI)
        {
            if( !active[bpI] )
                continue;

            vector acc =
                scalar(4.0)
               *neutralNormal[bpI];


            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ =
                    pointPoints(bpI, ppI);

                if
                (
                    bpJ < 0
                 || bpJ >= nBP
                 || !active[bpJ]
                )
                    continue;

                if
                (
                    wallPatch[bpJ]
                        != wallPatch[bpI]
                 || neutralPatch[bpJ]
                        != neutralPatch[bpI]
                )
                    continue;


                vector nJ =
                    neutralNormal[bpJ];

                if
                (
                    (nJ & neutralNormal[bpI])
                    < scalar(0)
                )
                {
                    nJ = -nJ;
                }

                acc += nJ;
            }


            const scalar aMag =
                mag(acc);

            if( aMag > VSMALL )
            {
                nextNormal[bpI] =
                    acc/aMag;
            }
        }


        neutralNormal =
            nextNormal;
    }


    // ------------------------------------------------------------
    // Constraint-plane coherence audit.
    // ------------------------------------------------------------

    scalar maxNormalUnitDiffBefore =
        scalar(0);

    scalar maxNormalUnitDiffAfter =
        scalar(0);

    scalar minNormalDotBefore =
        scalar(1);

    scalar minNormalDotAfter =
        scalar(1);

    label nNormalEdges = 0;


    forAll(active, bpI)
    {
        if( !active[bpI] )
            continue;


        vector oldNI =
            neutralBeforeStencil[bpI];

        vector newNI =
            neutralNormal[bpI];

        if
        (
            mag(oldNI) < VSMALL
         || mag(newNI) < VSMALL
        )
            continue;

        oldNI /= mag(oldNI);
        newNI /= mag(newNI);


        forAllRow(pointPoints, bpI, ppI)
        {
            const label bpJ =
                pointPoints(bpI, ppI);

            if
            (
                bpJ <= bpI
             || bpJ < 0
             || bpJ >= nBP
             || !active[bpJ]
            )
                continue;

            if
            (
                wallPatch[bpJ]
                    != wallPatch[bpI]
             || neutralPatch[bpJ]
                    != neutralPatch[bpI]
            )
                continue;


            vector oldNJ =
                neutralBeforeStencil[bpJ];

            vector newNJ =
                neutralNormal[bpJ];

            if
            (
                mag(oldNJ) < VSMALL
             || mag(newNJ) < VSMALL
            )
                continue;

            oldNJ /= mag(oldNJ);
            newNJ /= mag(newNJ);


            // Constraint planes are orientation-insensitive.
            if( (oldNI & oldNJ) < scalar(0) )
                oldNJ = -oldNJ;

            if( (newNI & newNJ) < scalar(0) )
                newNJ = -newNJ;


            const scalar oldDot =
                oldNI & oldNJ;

            const scalar newDot =
                newNI & newNJ;


            maxNormalUnitDiffBefore =
                Foam::max
                (
                    maxNormalUnitDiffBefore,
                    mag(oldNI-oldNJ)
                );

            maxNormalUnitDiffAfter =
                Foam::max
                (
                    maxNormalUnitDiffAfter,
                    mag(newNI-newNJ)
                );


            minNormalDotBefore =
                Foam::min
                (
                    minNormalDotBefore,
                    oldDot
                );

            minNormalDotAfter =
                Foam::min
                (
                    minNormalDotAfter,
                    newDot
                );

            ++nNormalEdges;
        }
    }


    const scalar avgStencilFaces =
    (
        nStencilSolved > 0
      ? scalar(nStencilFacesTotal)
       /scalar(nStencilSolved)
      : scalar(0)
    );


    if( minStencilFaces == 1000000000 )
        minStencilFaces = 0;


    Info
        << "CFMITCH NEUTRAL MANIFOLD V6:"
        << " depth=" << neutralStencilDepth
        << " solved=" << nStencilSolved
        << " fallback=" << nStencilFallback
        << " minFaces=" << minStencilFaces
        << " maxFaces=" << maxStencilFaces
        << " avgFaces=" << avgStencilFaces
        << " smoothPasses=" << manifoldSmoothPasses
        << " normalEdges=" << nNormalEdges
        << " maxNormalDiffBefore="
        << maxNormalUnitDiffBefore
        << " maxNormalDiffAfter="
        << maxNormalUnitDiffAfter
        << " minNormalDotBefore="
        << minNormalDotBefore
        << " minNormalDotAfter="
        << minNormalDotAfter
        << endl;


    // ============================================================
    // CONSTRUCT CONNECTED SEAM VECTOR FIELD -- V1B
    //
    // V1 failure:
    //
    //     d = n ^ tau
    //     if (d & oldHair) < 0 => flip d
    //
    // was done independently per vertex.
    //
    // Because d is an unoriented axis until the seam component is
    // oriented, that local sign choice produced adjacent ~180-degree
    // reversals.
    //
    // V1B:
    //
    //   A. build raw seam tangent axes
    //   B. orient tangent axes component-by-component
    //   C. smooth the connected tangent field
    //   D. construct preferred cross-seam directions
    //   E. orient growth directions component-by-component
    //   F. projected graph smoothing inside neutral tangent planes
    //
    // Hard invariant:
    //
    //     neutralNormal & solvedDir == 0
    //
    // approximately to floating-point precision.
    // ============================================================

    label nBadTangent = 0;
    label nBadDirection = 0;


    // ------------------------------------------------------------
    // A. RAW SEAM TANGENT AXES
    // ------------------------------------------------------------

    forAll(active, bpI)
    {
        if( !active[bpI] )
            continue;

        label nb0 = -1;
        label nb1 = -1;
        label degree = 0;

        forAllRow(pointPoints, bpI, ppI)
        {
            const label bpJ =
                pointPoints(bpI, ppI);

            if
            (
                bpJ < 0
             || bpJ >= nBP
             || !active[bpJ]
            )
                continue;

            if
            (
                wallPatch[bpJ]
                    != wallPatch[bpI]
             || neutralPatch[bpJ]
                    != neutralPatch[bpI]
            )
                continue;

            if( degree == 0 )
                nb0 = bpJ;
            else if( degree == 1 )
                nb1 = bpJ;

            ++degree;
        }

        if
        (
            degree < 1
         || degree > 2
        )
        {
            ++nBadTangent;
            continue;
        }

        const point& root =
            points[bPoints[bpI]];

        vector tau(vector::zero);

        if( degree == 1 )
        {
            tau =
                points[bPoints[nb0]]
              - root;
        }
        else
        {
            vector v0 =
                points[bPoints[nb0]]
              - root;

            vector v1 =
                points[bPoints[nb1]]
              - root;

            const scalar m0 = mag(v0);
            const scalar m1 = mag(v1);

            if
            (
                m0 < VSMALL
             || m1 < VSMALL
            )
            {
                ++nBadTangent;
                continue;
            }

            v0 /= m0;
            v1 /= m1;

            // Neighbours lie on opposite sides of a regular chain
            // vertex.  Their difference is the centred tangent axis.
            tau = v0 - v1;

            if( mag(tau) < VSMALL )
                tau = v0;
        }

        vector n =
            neutralNormal[bpI];

        const scalar nMag =
            mag(n);

        if( nMag < VSMALL )
        {
            ++nBadTangent;
            continue;
        }

        n /= nMag;

        // Tangent must lie in the neutral tangent plane.
        tau -= (tau & n)*n;

        const scalar tauMag =
            mag(tau);

        if( tauMag < VSMALL )
        {
            ++nBadTangent;
            continue;
        }

        seamTangent[bpI] =
            tau/tauMag;

        tangentValid[bpI] = true;
    }


    // ------------------------------------------------------------
    // B. ORIENT TANGENT AXES PER CONNECTED COMPONENT
    //
    // Tangent sign is geometrically arbitrary.  Pick a seed and
    // propagate orientation through the graph so neighbours point
    // consistently along the same component.
    // ------------------------------------------------------------

    boolList tangentOriented
    (
        nBP,
        false
    );

    label nTangentComponents = 0;

    for(label seed=0; seed<nBP; ++seed)
    {
        if
        (
            !active[seed]
         || !tangentValid[seed]
         || tangentOriented[seed]
        )
            continue;

        ++nTangentComponents;

        labelLongList queue;
        queue.append(seed);

        tangentOriented[seed] = true;

        for(label qI=0; qI<label(queue.size()); ++qI)
        {
            const label bpI =
                queue[qI];

            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ =
                    pointPoints(bpI, ppI);

                if
                (
                    bpJ < 0
                 || bpJ >= nBP
                 || !active[bpJ]
                 || !tangentValid[bpJ]
                )
                    continue;

                if
                (
                    wallPatch[bpJ]
                        != wallPatch[bpI]
                 || neutralPatch[bpJ]
                        != neutralPatch[bpI]
                )
                    continue;

                if( tangentOriented[bpJ] )
                    continue;

                if
                (
                    (seamTangent[bpJ]
                    & seamTangent[bpI])
                    < scalar(0)
                )
                {
                    seamTangent[bpJ] =
                        -seamTangent[bpJ];
                }

                tangentOriented[bpJ] = true;
                queue.append(bpJ);
            }
        }
    }


    // ------------------------------------------------------------
    // C. SMOOTH CONNECTED TANGENT FIELD
    //
    // Reproject after every Jacobi pass so neutral tangency remains
    // a hard constraint.
    // ------------------------------------------------------------

    const label tangentSmoothPasses = 8;

    for
    (
        label pass=0;
        pass<tangentSmoothPasses;
        ++pass
    )
    {
        vectorField nextTangent
        (
            seamTangent
        );

        forAll(active, bpI)
        {
            if
            (
                !active[bpI]
             || !tangentValid[bpI]
            )
                continue;

            vector acc =
                scalar(2.0)
               *seamTangent[bpI];

            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ =
                    pointPoints(bpI, ppI);

                if
                (
                    bpJ < 0
                 || bpJ >= nBP
                 || !active[bpJ]
                 || !tangentValid[bpJ]
                )
                    continue;

                if
                (
                    wallPatch[bpJ]
                        != wallPatch[bpI]
                 || neutralPatch[bpJ]
                        != neutralPatch[bpI]
                )
                    continue;

                vector tJ =
                    seamTangent[bpJ];

                if
                (
                    (tJ & seamTangent[bpI])
                    < scalar(0)
                )
                {
                    tJ = -tJ;
                }

                acc += tJ;
            }

            vector n =
                neutralNormal[bpI];

            const scalar nMag = mag(n);

            if( nMag < VSMALL )
                continue;

            n /= nMag;

            acc -= (acc & n)*n;

            const scalar aMag =
                mag(acc);

            if( aMag < VSMALL )
                continue;

            vector candidate =
                acc/aMag;

            if
            (
                (candidate & seamTangent[bpI])
                < scalar(0)
            )
            {
                candidate = -candidate;
            }

            nextTangent[bpI] =
                candidate;
        }

        seamTangent =
            nextTangent;
    }


    // ------------------------------------------------------------
    // D. PREFERRED CROSS-SEAM DIRECTION
    //
    // n ^ tau is tangent to the neutral manifold and transverse to
    // the connected seam.
    // ------------------------------------------------------------

    forAll(active, bpI)
    {
        if
        (
            !active[bpI]
         || !tangentValid[bpI]
        )
            continue;

        vector n =
            neutralNormal[bpI];

        const scalar nMag =
            mag(n);

        if( nMag < VSMALL )
        {
            ++nBadDirection;
            continue;
        }

        n /= nMag;

        vector d =
            n ^ seamTangent[bpI];

        const scalar dMag =
            mag(d);

        if( dMag < VSMALL )
        {
            ++nBadDirection;
            continue;
        }

        preferredDir[bpI] =
            d/dMag;

        directionValid[bpI] = true;
    }


    // ============================================================
    // CFMITCH PROJECTED-NATIVE SEAM COLLAR V5
    //
    // The previous architectures treated
    //
    //     nNeutral ^ seamTangent
    //
    // as an unoriented geometric axis and then attempted to solve
    // its +/- sign.
    //
    // Rotor37 demonstrated that neither:
    //
    //     independent sign selection
    //
    // nor:
    //
    //     one sign per connected component
    //
    // is generally sufficient.
    //
    // V5 removes the sign problem entirely.
    //
    // The original createNewVertex displacement is a PHYSICAL,
    // oriented vector.  At a clean BL/neutral seam we project that
    // vector into the neutral tangent plane:
    //
    //     p = oldU - (oldU & nNeutral) nNeutral
    //
    // This satisfies the hard neutral-manifold constraint while
    // preserving the physical extrusion side.
    //
    // The projected field is then smoothed as an ORIENTED connected
    // vector field along the seam.
    //
    // n^tau is retained only as a weak/fallback geometric reference
    // when the native projection is nearly singular.
    // ============================================================


    vectorField projectedNative
    (
        nBP,
        vector::zero
    );

    scalarField nativeProjectionStrength
    (
        nBP,
        scalar(0)
    );

    boolList nativeProjectionValid
    (
        nBP,
        false
    );

    label nProjectedNative = 0;
    label nWeakProjectedNative = 0;
    label nFallbackAxis = 0;

    scalar minProjectionStrength =
        GREAT;

    scalar sumProjectionStrength =
        scalar(0);


    // ------------------------------------------------------------
    // A. PHYSICAL NATIVE PROJECTION
    // ------------------------------------------------------------

    forAll(active, bpI)
    {
        if
        (
            !active[bpI]
         || !directionValid[bpI]
        )
            continue;

        vector n =
            neutralNormal[bpI];

        const scalar nMag =
            mag(n);

        if( nMag < VSMALL )
        {
            directionValid[bpI] = false;
            continue;
        }

        n /= nMag;


        vector oldU =
            originalDisp[bpI];

        const scalar oldMag =
            mag(oldU);

        if( oldMag < VSMALL )
        {
            directionValid[bpI] = false;
            continue;
        }

        oldU /= oldMag;


        vector p =
            oldU
          - (oldU & n)*n;

        const scalar pMag =
            mag(p);

        nativeProjectionStrength[bpI] =
            pMag;

        minProjectionStrength =
            Foam::min
            (
                minProjectionStrength,
                pMag
            );

        sumProjectionStrength +=
            pMag;

        ++nProjectedNative;


        if( pMag > scalar(0.05) )
        {
            p /= pMag;

            projectedNative[bpI] =
                p;

            preferredDir[bpI] =
                p;

            solvedDir[bpI] =
                p;

            nativeProjectionValid[bpI] =
                true;
        }
        else
        {
            ++nWeakProjectedNative;


            // ----------------------------------------------------
            // Weakly constrained physical projection.
            //
            // Keep the geometric n^tau direction as a temporary
            // fallback, but orient it to the original physical hair.
            //
            // Connected smoothing below is allowed to replace this
            // with neighbour information.
            // ----------------------------------------------------

            vector axis =
                preferredDir[bpI];

            const scalar axisMag =
                mag(axis);

            if( axisMag < VSMALL )
            {
                directionValid[bpI] = false;
                continue;
            }

            axis /= axisMag;

            if( (axis & oldU) < scalar(0) )
                axis = -axis;

            projectedNative[bpI] =
                axis;

            preferredDir[bpI] =
                axis;

            solvedDir[bpI] =
                axis;

            ++nFallbackAxis;
        }
    }


    // ------------------------------------------------------------
    // B. CONNECTED ORIENTED SEAM SMOOTHING
    //
    // There are deliberately NO vector sign flips.
    //
    // A negative neighbour contribution is real geometric
    // disagreement and participates in the average as such.
    //
    // Reliable native projections receive stronger self-weight.
    // Weak projections are driven primarily by their neighbours.
    // ------------------------------------------------------------

    const label seamPhysicalSmoothPasses =
        12;

    for
    (
        label pass=0;
        pass<seamPhysicalSmoothPasses;
        ++pass
    )
    {
        vectorField nextDir
        (
            solvedDir
        );

        forAll(active, bpI)
        {
            if
            (
                !active[bpI]
             || !directionValid[bpI]
            )
                continue;


            scalar selfWeight =
                nativeProjectionValid[bpI]
              ? scalar(2.0)
              : scalar(0.25);


            vector acc =
                selfWeight
               *projectedNative[bpI];

            scalar totalWeight =
                selfWeight;

            label nNbr = 0;


            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ =
                    pointPoints(bpI, ppI);

                if
                (
                    bpJ < 0
                 || bpJ >= nBP
                 || !active[bpJ]
                 || !directionValid[bpJ]
                )
                    continue;

                if
                (
                    wallPatch[bpJ]
                        != wallPatch[bpI]
                 || neutralPatch[bpJ]
                        != neutralPatch[bpI]
                )
                    continue;


                // ORIENTED vector: never flip dJ.
                acc +=
                    solvedDir[bpJ];

                totalWeight +=
                    scalar(1.0);

                ++nNbr;
            }


            if
            (
                nNbr == 0
             || totalWeight < VSMALL
            )
                continue;


            acc /=
                totalWeight;


            // ----------------------------------------------------
            // Hard neutral-manifold tangency.
            // ----------------------------------------------------

            vector n =
                neutralNormal[bpI];

            const scalar nMag =
                mag(n);

            if( nMag < VSMALL )
                continue;

            n /= nMag;

            acc -=
                (acc & n)*n;


            const scalar aMag =
                mag(acc);

            if( aMag < scalar(1e-8) )
            {
                // Near cancellation means local orientation is
                // genuinely unresolved by this iteration.
                //
                // Retain previous value rather than invent a sign.
                continue;
            }


            vector candidate =
                acc/aMag;


            // ----------------------------------------------------
            // Physical hemisphere guard.
            //
            // For a reliable projected-native direction, never
            // reverse to the opposite extrusion side.
            //
            // Weak points are allowed to inherit their orientation
            // from the connected field.
            // ----------------------------------------------------

            if( nativeProjectionValid[bpI] )
            {
                if
                (
                    (candidate
                    & projectedNative[bpI])
                    <= scalar(0)
                )
                {
                    candidate =
                        projectedNative[bpI];
                }
            }


            nextDir[bpI] =
                candidate;
        }


        solvedDir =
            nextDir;
    }


    // ------------------------------------------------------------
    // C. FINAL EXACT TANGENCY PROJECTION
    // ------------------------------------------------------------

    scalar maxProjectedNativeLeak =
        scalar(0);

    forAll(active, bpI)
    {
        if
        (
            !active[bpI]
         || !directionValid[bpI]
        )
            continue;

        vector n =
            neutralNormal[bpI];

        const scalar nMag =
            mag(n);

        if( nMag < VSMALL )
            continue;

        n /= nMag;


        vector d =
            solvedDir[bpI];

        d -=
            (d & n)*n;

        const scalar dMag =
            mag(d);

        if( dMag < VSMALL )
            continue;

        d /= dMag;


        solvedDir[bpI] =
            d;

        preferredDir[bpI] =
            d;


        maxProjectedNativeLeak =
            Foam::max
            (
                maxProjectedNativeLeak,
                Foam::mag(d & n)
            );
    }


    const scalar avgProjectionStrength =
    (
        nProjectedNative > 0
      ? sumProjectionStrength
       /scalar(nProjectedNative)
      : scalar(0)
    );


    Info
        << "CFMITCH SEAMCOLLAR V6 PROJECTION:"
        << " projected=" << nProjectedNative
        << " weak=" << nWeakProjectedNative
        << " fallbackAxis=" << nFallbackAxis
        << " smoothPasses="
        << seamPhysicalSmoothPasses
        << " minStrength="
        << minProjectionStrength
        << " avgStrength="
        << avgProjectionStrength
        << " maxNeutralLeak="
        << maxProjectedNativeLeak
        << endl;


    // ============================================================
    // CFMITCH CONNECTED SEAM COLLAR V2
    //
    // V1B successfully made the BL/neutral seam itself coherent and
    // exactly tangent to the neutral manifold, but left the adjacent
    // wall-only BL field frozen.
    //
    // Rotor37 V1B result:
    //
    //   seam0/seam1 improved ~45 deg -> ~16 deg
    //
    // but:
    //
    //   seam1/interior remained ~45 deg
    //
    // and POST_BL_CREATE negative volumes increased 39 -> 130.
    //
    // V2 treats the constrained seam as the boundary condition of a
    // small 2-D wall collar.
    //
    // ring 0 : BL/neutral seam, hard neutral tangency
    // ring 1 : first ordinary wall collar
    // ring 2 : second ordinary wall collar
    // ring 3 : third ordinary wall collar
    //
    // Ring 1 is allowed to respond strongly to the seam.
    // Rings 2/3 increasingly prefer their original wall extrusion.
    //
    // Macro hair LENGTH is still preserved exactly in V2.
    // Height/clearance compatibility is a later architecture layer.
    // ============================================================


    // ------------------------------------------------------------
    // F1. BUILD WALL COLLAR
    // ------------------------------------------------------------

    const label maxCollarRing = 3;

    labelList collarRing
    (
        nBP,
        -1
    );

    boolList solvePoint
    (
        nBP,
        false
    );

    label nRing0 = 0;
    label nRing1 = 0;
    label nRing2 = 0;
    label nRing3 = 0;

    // Existing valid seam field becomes ring 0.
    forAll(active, bpI)
    {
        if
        (
            !active[bpI]
         || !directionValid[bpI]
        )
            continue;

        collarRing[bpI] = 0;
        solvePoint[bpI] = true;

        ++nRing0;
    }


    // ------------------------------------------------------------
    // Expand only through ordinary wall-only points belonging to
    // the SAME BL patch.
    //
    // Do not cross:
    //   - another neutral constraint
    //   - termination boundary
    //   - BL/BL edge
    //   - BL/no-BL edge
    //   - corner/junction
    // ------------------------------------------------------------

    for
    (
        label ring=1;
        ring<=maxCollarRing;
        ++ring
    )
    {
        forAll(collarRing, bpI)
        {
            if( collarRing[bpI] != ring-1 )
                continue;

            const label sourceWallPatch =
                wallPatch[bpI];

            if( sourceWallPatch < 0 )
                continue;

            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ =
                    pointPoints(bpI, ppI);

                if
                (
                    bpJ < 0
                 || bpJ >= nBP
                 || collarRing[bpJ] >= 0
                )
                    continue;

                if
                (
                    blblCornerPoints_.found(bpJ)
                 || blblJunctionPoints_.found(bpJ)
                 || blNoBlEdgePoints_.found(bpJ)
                 || blNeutralEdgePoints_.found(bpJ)
                )
                    continue;


                // -----------------------------------------------
                // Require an ordinary wall-only incidence class.
                // -----------------------------------------------

                label nWallRoles = 0;
                label nNonWallRoles = 0;
                bool hasSourceWall = false;

                forAllRow(pPatches, bpJ, ppi)
                {
                    const label patchI =
                        pPatches(bpJ, ppi);

                    if
                    (
                        patchI < 0
                     || patchI >= label(patchRole_.size())
                    )
                        continue;

                    const label role =
                        patchRole_[patchI];

                    if( role == 0 )
                    {
                        ++nWallRoles;

                        if( patchI == sourceWallPatch )
                            hasSourceWall = true;
                    }
                    else
                    {
                        ++nNonWallRoles;
                    }
                }

                if
                (
                    !hasSourceWall
                 || nWallRoles != 1
                 || nNonWallRoles != 0
                )
                    continue;


                const label pointI =
                    bPoints[bpJ];

                if
                (
                    pointI < 0
                 || pointI >= label(newLabelForVertex_.size())
                )
                    continue;

                const label topI =
                    newLabelForVertex_[pointI];

                if
                (
                    topI < 0
                 || topI >= label(points.size())
                )
                    continue;

                const vector E0 =
                    points[topI]
                  - points[pointI];

                const scalar L0 =
                    mag(E0);

                if( L0 < VSMALL )
                    continue;


                collarRing[bpJ] =
                    ring;

                solvePoint[bpJ] =
                    true;

                wallPatch[bpJ] =
                    sourceWallPatch;

                // No neutral equality constraint on wall-only collar.
                neutralPatch[bpJ] =
                    -1;

                topLabel[bpJ] =
                    topI;

                originalDisp[bpJ] =
                    E0;

                originalLength[bpJ] =
                    L0;

                preferredDir[bpJ] =
                    E0/L0;

                solvedDir[bpJ] =
                    preferredDir[bpJ];

                directionValid[bpJ] =
                    true;

                if( ring == 1 )
                    ++nRing1;
                else if( ring == 2 )
                    ++nRing2;
                else if( ring == 3 )
                    ++nRing3;
            }
        }
    }


    // ------------------------------------------------------------
    // F2. INITIALISE SEAM FIELD FROM V1B
    // ------------------------------------------------------------

    forAll(active, bpI)
    {
        if
        (
            !active[bpI]
         || !directionValid[bpI]
        )
            continue;

        solvedDir[bpI] =
            preferredDir[bpI];
    }


    // ------------------------------------------------------------
    // F3. V3 ORIENTATION INITIALISATION
    //
    // preferredDir is now an ORIENTED vector field.
    //
    // There is deliberately NO component-level sign propagation
    // here.  A negative dot product between two neighbouring
    // vectors is diagnostic information, not permission to turn
    // one vector around.
    // ------------------------------------------------------------

    // V5 uses an intrinsically oriented projected-native field.
    // No binary component-sign solve exists anymore.
    label nV2Components = 0;
    label nV2SignFlips = 0;

    forAll(solvePoint, bpI)
    {
        if
        (
            !solvePoint[bpI]
         || !directionValid[bpI]
        )
            continue;

        solvedDir[bpI] =
            preferredDir[bpI];
    }


    // ------------------------------------------------------------
    // F4. COUPLED 2-D COLLAR FIELD SOLVE
    //
    // Objective, conceptually:
    //
    //   sum_edges |d_i-d_j|^2
    //     +
    //   w(ring) |d_i-dPreferred_i|^2
    //
    // with hard seam constraint:
    //
    //   nNeutral_i . d_i = 0
    //
    // Ring weights:
    //
    //   0 -> 1.0    constrained seam
    //   1 -> 1.0    strong collar response
    //   2 -> 2.0
    //   3 -> 4.0    transition back to native wall field
    // ------------------------------------------------------------

    const label collarSmoothPasses = 20;

    for
    (
        label pass=0;
        pass<collarSmoothPasses;
        ++pass
    )
    {
        vectorField nextDir
        (
            solvedDir
        );

        forAll(solvePoint, bpI)
        {
            if
            (
                !solvePoint[bpI]
             || !directionValid[bpI]
            )
                continue;

            const label ring =
                collarRing[bpI];

            scalar prefWeight =
                scalar(1.0);

            if( ring == 2 )
                prefWeight = scalar(2.0);
            else if( ring >= 3 )
                prefWeight = scalar(4.0);


            // V3: preferredDir is an oriented physical vector.
            // Never turn it into an unoriented axis here.
            vector pref =
                preferredDir[bpI];

            vector acc =
                prefWeight*pref;

            scalar totalWeight =
                prefWeight;

            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ =
                    pointPoints(bpI, ppI);

                if
                (
                    bpJ < 0
                 || bpJ >= nBP
                 || !solvePoint[bpJ]
                 || !directionValid[bpJ]
                )
                    continue;

                // Collar solve remains on the same physical wall patch.
                if
                (
                    wallPatch[bpJ]
                    != wallPatch[bpI]
                )
                    continue;

                // V3: displacement is an ORIENTED vector.
                //
                // Do not flip a neighbour to make the local average
                // numerically convenient.  Opposing vectors must
                // cancel / expose a discontinuity rather than being
                // silently identified as the same axis.
                const vector dJ =
                    solvedDir[bpJ];

                const scalar nbrWeight =
                    scalar(1.0);

                acc +=
                    nbrWeight*dJ;

                totalWeight +=
                    nbrWeight;
            }

            if( totalWeight < VSMALL )
                continue;

            acc /= totalWeight;


            // -----------------------------------------------
            // Ring 0: HARD neutral-manifold tangency.
            // -----------------------------------------------

            if( ring == 0 )
            {
                vector n =
                    neutralNormal[bpI];

                const scalar nMag =
                    mag(n);

                if( nMag < VSMALL )
                    continue;

                n /= nMag;

                acc -=
                    (acc & n)*n;
            }


            const scalar aMag =
                mag(acc);

            if( aMag < VSMALL )
                continue;

            vector candidate =
                acc/aMag;


            // -----------------------------------------------
            // V3 HARD PHYSICAL-HEMISPHERE INVARIANT
            //
            // preferredDir is the oriented authoritative reference:
            //
            //   seam   -> connected constrained axis, sign selected
            //             by projected native extrusion
            //
            //   collar -> original wall extrusion
            //
            // The smoothed vector may rotate within that hemisphere
            // but may never cross 90 degrees and reverse physical
            // extrusion side.
            // -----------------------------------------------

            vector physicalRef =
                preferredDir[bpI];

            const scalar refMag =
                mag(physicalRef);

            if( refMag < VSMALL )
                continue;

            physicalRef /=
                refMag;

            if
            (
                (candidate & physicalRef)
                <= scalar(0)
            )
            {
                candidate =
                    physicalRef;
            }

            nextDir[bpI] =
                candidate;
        }

        solvedDir =
            nextDir;
    }


    // ------------------------------------------------------------
    // F5. V2 FIELD QUALITY AUDIT
    // ------------------------------------------------------------

    scalar maxNeutralLeak =
        scalar(0);

    scalar maxAllUnitDiffBefore =
        scalar(0);

    scalar maxAllUnitDiffAfter =
        scalar(0);

    scalar maxSeamCollarDiffAfter =
        scalar(0);

    // CFMITCH SEAMCOLLAR V2 WORST EDGE AUDIT
    label worstAllI = -1;
    label worstAllJ = -1;

    label worstSeamCollarI = -1;
    label worstSeamCollarJ = -1;

    label nAllComparedEdges = 0;
    label nSeamCollarEdges = 0;

    forAll(solvePoint, bpI)
    {
        if
        (
            !solvePoint[bpI]
         || !directionValid[bpI]
        )
            continue;

        if( collarRing[bpI] == 0 )
        {
            vector n =
                neutralNormal[bpI];

            const scalar nMag =
                mag(n);

            if( nMag > VSMALL )
            {
                n /= nMag;

                maxNeutralLeak =
                    Foam::max
                    (
                        maxNeutralLeak,
                        Foam::mag
                        (
                            solvedDir[bpI] & n
                        )
                    );
            }
        }


        const vector oldU0 =
            originalDisp[bpI]
           /(mag(originalDisp[bpI]) + VSMALL);

        const vector newU0 =
            solvedDir[bpI]
           /(mag(solvedDir[bpI]) + VSMALL);

        forAllRow(pointPoints, bpI, ppI)
        {
            const label bpJ =
                pointPoints(bpI, ppI);

            if
            (
                bpJ <= bpI
             || bpJ < 0
             || bpJ >= nBP
             || !solvePoint[bpJ]
             || !directionValid[bpJ]
            )
                continue;

            if
            (
                wallPatch[bpJ]
                != wallPatch[bpI]
            )
                continue;

            const vector oldU1 =
                originalDisp[bpJ]
               /(mag(originalDisp[bpJ]) + VSMALL);

            const vector newU1 =
                solvedDir[bpJ]
               /(mag(solvedDir[bpJ]) + VSMALL);

            const scalar before =
                mag(oldU0-oldU1);

            const scalar after =
                mag(newU0-newU1);

            maxAllUnitDiffBefore =
                Foam::max
                (
                    maxAllUnitDiffBefore,
                    before
                );

            if( after > maxAllUnitDiffAfter )
            {
                maxAllUnitDiffAfter = after;
                worstAllI = bpI;
                worstAllJ = bpJ;
            }

            if
            (
                (
                    collarRing[bpI] == 0
                 && collarRing[bpJ] == 1
                )
             ||
                (
                    collarRing[bpI] == 1
                 && collarRing[bpJ] == 0
                )
            )
            {
                if( after > maxSeamCollarDiffAfter )
                {
                    maxSeamCollarDiffAfter = after;
                    worstSeamCollarI = bpI;
                    worstSeamCollarJ = bpJ;
                }

                ++nSeamCollarEdges;
            }

            ++nAllComparedEdges;
        }
    }


    // ------------------------------------------------------------
    // F6. PRESERVE MACRO LENGTH, CHANGE DIRECTION ONLY
    // ------------------------------------------------------------

    label nSolved = 0;

    forAll(solvePoint, bpI)
    {
        if
        (
            !solvePoint[bpI]
         || !directionValid[bpI]
        )
            continue;

        vector d =
            solvedDir[bpI];

        const scalar dMag =
            mag(d);

        if( dMag < VSMALL )
            continue;

        d /= dMag;

        proposedDisp[bpI] =
            originalLength[bpI]*d;

        ++nSolved;
    }


    if
    (
        worstAllI >= 0
     && worstAllJ >= 0
    )
    {
        const label pI = bPoints[worstAllI];
        const label pJ = bPoints[worstAllJ];

        const vector oldI =
            originalDisp[worstAllI]
           /(mag(originalDisp[worstAllI]) + VSMALL);

        const vector oldJ =
            originalDisp[worstAllJ]
           /(mag(originalDisp[worstAllJ]) + VSMALL);

        const vector newI =
            solvedDir[worstAllI]
           /(mag(solvedDir[worstAllI]) + VSMALL);

        const vector newJ =
            solvedDir[worstAllJ]
           /(mag(solvedDir[worstAllJ]) + VSMALL);

        Info
            << "CFMITCH SEAMCOLLAR_V6_WORST_ALL"
            << " bpI=" << worstAllI
            << " pointI=" << pI
            << " ringI=" << collarRing[worstAllI]
            << " wallPatchI=" << wallPatch[worstAllI]
            << " neutralPatchI=" << neutralPatch[worstAllI]
            << " rootI=" << points[pI]
            << " oldI=" << oldI
            << " prefI=" << preferredDir[worstAllI]
            << " newI=" << newI

            << " bpJ=" << worstAllJ
            << " pointJ=" << pJ
            << " ringJ=" << collarRing[worstAllJ]
            << " wallPatchJ=" << wallPatch[worstAllJ]
            << " neutralPatchJ=" << neutralPatch[worstAllJ]
            << " rootJ=" << points[pJ]
            << " oldJ=" << oldJ
            << " prefJ=" << preferredDir[worstAllJ]
            << " newJ=" << newJ

            << " oldDot=" << (oldI & oldJ)
            << " newDot=" << (newI & newJ)
            << " oldUnitDiff=" << mag(oldI-oldJ)
            << " newUnitDiff=" << mag(newI-newJ)
            << " rootDistance="
            << mag(points[pI]-points[pJ])
            << endl;
    }

    if
    (
        worstSeamCollarI >= 0
     && worstSeamCollarJ >= 0
    )
    {
        const label pI =
            bPoints[worstSeamCollarI];

        const label pJ =
            bPoints[worstSeamCollarJ];

        const vector oldI =
            originalDisp[worstSeamCollarI]
           /(mag(originalDisp[worstSeamCollarI]) + VSMALL);

        const vector oldJ =
            originalDisp[worstSeamCollarJ]
           /(mag(originalDisp[worstSeamCollarJ]) + VSMALL);

        const vector newI =
            solvedDir[worstSeamCollarI]
           /(mag(solvedDir[worstSeamCollarI]) + VSMALL);

        const vector newJ =
            solvedDir[worstSeamCollarJ]
           /(mag(solvedDir[worstSeamCollarJ]) + VSMALL);

        Info
            << "CFMITCH SEAMCOLLAR_V6_WORST_SEAMCOLLAR"
            << " bpI=" << worstSeamCollarI
            << " pointI=" << pI
            << " ringI=" << collarRing[worstSeamCollarI]
            << " wallPatchI=" << wallPatch[worstSeamCollarI]
            << " neutralPatchI="
            << neutralPatch[worstSeamCollarI]
            << " rootI=" << points[pI]
            << " oldI=" << oldI
            << " prefI="
            << preferredDir[worstSeamCollarI]
            << " newI=" << newI

            << " bpJ=" << worstSeamCollarJ
            << " pointJ=" << pJ
            << " ringJ=" << collarRing[worstSeamCollarJ]
            << " wallPatchJ=" << wallPatch[worstSeamCollarJ]
            << " neutralPatchJ="
            << neutralPatch[worstSeamCollarJ]
            << " rootJ=" << points[pJ]
            << " oldJ=" << oldJ
            << " prefJ="
            << preferredDir[worstSeamCollarJ]
            << " newJ=" << newJ

            << " oldDot=" << (oldI & oldJ)
            << " newDot=" << (newI & newJ)
            << " oldUnitDiff=" << mag(oldI-oldJ)
            << " newUnitDiff=" << mag(newI-newJ)
            << " rootDistance="
            << mag(points[pI]-points[pJ])
            << endl;
    }

    Info
        << "CFMITCH SEAMCOLLAR V6:"
        << " ring0=" << nRing0
        << " ring1=" << nRing1
        << " ring2=" << nRing2
        << " ring3=" << nRing3
        << " components=" << nV2Components
        << " signFlips=" << nV2SignFlips
        << " smoothPasses=" << collarSmoothPasses
        << " comparedEdges=" << nAllComparedEdges
        << " seamCollarEdges=" << nSeamCollarEdges
        << " maxNeutralLeak=" << maxNeutralLeak
        << " maxAllUnitDiffBefore=" << maxAllUnitDiffBefore
        << " maxAllUnitDiffAfter=" << maxAllUnitDiffAfter
        << " maxSeamCollarDiffAfter=" << maxSeamCollarDiffAfter
        << endl;


    // ============================================================
    // FIELD COHERENCE AUDIT
    // ============================================================

    scalar maxUnitDiffBefore =
        scalar(0);

    scalar maxUnitDiffAfter =
        scalar(0);

    label nComparedEdges = 0;

    forAll(active, bpI)
    {
        if
        (
            !active[bpI]
         || mag(proposedDisp[bpI]) < VSMALL
        )
            continue;

        const vector oldU0 =
            originalDisp[bpI]
           /(mag(originalDisp[bpI]) + VSMALL);

        const vector newU0 =
            proposedDisp[bpI]
           /(mag(proposedDisp[bpI]) + VSMALL);

        forAllRow(pointPoints, bpI, ppI)
        {
            const label bpJ =
                pointPoints(bpI, ppI);

            if
            (
                bpJ <= bpI
             || bpJ < 0
             || bpJ >= nBP
             || !active[bpJ]
             || mag(proposedDisp[bpJ]) < VSMALL
            )
                continue;

            if
            (
                wallPatch[bpJ]
                    != wallPatch[bpI]
             || neutralPatch[bpJ]
                    != neutralPatch[bpI]
            )
                continue;

            const vector oldU1 =
                originalDisp[bpJ]
               /(mag(originalDisp[bpJ]) + VSMALL);

            const vector newU1 =
                proposedDisp[bpJ]
               /(mag(proposedDisp[bpJ]) + VSMALL);

            maxUnitDiffBefore =
                Foam::max
                (
                    maxUnitDiffBefore,
                    mag(oldU0-oldU1)
                );

            maxUnitDiffAfter =
                Foam::max
                (
                    maxUnitDiffAfter,
                    mag(newU0-newU1)
                );

            ++nComparedEdges;
        }
    }


    // ============================================================
    // COMMIT CANDIDATE TOP COORDINATES
    //
    // Root/original boundary points are untouched here.
    // ============================================================

    label nCommitted = 0;

    forAll(solvePoint, bpI)
    {
        if
        (
            !solvePoint[bpI]
         || mag(proposedDisp[bpI]) < VSMALL
        )
            continue;

        const label topI =
            topLabel[bpI];

        if
        (
            topI < 0
         || topI >= label(points.size())
        )
            continue;

        const label pointI =
            bPoints[bpI];

        const point root =
            points[pointI];

        const vector oldE =
            originalDisp[bpI];

        const vector newE =
            proposedDisp[bpI];

        points[topI] =
            root + newE;


        // Existing Rotor37 forensic pair.
        // Diagnostic only; no logic depends on these labels.
        if
        (
            pointI == 337809
         || pointI == 337818
         || pointI == 337820
        )
        {
            Info
                << "CFMITCH SEAMSOLVE_TARGET"
                << " bpI=" << bpI
                << " pointI=" << pointI
                << " ring=" << collarRing[bpI]
                << " wallPatch=" << wallPatch[bpI]
                << " neutralPatch="
                << neutralPatch[bpI]
                << " oldE=" << oldE
                << " oldL=" << mag(oldE)
                << " rawNeutralN="
                << rawNeutralNormal[bpI]
                << " smoothNeutralN="
                << neutralNormal[bpI]
                << " newE=" << newE
                << " newL=" << mag(newE)
                << " oldNewUnitDiff="
                << mag
                   (
                       oldE
                      /(mag(oldE)+VSMALL)
                     -
                       newE
                      /(mag(newE)+VSMALL)
                   )
                << endl;
        }

        ++nCommitted;
    }


    Info
        << "CFMITCH SEAMSOLVE V6:"
        << " neutralSeeds="
        << blNeutralEdgePoints_.size()
        << " initiallyEligible="
        << nInitiallyEligible
        << " solved=" << nSolved
        << " committed=" << nCommitted
        << " graphPruned=" << nGraphPruned
        << " badRole=" << nBadRole
        << " specialSkipped=" << nSpecialSkipped
        << " noTop=" << nNoTop
        << " badNormal=" << nBadNormal
        << " badTangent=" << nBadTangent
        << " badDirection=" << nBadDirection
        << " comparedEdges=" << nComparedEdges
        << " maxUnitDiffBefore="
        << maxUnitDiffBefore
        << " maxUnitDiffAfter="
        << maxUnitDiffAfter
        << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::createNewVertices(const boolList& treatPatches)
{
    Info << "Creating vertices for layer cells" << endl;

    List<direction> patchVertex;
    findPatchVertices(treatPatches, patchVertex);

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();

    //- the following is needed for parallel runs
    //- it is ugly, but must stay for now :(
    if( Pstream::parRun() )
    {
        mse.pointNormals();
        mse.pointPoints();
    }

    pointFieldPMG& points = mesh_.points();

    label nExtrudedVertices(0);
    forAll(patchVertex, bpI)
        if( patchVertex[bpI] )
            ++nExtrudedVertices;

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const label nPatches = boundaries.size();
    const meshSurfacePartitioner& mPartDiag = surfacePartitioner();
    const VRWGraph& pPatchesDiag = mPartDiag.pointPatches();
    List<label> nDiagPts(nPatches, 0);
    List<label> nDiagMultiPts(nPatches, 0);
    List<label> nDiagLeVSmall(nPatches, 0);
    List<label> nDiagLe1e10(nPatches, 0);
    List<label> nDiagLe1e8(nPatches, 0);
    List<label> nDiagLe1e6(nPatches, 0);
    List<scalar> minDiagDist(nPatches, GREAT);
    List<scalar> sumDiagDist(nPatches, 0.0);

    points.setSize(points.size() + nExtrudedVertices);

    labelLongList procPoints;
    forAll(bPoints, bpI)
        if( patchVertex[bpI] )
        {
            if( patchVertex[bpI] & PARALLELBOUNDARY )
            {
                procPoints.append(bpI);
                continue;
            }

            const point oldPt = points[bPoints[bpI]];
            const point newPt = createNewVertex(bpI, treatPatches, patchVertex);
            const scalar extrudeDist = mag(newPt - oldPt);
            const bool isMultiPatch =
                bpI >= 0 && bpI < label(pPatchesDiag.size())
             && pPatchesDiag.sizeOfRow(bpI) > 1;
            if( bpI >= 0 && bpI < label(pPatchesDiag.size()) )
            {
                forAllRow(pPatchesDiag, bpI, ppi)
                {
                    const label patchI = pPatchesDiag(bpI, ppi);
                    if( patchI < 0 || patchI >= nPatches ) continue;
                    ++nDiagPts[patchI];
                    if( isMultiPatch ) ++nDiagMultiPts[patchI];
                    minDiagDist[patchI] = Foam::min(minDiagDist[patchI], extrudeDist);
                    sumDiagDist[patchI] += extrudeDist;
                    if( extrudeDist <= VSMALL )        ++nDiagLeVSmall[patchI];
                    if( extrudeDist <= scalar(1e-10) ) ++nDiagLe1e10[patchI];
                    if( extrudeDist <= scalar(1e-8) )  ++nDiagLe1e8[patchI];
                    if( extrudeDist <= scalar(1e-6) )  ++nDiagLe1e6[patchI];
                }
            }
            points[nPoints_] = newPt;
            newLabelForVertex_[bPoints[bpI]] = nPoints_;
            ++nPoints_;
        }

    Info << "BL extrusion-distance audit per patch:" << endl;
    forAll(boundaries, patchI)
    {
        if( nDiagPts[patchI] == 0 ) continue;
        const scalar avgDist =
            sumDiagDist[patchI] / scalar(Foam::max(label(1), nDiagPts[patchI]));
        Info << "  " << boundaries[patchI].patchName() << ":"
             << " pts=" << nDiagPts[patchI]
             << " multiPts=" << nDiagMultiPts[patchI]
             << " minDist=" << minDiagDist[patchI]
             << " avgDist=" << avgDist
             << " <=VSMALL=" << nDiagLeVSmall[patchI]
             << " <=1e-10=" << nDiagLe1e10[patchI]
             << " <=1e-8=" << nDiagLe1e8[patchI]
             << " <=1e-6=" << nDiagLe1e6[patchI]
             << endl;
    }

    if( Pstream::parRun() )
    {
        createNewPartitionVerticesParallel
        (
            procPoints,
            patchVertex,
            treatPatches
        );

        createNewEdgeVerticesParallel
        (
            procPoints,
            patchVertex,
            treatPatches
        );
    }

    // CFMitch connected neutral-seam constraint.
    //
    // Candidate tops exist, but original boundary geometry remains
    // untouched and authoritative.
    constrainNeutralSeamCandidatesBeforeSwap();

    //- swap coordinates of new and old points
    forAll(bPoints, bpI)
    {
        const label pLabel = newLabelForVertex_[bPoints[bpI]];
        if( pLabel != -1 )
        {
            const point p = points[pLabel];
            points[pLabel] = points[bPoints[bpI]];
            points[bPoints[bpI]] = p;
        }
    }

    if( nPoints_ != points.size() )
        FatalErrorIn
        (
            "void boundaryLayers::createNewVertices("
            "const meshSurfaceEngine& mse,"
            "const boolList& treatPatches,"
            "labelList& newLabelForVertex)"
        ) << "Number of vertices " << nPoints_
            << " does not match the list size "
            << abort(FatalError);

    // Local topology-aware layer rollback.
    // Restricted to topology-sensitive points only:
    // BL/neutral edges, BL/no-BL edges, BL/BL junction points.
    // Uses correct sign logic and area-scaled tolerance.
    {
        const label maxRollbackIter = 5;
        const scalar dampFactor = 0.5;
        const VRWGraph& ptFacesRB = mse.pointFaces();
        const faceList::subList& bFacesRB = mse.boundaryFaces();
        label nRolledBack = 0;

        // Build restricted rollback set: topology-sensitive points only
        labelHashSet rollbackSet;
        forAllConstIter(labelHashSet, blNeutralEdgePoints_, it)
            rollbackSet.insert(it.key());
        forAllConstIter(labelHashSet, blNoBlEdgePoints_, it)
            rollbackSet.insert(it.key());
        forAllConstIter(labelHashSet, blblJunctionPoints_, it)
            rollbackSet.insert(it.key());

        for(label iter=0; iter<maxRollbackIter; ++iter)
        {
            label nBad = 0;
            forAll(bPoints, bpI)  // Fix: was procPoints (serial no-op)
            {
                if( !rollbackSet.found(bpI) ) continue;
                const label meshPtI = bPoints[bpI];
                const label origPtI = newLabelForVertex_[meshPtI];
                if( origPtI < 0 ) continue;

                // Value copies -- avoid reference aliasing during mutation
                const point layerPt = points[meshPtI];
                const point basePt  = points[origPtI];

                bool bad = false;
                forAllRow(ptFacesRB, bpI, pfI)
                {
                    const face& f = bFacesRB[ptFacesRB(bpI, pfI)];
                    point fc = point::zero;
                    forAll(f, fi) fc += points[f[fi]];
                    fc /= scalar(f.size());
                    vector fn = vector::zero;
                    const point& fp0 = points[f[0]];
                    for(label pi=1; pi<f.size()-1; ++pi)
                        fn += (points[f[pi]]-fp0)^(points[f[pi+1]]-fp0);
                    const scalar areaMag = mag(fn);
                    if( areaMag < VSMALL ) continue;
                    const scalar tol = 1e-12 * areaMag;
                    const scalar volLayer = (layerPt - fc) & fn;
                    const scalar volBase  = (basePt  - fc) & fn;
                    // Bad: base and layer on opposite sides of face plane
                    if( mag(volBase) > tol && volBase*volLayer < -tol )
                    { bad = true; break; }
                }

                if( bad )
                {
                    points[meshPtI] = dampFactor*layerPt
                                   + (1.0-dampFactor)*basePt;
                    ++nBad;
                    ++nRolledBack;
                }
            }
            if( nBad == 0 ) break;
        }

        if( nRolledBack > 0 )
            Info << "Layer rollback: " << nRolledBack
                 << " topology-sensitive vertices relaxed" << endl;
    }

    Info << "Finished creating layer vertices" << endl;
}

void boundaryLayers::createNewVertices(const labelList& patchLabels)
{
    otherVrts_.clear();

    patchKey_.setSize(mesh_.boundaries().size());
    patchKey_ = -1;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    //- populate zeroDistPoints_ for BL-transition vertices
    if( terminateLayersAtConcaveEdges_ )
    {
        boolList _skipPoint(bPoints.size(), false);
        markConcaveEdgePoints(_skipPoint);
    }

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    //- the following is needed for parallel runs
    //- it is ugly, but must stay for now :(
    mse.boundaryFaces();
    mse.pointNormals();
    mse.pointFaces();
    const VRWGraph& pointPoints = mse.pointPoints();

    pointFieldPMG& points = mesh_.points();
    boolList treatPatches(mesh_.boundaries().size());
    List<direction> patchVertex(bPoints.size());

    //- make sure than the points are never re-allocated during the process
    points.reserve(points.size() + 2 * bPoints.size());

    const PtrList<boundaryPatch>& boundaries2 = mesh_.boundaries();
    const label nPatches2 = boundaries2.size();
    const VRWGraph& pPatchesDiag2 = mPart.pointPatches();
    List<label> nDiagPts2(nPatches2, 0);
    List<label> nDiagLeVSmall2(nPatches2, 0);
    List<label> nDiagLe1e102(nPatches2, 0);
    List<label> nDiagLe1e82(nPatches2, 0);
    List<label> nDiagLe1e62(nPatches2, 0);
    List<scalar> minDiagDist2(nPatches2, GREAT);
    List<scalar> sumDiagDist2(nPatches2, 0.0);

    //- generate new layer vertices for each patch
    forAll(patchLabels, patchI)
    {
        const label pLabel = patchLabels[patchI];
        treatPatches = false;

        bool treat(true);
        forAll(treatPatchesWithPatch_[pLabel], pI)
        {
            const label otherPatch = treatPatchesWithPatch_[pLabel][pI];
            treatPatches[otherPatch] = true;

            if( patchKey_[otherPatch] == -1 )
            {
                patchKey_[otherPatch] = patchI;
            }
            else
            {
                Info << "BLGROUPCOLLISION loopI=" << patchI
                     << " pName=" << boundaries2[pLabel].patchName()
                     << " member=" << boundaries2[otherPatch].patchName()
                     << " existingKey=" << patchKey_[otherPatch];
                if( patchKey_[otherPatch] >= 0
                 && patchKey_[otherPatch] < label(patchLabels.size()) )
                    Info << " existingRepName="
                         << boundaries2[patchLabels[patchKey_[otherPatch]]].patchName();
                Info << endl;
                treat = false;
            }
        }

        if( !treat )
        {
            Info << "BLGROUPSKIP pName=" << boundaries2[pLabel].patchName()
                 << " group=";
            forAll(treatPatchesWithPatch_[pLabel], gI)
            {
                const label m = treatPatchesWithPatch_[pLabel][gI];
                Info << " [" << boundaries2[m].patchName()
                     << " key=" << patchKey_[m] << "]";
            }
            Info << endl;
            continue;
        }

        const label pKey = patchKey_[pLabel];

        //- classify vertices belonging to this patch
        findPatchVertices(treatPatches, patchVertex);

        //- create indices and allocate maps for new points
        labelLongList procPoints, patchPoints;
        forAll(bPoints, bpI)
        {
            if( !patchVertex[bpI] )
                continue;

            //- skip vertices at parallel boundaries
            if( patchVertex[bpI] & PARALLELBOUNDARY )
            {
                procPoints.append(bpI);

                continue;
            }

            patchPoints.append(bpI);
            const label pointI = bPoints[bpI];

            if( patchVertex[bpI] & EDGENODE )
            {
                if( otherVrts_.find(pointI) == otherVrts_.end() )
                {
                    std::map<std::pair<label, label>, label> m;

                    otherVrts_.insert(std::make_pair(pointI, m));
                }

                std::pair<label, label> pr(pKey, pKey);
                otherVrts_[pointI].insert(std::make_pair(pr, nPoints_++));
            }
            else
            {
                //- this the only new point
                newLabelForVertex_[pointI] = nPoints_++;
            }
        }

        //- set the size of points
        points.setSize(nPoints_);

        //- calculate coordinates of new points
        // Serial: otherVrts_ is shared std::map -- concurrent find/insert/[]
        // is not thread-safe. Races corrupt BL edge/corner vertex positions.
        forAll(patchPoints, i)
        {
            const label bpI = patchPoints[i];

            const label pointI = bPoints[bpI];

            //- create new point
            const point pBefore = points[bPoints[bpI]];
            const point p = createNewVertex(bpI, treatPatches, patchVertex);
            const scalar eDist2 = mag(p - pBefore);
            if( bpI >= 0 && bpI < label(pPatchesDiag2.size()) )
            {
                forAllRow(pPatchesDiag2, bpI, ppi)
                {
                    const label pIdx2 = pPatchesDiag2(bpI, ppi);
                    if( pIdx2 < 0 || pIdx2 >= nPatches2 ) continue;
                    ++nDiagPts2[pIdx2];
                    minDiagDist2[pIdx2] = Foam::min(minDiagDist2[pIdx2], eDist2);
                    sumDiagDist2[pIdx2] += eDist2;
                    if( eDist2 <= VSMALL )        ++nDiagLeVSmall2[pIdx2];
                    if( eDist2 <= scalar(1e-10) ) ++nDiagLe1e102[pIdx2];
                    if( eDist2 <= scalar(1e-8) )  ++nDiagLe1e82[pIdx2];
                    if( eDist2 <= scalar(1e-6) )  ++nDiagLe1e62[pIdx2];
                }
            }

            if( patchVertex[bpI] & EDGENODE )
            {
                //- set the new point or an edge point
                if( otherVrts_.find(pointI) == otherVrts_.end() )
                {
                    std::map<std::pair<label, label>, label> m;

                    otherVrts_.insert(std::make_pair(pointI, m));
                }

                std::pair<label, label> pr(pKey, pKey);
                std::map<std::pair<label,label>,label>::const_iterator npIter =
                    otherVrts_[pointI].find(pr);
                if( npIter == otherVrts_[pointI].end() )
                {
                    FatalErrorIn("boundaryLayers::createNewVertices")
                        << "Missing otherVrts_ entry for point " << pointI
                        << abort(FatalError);
                }
                const label npI = npIter->second;
                points[npI] = p;
            }
            else
            {
                //- set the new point
                points[newLabelForVertex_[pointI]] = p;
            }
        }

        if( Pstream::parRun() )
        {
            points.setSize(nPoints_+procPoints.size());

            createNewPartitionVerticesParallel
            (
                procPoints,
                patchVertex,
                treatPatches
            );

            createNewEdgeVerticesParallel
            (
                procPoints,
                patchVertex,
                treatPatches
            );
        }
    }

    //- create missing vertices for edge and corner vertices
    //- they should be stored in the otherNodes map
    forAll(bPoints, bpI)
    {
        const label pointI = bPoints[bpI];

        if( otherVrts_.find(pointI) == otherVrts_.end() )
            continue;

        const point& p = points[pointI];

        std::map<std::pair<label, label>, label>& m = otherVrts_[pointI];
        DynList<label> usedPatches;
        DynList<label> newNodeLabel;
        DynList<vector> newPatchPenetrationVector;
        forAllRow(pPatches, bpI, patchI)
        {
            const label pKey = patchKey_[pPatches(bpI, patchI)];
            const std::pair<label, label> pr(pKey, pKey);
            const std::map<std::pair<label, label>, label>::const_iterator it =
                m.find(pr);
            if( (it != m.end()) && !usedPatches.contains(pKey) )
            {
                usedPatches.append(pKey);
                newNodeLabel.append(it->second);
                newPatchPenetrationVector.append(points[it->second] - p);
            }
        }

        if( newNodeLabel.size() == 1 )
        {
            //- only one patch is treated
            //- PRE-CLEANUP AUDIT: cap candidate erased from otherVrts_
            //- because only one patch survived edge/corner vertex cleanup.
            if
            (
                useHardBLBLCapCells_
             && !blCapCellEndwallPatch_.empty()
             && blCapCellEndwallPatch_.found(bpI)
            )
            {
                static label nCapErased = 0;
                ++nCapErased;
                if( nCapErased <= 20 )
                {
                    const label ewPatch = blCapCellEndwallPatch_[bpI];
                    const label blPatch =
                        blCapCellBladePatch_.found(bpI) ?
                        blCapCellBladePatch_[bpI] : -1;
                    Info << "CAP_PRECLEANUP_ERASE:"
                         << " bpI=" << bpI
                         << " pointI=" << pointI
                         << " endwallPatch="
                         << ((ewPatch >= 0 && ewPatch < label(patchNames_.size()))
                                ? patchNames_[ewPatch] : word("?"))
                         << " bladePatch="
                         << ((blPatch >= 0 && blPatch < label(patchNames_.size()))
                                ? patchNames_[blPatch] : word("?"))
                         << " otherVrtsSize=" << label(otherVrts_[pointI].size())
                         << " newNodeLabel.size=" << newNodeLabel.size()
                         << endl;
                }
                if( nCapErased == 20 )
                    Info << "CAP_PRECLEANUP_ERASE: (further suppressed)" << endl;
            }
            newLabelForVertex_[pointI] = newNodeLabel[0];
            otherVrts_.erase(pointI);
        }
        else if( newNodeLabel.size() == 2 )
        {
            //- point is located at an extrusion edge
            //- create the new position for the existing point
            point newP(p);
            newP += newPatchPenetrationVector[0];
            newP += newPatchPenetrationVector[1];

            if( !help::isnan(newP) && !help::isinf(newP) )
            {
                points.append(newP);
            }
            else
            {
                points.append(p);
            }
            newLabelForVertex_[pointI] = nPoints_;
            ++nPoints_;
        }
        else if( newNodeLabel.size() == 3 )
        {
            //- point is located at an extrusion corner
            //- create 3 points and the new position for the existing point
            point newP(p);
            for(label i=0;i<3;++i)
            {
                newP += newPatchPenetrationVector[i];
                for(label j=i+1;j<3;++j)
                {
                    const point np =
                        p + newPatchPenetrationVector[i] +
                        newPatchPenetrationVector[j];

                    if( !help::isnan(np) && !help::isinf(np) )
                    {
                        points.append(np);
                    }
                    else
                    {
                        points.append(p);
                    }

                    m.insert
                    (
                        std::make_pair
                        (
                            std::make_pair(usedPatches[i], usedPatches[j]),
                            nPoints_
                        )
                    );
                    ++nPoints_;
                }
            }

            //- create new position for the existing point
            if( !help::isnan(newP) && !help::isinf(newP) )
            {
                points.append(newP);
            }
            else
            {
                points.append(p);
            }

            newLabelForVertex_[pointI] = nPoints_;
            ++nPoints_;
        }
        else
        {
            FatalErrorIn
            (
                "void boundaryLayers::createNewVertices("
                "const labelList& patchLabels, labelLongList& newLabelForVertex,"
                "std::map<label, std::map<std::pair<label, label>, label> >&)"
            ) << "Boundary node " << bpI << " is not at an edge!"
                << abort(FatalError);
        }
    }

    Info << "BL extrusion-distance audit (labelList path) per patch:" << endl;
    forAll(boundaries2, patchI)
    {
        if( nDiagPts2[patchI] == 0 ) continue;
        const scalar avgD =
            sumDiagDist2[patchI] / scalar(Foam::max(label(1), nDiagPts2[patchI]));
        Info << "  " << boundaries2[patchI].patchName() << ":"
             << " pts=" << nDiagPts2[patchI]
             << " minDist=" << minDiagDist2[patchI]
             << " avgDist=" << avgD
             << " <=VSMALL=" << nDiagLeVSmall2[patchI]
             << " <=1e-10=" << nDiagLe1e102[patchI]
             << " <=1e-8=" << nDiagLe1e82[patchI]
             << " <=1e-6=" << nDiagLe1e62[patchI]
             << endl;
    }

    // CFMitch connected neutral-seam constraint.
    //
    // Candidate tops exist, but original boundary geometry remains
    // untouched and authoritative.
    constrainNeutralSeamCandidatesBeforeSwap();

    //- swap coordinates of new and old points
    // Serial: OMP parallelism here causes coordinate corruption at BL/BL
    // junctions -- primary source of 67-114 bad pyramid faces per run.
    // Build the exact union treatment mask used by createLayerCells().
    // Do not use treatedPatch_ here: patchwise layer generation constructs
    // its active treatment set from patchLabels/treatPatchesWithPatch_.
    boolList contactSweepTreatPatches(mesh_.boundaries().size(), false);

    forAll(patchLabels, cspI)
    {
        const label pLabel = patchLabels[cspI];

        if
        (
            pLabel < 0
         || pLabel >= label(treatPatchesWithPatch_.size())
        )
            continue;

        forAll(treatPatchesWithPatch_[pLabel], gI)
        {
            const label patchI =
                treatPatchesWithPatch_[pLabel][gI];

            if
            (
                patchI >= 0
             && patchI < label(contactSweepTreatPatches.size())
            )
                contactSweepTreatPatches[patchI] = true;
        }
    }

    label nContactSweepTreatPatches = 0;
    forAll(contactSweepTreatPatches, patchI)
        if( contactSweepTreatPatches[patchI] )
            ++nContactSweepTreatPatches;

    Info << "CONTACT_SWEEP_PATCH_MASK activePatches="
         << nContactSweepTreatPatches << endl;

    // Snapshot topology references while surfaceEngine still corresponds
    // to the pre-swap mesh state.  The contact smoother is deliberately
    // forbidden from rebuilding surfaceEngine after coordinates are swapped.
    const meshSurfaceEngine& mseContactSweep = surfaceEngine();
    const faceList::subList& bFacesContactSweep =
        mseContactSweep.boundaryFaces();
    const VRWGraph& pointFacesContactSweep =
        mseContactSweep.pointFaces();
    const labelList& boundaryFacePatchesContactSweep =
        mseContactSweep.boundaryFacePatches();

    forAll(bPoints, bpI)
    {
        const label pLabel = newLabelForVertex_[bPoints[bpI]];

        if( pLabel != -1 )
        {
            const point p = points[pLabel];
            points[pLabel] = points[bPoints[bpI]];
            points[bPoints[bpI]] = p;
        }
    }

    // ------------------------------------------------------------
    // CFMitch v2: capture the raw as-born BL height field.
    //
    // This intentionally runs AFTER coordinate swap, so meshPtI is
    // the extruded point and newLabelForVertex_[meshPtI] is its base,
    // but BEFORE the historical contact-line smoother or rollback
    // changes the geometry.
    //
    // Diagnostic only.
    // ------------------------------------------------------------
    if( cfmitchHeightDiagnostics_ )
    {
        writeCFMitchHeightAtlas
        (
            "preSmoother",
            bPoints,
            pointPoints,
            contactSweepTreatPatches,
            boundaryFacePatchesContactSweep,
            pointFacesContactSweep
        );
    }

    if( cfmitchHeightSmoothing_ )
    {
        applyCFMitchHeightCompatibility
        (
            "graphProjection",
            bPoints,
            pointPoints,
            contactSweepTreatPatches,
            boundaryFacePatchesContactSweep,
            pointFacesContactSweep
        );

        if( cfmitchHeightDiagnostics_ )
        {
            writeCFMitchHeightAtlas
            (
                "postCFMitch",
                bPoints,
                pointPoints,
                contactSweepTreatPatches,
                boundaryFacePatchesContactSweep,
                pointFacesContactSweep
            );
        }
    }

    // Contact-line height limiter: runs after coordinate swap,
    // before rollback. Detects and limits local height spikes/collapses
    // along contact-line points. Diagnostic atlas always written;
    // point movement gated by useContactLineHeightSmoother.
    if( useContactLineHeightSmoother_ || writeContactLineHeightSmootherAtlas_ )
    {
        for( label iter=0; iter<contactLineSmootherIterations_; ++iter )
            smoothContactLineHeights
            (
                "pass2_iter" + Foam::name(iter),
                bPoints,
                pointPoints,
                contactSweepTreatPatches,
                bFacesContactSweep,
                pointFacesContactSweep,
                boundaryFacePatchesContactSweep
            );
    }

    // Local topology-aware layer rollback with minimum-height floor.
    // Relaxes inverted layer points toward base, but clamps so a point never
    // drops below rollbackMinHeightFraction of its intended extrusion height,
    // preventing collapse-to-wall slivers at neutral/periodic/junction edges.
    {
        const meshSurfaceEngine& mseRB = surfaceEngine();
        const VRWGraph& ptFacesRB = mseRB.pointFaces();
        const faceList::subList& bFacesRB = mseRB.boundaryFaces();
        const label maxRollbackIter = 5;
        const scalar dampFactor = 0.5;
        const scalar rollbackMinHeightFraction = 0.30;
        label nRolledBack = 0;
        Info << "Layer rollback (pass2): minHeightFraction="
             << rollbackMinHeightFraction << endl;

        labelHashSet rollbackSet;
        forAllConstIter(labelHashSet, blNeutralEdgePoints_, it)
            rollbackSet.insert(it.key());
        forAllConstIter(labelHashSet, blNoBlEdgePoints_, it)
            rollbackSet.insert(it.key());
        forAllConstIter(labelHashSet, blblJunctionPoints_, it)
            rollbackSet.insert(it.key());

        for(label iter=0; iter<maxRollbackIter; ++iter)
        {
            label nBad = 0;
            forAll(bPoints, bpI)
            {
                if( !rollbackSet.found(bpI) ) continue;
                const label meshPtI = bPoints[bpI];
                const label origPtI = newLabelForVertex_[meshPtI];
                if( origPtI < 0 ) continue;
                const point layerPt = points[meshPtI];
                const point basePt  = points[origPtI];
                bool bad = false;
                forAllRow(ptFacesRB, bpI, pfI)
                {
                    const face& f = bFacesRB[ptFacesRB(bpI, pfI)];
                    point fc = point::zero;
                    forAll(f, fi) fc += points[f[fi]];
                    fc /= scalar(f.size());
                    vector fn = vector::zero;
                    const point& fp0 = points[f[0]];
                    for(label pi=1; pi<f.size()-1; ++pi)
                        fn += (points[f[pi]]-fp0)^(points[f[pi+1]]-fp0);
                    const scalar areaMag = mag(fn);
                    if( areaMag < VSMALL ) continue;
                    const scalar tol = 1e-12 * areaMag;
                    const scalar volLayer = (layerPt - fc) & fn;
                    const scalar volBase  = (basePt  - fc) & fn;
                    if( mag(volBase) > tol && volBase*volLayer < -tol )
                    { bad = true; break; }
                }
                if( bad )
                {
                    const vector h = layerPt - basePt;
                    const scalar hMag = mag(h);
                    point relaxed = dampFactor*layerPt
                                  + (1.0-dampFactor)*basePt;
                    if( hMag > VSMALL )
                    {
                        const vector r = relaxed - basePt;
                        const scalar rMag = mag(r);
                        const scalar minMag = rollbackMinHeightFraction*hMag;
                        if( rMag < minMag )
                            relaxed = basePt + (minMag/hMag)*h;
                    }
                    points[meshPtI] = relaxed;
                    ++nBad;
                    ++nRolledBack;
                }
            }
            if( nBad == 0 ) break;
        }
        if( nRolledBack > 0 )
            Info << "Layer rollback (pass2): " << nRolledBack
                 << " topology-sensitive vertices relaxed (floored)" << endl;
    }

    // Resolve patchKey_ for cap cell candidates now that patchKey_ is populated.
    if( useHardBLBLCapCells_ && !blCapCellEndwallPatch_.empty() )
    {
        forAllConstIter(Map<label>, blCapCellEndwallPatch_, it)
        {
            const label bpI = it.key();
            const label ewPatch = it();
            const label blPatch =
                blCapCellBladePatch_.found(bpI) ?
                blCapCellBladePatch_[bpI] : -1;
            if( ewPatch >= 0 && ewPatch < label(patchKey_.size()) )
                blCapCellEndwallPKey_[bpI] = patchKey_[ewPatch];
            if( blPatch >= 0 && blPatch < label(patchKey_.size()) )
                blCapCellBladePKey_[bpI] = patchKey_[blPatch];
        }

        //- Step 2/4: populate capSideVrtMap_.
        //- Step 2 (useHardBLBLCapVertexInsertion=false): no-op, same labels.
        //- Step 4 (useHardBLBLCapVertexInsertion=true): real asymmetric vertices.
        //- Safety interlock in BL.C prevents step4 without reducedCells.
        capSideVrtMap_.clear();
        label nCapSideInserted = 0;
        label nCapSideMissing = 0;
        label nCapEWCreated = 0;
        const meshSurfaceEngine& mseCap = surfaceEngine();
        const labelList& bPointsCap = mseCap.boundaryPoints();
        const faceList::subList& bFacesCap = mseCap.boundaryFaces();
        const labelList& bFacePatchesCap = mseCap.boundaryFacePatches();
        const VRWGraph& pFacesCap = mseCap.pointFaces();
        const VRWGraph& ppCap = mseCap.pointPoints();
        OFstream* capAtlasOsPtr = nullptr;
        if( useHardBLBLCapVertexInsertion_ )
        {
            capAtlasOsPtr = new OFstream("blCapCellVertexAtlas.csv");
            *capAtlasOsPtr << "bpI,pointI,endwallPatch,bladePatch,ewPKey,blPKey,"
                          << "ewVertLabel,blVertLabel,rawDist,capDist,"
                          << "ewX,ewY,ewZ" << nl;
        }
        forAllConstIter(Map<label>, blCapCellEndwallPatch_, it2)
        {
            const label bpI = it2.key();
            if( bpI < 0 || bpI >= label(bPointsCap.size()) )
            { ++nCapSideMissing; continue; }
            const label pointI = bPointsCap[bpI];
            const label ewPatch = it2();
            const label bladePatch =
                blCapCellBladePatch_.found(bpI) ?
                blCapCellBladePatch_[bpI] : -1;
            const label ewPKey =
                blCapCellEndwallPKey_.found(bpI) ?
                blCapCellEndwallPKey_[bpI] : -1;
            const label blPKey =
                blCapCellBladePKey_.found(bpI) ?
                blCapCellBladePKey_[bpI] : -1;
            if( ewPKey < 0 || blPKey < 0 )
            { ++nCapSideMissing; continue; }

            if( !useHardBLBLCapVertexInsertion_ )
            {
                //- Step 2: no-op -- key by patchI not pKey
                //- hub+blade share pKey=0, so pKey cannot distinguish sides
                const label ewOld = newLabelForVertex_[pointI];
                const label blOld = newLabelForVertex_[pointI];
                if( ewOld >= 0 )
                { capSideVrtMap_[std::make_pair(pointI,ewPatch)] = ewOld; ++nCapSideInserted; }
                else { ++nCapSideMissing; }
                if( blOld >= 0 )
                { capSideVrtMap_[std::make_pair(pointI,bladePatch)] = blOld; ++nCapSideInserted; }
                else { ++nCapSideMissing; }
                continue;
            }

            //- Step 4: real asymmetric cap vertices
            const label blVertLabel = newLabelForVertex_[pointI];
            if( blVertLabel < 0 || blVertLabel >= label(points.size()) )
            { ++nCapSideMissing; continue; }
            const point& p0 = points[pointI];
            vector ewNormal = vector::zero;
            forAllRow(pFacesCap, bpI, pfI)
            {
                const label bfI2 = pFacesCap(bpI, pfI);
                if( bfI2 < 0 || bfI2 >= label(bFacePatchesCap.size()) ) continue;
                if( bFacePatchesCap[bfI2] != ewPatch ) continue;
                const face& f2 = bFacesCap[bfI2];
                if( f2.size() < 3 ) continue;
                vector n = vector::zero;
                const point& fp0 = points[f2[0]];
                for(label i=1; i<f2.size()-1; ++i)
                    n += (points[f2[i]]-fp0)^(points[f2[i+1]]-fp0);
                ewNormal += n;
            }
            const scalar magN = mag(ewNormal);
            if( magN < VSMALL ) { ++nCapSideMissing; continue; }
            ewNormal /= magN;
            scalar rawDist = GREAT;
            forAllRow(ppCap, bpI, ppI)
            {
                const label bpJ = ppCap(bpI, ppI);
                if( bpJ < 0 || bpJ >= label(bPointsCap.size()) ) continue;
                const scalar d = 0.5*mag(points[bPointsCap[bpJ]] - p0);
                rawDist = Foam::min(rawDist, d);
            }
            if( rawDist >= GREAT || rawDist < VSMALL ) { ++nCapSideMissing; continue; }
            const scalar capDist = hardBLBLCapScale_ * rawDist;
            const point ewPt = p0 - capDist * ewNormal;
            if( help::isnan(ewPt) || help::isinf(ewPt) ) { ++nCapSideMissing; continue; }
            const label ewVertLabel = nPoints_++;
            if( ewVertLabel >= label(points.size()) )
                points.setSize(ewVertLabel + 1);
            points[ewVertLabel] = ewPt;
            //- Key by patchI not pKey (hub+blade share pKey=0)
            capSideVrtMap_[std::make_pair(pointI, ewPatch)] = ewVertLabel;
            capSideVrtMap_[std::make_pair(pointI, bladePatch)] = blVertLabel;
            nCapSideInserted += 2;
            ++nCapEWCreated;
            if( capAtlasOsPtr )
            {
                const word ewName = (ewPatch>=0&&ewPatch<label(patchNames_.size())) ?
                    patchNames_[ewPatch] : word("?");
                const word blName = (bladePatch>=0&&bladePatch<label(patchNames_.size())) ?
                    patchNames_[bladePatch] : word("?");
                *capAtlasOsPtr << bpI << "," << pointI << ","
                    << ewName << "," << blName << ","
                    << ewPKey << "," << blPKey << ","
                    << ewVertLabel << "," << blVertLabel << ","
                    << rawDist << "," << capDist << ","
                    << ewPt.x() << "," << ewPt.y() << "," << ewPt.z() << nl;
            }
        }
        delete capAtlasOsPtr; capAtlasOsPtr = nullptr;
        if( useHardBLBLCapVertexInsertion_ )
            Info << "BL cap side vertex map step4: inserted=" << nCapSideInserted
                 << " missing=" << nCapSideMissing
                 << " ewCreated=" << nCapEWCreated
                 << " (asymmetric cap vertices)" << endl;
        else
            Info << "BL cap side vertex map step2: inserted=" << nCapSideInserted
                 << " missing=" << nCapSideMissing
                 << " (provably no-op labels)" << endl;

        writeCapCellGeometryDryRun();
        writeCapCellRoutingAtlas();
    }
}

void boundaryLayers::createNewPartitionVerticesParallel
(
    const labelLongList& procPoints,
    const List<direction>& pVertices,
    const boolList& /*treatPatches*/
)
{
    if( !Pstream::parRun() )
        return;

    if( returnReduce(procPoints.size(), sumOp<label>()) == 0 )
        return;

    const meshSurfaceEngine& mse = surfaceEngine();
    pointFieldPMG& points = mesh_.points();
    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointPoints = mse.pointPoints();
    const VRWGraph& bpAtProcs = mse.bpAtProcs();
    const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
    const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

    scalarField penetrationDistances(bPoints.size(), VGREAT);

    std::map<label, LongList<labelledScalar> > exchangeDistances;

    forAll(procPoints, pointI)
    {
        const label bpI = procPoints[pointI];
        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            if( exchangeDistances.find(neiProc) == exchangeDistances.end() )
            {
                exchangeDistances.insert
                (
                    std::make_pair(neiProc, LongList<labelledScalar>())
                );
            }
        }

        if( pVertices[bpI] & EDGENODE )
            continue;

        scalar dist(VGREAT);
        const point& p = points[bPoints[bpI]];
        forAllRow(pointPoints, bpI, ppI)
        {
            const scalar d =
                0.5 * mag(points[bPoints[pointPoints(bpI, ppI)]] - p);

            if( d < dist )
                dist = d;
        }

        penetrationDistances[bpI] = dist;

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeDistances[neiProc].append
            (
                labelledScalar(globalPointLabel[bpI], dist)
            );
        }
    }

    //- exchange distances with other processors
    LongList<labelledScalar> receivedData;
    help::exchangeMap(exchangeDistances, receivedData);
    forAll(receivedData, i)
    {
        if( !globalToLocal.found(receivedData[i].scalarLabel()) ) continue;
        const label bpI = globalToLocal[receivedData[i].scalarLabel()];

        if( penetrationDistances[bpI] > receivedData[i].value() )
            penetrationDistances[bpI] = receivedData[i].value();
    }

    //- Finally, create the points
    // Override point normals at BL/no-BL interface points to use
    // BL-side faces only. The global pointNormals() averages all
    // adjacent faces including no-BL patch faces, which tilts the
    // extrusion normal toward the no-BL surface and causes protrusions.
    vectorField pNormals = mse.pointNormals();
    if( !blNoBlEdgePoints_.empty() )
    {
        const VRWGraph& pFaces = mse.pointFaces();
        const labelList& boundaryFacePatches = mse.boundaryFacePatches();
        const faceList::subList& bFaces = mse.boundaryFaces();
        const pointFieldPMG& pts = mesh_.points();

        forAllConstIter(labelHashSet, blNoBlEdgePoints_, iter)
        {
            const label bpI = iter.key();

            // Find the BL-side patch for this point
            Map<label>::const_iterator patchIt = blNoBlPointPatch_.find(bpI);
            if( patchIt == blNoBlPointPatch_.end() || patchIt() < 0 )
                continue;
            const label blPatch = patchIt();

            // Recompute normal using only BL-side faces
            vector blNormal(vector::zero);
            forAllRow(pFaces, bpI, pfI)
            {
                if( boundaryFacePatches[pFaces(bpI, pfI)] != blPatch )
                    continue;
                const face& f = bFaces[pFaces(bpI, pfI)];
                vector fn = vector::zero;
                const point& p0 = pts[f[0]];
                for(label pi=1; pi<f.size()-1; ++pi)
                    fn += (pts[f[pi]]-p0)^(pts[f[pi+1]]-p0);
                blNormal += fn;
            }
            const scalar magN = mag(blNormal);
            if( magN > VSMALL )
                pNormals[bpI] = blNormal / magN;
        }
    }

    // Override point normals at BL/neutral interface points (blade/periodic)
    // Same logic as BL/no-BL: use only BL-side faces to compute normal.
    // Prevents extrusion direction tilting toward periodic plane.
    if( !blNeutralEdgePoints_.empty() )
    {
        const VRWGraph& pFaces2 = mse.pointFaces();
        const labelList& bFacePatches2 = mse.boundaryFacePatches();
        const faceList::subList& bFaces2 = mse.boundaryFaces();
        const pointFieldPMG& pts2 = mesh_.points();

        forAllConstIter(labelHashSet, blNeutralEdgePoints_, iter)
        {
            const label bpI = iter.key();

            Map<label>::const_iterator patchIt = blNeutralPointPatch_.find(bpI);
            if( patchIt == blNeutralPointPatch_.end() || patchIt() < 0 )
                continue;
            const label blPatch = patchIt();

            vector blNormal(vector::zero);
            forAllRow(pFaces2, bpI, pfI)
            {
                if( bFacePatches2[pFaces2(bpI, pfI)] != blPatch )
                    continue;
                const face& f = bFaces2[pFaces2(bpI, pfI)];
                vector fn = vector::zero;
                const point& p0 = pts2[f[0]];
                for(label pi=1; pi<f.size()-1; ++pi)
                    fn += (pts2[f[pi]]-p0)^(pts2[f[pi+1]]-p0);
                blNormal += fn;
            }
            const scalar magN = mag(blNormal);
            if( magN > VSMALL )
                pNormals[bpI] = blNormal / magN;
        }
    }

    forAll(procPoints, pointI)
    {
        const label bpI = procPoints[pointI];

        if( pVertices[bpI] & EDGENODE )
            continue;

        const point& p = points[bPoints[bpI]];
        scalar layerDist = penetrationDistances[bpI];
        if( terminateLayersAtConcaveEdges_
         && layerScale_.size() > bpI )
            layerDist *= layerScale_[bpI];
        const point np = p - pNormals[bpI] * layerDist;
        if( !help::isnan(np) && !help::isinf(np) )
        {
            points[nPoints_] = np;
        }
        else
        {
            points[nPoints_] = p;
        }
        newLabelForVertex_[bPoints[bpI]] = nPoints_;
        ++nPoints_;
    }
}

void boundaryLayers::createNewEdgeVerticesParallel
(
    const labelLongList& procPoints,
    const List<direction>& pVertices,
    const boolList& treatPatches
)
{
    if( !Pstream::parRun() )
        return;

    if( returnReduce(procPoints.size(), sumOp<label>()) == 0 )
        return;

    const meshSurfaceEngine& mse = surfaceEngine();
    pointFieldPMG& points = mesh_.points();
    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointPoints = mse.pointPoints();
    const VRWGraph& bpAtProcs = mse.bpAtProcs();
    const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
    const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

    DynList<label> neiProcs;
    labelLongList edgePoints;
    Map<label> bpToEdgePoint;
    forAll(procPoints, pointI)
    {
        const label bpI = procPoints[pointI];
        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            neiProcs.appendIfNotIn(neiProc);
        }

        if( pVertices[bpI] & EDGENODE )
        {
            bpToEdgePoint.insert(bpI, edgePoints.size());
            edgePoints.append(bpI);
        }
    }

    if( returnReduce(edgePoints.size(), sumOp<label>()) == 0 )
        return;

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    const VRWGraph& pFaces = mse.pointFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();

    scalarField dist(edgePoints.size(), VGREAT);
    vectorField normal(edgePoints.size(), vector::zero);
    vectorField v(edgePoints.size(), vector::zero);

    label pKey(-1);
    if( patchKey_.size() )
    {
        forAll(treatPatches, patchI)
            if( treatPatches[patchI] )
            {
                pKey = patchKey_[patchI];
                break;
            }
    }

    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];
        const point& p = points[bPoints[bpI]];

        //- find patches for the given point
        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn(pPatches(bpI, patchI));

        //- find local values of normals and v
        if( otherPatches.size() == 1 )
        {
            forAllRow(pFaces, bpI, pfI)
            {
                const face& f = bFaces[pFaces(bpI, pfI)];
                const label patchLabel =
                    boundaryFacePatches[pFaces(bpI, pfI)];

                if( treatPatches[patchLabel] )
                {
                    { vector _n=vector::zero; const point& _p0=points[f[0]]; for(label _pi=1;_pi<f.size()-1;++_pi) _n+=(points[f[_pi]]-_p0)^(points[f[_pi+1]]-_p0); normal[epI] += _n; }
                }
                else
                {
                    { vector _n=vector::zero; const point& _p0=points[f[0]]; for(label _pi=1;_pi<f.size()-1;++_pi) _n+=(points[f[_pi]]-_p0)^(points[f[_pi+1]]-_p0); v[epI] += _n; }
                }
            }
        }
        else if( otherPatches.size() == 2 )
        {
            label otherVertex(-1);
            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ = pointPoints(bpI, ppI);

                bool found(true);
                forAll(otherPatches, opI)
                    if( !pPatches.contains(bpJ, otherPatches[opI]) )
                    {
                        found = false;
                        break;
                    }

                if( found )
                {
                    otherVertex = bpJ;
                    break;
                }
            }

            if( otherVertex == -1 )
                continue;

            //- normal vector is co-linear with that edge
            normal[epI] = p - points[bPoints[otherVertex]];
            dist[epI] = mag(normal[epI]);
        }
        else
        {
            // Multi-patch singularity in parallel edge vertex creation.
            // Force local termination -- zero extrusion at this point.
            normal[epI] = mse.pointNormals()[bpI];
            dist[epI] = 0.0;
        }
    }

    //- prepare normals and v for sending to other procs
    std::map<label, LongList<labelledPoint> > exchangeNormals;
    forAll(neiProcs, procI)
        exchangeNormals.insert
        (
            std::make_pair(neiProcs[procI], LongList<labelledPoint>())
        );

    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            //- store values in the list for sending
            LongList<labelledPoint>& dataToSend = exchangeNormals[neiProc];
            dataToSend.append
            (
                labelledPoint(globalPointLabel[bpI], normal[epI])
            );
            dataToSend.append(labelledPoint(globalPointLabel[bpI], v[epI]));
        }
    }

    //- exchange data with other processors
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeNormals, receivedData);
    exchangeNormals.clear();

    label counter(0);
    while( counter < receivedData.size() )
    {
        const labelledPoint& otherNormal = receivedData[counter++];
        const labelledPoint& otherV = receivedData[counter++];

        if( !globalToLocal.found(otherNormal.pointLabel()) ) continue;
        const label bpI = globalToLocal[otherNormal.pointLabel()];
        normal[bpToEdgePoint[bpI]] += otherNormal.coordinates();
        v[bpToEdgePoint[bpI]] += otherV.coordinates();
    }

    //- calculate normals
    forAll(normal, epI)
    {
        const label bpI = edgePoints[epI];

        //- find patches for the given point
        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn(pPatches(bpI, patchI));

        if( otherPatches.size() == 1 )
        {
            const scalar magV = mag(v[epI]) + VSMALL;
            v[epI] /= magV;
            normal[epI] -= (normal[epI] & v[epI]) * v[epI];
        }

        const scalar magN = mag(normal[epI]) + VSMALL;
        normal[epI] /= magN;
    }

    //- calculate distances
    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];
        const point& p = points[bPoints[bpI]];

        //- find patches for the given point
        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
            if( !treatPatches[pPatches(bpI, patchI)] )
                otherPatches.appendIfNotIn(pPatches(bpI, patchI));

        if( otherPatches.size() == 1 )
        {
            forAllRow(pointPoints, bpI, ppI)
            {
                if( pVertices[pointPoints(bpI, ppI)] )
                    continue;

                const vector vec = points[bPoints[pointPoints(bpI, ppI)]] - p;
                const scalar prod = 0.5 * mag(vec & normal[epI]);

                if( prod < dist[epI] )
                    dist[epI] = prod;
            }
        }

        //- limit distances
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceLabel = pFaces(bpI, pfI);
            if( otherPatches.contains(boundaryFacePatches[faceLabel]) )
            {
                const face& f = bFaces[faceLabel];
                const label pos = f.which(bPoints[bpI]);

                if( pos != -1 )
                {
                    const point& ep1 = points[f.prevLabel(pos)];
                    const point& ep2 = points[f.nextLabel(pos)];

                    const scalar dst =
                        help::distanceOfPointFromTheEdge(ep1, ep2, p);

                    if( dst < dist[epI] )
                        dist[epI] = 0.9 * dst;
                }
                else
                {
                    FatalErrorIn
                    (
                        "void boundaryLayers::createNewEdgeVerticesParallel"
                        "("
                            "const labelLongList& procPoints,"
                            "const List<direction>& pVertices,"
                            "const boolList& treatPatches,"
                            "labelList& newLabelForVertex"
                        ") const"
                    ) << "Face does not contains this vertex!"
                        << abort(FatalError);
                }
            }
        }
    }

    //- exchange distances with other processors
    std::map<label, LongList<labelledScalar> > exchangeDistances;
    forAll(neiProcs, procI)
        exchangeDistances.insert
        (
            std::make_pair(neiProcs[procI], LongList<labelledScalar>())
        );

    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];
        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            LongList<labelledScalar>& ls = exchangeDistances[neiProc];
            ls.append(labelledScalar(globalPointLabel[bpI], dist[epI]));
        }
    }

    //- exchange distances with other processors
    LongList<labelledScalar> receivedDistances;
    help::exchangeMap(exchangeDistances, receivedDistances);
    exchangeDistances.clear();

    forAll(receivedDistances, i)
    {
        if( !globalToLocal.found(receivedDistances[i].scalarLabel()) ) continue;
        const label bpI = globalToLocal[receivedDistances[i].scalarLabel()];
        const label epI = bpToEdgePoint[bpI];
        if( dist[epI] > receivedDistances[i].value() )
            dist[epI] = receivedDistances[i].value();
    }

    //- Finally, create new points
    forAll(edgePoints, epI)
    {
        const label bpI = edgePoints[epI];

        const point& p = points[bPoints[bpI]];
        const point np = p - normal[epI] * dist[epI];
        if( !help::isnan(np) && !help::isinf(np) )
        {
            points[nPoints_] = np;
        }
        else
        {
            points[nPoints_] = p;
        }

        if( pKey == -1 )
        {
            //- extrusion for one patch in a single go
            newLabelForVertex_[bPoints[bpI]] = nPoints_;
        }
        else
        {
            const label pointI = bPoints[bpI];

            if( otherVrts_.find(pointI) == otherVrts_.end() )
            {
                std::map<std::pair<label, label>, label> m;
                otherVrts_.insert(std::make_pair(pointI, m));
            }

            std::pair<label, label> pr(pKey, pKey);
            otherVrts_[pointI].insert(std::make_pair(pr, nPoints_));
        }
        ++nPoints_;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
