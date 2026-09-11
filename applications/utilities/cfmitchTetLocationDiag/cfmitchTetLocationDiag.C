#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyMeshTetDecomposition.H"
#include "HashSet.H"
#include "OFstream.H"

#include <algorithm>
#include <deque>
#include <map>

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "addMeshOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createSpecifiedPolyMesh.H"

    Info<< nl
        << "CFMITCH TET LOCATION DIAG"
        << nl << endl;

    const word wallPatchName("propeller");

    const label wallPatchI =
        mesh.boundaryMesh().findIndex(wallPatchName);

    if (wallPatchI < 0)
    {
        FatalErrorInFunction
            << "Cannot find patch " << wallPatchName
            << exit(FatalError);
    }

    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();

    const vectorField& cellCentres = mesh.cellCentres();
    const vectorField& faceCentres = mesh.faceCentres();
    const pointField& points = mesh.points();

    const scalar tetTol =
        polyMeshTetDecomposition::minTetQuality;

    // -------------------------------------------------------------
    // Reproduce OpenFOAM's unique lowQualityTetFaces population.
    // -------------------------------------------------------------

    labelHashSet badFaceSet(mesh.nFaces()/100 + 1);

    polyMeshTetDecomposition::checkFaceTets
    (
        mesh,
        polyMeshTetDecomposition::minTetQuality,
        false,
        &badFaceSet
    );

    const labelList badFaces(badFaceSet.toc());

    Info<< "CFMITCH TET LOCATION: uniqueFailedFaces="
        << badFaces.size() << nl;

    // -------------------------------------------------------------
    // Cell-graph distance from propeller wall.
    //
    // Wall-adjacent cells are distance 0.
    // The next cell across an internal face is distance 1, etc.
    // -------------------------------------------------------------

    labelList wallDistance(mesh.nCells(), -1);

    std::deque<label> queue;

    const polyPatch& wallPatch = mesh.boundaryMesh()[wallPatchI];

    forAll(wallPatch, patchFaceI)
    {
        const label faceI = wallPatch.start() + patchFaceI;
        const label cellI = owner[faceI];

        if (wallDistance[cellI] == -1)
        {
            wallDistance[cellI] = 0;
            queue.push_back(cellI);
        }
    }

    const labelListList& cellCells = mesh.cellCells();

    while (!queue.empty())
    {
        const label cellI = queue.front();
        queue.pop_front();

        const label nextDistance = wallDistance[cellI] + 1;

        const labelList& nbrs = cellCells[cellI];

        forAll(nbrs, nbrI)
        {
            const label nbrCellI = nbrs[nbrI];

            if (wallDistance[nbrCellI] == -1)
            {
                wallDistance[nbrCellI] = nextDistance;
                queue.push_back(nbrCellI);
            }
        }
    }

    // -------------------------------------------------------------
    // Classification.
    // -------------------------------------------------------------

    std::map<label, label> distanceHist;
    std::map<label, label> ownerFaceCountHist;
    std::map<label, label> neighbourFaceCountHist;

    label internalFailed = 0;
    label boundaryFailed = 0;
    label propellerBoundaryFailed = 0;
    label farfieldBoundaryFailed = 0;

    label bothSixFaceCells = 0;
    label eitherNonSixFaceCell = 0;

    // Exact checkFaceTets() failure-mode classification.
    label ownerTetSignFail = 0;
    label neighbourTetSignFail = 0;
    label sharedBasePointFail = 0;
    label boundaryBasePointFail = 0;
    label coupledBoundaryUnclassified = 0;

    // Shared-base orientation diagnostic.
    //
    // In a regular layered six-face BL:
    //   same owner/neighbour wall distance ~= column side face
    //   distance jump of one             ~= layer-to-layer face
    label sharedSameDistance = 0;
    label sharedCrossDistance = 0;
    label sharedDistanceJump1 = 0;
    label sharedDistanceJumpGt1 = 0;
    label sharedDistanceUnknown = 0;

    label sharedSameDistanceBothSix = 0;
    label sharedCrossDistanceBothSix = 0;

    label sharedTriFaces = 0;
    label sharedQuadFaces = 0;
    label sharedOtherFaces = 0;

    label ownerOnlySignFail = 0;
    label neighbourOnlySignFail = 0;
    label ownerAndNeighbourSignFail = 0;
    label neitherInternalSignFail = 0;

    std::map<label, label> ownerSignDistanceHist;
    std::map<label, label> neighbourSignDistanceHist;
    std::map<label, label> sharedBaseDistanceHist;
    std::map<label, label> boundaryBaseDistanceHist;

    label dist0to14 = 0;
    label dist15 = 0;
    label dist16to20 = 0;
    label dist21to30 = 0;
    label dist31to50 = 0;
    label dist51plus = 0;
    label unreachable = 0;

    OFstream csv("cfmitchTetLocationDiag.csv");

    csv
        << "face,internal,patch,owner,neighbour,"
        << "wallDistance,ownerNFaces,neighbourNFaces"
        << nl;

    forAll(badFaces, badI)
    {
        const label faceI = badFaces[badI];

        const label ownCell = owner[faceI];
        const label ownNFaces = mesh.cells()[ownCell].size();

        label neiCell = -1;
        label neiNFaces = -1;

        label faceDistance = wallDistance[ownCell];

        word patchName("internal");

        const bool internal = mesh.isInternalFace(faceI);

        if (internal)
        {
            ++internalFailed;

            neiCell = neighbour[faceI];
            neiNFaces = mesh.cells()[neiCell].size();

            if
            (
                faceDistance < 0
             || (
                    wallDistance[neiCell] >= 0
                 && wallDistance[neiCell] < faceDistance
                )
            )
            {
                faceDistance = wallDistance[neiCell];
            }

            if (ownNFaces == 6 && neiNFaces == 6)
            {
                ++bothSixFaceCells;
            }
            else
            {
                ++eitherNonSixFaceCell;
            }

            ++neighbourFaceCountHist[neiNFaces];
        }
        else
        {
            ++boundaryFailed;

            const label patchI =
                mesh.boundaryMesh().whichPatch(faceI);

            patchName = mesh.boundaryMesh()[patchI].name();

            if (patchName == "propeller")
            {
                ++propellerBoundaryFailed;
            }
            else if (patchName == "farfield")
            {
                ++farfieldBoundaryFailed;
            }

            if (ownNFaces == 6)
            {
                ++bothSixFaceCells;
            }
            else
            {
                ++eitherNonSixFaceCell;
            }
        }

        // ---------------------------------------------------------
        // Reproduce the individual failure branches from
        // OpenFOAM-13 polyMeshTetDecomposition::checkFaceTets().
        //
        // The original check can increment its event counter more than
        // once for a single face:
        //
        //   owner tet sign/degeneracy
        //   neighbour tet sign/degeneracy
        //   shared/base-point decomposition failure
        //
        // badFaceSet collapses those events to one unique face.  These
        // counters recover the event-level classification.
        // ---------------------------------------------------------

        const face& failedFace = mesh.faces()[faceI];

        bool ownerSignBad(false);
        bool neighbourSignBad(false);
        bool sharedBaseBad(false);
        bool boundaryBaseBad(false);

        forAll(failedFace, fPtI)
        {
            const scalar tetQual =
                tetPointRef
                (
                    points[failedFace[fPtI]],
                    points[failedFace.nextLabel(fPtI)],
                    faceCentres[faceI],
                    cellCentres[ownCell]
                ).quality();

            // checkFaceTets(): owner tet is expected negative.
            if( tetQual > -tetTol )
            {
                ownerSignBad = true;
                break;
            }
        }

        if( ownerSignBad )
        {
            ++ownerTetSignFail;
            ++ownerSignDistanceHist[faceDistance];
        }

        if( internal )
        {
            forAll(failedFace, fPtI)
            {
                const scalar tetQual =
                    tetPointRef
                    (
                        points[failedFace[fPtI]],
                        points[failedFace.nextLabel(fPtI)],
                        faceCentres[faceI],
                        cellCentres[neiCell]
                    ).quality();

                // checkFaceTets(): neighbour tet is expected positive.
                if( tetQual < tetTol )
                {
                    neighbourSignBad = true;
                    break;
                }
            }

            if( neighbourSignBad )
            {
                ++neighbourTetSignFail;
                ++neighbourSignDistanceHist[faceDistance];
            }

            sharedBaseBad =
            (
                polyMeshTetDecomposition::findSharedBasePoint
                (
                    mesh,
                    faceI,
                    tetTol,
                    false
                ) == -1
            );

            if( sharedBaseBad )
            {
                ++sharedBasePointFail;
                ++sharedBaseDistanceHist[faceDistance];

                // Face-size classification.
                if( failedFace.size() == 3 )
                {
                    ++sharedTriFaces;
                }
                else if( failedFace.size() == 4 )
                {
                    ++sharedQuadFaces;
                }
                else
                {
                    ++sharedOtherFaces;
                }

                // Wall-distance classification.
                const label ownWallDistance =
                    wallDistance[ownCell];

                const label neiWallDistance =
                    wallDistance[neiCell];

                const bool sharedBothSix =
                    ownNFaces == 6 && neiNFaces == 6;

                if
                (
                    ownWallDistance < 0
                 || neiWallDistance < 0
                )
                {
                    ++sharedDistanceUnknown;
                }
                else
                {
                    label distanceJump =
                        ownWallDistance - neiWallDistance;

                    if( distanceJump < 0 )
                        distanceJump = -distanceJump;

                    if( distanceJump == 0 )
                    {
                        ++sharedSameDistance;

                        if( sharedBothSix )
                            ++sharedSameDistanceBothSix;
                    }
                    else
                    {
                        ++sharedCrossDistance;

                        if( sharedBothSix )
                            ++sharedCrossDistanceBothSix;

                        if( distanceJump == 1 )
                            ++sharedDistanceJump1;
                        else
                            ++sharedDistanceJumpGt1;
                    }
                }
            }

            if( ownerSignBad && neighbourSignBad )
                ++ownerAndNeighbourSignFail;
            else if( ownerSignBad )
                ++ownerOnlySignFail;
            else if( neighbourSignBad )
                ++neighbourOnlySignFail;
            else
                ++neitherInternalSignFail;
        }
        else
        {
            const label patchI =
                mesh.boundaryMesh().whichPatch(faceI);

            if
            (
                patchI >= 0
             && !mesh.boundaryMesh()[patchI].coupled()
            )
            {
                boundaryBaseBad =
                (
                    polyMeshTetDecomposition::findBasePoint
                    (
                        mesh,
                        faceI,
                        tetTol,
                        false
                    ) == -1
                );

                if( boundaryBaseBad )
                {
                    ++boundaryBasePointFail;
                    ++boundaryBaseDistanceHist[faceDistance];
                }
            }
            else
            {
                // Current propeller/farfield case has no failed coupled
                // boundary faces.  Keep this explicit rather than silently
                // misclassifying a future coupled case.
                ++coupledBoundaryUnclassified;
            }
        }

        ++distanceHist[faceDistance];
        ++ownerFaceCountHist[ownNFaces];

        if (faceDistance < 0)
        {
            ++unreachable;
        }
        else if (faceDistance <= 14)
        {
            ++dist0to14;
        }
        else if (faceDistance == 15)
        {
            ++dist15;
        }
        else if (faceDistance <= 20)
        {
            ++dist16to20;
        }
        else if (faceDistance <= 30)
        {
            ++dist21to30;
        }
        else if (faceDistance <= 50)
        {
            ++dist31to50;
        }
        else
        {
            ++dist51plus;
        }

        csv
            << faceI << ','
            << label(internal) << ','
            << patchName << ','
            << ownCell << ','
            << neiCell << ','
            << faceDistance << ','
            << ownNFaces << ','
            << neiNFaces
            << nl;
    }

    Info<< nl
        << "CFMITCH TET LOCATION SUMMARY:"
        << " uniqueFailedFaces=" << badFaces.size()
        << " internal=" << internalFailed
        << " boundary=" << boundaryFailed
        << " propellerBoundary=" << propellerBoundaryFailed
        << " farfieldBoundary=" << farfieldBoundaryFailed
        << nl;

    Info<<
        "CFMITCH TET LOCATION CELL TYPES:"
        << " bothSixFace=" << bothSixFaceCells
        << " eitherNonSixFace=" << eitherNonSixFaceCell
        << nl;

    Info<<
        "CFMITCH TET LOCATION DISTANCE BANDS:"
        << " d0to14=" << dist0to14
        << " d15=" << dist15
        << " d16to20=" << dist16to20
        << " d21to30=" << dist21to30
        << " d31to50=" << dist31to50
        << " d51plus=" << dist51plus
        << " unreachable=" << unreachable
        << nl;

    // ---------------------------------------------------------
    // CFMitch cross-layer quad population diagnostic.
    //
    // Compare FAILING and PASSING ordinary six-face, quad,
    // wall-distance-jump-1 internal faces in the near-wall band.
    //
    // This uses OpenFOAM's exact:
    //   * face-flatness definition
    //   * polyMeshTetDecomposition::minQuality()
    //   * findSharedBasePoint()
    //
    // No mesh modification.
    // ---------------------------------------------------------

    const vectorField& faceAreas = mesh.faceAreas();

    label crossQuadTotal(0);
    label crossQuadFail(0);
    label crossQuadPass(0);

    scalar failFlatSum(0.0);
    scalar passFlatSum(0.0);
    scalar failFlatMin(GREAT);
    scalar passFlatMin(GREAT);
    scalar failFlatMax(-GREAT);
    scalar passFlatMax(-GREAT);

    scalar failBestSharedSum(0.0);
    scalar passBestSharedSum(0.0);
    scalar failBestSharedMin(GREAT);
    scalar passBestSharedMin(GREAT);
    scalar failBestSharedMax(-GREAT);
    scalar passBestSharedMax(-GREAT);

    label failDiagonalConflict(0);
    label passDiagonalConflict(0);

    // ---------------------------------------------------------
    // CFMitch physical warp / centre-clearance diagnostic.
    //
    // OpenFOAM's normalized faceFlatness can hide physically
    // important micron-scale saddling on extremely thin BL cells.
    //
    // Measure:
    //   * absolute quad non-planarity in metres;
    //   * diagonal-specific owner/neighbour plane clearances;
    //   * face-normal cell-centre offsets;
    //   * warp / centre-offset ratio.
    //
    // This is strictly read-only.
    // ---------------------------------------------------------

    label physicalValidFaces(0);
    label physicalDegenerateFaces(0);

    label physicalFailFaces(0);
    label physicalPassFaces(0);

    scalar failWarpSum(0.0);
    scalar passWarpSum(0.0);
    scalar failWarpMin(GREAT);
    scalar passWarpMin(GREAT);
    scalar failWarpMax(-GREAT);
    scalar passWarpMax(-GREAT);

    scalar failWarpRatioSum(0.0);
    scalar passWarpRatioSum(0.0);
    scalar failWarpRatioMin(GREAT);
    scalar passWarpRatioMin(GREAT);
    scalar failWarpRatioMax(-GREAT);
    scalar passWarpRatioMax(-GREAT);

    scalar failPhysicalClearanceSum(0.0);
    scalar passPhysicalClearanceSum(0.0);
    scalar failPhysicalClearanceMin(GREAT);
    scalar passPhysicalClearanceMin(GREAT);
    scalar failPhysicalClearanceMax(-GREAT);
    scalar passPhysicalClearanceMax(-GREAT);

    label physicalSharedFail(0);
    label physicalSharedPass(0);

    label sharedFailPhysicalPass(0);
    label sharedPassPhysicalFail(0);

    // ---------------------------------------------------------
    // CFMitch centre-fan / planarization feasibility.
    //
    // Read-only virtual tests:
    //
    // 1. Replace one quad by four triangles using the existing
    //    OpenFOAM face centre.  No point movement.
    //
    // 2. Measure the maximum vertex displacement required to
    //    project the quad onto its current face-area plane.
    // ---------------------------------------------------------

    label centreFanValid(0);
    label centreFanDegenerate(0);

    label centreFanRecoverFail(0);
    label centreFanStillFail(0);

    label centreFanKeepPass(0);
    label centreFanBreakPass(0);

    scalar failCentreFanClearSum(0.0);
    scalar failCentreFanClearMin(GREAT);
    scalar failCentreFanClearMax(-GREAT);

    scalar failPlanarMoveSum(0.0);
    scalar failPlanarMoveMin(GREAT);
    scalar failPlanarMoveMax(-GREAT);

    scalar failPlanarMoveRatioSum(0.0);
    scalar failPlanarMoveRatioMin(GREAT);
    scalar failPlanarMoveRatioMax(-GREAT);

    label failFacePlaneSeparates(0);
    label failFacePlaneDoesNotSeparate(0);

    // Histogram of warp / min(ownerOffset, neighbourOffset)
    //
    // 0: <0.05
    // 1: 0.05-0.10
    // 2: 0.10-0.20
    // 3: 0.20-0.50
    // 4: 0.50-1.00
    // 5: 1.00-2.00
    // 6: 2.00-5.00
    // 7: >=5.00
    labelList warpRatioTotalBins(8, 0);
    labelList warpRatioFailBins(8, 0);

    OFstream warpCsv
    (
        "cfmitchCrossLayerWarpClearance.csv"
    );

    warpCsv
        << "face,failed,nearDistance,flatness,bestShared,"
        << "bestOwnerBase,bestNeighbourBase,diagonalConflict,"
        << "physicalValid,"
        << "warpDiag02,warpDiag13,warpMax,"
        << "ownerFaceOffset,neighbourFaceOffset,"
        << "centreNormalSeparation,"
        << "warpOverOwnerOffset,"
        << "warpOverNeighbourOffset,"
        << "warpOverMinOffset,"
        << "diag02OwnerClearance,"
        << "diag02NeighbourClearance,"
        << "diag02CommonClearance,"
        << "diag13OwnerClearance,"
        << "diag13NeighbourClearance,"
        << "diag13CommonClearance,"
        << "bestPhysicalCommonClearance,"
        << "physicalSharedFail"
        << nl;

    labelList crossQuadTotalByDistance(14, 0);
    labelList crossQuadFailByDistance(14, 0);

    for
    (
        label faceI = 0;
        faceI < mesh.nInternalFaces();
        ++faceI
    )
    {
        const label ownCell = owner[faceI];
        const label neiCell = neighbour[faceI];

        if
        (
            ownCell < 0
         || neiCell < 0
         || ownCell >= label(mesh.cells().size())
         || neiCell >= label(mesh.cells().size())
        )
            continue;

        // Restrict to the ordinary six-face population identified by
        // the previous diagnostic.
        if
        (
            mesh.cells()[ownCell].size() != 6
         || mesh.cells()[neiCell].size() != 6
        )
            continue;

        const face& f = mesh.faces()[faceI];

        if( f.size() != 4 )
            continue;

        const label ownD = wallDistance[ownCell];
        const label neiD = wallDistance[neiCell];

        if( ownD < 0 || neiD < 0 )
            continue;

        label jump = ownD - neiD;
        if( jump < 0 )
            jump = -jump;

        if( jump != 1 )
            continue;

        const label nearD = Foam::min(ownD, neiD);

        // Our failure population is confined to the BL band.
        // Include interfaces 0-1 through 13-14.
        if( nearD < 0 || nearD > 13 )
            continue;

        ++crossQuadTotal;
        ++crossQuadTotalByDistance[nearD];

        // OpenFOAM faceFlatness() definition.
        const point& fc = faceCentres[faceI];

        scalar sumA(0.0);

        forAll(f, fp)
        {
            const point& thisPoint = points[f[fp]];
            const point& nextPoint =
                points[f.nextLabel(fp)];

            const vector triArea =
                scalar(0.5)
               *((nextPoint - thisPoint)
                ^ (fc - thisPoint));

            sumA += mag(triArea);
        }

        const scalar flatness =
            mag(faceAreas[faceI])
           /(sumA + rootVSmall);

        // Reproduce the base-point search, but retain the best values.
        scalar bestShared(-GREAT);
        scalar bestOwner(-GREAT);
        scalar bestNeighbour(-GREAT);

        label bestOwnerBase(-1);
        label bestNeighbourBase(-1);

        forAll(f, baseI)
        {
            const scalar ownerQ =
                polyMeshTetDecomposition::minQuality
                (
                    mesh,
                    cellCentres[ownCell],
                    faceI,
                    true,
                    baseI
                );

            const scalar neighbourQ =
                polyMeshTetDecomposition::minQuality
                (
                    mesh,
                    cellCentres[neiCell],
                    faceI,
                    false,
                    baseI
                );

            const scalar sharedQ =
                Foam::min(ownerQ, neighbourQ);

            if( sharedQ > bestShared )
                bestShared = sharedQ;

            if( ownerQ > bestOwner )
            {
                bestOwner = ownerQ;
                bestOwnerBase = baseI;
            }

            if( neighbourQ > bestNeighbour )
            {
                bestNeighbour = neighbourQ;
                bestNeighbourBase = baseI;
            }
        }

        const bool sharedFail =
        (
            polyMeshTetDecomposition::findSharedBasePoint
            (
                mesh,
                faceI,
                tetTol,
                false
            ) == -1
        );

        // For a quad, base vertices 0/2 imply one diagonal and
        // base vertices 1/3 imply the other.
        const bool diagonalConflict =
        (
            bestOwnerBase >= 0
         && bestNeighbourBase >= 0
         && (bestOwnerBase % 2) != (bestNeighbourBase % 2)
        );

        // -----------------------------------------------------
        // Physical quad warp / centre-clearance measurement.
        //
        // Quad ordering:
        //
        //     p0 ----- p1
        //      |       |
        //      |       |
        //     p3 ----- p2
        //
        // Diagonal 0-2:
        //     triangles (0,1,2) and (0,2,3)
        //
        // Diagonal 1-3:
        //     triangles (0,1,3) and (1,2,3)
        //
        // Triangle normals use the face's vertex orientation.
        // For a correctly oriented internal face:
        //
        //     owner centre     -> negative side
        //     neighbour centre -> positive side
        //
        // Thus positive clearance below means the cell centre
        // lies on the expected side of that triangle plane.
        // -----------------------------------------------------

        bool physicalValid(false);

        scalar warpDiag02(-1.0);
        scalar warpDiag13(-1.0);
        scalar warpMax(-1.0);

        scalar ownerFaceOffset(-1.0);
        scalar neighbourFaceOffset(-1.0);
        scalar centreNormalSeparation(-1.0);

        scalar warpOverOwnerOffset(-1.0);
        scalar warpOverNeighbourOffset(-1.0);
        scalar warpOverMinOffset(-1.0);

        scalar diag02OwnerClearance(-GREAT);
        scalar diag02NeighbourClearance(-GREAT);
        scalar diag02CommonClearance(-GREAT);

        scalar diag13OwnerClearance(-GREAT);
        scalar diag13NeighbourClearance(-GREAT);
        scalar diag13CommonClearance(-GREAT);

        scalar bestPhysicalCommonClearance(-GREAT);

        bool physicalFail(false);

        const point& p0 = points[f[0]];
        const point& p1 = points[f[1]];
        const point& p2 = points[f[2]];
        const point& p3 = points[f[3]];

        const vector n012 = (p1 - p0) ^ (p2 - p0);
        const vector n023 = (p2 - p0) ^ (p3 - p0);

        const vector n013 = (p1 - p0) ^ (p3 - p0);
        const vector n123 = (p2 - p1) ^ (p3 - p1);

        const scalar m012 = mag(n012);
        const scalar m023 = mag(n023);
        const scalar m013 = mag(n013);
        const scalar m123 = mag(n123);

        const scalar faceAreaMag =
            mag(faceAreas[faceI]);

        if
        (
            m012 > rootVSmall
         && m023 > rootVSmall
         && m013 > rootVSmall
         && m123 > rootVSmall
         && faceAreaMag > rootVSmall
        )
        {
            physicalValid = true;
            ++physicalValidFaces;

            // -------------------------------------------------
            // Absolute quad non-planarity.
            //
            // Each diagonal gets a symmetric two-sided measure:
            // maximum distance of either opposite vertex from
            // the plane of the complementary triangle.
            // Units: metres.
            // -------------------------------------------------

            const scalar p3From012 =
                ((p3 - p0) & n012) / m012;

            const scalar p1From023 =
                ((p1 - p0) & n023) / m023;

            const scalar p2From013 =
                ((p2 - p0) & n013) / m013;

            const scalar p0From123 =
                ((p0 - p1) & n123) / m123;

            warpDiag02 =
                Foam::max
                (
                    mag(p3From012),
                    mag(p1From023)
                );

            warpDiag13 =
                Foam::max
                (
                    mag(p2From013),
                    mag(p0From123)
                );

            warpMax =
                Foam::max(warpDiag02, warpDiag13);

            // -------------------------------------------------
            // Face-normal centre offsets.
            //
            // Use the OpenFOAM face-area direction merely as the
            // reference normal.  Magnitudes are reported here,
            // while the diagonal-specific quantities below retain
            // orientation/sign information.
            // -------------------------------------------------

            const vector faceNormal =
                faceAreas[faceI] / faceAreaMag;

            ownerFaceOffset =
                mag
                (
                    (cellCentres[ownCell] - fc)
                  & faceNormal
                );

            neighbourFaceOffset =
                mag
                (
                    (cellCentres[neiCell] - fc)
                  & faceNormal
                );

            centreNormalSeparation =
                mag
                (
                    (cellCentres[neiCell]
                   - cellCentres[ownCell])
                  & faceNormal
                );

            warpOverOwnerOffset =
                warpMax
               /Foam::max(ownerFaceOffset, rootVSmall);

            warpOverNeighbourOffset =
                warpMax
               /Foam::max(neighbourFaceOffset, rootVSmall);

            const scalar minCentreOffset =
                Foam::min
                (
                    ownerFaceOffset,
                    neighbourFaceOffset
                );

            warpOverMinOffset =
                warpMax
               /Foam::max(minCentreOffset, rootVSmall);

            // -------------------------------------------------
            // Signed physical clearance for diagonal 0-2.
            //
            // Face-oriented triangle normals point from owner
            // toward neighbour for a normal internal face.
            //
            // owner clearance = -signedDistance
            // neighbour       = +signedDistance
            //
            // A positive value means correct side of the plane.
            // -------------------------------------------------

            const scalar owner012 =
                (
                    (cellCentres[ownCell] - p0)
                  & n012
                ) / m012;

            const scalar owner023 =
                (
                    (cellCentres[ownCell] - p0)
                  & n023
                ) / m023;

            const scalar neighbour012 =
                (
                    (cellCentres[neiCell] - p0)
                  & n012
                ) / m012;

            const scalar neighbour023 =
                (
                    (cellCentres[neiCell] - p0)
                  & n023
                ) / m023;

            diag02OwnerClearance =
                Foam::min(-owner012, -owner023);

            diag02NeighbourClearance =
                Foam::min(neighbour012, neighbour023);

            diag02CommonClearance =
                Foam::min
                (
                    diag02OwnerClearance,
                    diag02NeighbourClearance
                );

            // -------------------------------------------------
            // Signed physical clearance for diagonal 1-3.
            // -------------------------------------------------

            const scalar owner013 =
                (
                    (cellCentres[ownCell] - p0)
                  & n013
                ) / m013;

            const scalar owner123 =
                (
                    (cellCentres[ownCell] - p1)
                  & n123
                ) / m123;

            const scalar neighbour013 =
                (
                    (cellCentres[neiCell] - p0)
                  & n013
                ) / m013;

            const scalar neighbour123 =
                (
                    (cellCentres[neiCell] - p1)
                  & n123
                ) / m123;

            diag13OwnerClearance =
                Foam::min(-owner013, -owner123);

            diag13NeighbourClearance =
                Foam::min(neighbour013, neighbour123);

            diag13CommonClearance =
                Foam::min
                (
                    diag13OwnerClearance,
                    diag13NeighbourClearance
                );

            bestPhysicalCommonClearance =
                Foam::max
                (
                    diag02CommonClearance,
                    diag13CommonClearance
                );

            // Sign-only physical analogue of findSharedBasePoint.
            // OpenFOAM additionally normalizes tet quality, but at
            // its near-zero tolerance the sign boundary should be
            // strongly correlated with this quantity.
            physicalFail =
                !(bestPhysicalCommonClearance > scalar(0.0));

            // -------------------------------------------------
            // Virtual centre-fan topology test.
            //
            // Replace the quad interface conceptually by:
            //
            //   (p0,p1,fc)
            //   (p1,p2,fc)
            //   (p2,p3,fc)
            //   (p3,p0,fc)
            //
            // Each triangle uses the original face orientation.
            // A positive minimum common clearance means every
            // triangle separates owner and neighbour correctly.
            //
            // Geometry and mesh remain completely unchanged.
            // -------------------------------------------------

            bool centreFanIsValid(true);
            scalar centreFanCommonClearance(GREAT);

            forAll(f, fp)
            {
                const point& fanP0 =
                    points[f[fp]];

                const point& fanP1 =
                    points[f.nextLabel(fp)];

                const vector fanNormal =
                    (fanP1 - fanP0)
                  ^ (fc - fanP0);

                const scalar fanMag =
                    mag(fanNormal);

                if( fanMag <= rootVSmall )
                {
                    centreFanIsValid = false;
                    break;
                }

                const scalar fanOwnerSigned =
                    (
                        (cellCentres[ownCell] - fanP0)
                      & fanNormal
                    ) / fanMag;

                const scalar fanNeighbourSigned =
                    (
                        (cellCentres[neiCell] - fanP0)
                      & fanNormal
                    ) / fanMag;

                const scalar fanCommon =
                    Foam::min
                    (
                        -fanOwnerSigned,
                        fanNeighbourSigned
                    );

                centreFanCommonClearance =
                    Foam::min
                    (
                        centreFanCommonClearance,
                        fanCommon
                    );
            }

            if( centreFanIsValid )
            {
                ++centreFanValid;

                const bool centreFanPass =
                    centreFanCommonClearance > scalar(0.0);

                if( sharedFail )
                {
                    if( centreFanPass )
                        ++centreFanRecoverFail;
                    else
                        ++centreFanStillFail;

                    failCentreFanClearSum +=
                        centreFanCommonClearance;

                    failCentreFanClearMin =
                        Foam::min
                        (
                            failCentreFanClearMin,
                            centreFanCommonClearance
                        );

                    failCentreFanClearMax =
                        Foam::max
                        (
                            failCentreFanClearMax,
                            centreFanCommonClearance
                        );
                }
                else
                {
                    if( centreFanPass )
                        ++centreFanKeepPass;
                    else
                        ++centreFanBreakPass;
                }
            }
            else
            {
                ++centreFanDegenerate;
            }

            // -------------------------------------------------
            // Virtual planarization displacement.
            //
            // Projecting every quad vertex to the plane through
            // fc normal to the current OpenFOAM face-area vector
            // would make the face exactly planar.
            //
            // We do NOT perform the move.  We merely measure
            // the maximum required normal displacement.
            //
            // Also determine whether that plane itself separates
            // owner and neighbour correctly.
            // -------------------------------------------------

            scalar planarMaxMove(0.0);

            forAll(f, fp)
            {
                const scalar vertexPlaneDistance =
                    mag
                    (
                        (points[f[fp]] - fc)
                      & faceNormal
                    );

                planarMaxMove =
                    Foam::max
                    (
                        planarMaxMove,
                        vertexPlaneDistance
                    );
            }

            const scalar ownerPlaneSigned =
                (
                    cellCentres[ownCell] - fc
                ) & faceNormal;

            const scalar neighbourPlaneSigned =
                (
                    cellCentres[neiCell] - fc
                ) & faceNormal;

            const bool facePlaneSeparates =
            (
                ownerPlaneSigned < scalar(0.0)
             && neighbourPlaneSigned > scalar(0.0)
            );

            const scalar planarMoveRatio =
                planarMaxMove
               /Foam::max
                (
                    minCentreOffset,
                    rootVSmall
                );

            if( sharedFail )
            {
                failPlanarMoveSum += planarMaxMove;

                failPlanarMoveMin =
                    Foam::min
                    (
                        failPlanarMoveMin,
                        planarMaxMove
                    );

                failPlanarMoveMax =
                    Foam::max
                    (
                        failPlanarMoveMax,
                        planarMaxMove
                    );

                failPlanarMoveRatioSum +=
                    planarMoveRatio;

                failPlanarMoveRatioMin =
                    Foam::min
                    (
                        failPlanarMoveRatioMin,
                        planarMoveRatio
                    );

                failPlanarMoveRatioMax =
                    Foam::max
                    (
                        failPlanarMoveRatioMax,
                        planarMoveRatio
                    );

                if( facePlaneSeparates )
                    ++failFacePlaneSeparates;
                else
                    ++failFacePlaneDoesNotSeparate;
            }

            if( physicalFail )
                ++physicalSharedFail;
            else
                ++physicalSharedPass;

            if( sharedFail && !physicalFail )
                ++sharedFailPhysicalPass;

            if( !sharedFail && physicalFail )
                ++sharedPassPhysicalFail;

            // -------------------------------------------------
            // Fail/pass population statistics.
            // -------------------------------------------------

            if( sharedFail )
            {
                ++physicalFailFaces;

                failWarpSum += warpMax;
                failWarpMin =
                    Foam::min(failWarpMin, warpMax);
                failWarpMax =
                    Foam::max(failWarpMax, warpMax);

                failWarpRatioSum += warpOverMinOffset;
                failWarpRatioMin =
                    Foam::min
                    (
                        failWarpRatioMin,
                        warpOverMinOffset
                    );
                failWarpRatioMax =
                    Foam::max
                    (
                        failWarpRatioMax,
                        warpOverMinOffset
                    );

                failPhysicalClearanceSum +=
                    bestPhysicalCommonClearance;

                failPhysicalClearanceMin =
                    Foam::min
                    (
                        failPhysicalClearanceMin,
                        bestPhysicalCommonClearance
                    );

                failPhysicalClearanceMax =
                    Foam::max
                    (
                        failPhysicalClearanceMax,
                        bestPhysicalCommonClearance
                    );
            }
            else
            {
                ++physicalPassFaces;

                passWarpSum += warpMax;
                passWarpMin =
                    Foam::min(passWarpMin, warpMax);
                passWarpMax =
                    Foam::max(passWarpMax, warpMax);

                passWarpRatioSum += warpOverMinOffset;
                passWarpRatioMin =
                    Foam::min
                    (
                        passWarpRatioMin,
                        warpOverMinOffset
                    );
                passWarpRatioMax =
                    Foam::max
                    (
                        passWarpRatioMax,
                        warpOverMinOffset
                    );

                passPhysicalClearanceSum +=
                    bestPhysicalCommonClearance;

                passPhysicalClearanceMin =
                    Foam::min
                    (
                        passPhysicalClearanceMin,
                        bestPhysicalCommonClearance
                    );

                passPhysicalClearanceMax =
                    Foam::max
                    (
                        passPhysicalClearanceMax,
                        bestPhysicalCommonClearance
                    );
            }

            // -------------------------------------------------
            // warp / min-centre-offset histogram.
            // -------------------------------------------------

            label ratioBin = 7;

            if( warpOverMinOffset < scalar(0.05) )
                ratioBin = 0;
            else if( warpOverMinOffset < scalar(0.10) )
                ratioBin = 1;
            else if( warpOverMinOffset < scalar(0.20) )
                ratioBin = 2;
            else if( warpOverMinOffset < scalar(0.50) )
                ratioBin = 3;
            else if( warpOverMinOffset < scalar(1.00) )
                ratioBin = 4;
            else if( warpOverMinOffset < scalar(2.00) )
                ratioBin = 5;
            else if( warpOverMinOffset < scalar(5.00) )
                ratioBin = 6;

            ++warpRatioTotalBins[ratioBin];

            if( sharedFail )
                ++warpRatioFailBins[ratioBin];
        }
        else
        {
            ++physicalDegenerateFaces;
        }

        // One row per ordinary cross-layer quad.
        warpCsv
            << faceI << ','
            << label(sharedFail) << ','
            << nearD << ','
            << flatness << ','
            << bestShared << ','
            << bestOwnerBase << ','
            << bestNeighbourBase << ','
            << label(diagonalConflict) << ','
            << label(physicalValid) << ','
            << warpDiag02 << ','
            << warpDiag13 << ','
            << warpMax << ','
            << ownerFaceOffset << ','
            << neighbourFaceOffset << ','
            << centreNormalSeparation << ','
            << warpOverOwnerOffset << ','
            << warpOverNeighbourOffset << ','
            << warpOverMinOffset << ','
            << diag02OwnerClearance << ','
            << diag02NeighbourClearance << ','
            << diag02CommonClearance << ','
            << diag13OwnerClearance << ','
            << diag13NeighbourClearance << ','
            << diag13CommonClearance << ','
            << bestPhysicalCommonClearance << ','
            << label(physicalFail)
            << nl;

        if( sharedFail )
        {
            ++crossQuadFail;
            ++crossQuadFailByDistance[nearD];

            failFlatSum += flatness;
            failFlatMin = Foam::min(failFlatMin, flatness);
            failFlatMax = Foam::max(failFlatMax, flatness);

            failBestSharedSum += bestShared;
            failBestSharedMin =
                Foam::min(failBestSharedMin, bestShared);
            failBestSharedMax =
                Foam::max(failBestSharedMax, bestShared);

            if( diagonalConflict )
                ++failDiagonalConflict;
        }
        else
        {
            ++crossQuadPass;

            passFlatSum += flatness;
            passFlatMin = Foam::min(passFlatMin, flatness);
            passFlatMax = Foam::max(passFlatMax, flatness);

            passBestSharedSum += bestShared;
            passBestSharedMin =
                Foam::min(passBestSharedMin, bestShared);
            passBestSharedMax =
                Foam::max(passBestSharedMax, bestShared);

            if( diagonalConflict )
                ++passDiagonalConflict;
        }
    }

    Info<< nl
        << "CFMITCH CROSS-LAYER QUAD POPULATION:"
        << " total=" << crossQuadTotal
        << " fail=" << crossQuadFail
        << " pass=" << crossQuadPass;

    if( crossQuadTotal > 0 )
    {
        Info<< " failPct="
            << scalar(100.0)
              *scalar(crossQuadFail)
              /scalar(crossQuadTotal);
    }

    Info<< nl;

    Info<<
        "CFMITCH CROSS-LAYER QUAD FLATNESS:"
        << " failMin=" << failFlatMin
        << " failAvg="
        << (
               crossQuadFail > 0
             ? failFlatSum/scalar(crossQuadFail)
             : scalar(0)
           )
        << " failMax=" << failFlatMax
        << " passMin=" << passFlatMin
        << " passAvg="
        << (
               crossQuadPass > 0
             ? passFlatSum/scalar(crossQuadPass)
             : scalar(0)
           )
        << " passMax=" << passFlatMax
        << nl;

    Info<<
        "CFMITCH CROSS-LAYER QUAD BEST_SHARED:"
        << " failMin=" << failBestSharedMin
        << " failAvg="
        << (
               crossQuadFail > 0
             ? failBestSharedSum/scalar(crossQuadFail)
             : scalar(0)
           )
        << " failMax=" << failBestSharedMax
        << " passMin=" << passBestSharedMin
        << " passAvg="
        << (
               crossQuadPass > 0
             ? passBestSharedSum/scalar(crossQuadPass)
             : scalar(0)
           )
        << " passMax=" << passBestSharedMax
        << nl;

    Info<<
        "CFMITCH CROSS-LAYER QUAD DIAGONAL:"
        << " failConflict=" << failDiagonalConflict
        << " failConflictPct="
        << (
               crossQuadFail > 0
             ? scalar(100.0)
               *scalar(failDiagonalConflict)
               /scalar(crossQuadFail)
             : scalar(0)
           )
        << " passConflict=" << passDiagonalConflict
        << " passConflictPct="
        << (
               crossQuadPass > 0
             ? scalar(100.0)
               *scalar(passDiagonalConflict)
               /scalar(crossQuadPass)
             : scalar(0)
           )
        << nl;

    Info<<
        "CFMITCH CROSS-LAYER PHYSICAL WARP:"
        << " valid=" << physicalValidFaces
        << " degenerate=" << physicalDegenerateFaces
        << " failCount=" << physicalFailFaces
        << " failMinM="
        << (
               physicalFailFaces > 0
             ? failWarpMin
             : scalar(0)
           )
        << " failAvgM="
        << (
               physicalFailFaces > 0
             ? failWarpSum/scalar(physicalFailFaces)
             : scalar(0)
           )
        << " failMaxM="
        << (
               physicalFailFaces > 0
             ? failWarpMax
             : scalar(0)
           )
        << " passCount=" << physicalPassFaces
        << " passMinM="
        << (
               physicalPassFaces > 0
             ? passWarpMin
             : scalar(0)
           )
        << " passAvgM="
        << (
               physicalPassFaces > 0
             ? passWarpSum/scalar(physicalPassFaces)
             : scalar(0)
           )
        << " passMaxM="
        << (
               physicalPassFaces > 0
             ? passWarpMax
             : scalar(0)
           )
        << nl;

    Info<<
        "CFMITCH CROSS-LAYER WARP/MIN-CENTRE-OFFSET:"
        << " failMin="
        << (
               physicalFailFaces > 0
             ? failWarpRatioMin
             : scalar(0)
           )
        << " failAvg="
        << (
               physicalFailFaces > 0
             ? failWarpRatioSum/scalar(physicalFailFaces)
             : scalar(0)
           )
        << " failMax="
        << (
               physicalFailFaces > 0
             ? failWarpRatioMax
             : scalar(0)
           )
        << " passMin="
        << (
               physicalPassFaces > 0
             ? passWarpRatioMin
             : scalar(0)
           )
        << " passAvg="
        << (
               physicalPassFaces > 0
             ? passWarpRatioSum/scalar(physicalPassFaces)
             : scalar(0)
           )
        << " passMax="
        << (
               physicalPassFaces > 0
             ? passWarpRatioMax
             : scalar(0)
           )
        << nl;

    Info<<
        "CFMITCH CROSS-LAYER PHYSICAL SHARED CLEARANCE:"
        << " failMinM="
        << (
               physicalFailFaces > 0
             ? failPhysicalClearanceMin
             : scalar(0)
           )
        << " failAvgM="
        << (
               physicalFailFaces > 0
             ? failPhysicalClearanceSum
              /scalar(physicalFailFaces)
             : scalar(0)
           )
        << " failMaxM="
        << (
               physicalFailFaces > 0
             ? failPhysicalClearanceMax
             : scalar(0)
           )
        << " passMinM="
        << (
               physicalPassFaces > 0
             ? passPhysicalClearanceMin
             : scalar(0)
           )
        << " passAvgM="
        << (
               physicalPassFaces > 0
             ? passPhysicalClearanceSum
              /scalar(physicalPassFaces)
             : scalar(0)
           )
        << " passMaxM="
        << (
               physicalPassFaces > 0
             ? passPhysicalClearanceMax
             : scalar(0)
           )
        << nl;

    Info<<
        "CFMITCH CROSS-LAYER PHYSICAL SIGN PARITY:"
        << " physicalFail=" << physicalSharedFail
        << " physicalPass=" << physicalSharedPass
        << " sharedFailPhysicalPass="
        << sharedFailPhysicalPass
        << " sharedPassPhysicalFail="
        << sharedPassPhysicalFail
        << nl;

    Info<<
        "CFMITCH CENTRE-FAN FEASIBILITY:"
        << " valid=" << centreFanValid
        << " degenerate=" << centreFanDegenerate
        << " sharedFails=" << crossQuadFail
        << " recovered=" << centreFanRecoverFail
        << " stillFail=" << centreFanStillFail
        << " recoveryPct="
        << (
               crossQuadFail > 0
             ? scalar(100.0)
              *scalar(centreFanRecoverFail)
              /scalar(crossQuadFail)
             : scalar(0)
           )
        << " existingPass=" << crossQuadPass
        << " keepPass=" << centreFanKeepPass
        << " breakPass=" << centreFanBreakPass
        << nl;

    Info<<
        "CFMITCH CENTRE-FAN FAIL CLEARANCE:"
        << " minM="
        << (
               crossQuadFail > 0
             ? failCentreFanClearMin
             : scalar(0)
           )
        << " avgM="
        << (
               crossQuadFail > 0
             ? failCentreFanClearSum
              /scalar(crossQuadFail)
             : scalar(0)
           )
        << " maxM="
        << (
               crossQuadFail > 0
             ? failCentreFanClearMax
             : scalar(0)
           )
        << nl;

    Info<<
        "CFMITCH PLANARIZATION FEASIBILITY:"
        << " sharedFails=" << crossQuadFail
        << " facePlaneSeparates="
        << failFacePlaneSeparates
        << " facePlaneDoesNotSeparate="
        << failFacePlaneDoesNotSeparate
        << " moveMinM="
        << (
               crossQuadFail > 0
             ? failPlanarMoveMin
             : scalar(0)
           )
        << " moveAvgM="
        << (
               crossQuadFail > 0
             ? failPlanarMoveSum
              /scalar(crossQuadFail)
             : scalar(0)
           )
        << " moveMaxM="
        << (
               crossQuadFail > 0
             ? failPlanarMoveMax
             : scalar(0)
           )
        << " moveRatioMin="
        << (
               crossQuadFail > 0
             ? failPlanarMoveRatioMin
             : scalar(0)
           )
        << " moveRatioAvg="
        << (
               crossQuadFail > 0
             ? failPlanarMoveRatioSum
              /scalar(crossQuadFail)
             : scalar(0)
           )
        << " moveRatioMax="
        << (
               crossQuadFail > 0
             ? failPlanarMoveRatioMax
             : scalar(0)
           )
        << nl;

    Info<<
        "CFMITCH WARP/MIN-CENTRE-OFFSET HISTOGRAM:"
        << nl;

    Info<<
        "  bin=<0.05"
        << " total=" << warpRatioTotalBins[0]
        << " fail=" << warpRatioFailBins[0]
        << " failPct="
        << (
               warpRatioTotalBins[0] > 0
             ? scalar(100.0)
              *scalar(warpRatioFailBins[0])
              /scalar(warpRatioTotalBins[0])
             : scalar(0)
           )
        << nl;

    Info<<
        "  bin=0.05-0.10"
        << " total=" << warpRatioTotalBins[1]
        << " fail=" << warpRatioFailBins[1]
        << " failPct="
        << (
               warpRatioTotalBins[1] > 0
             ? scalar(100.0)
              *scalar(warpRatioFailBins[1])
              /scalar(warpRatioTotalBins[1])
             : scalar(0)
           )
        << nl;

    Info<<
        "  bin=0.10-0.20"
        << " total=" << warpRatioTotalBins[2]
        << " fail=" << warpRatioFailBins[2]
        << " failPct="
        << (
               warpRatioTotalBins[2] > 0
             ? scalar(100.0)
              *scalar(warpRatioFailBins[2])
              /scalar(warpRatioTotalBins[2])
             : scalar(0)
           )
        << nl;

    Info<<
        "  bin=0.20-0.50"
        << " total=" << warpRatioTotalBins[3]
        << " fail=" << warpRatioFailBins[3]
        << " failPct="
        << (
               warpRatioTotalBins[3] > 0
             ? scalar(100.0)
              *scalar(warpRatioFailBins[3])
              /scalar(warpRatioTotalBins[3])
             : scalar(0)
           )
        << nl;

    Info<<
        "  bin=0.50-1.00"
        << " total=" << warpRatioTotalBins[4]
        << " fail=" << warpRatioFailBins[4]
        << " failPct="
        << (
               warpRatioTotalBins[4] > 0
             ? scalar(100.0)
              *scalar(warpRatioFailBins[4])
              /scalar(warpRatioTotalBins[4])
             : scalar(0)
           )
        << nl;

    Info<<
        "  bin=1.00-2.00"
        << " total=" << warpRatioTotalBins[5]
        << " fail=" << warpRatioFailBins[5]
        << " failPct="
        << (
               warpRatioTotalBins[5] > 0
             ? scalar(100.0)
              *scalar(warpRatioFailBins[5])
              /scalar(warpRatioTotalBins[5])
             : scalar(0)
           )
        << nl;

    Info<<
        "  bin=2.00-5.00"
        << " total=" << warpRatioTotalBins[6]
        << " fail=" << warpRatioFailBins[6]
        << " failPct="
        << (
               warpRatioTotalBins[6] > 0
             ? scalar(100.0)
              *scalar(warpRatioFailBins[6])
              /scalar(warpRatioTotalBins[6])
             : scalar(0)
           )
        << nl;

    Info<<
        "  bin=>=5.00"
        << " total=" << warpRatioTotalBins[7]
        << " fail=" << warpRatioFailBins[7]
        << " failPct="
        << (
               warpRatioTotalBins[7] > 0
             ? scalar(100.0)
              *scalar(warpRatioFailBins[7])
              /scalar(warpRatioTotalBins[7])
             : scalar(0)
           )
        << nl;

    Info<<
        "CFMITCH CROSS-LAYER WARP CSV:"
        << " file=cfmitchCrossLayerWarpClearance.csv"
        << " good=" << warpCsv.good()
        << nl;

    Info<< "CFMITCH CROSS-LAYER QUAD BY DISTANCE:" << nl;

    forAll(crossQuadTotalByDistance, d)
    {
        const label total = crossQuadTotalByDistance[d];
        const label fail = crossQuadFailByDistance[d];

        if( total == 0 )
            continue;

        Info<< "  distance=" << d
            << " total=" << total
            << " fail=" << fail
            << " pass=" << (total - fail)
            << " failPct="
            << scalar(100.0)*scalar(fail)/scalar(total)
            << nl;
    }

    const label classifiedFailureEvents =
        ownerTetSignFail
      + neighbourTetSignFail
      + sharedBasePointFail
      + boundaryBasePointFail;

    Info<< nl
        << "CFMITCH TET FAILURE MODES:"
        << " ownerSign=" << ownerTetSignFail
        << " neighbourSign=" << neighbourTetSignFail
        << " sharedBasePoint=" << sharedBasePointFail
        << " boundaryBasePoint=" << boundaryBasePointFail
        << " coupledBoundaryUnclassified="
        << coupledBoundaryUnclassified
        << " classifiedEvents=" << classifiedFailureEvents
        << nl;

    Info<<
        "CFMITCH TET SHARED BASE ORIENTATION:"
        << " sameDistance=" << sharedSameDistance
        << " crossDistance=" << sharedCrossDistance
        << " jump1=" << sharedDistanceJump1
        << " jumpGt1=" << sharedDistanceJumpGt1
        << " unknown=" << sharedDistanceUnknown
        << " sameBothSix=" << sharedSameDistanceBothSix
        << " crossBothSix=" << sharedCrossDistanceBothSix
        << nl;

    Info<<
        "CFMITCH TET SHARED BASE FACE SIZE:"
        << " tri=" << sharedTriFaces
        << " quad=" << sharedQuadFaces
        << " other=" << sharedOtherFaces
        << nl;

    Info<<
        "CFMITCH TET INTERNAL SIGN OVERLAP:"
        << " ownerOnly=" << ownerOnlySignFail
        << " neighbourOnly=" << neighbourOnlySignFail
        << " both=" << ownerAndNeighbourSignFail
        << " neither=" << neitherInternalSignFail
        << nl;

    auto printFailureDistanceHistogram =
    [&](const char* name, const std::map<label,label>& hist)
    {
        Info<< "CFMITCH TET FAILURE DISTANCE " << name << ":" << nl;

        for
        (
            std::map<label,label>::const_iterator iter = hist.begin();
            iter != hist.end();
            ++iter
        )
        {
            Info<< "  distance=" << iter->first
                << " faces=" << iter->second << nl;
        }
    };

    printFailureDistanceHistogram
    (
        "OWNER_SIGN",
        ownerSignDistanceHist
    );

    printFailureDistanceHistogram
    (
        "NEIGHBOUR_SIGN",
        neighbourSignDistanceHist
    );

    printFailureDistanceHistogram
    (
        "SHARED_BASE",
        sharedBaseDistanceHist
    );

    printFailureDistanceHistogram
    (
        "BOUNDARY_BASE",
        boundaryBaseDistanceHist
    );

    Info<< nl
        << "CFMITCH TET LOCATION EXACT DISTANCE HISTOGRAM:"
        << nl;

    for
    (
        std::map<label,label>::const_iterator iter = distanceHist.begin();
        iter != distanceHist.end();
        ++iter
    )
    {
        Info<< "  distance=" << iter->first
            << " faces=" << iter->second << nl;
    }

    Info<< nl
        << "CFMITCH TET LOCATION OWNER FACE-COUNT HISTOGRAM:"
        << nl;

    for
    (
        std::map<label,label>::const_iterator iter =
            ownerFaceCountHist.begin();
        iter != ownerFaceCountHist.end();
        ++iter
    )
    {
        Info<< "  nFaces=" << iter->first
            << " count=" << iter->second << nl;
    }

    Info<< nl
        << "CFMITCH TET LOCATION NEIGHBOUR FACE-COUNT HISTOGRAM:"
        << nl;

    for
    (
        std::map<label,label>::const_iterator iter =
            neighbourFaceCountHist.begin();
        iter != neighbourFaceCountHist.end();
        ++iter
    )
    {
        Info<< "  nFaces=" << iter->first
            << " count=" << iter->second << nl;
    }

    // ============================================================
    // CFMitch CAUSAL AUDIT V1
    //
    // Read-only investigation of the three architectural hypotheses:
    //
    //  A. Tangential resolution / anisotropy is incompatible with h1.
    //  B. Significant layer-face warp is inherited from the wall footprint.
    //  C. The V5.7 same-owner/same-neighbour fan exposes folded interface
    //     planes and therefore creates OpenFOAM face-plane concavity.
    //
    // No mesh modification.
    //
    // The determinant calculation below intentionally reproduces
    // OpenFOAM Foundation 13 meshCheck::cellDeterminant().
    // ============================================================

    Info<< nl
        << "CFMITCH CAUSAL AUDIT V1:"
        << " begin"
        << nl;


    // ------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------

    auto auditFaceCharacteristicLength =
    [&](const face& f) -> scalar
    {
        if( f.size() < 2 )
            return scalar(0);

        scalar sumLength = scalar(0);

        forAll(f, fp)
        {
            const point& a = points[f[fp]];
            const point& b = points[f.nextLabel(fp)];

            sumLength += mag(b - a);
        }

        return sumLength/scalar(f.size());
    };


    auto auditQuadWarp =
    [&](const face& f) -> scalar
    {
        if( f.size() != 4 )
            return scalar(-1);

        const point& p0 = points[f[0]];
        const point& p1 = points[f[1]];
        const point& p2 = points[f[2]];
        const point& p3 = points[f[3]];

        const vector n012 = (p1 - p0) ^ (p2 - p0);
        const vector n023 = (p2 - p0) ^ (p3 - p0);

        const vector n013 = (p1 - p0) ^ (p3 - p0);
        const vector n123 = (p2 - p1) ^ (p3 - p1);

        const scalar m012 = mag(n012);
        const scalar m023 = mag(n023);
        const scalar m013 = mag(n013);
        const scalar m123 = mag(n123);

        if
        (
            m012 <= rootVSmall
         || m023 <= rootVSmall
         || m013 <= rootVSmall
         || m123 <= rootVSmall
        )
            return scalar(-1);

        const scalar p3From012 =
            mag(((p3 - p0) & n012)/m012);

        const scalar p1From023 =
            mag(((p1 - p0) & n023)/m023);

        const scalar p2From013 =
            mag(((p2 - p0) & n013)/m013);

        const scalar p0From123 =
            mag(((p0 - p1) & n123)/m123);

        const scalar warp02 =
            Foam::max(p3From012, p1From023);

        const scalar warp13 =
            Foam::max(p2From013, p0From123);

        return Foam::max(warp02, warp13);
    };


    auto auditOtherCell =
    [&](const label faceI, const label cellI) -> label
    {
        if
        (
            faceI < 0
         || faceI >= mesh.nInternalFaces()
        )
            return -1;

        if( owner[faceI] == cellI )
            return neighbour[faceI];

        if( neighbour[faceI] == cellI )
            return owner[faceI];

        return -1;
    };


    // Idealized six-internal-face square-prism reference.
    //
    // For dimensions L,L,h and A=L/h, the exact OF13 normalized
    // face-area tensor expression reduces to:
    //
    //   D = 5832 A^2 / (A+2)^6
    //
    // This is only a reference curve. Boundary cells do not have all
    // physical boundary faces included in the OF determinant stencil.
    auto auditIdealInteriorDet =
    [&](const scalar aspect) -> scalar
    {
        if( aspect <= rootVSmall )
            return scalar(1);

        const scalar q = aspect + scalar(2);

        const scalar q2 = q*q;
        const scalar q4 = q2*q2;
        const scalar q6 = q4*q2;

        return
            scalar(5832)
           *aspect
           *aspect
           /Foam::max(q6, rootVSmall);
    };


    // ------------------------------------------------------------
    // Exact OpenFOAM-13 cell determinant.
    // ------------------------------------------------------------

    boolList auditInternalOrCoupled
    (
        mesh.nFaces(),
        false
    );

    for
    (
        label faceI=0;
        faceI<mesh.nInternalFaces();
        ++faceI
    )
    {
        auditInternalOrCoupled[faceI] = true;
    }

    const polyBoundaryMesh& auditPatches =
        mesh.boundaryMesh();

    forAll(auditPatches, patchI)
    {
        const polyPatch& pp = auditPatches[patchI];

        if( !pp.coupled() )
            continue;

        forAll(pp, patchFaceI)
        {
            auditInternalOrCoupled
            [
                pp.start() + patchFaceI
            ] = true;
        }
    }


    scalarField auditCellDet
    (
        mesh.nCells(),
        scalar(0)
    );

    const Vector<label>& auditMeshD =
        mesh.geometricD();

    label auditNDims = 0;
    label auditTwoD = -1;

    for
    (
        direction dir=0;
        dir<vector::nComponents;
        ++dir
    )
    {
        if( auditMeshD[dir] == 1 )
        {
            ++auditNDims;
        }
        else
        {
            auditTwoD = dir;
        }
    }


    const cellList& auditCells =
        mesh.cells();

    if( auditNDims == 1 )
    {
        auditCellDet = scalar(1);
    }
    else
    {
        forAll(auditCells, cellI)
        {
            const cell& cFaces =
                auditCells[cellI];

            scalar avgArea = scalar(0);
            label nInternalFaces = 0;

            forAll(cFaces, cfI)
            {
                const label faceI =
                    cFaces[cfI];

                if
                (
                    faceI >= 0
                 && faceI < label
                    (
                        auditInternalOrCoupled.size()
                    )
                 && auditInternalOrCoupled[faceI]
                )
                {
                    avgArea +=
                        mag(faceAreas[faceI]);

                    ++nInternalFaces;
                }
            }

            if( nInternalFaces == 0 )
            {
                auditCellDet[cellI] =
                    scalar(0);

                continue;
            }

            avgArea /=
                scalar(nInternalFaces);

            symmTensor areaTensor(Zero);

            forAll(cFaces, cfI)
            {
                const label faceI =
                    cFaces[cfI];

                if
                (
                    faceI >= 0
                 && faceI < label
                    (
                        auditInternalOrCoupled.size()
                    )
                 && auditInternalOrCoupled[faceI]
                )
                {
                    areaTensor +=
                        sqr
                        (
                            faceAreas[faceI]
                           /Foam::max
                            (
                                avgArea,
                                rootVSmall
                            )
                        );
                }
            }

            if( auditNDims == 2 )
            {
                if( auditTwoD == 0 )
                    areaTensor.xx() = 1;
                else if( auditTwoD == 1 )
                    areaTensor.yy() = 1;
                else
                    areaTensor.zz() = 1;
            }

            auditCellDet[cellI] =
                mag(det(areaTensor));
        }
    }


    const scalar auditWarnDet =
        scalar(1.0e-3);

    label auditDetFailAll = 0;
    scalar auditDetMin = GREAT;
    scalar auditDetSum = scalar(0);

    forAll(auditCellDet, cellI)
    {
        const scalar d =
            auditCellDet[cellI];

        auditDetMin =
            Foam::min(auditDetMin, d);

        auditDetSum += d;

        if( d < auditWarnDet )
            ++auditDetFailAll;
    }


    Info<<
        "CFMITCH CAUSAL AUDIT DETERMINANT PARITY:"
        << " cells=" << auditCellDet.size()
        << " fail=" << auditDetFailAll
        << " threshold=" << auditWarnDet
        << " min=" << auditDetMin
        << " avg="
        << (
               auditCellDet.size() > 0
             ? auditDetSum
              /scalar(auditCellDet.size())
             : scalar(0)
           )
        << nl;


    // ------------------------------------------------------------
    // Exact determinant population by wall-distance layer.
    //
    // Restrict the detailed layer histogram to ordinary six-face
    // cells so it represents the canonical BL population.
    // ------------------------------------------------------------

    const label auditNLayerBins = 15;

    labelList auditDetCellsByD
    (
        auditNLayerBins,
        label(0)
    );

    labelList auditDetFailByD
    (
        auditNLayerBins,
        label(0)
    );

    scalarField auditDetSumByD
    (
        auditNLayerBins,
        scalar(0)
    );

    scalarField auditDetMinByD
    (
        auditNLayerBins,
        GREAT
    );


    forAll(auditCells, cellI)
    {
        const label d =
            wallDistance[cellI];

        if
        (
            d < 0
         || d >= auditNLayerBins
         || auditCells[cellI].size() != 6
        )
            continue;

        ++auditDetCellsByD[d];

        auditDetSumByD[d] +=
            auditCellDet[cellI];

        auditDetMinByD[d] =
            Foam::min
            (
                auditDetMinByD[d],
                auditCellDet[cellI]
            );

        if
        (
            auditCellDet[cellI]
          < auditWarnDet
        )
        {
            ++auditDetFailByD[d];
        }
    }


    Info<<
        "CFMITCH CAUSAL AUDIT DETERMINANT BY DISTANCE:"
        << nl;

    for
    (
        label d=0;
        d<auditNLayerBins;
        ++d
    )
    {
        if( auditDetCellsByD[d] == 0 )
            continue;

        Info<<
            "  distance=" << d
            << " sixFaceCells="
            << auditDetCellsByD[d]
            << " detFail="
            << auditDetFailByD[d]
            << " failPct="
            << scalar(100)
             *scalar(auditDetFailByD[d])
             /scalar(auditDetCellsByD[d])
            << " minDet="
            << auditDetMinByD[d]
            << " avgDet="
            << auditDetSumByD[d]
             /scalar(auditDetCellsByD[d])
            << nl;
    }


    // ------------------------------------------------------------
    // Cross-layer causal audit.
    //
    // L:
    //   mean edge length of layer-parallel quad.
    //
    // ownerH / neighbourH:
    //   2 * centre-to-interface normal offset.
    //
    // This is intentionally a geometric thickness proxy rather than
    // claiming reconstruction of the original split-edge schedule.
    // ------------------------------------------------------------

    labelList auditCrossTotal
    (
        14,
        label(0)
    );

    labelList auditCrossFail
    (
        14,
        label(0)
    );

    scalarField auditFailWarpSum
    (
        14,
        scalar(0)
    );

    scalarField auditPassWarpSum
    (
        14,
        scalar(0)
    );

    scalarField auditFailLSum
    (
        14,
        scalar(0)
    );

    scalarField auditPassLSum
    (
        14,
        scalar(0)
    );

    scalarField auditFailWarpOverLSum
    (
        14,
        scalar(0)
    );

    scalarField auditPassWarpOverLSum
    (
        14,
        scalar(0)
    );

    scalarField auditFailAspectSum
    (
        14,
        scalar(0)
    );

    scalarField auditPassAspectSum
    (
        14,
        scalar(0)
    );

    labelList auditDetSideSamples
    (
        14,
        label(0)
    );

    labelList auditDetSideFails
    (
        14,
        label(0)
    );


    // L/h bins for determinant incidence.
    //
    // 0 <20
    // 1 20-40
    // 2 40-46
    // 3 46-60
    // 4 60-100
    // 5 100-200
    // 6 200-500
    // 7 >=500

    labelList auditAspectBinTotal
    (
        8,
        label(0)
    );

    labelList auditAspectBinDetFail
    (
        8,
        label(0)
    );


    auto auditAspectBin =
    [&](const scalar a) -> label
    {
        if( a < scalar(20) )
            return 0;

        if( a < scalar(40) )
            return 1;

        if( a < scalar(46) )
            return 2;

        if( a < scalar(60) )
            return 3;

        if( a < scalar(100) )
            return 4;

        if( a < scalar(200) )
            return 5;

        if( a < scalar(500) )
            return 6;

        return 7;
    };


    OFstream auditCrossCsv
    (
        "cfmitchCausalAuditCrossLayer.csv"
    );

    auditCrossCsv
        << "face,failed,distance,"
        << "warp,L,warpOverL,"
        << "ownerHalfOffset,neighbourHalfOffset,"
        << "ownerH,neighbourH,"
        << "ownerLoverH,neighbourLoverH,"
        << "ownerDet,neighbourDet,"
        << "ownerDetFail,neighbourDetFail,"
        << "ownerIdealInteriorDet,"
        << "neighbourIdealInteriorDet"
        << nl;


    for
    (
        label faceI=0;
        faceI<mesh.nInternalFaces();
        ++faceI
    )
    {
        const label ownCell =
            owner[faceI];

        const label neiCell =
            neighbour[faceI];

        if
        (
            ownCell < 0
         || neiCell < 0
         || ownCell >= mesh.nCells()
         || neiCell >= mesh.nCells()
        )
            continue;

        if
        (
            auditCells[ownCell].size() != 6
         || auditCells[neiCell].size() != 6
        )
            continue;

        const face& f =
            mesh.faces()[faceI];

        if( f.size() != 4 )
            continue;

        const label ownD =
            wallDistance[ownCell];

        const label neiD =
            wallDistance[neiCell];

        if( ownD < 0 || neiD < 0 )
            continue;

        label jump = ownD - neiD;

        if( jump < 0 )
            jump = -jump;

        if( jump != 1 )
            continue;

        const label nearD =
            Foam::min(ownD, neiD);

        if( nearD < 0 || nearD > 13 )
            continue;


        const scalar warp =
            auditQuadWarp(f);

        if( warp < scalar(0) )
            continue;


        const scalar L =
            auditFaceCharacteristicLength(f);

        const scalar areaMag =
            mag(faceAreas[faceI]);

        if
        (
            L <= rootVSmall
         || areaMag <= rootVSmall
        )
            continue;


        const vector n =
            faceAreas[faceI]/areaMag;

        const point& fc =
            faceCentres[faceI];

        const scalar ownHalfOffset =
            mag
            (
                (cellCentres[ownCell] - fc)
              & n
            );

        const scalar neiHalfOffset =
            mag
            (
                (cellCentres[neiCell] - fc)
              & n
            );


        const scalar ownerH =
            scalar(2)*ownHalfOffset;

        const scalar neighbourH =
            scalar(2)*neiHalfOffset;


        const scalar ownerAspect =
            L/Foam::max
            (
                ownerH,
                rootVSmall
            );

        const scalar neighbourAspect =
            L/Foam::max
            (
                neighbourH,
                rootVSmall
            );


        const bool failed =
            badFaceSet.found(faceI);

        ++auditCrossTotal[nearD];

        if( failed )
            ++auditCrossFail[nearD];


        const scalar warpOverL =
            warp/Foam::max(L, rootVSmall);

        const scalar meanAspect =
            scalar(0.5)
           *(
                ownerAspect
              + neighbourAspect
            );


        if( failed )
        {
            auditFailWarpSum[nearD] += warp;
            auditFailLSum[nearD] += L;
            auditFailWarpOverLSum[nearD] +=
                warpOverL;
            auditFailAspectSum[nearD] +=
                meanAspect;
        }
        else
        {
            auditPassWarpSum[nearD] += warp;
            auditPassLSum[nearD] += L;
            auditPassWarpOverLSum[nearD] +=
                warpOverL;
            auditPassAspectSum[nearD] +=
                meanAspect;
        }


        ++auditDetSideSamples[nearD];
        ++auditDetSideSamples[nearD];

        if
        (
            auditCellDet[ownCell]
          < auditWarnDet
        )
            ++auditDetSideFails[nearD];

        if
        (
            auditCellDet[neiCell]
          < auditWarnDet
        )
            ++auditDetSideFails[nearD];


        const label ownAspectBin =
            auditAspectBin(ownerAspect);

        const label neiAspectBin =
            auditAspectBin(neighbourAspect);

        ++auditAspectBinTotal[ownAspectBin];
        ++auditAspectBinTotal[neiAspectBin];

        if
        (
            auditCellDet[ownCell]
          < auditWarnDet
        )
            ++auditAspectBinDetFail[ownAspectBin];

        if
        (
            auditCellDet[neiCell]
          < auditWarnDet
        )
            ++auditAspectBinDetFail[neiAspectBin];


        auditCrossCsv
            << faceI << ','
            << label(failed) << ','
            << nearD << ','
            << warp << ','
            << L << ','
            << warpOverL << ','
            << ownHalfOffset << ','
            << neiHalfOffset << ','
            << ownerH << ','
            << neighbourH << ','
            << ownerAspect << ','
            << neighbourAspect << ','
            << auditCellDet[ownCell] << ','
            << auditCellDet[neiCell] << ','
            << label
               (
                   auditCellDet[ownCell]
                 < auditWarnDet
               ) << ','
            << label
               (
                   auditCellDet[neiCell]
                 < auditWarnDet
               ) << ','
            << auditIdealInteriorDet
               (
                   ownerAspect
               ) << ','
            << auditIdealInteriorDet
               (
                   neighbourAspect
               )
            << nl;
    }


    Info<<
        "CFMITCH CAUSAL AUDIT CROSS-LAYER BY DISTANCE:"
        << nl;

    forAll(auditCrossTotal, d)
    {
        const label total =
            auditCrossTotal[d];

        if( total == 0 )
            continue;

        const label fail =
            auditCrossFail[d];

        const label pass =
            total - fail;

        Info<<
            "  distance=" << d
            << " total=" << total
            << " fail=" << fail
            << " failPct="
            << scalar(100)
             *scalar(fail)
             /scalar(total);

        if( fail > 0 )
        {
            Info<<
                " failWarpAvgM="
                << auditFailWarpSum[d]
                 /scalar(fail)
                << " failLAvgM="
                << auditFailLSum[d]
                 /scalar(fail)
                << " failWarpOverLAvg="
                << auditFailWarpOverLSum[d]
                 /scalar(fail)
                << " failLoverHAvg="
                << auditFailAspectSum[d]
                 /scalar(fail);
        }

        if( pass > 0 )
        {
            Info<<
                " passWarpAvgM="
                << auditPassWarpSum[d]
                 /scalar(pass)
                << " passLAvgM="
                << auditPassLSum[d]
                 /scalar(pass)
                << " passWarpOverLAvg="
                << auditPassWarpOverLSum[d]
                 /scalar(pass)
                << " passLoverHAvg="
                << auditPassAspectSum[d]
                 /scalar(pass);
        }

        Info<<
            " determinantSideSamples="
            << auditDetSideSamples[d]
            << " determinantSideFails="
            << auditDetSideFails[d]
            << " determinantSideFailPct="
            << (
                   auditDetSideSamples[d] > 0
                 ? scalar(100)
                  *scalar(auditDetSideFails[d])
                  /scalar(auditDetSideSamples[d])
                 : scalar(0)
               )
            << nl;
    }


    Info<<
        "CFMITCH CAUSAL AUDIT L/H DETERMINANT BINS:"
        << nl;

    const char* auditAspectBinNames[8] =
    {
        "<20",
        "20-40",
        "40-46",
        "46-60",
        "60-100",
        "100-200",
        "200-500",
        ">=500"
    };

    forAll(auditAspectBinTotal, binI)
    {
        Info<<
            "  bin=" << auditAspectBinNames[binI]
            << " samples="
            << auditAspectBinTotal[binI]
            << " detFail="
            << auditAspectBinDetFail[binI]
            << " failPct="
            << (
                   auditAspectBinTotal[binI] > 0
                 ? scalar(100)
                  *scalar
                   (
                       auditAspectBinDetFail[binI]
                   )
                  /scalar
                   (
                       auditAspectBinTotal[binI]
                   )
                 : scalar(0)
               )
            << nl;
    }


    // ------------------------------------------------------------
    // Wall -> first internal layer warp inheritance.
    //
    // For a regular wall-adjacent six-face cell, locate its quad
    // internal face leading to wallDistance==1.
    //
    // This directly tests whether non-planarity of the wall footprint
    // is propagated into the first internal layer interface.
    // ------------------------------------------------------------

    label auditWallQuad = 0;
    label auditWallQuadFail = 0;

    scalar auditWallFailWarpSum = scalar(0);
    scalar auditWallPassWarpSum = scalar(0);

    scalar auditWallFailLSum = scalar(0);
    scalar auditWallPassLSum = scalar(0);

    scalar auditWallFailAspectSum = scalar(0);
    scalar auditWallPassAspectSum = scalar(0);

    label auditWallPaired = 0;

    scalar auditWallWarpSumX = scalar(0);
    scalar auditWallWarpSumY = scalar(0);
    scalar auditWallWarpSumXX = scalar(0);
    scalar auditWallWarpSumYY = scalar(0);
    scalar auditWallWarpSumXY = scalar(0);


    OFstream auditWallCsv
    (
        "cfmitchCausalAuditWall.csv"
    );

    auditWallCsv
        << "face,failed,owner,"
        << "wallWarp,wallL,"
        << "h1Proxy,LoverH1Proxy,"
        << "warpOverH1Proxy,"
        << "ownerDet,ownerDetFail,"
        << "firstInternalFace,"
        << "firstInternalWarp,"
        << "firstInternalL"
        << nl;


    forAll(wallPatch, patchFaceI)
    {
        const label faceI =
            wallPatch.start()
          + patchFaceI;

        const label cellI =
            owner[faceI];

        if
        (
            cellI < 0
         || cellI >= mesh.nCells()
        )
            continue;

        const face& f =
            mesh.faces()[faceI];

        if( f.size() != 4 )
            continue;

        const scalar wallWarp =
            auditQuadWarp(f);

        const scalar wallL =
            auditFaceCharacteristicLength(f);

        if
        (
            wallWarp < scalar(0)
         || wallL <= rootVSmall
        )
            continue;

        ++auditWallQuad;

        const bool failed =
            badFaceSet.found(faceI);

        if( failed )
            ++auditWallQuadFail;


        const scalar wallAreaMag =
            mag(faceAreas[faceI]);

        scalar h1Proxy = scalar(-1);
        scalar wallAspect = scalar(-1);
        scalar wallWarpOverH = scalar(-1);

        if( wallAreaMag > rootVSmall )
        {
            vector wallN =
                faceAreas[faceI]
               /wallAreaMag;

            // Boundary face owner orientation points out of the cell.
            const scalar halfThickness =
                mag
                (
                    (cellCentres[cellI]
                   - faceCentres[faceI])
                  & wallN
                );

            h1Proxy =
                scalar(2)
               *halfThickness;

            if( h1Proxy > rootVSmall )
            {
                wallAspect =
                    wallL/h1Proxy;

                wallWarpOverH =
                    wallWarp/h1Proxy;
            }
        }


        if( failed )
        {
            auditWallFailWarpSum +=
                wallWarp;

            auditWallFailLSum +=
                wallL;

            if( wallAspect >= scalar(0) )
                auditWallFailAspectSum +=
                    wallAspect;
        }
        else
        {
            auditWallPassWarpSum +=
                wallWarp;

            auditWallPassLSum +=
                wallL;

            if( wallAspect >= scalar(0) )
                auditWallPassAspectSum +=
                    wallAspect;
        }


        label firstInternalFace = -1;
        scalar firstInternalWarp = scalar(-1);
        scalar firstInternalL = scalar(-1);

        const cell& cFaces =
            auditCells[cellI];

        forAll(cFaces, cfI)
        {
            const label candidateFaceI =
                cFaces[cfI];

            if
            (
                candidateFaceI < 0
             || candidateFaceI >=
                mesh.nInternalFaces()
            )
                continue;

            const label otherCell =
                auditOtherCell
                (
                    candidateFaceI,
                    cellI
                );

            if
            (
                otherCell < 0
             || otherCell >= mesh.nCells()
             || wallDistance[otherCell] != 1
            )
                continue;

            const face& candidateFace =
                mesh.faces()[candidateFaceI];

            if( candidateFace.size() != 4 )
                continue;

            const scalar candidateWarp =
                auditQuadWarp(candidateFace);

            if( candidateWarp < scalar(0) )
                continue;

            firstInternalFace =
                candidateFaceI;

            firstInternalWarp =
                candidateWarp;

            firstInternalL =
                auditFaceCharacteristicLength
                (
                    candidateFace
                );

            break;
        }


        if
        (
            firstInternalFace >= 0
         && firstInternalWarp >= scalar(0)
        )
        {
            ++auditWallPaired;

            auditWallWarpSumX +=
                wallWarp;

            auditWallWarpSumY +=
                firstInternalWarp;

            auditWallWarpSumXX +=
                wallWarp*wallWarp;

            auditWallWarpSumYY +=
                firstInternalWarp
               *firstInternalWarp;

            auditWallWarpSumXY +=
                wallWarp
               *firstInternalWarp;
        }


        auditWallCsv
            << faceI << ','
            << label(failed) << ','
            << cellI << ','
            << wallWarp << ','
            << wallL << ','
            << h1Proxy << ','
            << wallAspect << ','
            << wallWarpOverH << ','
            << auditCellDet[cellI] << ','
            << label
               (
                   auditCellDet[cellI]
                 < auditWarnDet
               ) << ','
            << firstInternalFace << ','
            << firstInternalWarp << ','
            << firstInternalL
            << nl;
    }


    scalar auditWallFirstCorr = scalar(0);

    if( auditWallPaired > 1 )
    {
        const scalar n =
            scalar(auditWallPaired);

        const scalar numerator =
            n*auditWallWarpSumXY
          - auditWallWarpSumX
           *auditWallWarpSumY;

        const scalar denomX =
            n*auditWallWarpSumXX
          - auditWallWarpSumX
           *auditWallWarpSumX;

        const scalar denomY =
            n*auditWallWarpSumYY
          - auditWallWarpSumY
           *auditWallWarpSumY;

        const scalar denominator =
            Foam::sqrt
            (
                Foam::max
                (
                    scalar(0),
                    denomX*denomY
                )
            );

        if( denominator > rootVSmall )
        {
            auditWallFirstCorr =
                numerator/denominator;
        }
    }


    const label auditWallQuadPass =
        auditWallQuad
      - auditWallQuadFail;

    Info<<
        "CFMITCH CAUSAL AUDIT WALL:"
        << " quadFaces=" << auditWallQuad
        << " failed=" << auditWallQuadFail
        << " failPct="
        << (
               auditWallQuad > 0
             ? scalar(100)
              *scalar(auditWallQuadFail)
              /scalar(auditWallQuad)
             : scalar(0)
           )
        << " failWarpAvgM="
        << (
               auditWallQuadFail > 0
             ? auditWallFailWarpSum
              /scalar(auditWallQuadFail)
             : scalar(0)
           )
        << " passWarpAvgM="
        << (
               auditWallQuadPass > 0
             ? auditWallPassWarpSum
              /scalar(auditWallQuadPass)
             : scalar(0)
           )
        << " failLAvgM="
        << (
               auditWallQuadFail > 0
             ? auditWallFailLSum
              /scalar(auditWallQuadFail)
             : scalar(0)
           )
        << " passLAvgM="
        << (
               auditWallQuadPass > 0
             ? auditWallPassLSum
              /scalar(auditWallQuadPass)
             : scalar(0)
           )
        << " failLoverHAvg="
        << (
               auditWallQuadFail > 0
             ? auditWallFailAspectSum
              /scalar(auditWallQuadFail)
             : scalar(0)
           )
        << " passLoverHAvg="
        << (
               auditWallQuadPass > 0
             ? auditWallPassAspectSum
              /scalar(auditWallQuadPass)
             : scalar(0)
           )
        << " pairedWallFirstInternal="
        << auditWallPaired
        << " warpCorrelation="
        << auditWallFirstCorr
        << nl;


    // ------------------------------------------------------------
    // V5.7 same-pair fan concavity witness.
    //
    // A fan-touched cell is inferred conservatively when it contains
    // >=4 internal TRIANGULAR faces to the same neighbouring cell.
    //
    // Reproduce OF13 checkConcaveCells exactly for those cells and
    // classify the first witness pair:
    //
    //   fan-fan
    //   fan-other
    //   other-other
    //
    // This tests whether the late fan exposes incompatible folded
    // interface planes, as suggested by the external reviews.
    // ------------------------------------------------------------

    const scalar auditPlanarCosAngle =
        scalar(1.0e-6);

    label auditFanCells = 0;
    label auditFanConcaveCells = 0;

    label auditFanWitnessFanFan = 0;
    label auditFanWitnessFanOther = 0;
    label auditFanWitnessOtherOther = 0;


    OFstream auditFanCsv
    (
        "cfmitchFanConcavityWitness.csv"
    );

    auditFanCsv
        << "cell,nFaces,fanNeighbour,"
        << "repeatedFanTriangles,"
        << "concave,"
        << "planeFace,witnessFace,"
        << "planeFaceSize,witnessFaceSize,"
        << "planeIsFanTri,witnessIsFanTri,"
        << "dotMetric"
        << nl;


    forAll(auditCells, cellI)
    {
        const cell& cFaces =
            auditCells[cellI];

        std::map<label,label>
            fanNeighbourCounts;

        forAll(cFaces, cfI)
        {
            const label faceI =
                cFaces[cfI];

            if
            (
                faceI < 0
             || faceI >= mesh.nInternalFaces()
            )
                continue;

            if
            (
                mesh.faces()[faceI].size()
             != 3
            )
                continue;

            const label otherCell =
                auditOtherCell
                (
                    faceI,
                    cellI
                );

            if( otherCell >= 0 )
            {
                ++fanNeighbourCounts
                [
                    otherCell
                ];
            }
        }


        label fanNeighbour = -1;
        label repeatedFanTriangles = 0;

        for
        (
            std::map<label,label>::
                const_iterator
                it=fanNeighbourCounts.begin();
            it!=fanNeighbourCounts.end();
            ++it
        )
        {
            if
            (
                it->second
              > repeatedFanTriangles
            )
            {
                fanNeighbour =
                    it->first;

                repeatedFanTriangles =
                    it->second;
            }
        }


        if( repeatedFanTriangles < 4 )
            continue;


        ++auditFanCells;


        bool concave = false;

        label planeFaceI = -1;
        label witnessFaceI = -1;

        scalar witnessMetric = -GREAT;


        forAll(cFaces, i)
        {
            if( concave )
                break;

            const label fI =
                cFaces[i];

            const point& fC =
                faceCentres[fI];

            vector fN =
                faceAreas[fI];

            fN /=
                Foam::max
                (
                    mag(fN),
                    vSmall
                );

            if( owner[fI] != cellI )
                fN *= scalar(-1);


            forAll(cFaces, j)
            {
                if( j == i )
                    continue;

                const label fJ =
                    cFaces[j];

                const point& pt =
                    faceCentres[fJ];

                vector pC =
                    pt - fC;

                pC /=
                    Foam::max
                    (
                        mag(pC),
                        vSmall
                    );

                const scalar metric =
                    pC & fN;

                if
                (
                    metric
                  > -auditPlanarCosAngle
                )
                {
                    concave = true;

                    planeFaceI = fI;
                    witnessFaceI = fJ;
                    witnessMetric = metric;

                    break;
                }
            }
        }


        bool planeIsFanTri = false;
        bool witnessIsFanTri = false;

        if( planeFaceI >= 0 )
        {
            planeIsFanTri =
            (
                planeFaceI
              < mesh.nInternalFaces()
             && mesh.faces()[planeFaceI].size()
                == 3
             && auditOtherCell
                (
                    planeFaceI,
                    cellI
                )
                == fanNeighbour
            );
        }

        if( witnessFaceI >= 0 )
        {
            witnessIsFanTri =
            (
                witnessFaceI
              < mesh.nInternalFaces()
             && mesh.faces()[witnessFaceI].size()
                == 3
             && auditOtherCell
                (
                    witnessFaceI,
                    cellI
                )
                == fanNeighbour
            );
        }


        if( concave )
        {
            ++auditFanConcaveCells;

            if
            (
                planeIsFanTri
             && witnessIsFanTri
            )
            {
                ++auditFanWitnessFanFan;
            }
            else if
            (
                planeIsFanTri
             || witnessIsFanTri
            )
            {
                ++auditFanWitnessFanOther;
            }
            else
            {
                ++auditFanWitnessOtherOther;
            }
        }


        auditFanCsv
            << cellI << ','
            << cFaces.size() << ','
            << fanNeighbour << ','
            << repeatedFanTriangles << ','
            << label(concave) << ','
            << planeFaceI << ','
            << witnessFaceI << ','
            << (
                   planeFaceI >= 0
                 ? mesh.faces()[planeFaceI].size()
                 : label(-1)
               ) << ','
            << (
                   witnessFaceI >= 0
                 ? mesh.faces()[witnessFaceI].size()
                 : label(-1)
               ) << ','
            << label(planeIsFanTri) << ','
            << label(witnessIsFanTri) << ','
            << witnessMetric
            << nl;
    }


    Info<<
        "CFMITCH CAUSAL AUDIT FAN CONCAVITY:"
        << " inferredFanCells="
        << auditFanCells
        << " concave="
        << auditFanConcaveCells
        << " concavePct="
        << (
               auditFanCells > 0
             ? scalar(100)
              *scalar(auditFanConcaveCells)
              /scalar(auditFanCells)
             : scalar(0)
           )
        << " fanFanWitness="
        << auditFanWitnessFanFan
        << " fanOtherWitness="
        << auditFanWitnessFanOther
        << " otherOtherWitness="
        << auditFanWitnessOtherOther
        << nl;


    Info<<
        "CFMITCH CAUSAL AUDIT OUTPUT:"
        << " crossLayerCsv=cfmitchCausalAuditCrossLayer.csv"
        << " wallCsv=cfmitchCausalAuditWall.csv"
        << " fanCsv=cfmitchFanConcavityWitness.csv"
        << nl;

    Info<<
        "CFMITCH CAUSAL AUDIT V1:"
        << " end"
        << nl << endl;


    Info<< nl
        << "CFMITCH TET LOCATION CSV: cfmitchTetLocationDiag.csv"
        << nl
        << "End" << nl << endl;

    return 0;
}
