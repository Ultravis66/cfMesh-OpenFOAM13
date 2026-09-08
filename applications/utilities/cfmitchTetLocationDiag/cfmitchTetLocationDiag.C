#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyMeshTetDecomposition.H"
#include "labelHashSet.H"
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
        mesh.boundaryMesh().findPatchID(wallPatchName);

    if (wallPatchI < 0)
    {
        FatalErrorInFunction
            << "Cannot find patch " << wallPatchName
            << exit(FatalError);
    }

    const labelList& owner = mesh.faceOwner();
    const labelList& neighbour = mesh.faceNeighbour();

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

    Info<< nl
        << "CFMITCH TET LOCATION CSV: cfmitchTetLocationDiag.csv"
        << nl
        << "End" << nl << endl;

    return 0;
}
