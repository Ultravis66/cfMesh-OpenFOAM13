# cfMesh for OpenFOAM Foundation 13

A heavily modified and actively developed fork of **cfMesh** for **OpenFOAM Foundation 13**, focused on robust Cartesian/polyhedral meshing and boundary-layer generation for difficult CFD geometry.

This fork is aimed particularly at geometries where conventional boundary-layer extrusion becomes fragile: turbomachinery, propellers, very thin near-wall layers, periodic-sector boundaries, sharp features, multi-patch junctions, and strong refinement transitions.

> **Current release: v1.0.0**
>
> **Target platform: OpenFOAM Foundation 13**
>
> **CFD validation of the v1.0.0 meshes is now in progress.**

The original cfMesh was developed by **Franjo Juretić / Creative Fields**. This fork preserves the original GPL licensing and authorship while adding extensive OpenFOAM Foundation 13 compatibility, robustness, diagnostics, topology repair, and boundary-layer development.

---

## v1.0.0 milestone

Release 1 marks the transition from primarily mesher development to CFD validation.

The current code has demonstrated robust meshing on two particularly difficult geometries:

### NASA Rotor37

The current Rotor37 regression contains approximately **4.75 million cells** with **15 requested boundary layers** and a first-layer height of approximately **3 micrometres**.

Final ordinary OpenFOAM Foundation 13 `checkMesh` result:

| Metric | Result |
|---|---:|
| Points | 5,291,690 |
| Faces | 14,786,584 |
| Cells | 4,747,717 |
| Hexahedra | 4,444,069 |
| Prisms | 54,165 |
| Polyhedra | 249,483 |
| Negative-volume cells | 0 |
| Incorrect face pyramids | 0 |
| Maximum non-orthogonality | 88.678 degrees |
| Faces >= 90 degrees | 0 |
| Severe non-orthogonal faces (>70 degrees) | 881 |
| Maximum skewness | 3.992 |
| Ordinary `checkMesh` | **Mesh OK.** |

Boundary-layer coverage is approximately:

| Patch | Coverage |
|---|---:|
| Hub | 99.63% |
| Shroud | 89.57% |
| Periodic 1 | 99.04% |
| Periodic 2 | 99.24% |
| Blade 1 | 98.78% |
| Blade 2 | 98.99% |
| Blade 3 | 99.47% |
| Blade 4 | 99.78% |

Rotor37 exercises blade/hub/shroud junctions, periodic boundaries, sharp trailing-edge geometry, BL/no-BL transitions, cross-patch BL seams, and very small first-layer heights.

---

## Propeller regression

The second major regression is a complex marine-propeller geometry with a very thin wall-normal boundary layer.

| Metric | Result |
|---|---:|
| Points | 8,061,652 |
| Faces | 23,804,465 |
| Cells | 7,869,273 |
| Hexahedra | 7,776,252 |
| Prisms | 3,597 |
| Polyhedra | 89,424 |
| Requested propeller BL layers | 14 |
| Minimum first-layer height | approximately 7.5 micrometres |
| Nominal first-layer height | 7.5-10 micrometres |
| Propeller BL coverage | 96.33% |
| Negative-volume cells | 0 |
| Incorrect face pyramids | 0 |
| Maximum aspect ratio | 440.46 |
| Maximum non-orthogonality | 88.21 degrees |
| Faces >= 90 degrees | 0 |
| Severe non-orthogonal faces (>70 degrees) | 36 |
| Maximum skewness | 3.847 |

The very thin first layer is intentional for wall-resolved marine CFD.

### Advanced geometry diagnostics

The propeller still reports several failures under:

~~~bash
checkMesh -allTopology -allGeometry
~~~

Current populations include:

| Diagnostic | Count |
|---|---:|
| Face-tet decomposition failures | 138,732 |
| Cells with determinant < 0.001 | 684,625 |
| Concave cells | 10,612 |
| Faces with interpolation weight < 0.05 | 1,378 |
| Faces with volume ratio < 0.01 | 1,176 |

These remain active development targets.

They are disclosed deliberately. A successful ordinary mesh-quality result does not by itself prove CFD accuracy or solver robustness.

---

## CFD validation status

**CFD validation is now underway.**

An earlier development mesh from this project successfully advanced an OpenFOAM CFD solution and produced velocity, pressure, turbulence, and `yPlus` fields.

The exact **v1.0.0 Rotor37 and propeller meshes** are now entering solver validation.

At this stage the project does **not** claim validated quantitative:

- propeller thrust or torque;
- cavitation extent;
- compressor pressure ratio;
- compressor efficiency;
- compressor mass flow;
- aerodynamic or hydrodynamic performance accuracy.

Meshing validation and CFD validation are separate milestones.

---

## Boundary-layer architecture

The primary development architecture is:

~~~text
boundaryLayerArchitecture constraintPlanner;
~~~

The `constraintPlanner` path treats boundary-layer creation as a constrained topology problem instead of independently extruding every surface.

Major capabilities include:

- patch-specific layer counts and first-layer heights;
- neighboring-face layer-count compatibility;
- configurable maximum layer-count steps;
- cross-patch BL seam synchronization;
- stable BL provenance;
- BL/no-BL transition handling;
- topology validation;
- transactional repair and rollback;
- proof-gated topology reconstruction;
- OpenFOAM-aligned geometry diagnostics.

The goal is to preserve useful near-wall resolution whenever possible instead of simply deleting difficult layers.

---

## Automatic boundary-junction discovery

v1.0.0 includes automatic detection of multi-patch boundary junctions from actual surface topology and resolved boundary-layer intent.

The mesher can identify relationships such as:

- BL / BL seams;
- BL / no-BL termination seams;
- multi-patch junctions;
- BL/BL/neutral junctions;
- BL/neutral/neutral junctions.

This reduces reliance on hard-coded patch names.

### Current clearance limitation

Real non-touching clearances, such as a compressor blade tip separated from a shroud, are a different geometric problem.

The current gap system still uses explicit `gapPatchPairs` and proximity detection.

Automatic opposing-surface clearance detection and BL-envelope negotiation are planned future work.

---

## Major v1.0.0 improvements

Release 1 includes extensive work in:

- OpenFOAM Foundation 13 compatibility;
- boundary-layer topology robustness;
- cross-patch layer synchronization;
- periodic face-fan repair;
- multi-patch junction handling;
- stable boundary-layer provenance;
- prospective quality checks;
- transactional topology repair;
- lossless 17-digit mesh serialization;
- large polyhedral topology handling;
- OpenMP race hardening;
- surface projection and feature handling;
- OpenFOAM geometry parity diagnostics.

Several general correctness issues were also fixed, including:

- `meshSubset` equality logic;
- `processorBoundaryPatch` type comparison;
- `quadricFitting` eigenvector normal extraction;
- initializer ordering;
- compiler-warning cleanup.

The current tested `libmeshLibrary.so` build completes with:

~~~text
WARNING COUNT: 0
~~~

---

## Build

This project targets **OpenFOAM Foundation 13**.

Load OpenFOAM:

~~~bash
source /opt/openfoam13/etc/bashrc
~~~

Clone:

~~~bash
git clone https://github.com/Ultravis66/cfMesh-OpenFOAM13.git cfMesh
cd cfMesh
~~~

For the v1.0.0 release branch:

~~~bash
git switch release/v1.0.0
~~~

Build:

~~~bash
./Allwclean
./Allwmake
~~~

---

## Running cartesianMesh

The primary supported mesher in this fork is `cartesianMesh`.

For deterministic development and regression testing:

~~~bash
export OMP_NUM_THREADS=1
cartesianMesh
checkMesh
~~~

For deeper geometry diagnostics:

~~~bash
checkMesh -allTopology -allGeometry
~~~

Highly anisotropic boundary-layer meshes should not be judged from a single quality metric.

Relevant checks include:

- positive cell volume;
- face-pyramid orientation;
- non-orthogonality;
- skewness;
- interpolation weight;
- volume ratio;
- determinant;
- face-tet quality;
- BL continuity and coverage;
- actual CFD solver behavior.

---

## Current development priorities

Following v1.0.0, the main priorities are:

1. CFD validation of the Rotor37 and propeller meshes.
2. High-thread-count/OpenMP regression testing.
3. Reduction of advanced propeller geometry failures without destroying near-wall resolution.
4. Further automatic boundary-relation discovery.
5. Automatic non-touching clearance detection for turbomachinery blade-tip/shroud gaps.
6. Continued reduction of manual patch-role configuration.

---

## Project history

This project targets **OpenFOAM Foundation 13**.

Earlier work involved bringing cfMesh forward through newer OpenFOAM versions, but active mesher development and validation are now focused on Foundation 13.

OpenFOAM 12 is **not an actively developed target** of this repository.

---

## License and original authorship

cfMesh is free software distributed under the GNU General Public License.

The original cfMesh project was developed by **Franjo Juretić and Creative Fields**.

This repository is an OpenFOAM Foundation 13 port and extensive engineering fork. Original copyright, licensing, and authorship are preserved in accordance with the existing source headers and GPL requirements.

---

## Release status

**v1.0.0 — mesher validation milestone reached; CFD validation in progress.**
