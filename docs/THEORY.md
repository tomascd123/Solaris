# SOLARIS: theory and methods

This document describes the geometry, the descriptors and the statistics behind `solaris.py`. Every formula and
threshold below is the one implemented in the code; the names in `code font` are the column names of the tables.
For how to run the program see [TUTORIAL.md](../TUTORIAL.md).

## Contents

1. [The problem](#1-the-problem)
2. [A reference frame attached to the dicopper site](#2-a-reference-frame-attached-to-the-dicopper-site)
3. [Cylindrical coordinates](#3-cylindrical-coordinates)
4. [Pocket landmarks](#4-pocket-landmarks)
5. [Substrate descriptors (L-Tyr)](#5-substrate-descriptors-l-tyr)
6. [The Cu2O2 core](#6-the-cu2o2-core)
7. [π-stacking of the substrate with HisB2](#7-π-stacking-of-the-substrate-with-hisb2)
8. [Trajectories](#8-trajectories)
9. [Statistics](#9-statistics)
10. [Quality checks](#10-quality-checks)
11. [Figures](#11-figures)
12. [Assumptions and limitations](#12-assumptions-and-limitations)
13. [References](#13-references)

---

## 1. The problem

Type-3 copper proteins (tyrosinases, catechol oxidases, hemocyanins) share a dicopper site: two copper ions, CuA
and CuB, each coordinated by three histidines (HisA1–A3 and HisB1–B3). Around this conserved core, the residues
that shape the substrate pocket differ between enzymes in sequence position, residue type and conformation.

Comparing pockets across a family usually relies on sequence alignments or on superposing whole proteins. Both
are indirect: a sequence alignment says nothing about where a side chain points, and a global superposition is
dominated by the protein fold, not by the active site.

SOLARIS compares the pockets in a frame attached to the site itself. Every structure, docking pose or trajectory
frame is superposed on the dicopper core, and every pocket atom is then described by where it sits around the
pseudo-two-fold axis of the site. Positions become directly comparable between enzymes, between states of one
enzyme (minima, docking poses) and along a simulation.

## 2. A reference frame attached to the dicopper site

### 2.1 The core

The core has 32 atoms: CuA, CuB and the five imidazole atoms (CG, ND1, CD2, CE1, NE2) of each of the six
histidines. These atoms are present in every type-3 copper site and are much more rigid than the rest of the
pocket.

### 2.2 The family template

SOLARIS stores a consensus position for each of the 32 core atoms (`TEMPLATE` in the code). The consensus was
obtained by generalized Procrustes analysis (GPA; Gower, 1975): all cores are superposed onto their mean, the
mean is recomputed, and the process is repeated until it converges. The result is the average shape of the site
in the family, free of the orientation of any single structure.

The template is stored already expressed in the SOLARIS frame (section 2.4), so superposing a structure onto the
template places it directly in that frame.

### 2.3 Superposition of each structure or frame (Kabsch)

For each structure or trajectory frame, the core atoms present (rows *p<sub>i</sub>*) are superposed onto the
corresponding template atoms (rows *q<sub>i</sub>*) by the rotation *R* and translation *t* that minimise

$$\sum_i \lVert R\,p_i + t - q_i \rVert^2 .$$

The solution is given by the Kabsch algorithm (Kabsch, 1976, 1978). With centred coordinates and the covariance
matrix

$$H = \sum_i (p_i - \bar p)(q_i - \bar q)^{\mathsf T} = U S V^{\mathsf T} \quad\text{(singular value decomposition)},$$

$$d = \mathrm{sign}\,\det(V U^{\mathsf T}), \qquad R = V \mathrm{diag}(1, 1, d)\, U^{\mathsf T}, \qquad t = \bar q - R\,\bar p .$$

The factor *d* forbids reflections. The whole structure is then transformed with the same *R* and *t*.

The root-mean-square deviation of the fitted core from the template, `Frame.core_RMSD`, measures how much this
site deviates from the family average. In the example structures it is 0.2–0.7 Å.

**Which copper is CuA.** The two copper atoms are labelled by their histidines: SOLARIS compares the sum of the
six Cu–NE2 distances for the two possible assignments and keeps the one in which HisA1–A3 bind CuA and HisB1–B3
bind CuB. The assignment is made on the first frame and kept for the whole trajectory.

### 2.4 Axes of the frame

The template defines the frame:

| Element | Definition |
|---|---|
| Origin | Midpoint of CuA and CuB |
| *z* axis | Pseudo-C2 axis of the site (below) |
| *x* axis | Perpendicular to *z*, pointing towards CuA |
| *y* axis | *z* × *x* (right-handed frame) |
| Sign of *z* | +*z* points to the side of the conserved Glu/Asn, the entrance of the substrate |

**The pseudo-C2 axis.** The CuA site (CuA + HisA1–A3) and the CuB site (CuB + HisB1–B3) are related by an
approximate two-fold rotation that maps A*k* onto B*k*. A rotation by 180° about a unit axis *n* through the
origin is the matrix 2*nn*<sup>T</sup> − *I*. The best axis maximises the overlap of the rotated A atoms
(*a<sub>k</sub>*) with the B atoms (*b<sub>k</sub>*):

$$\sum_k b_k \cdot (2nn^{\mathsf T} - I)\, a_k \;=\; n^{\mathsf T} M n - \sum_k a_k \cdot b_k,
\qquad M = \sum_k \left(a_k b_k^{\mathsf T} + b_k a_k^{\mathsf T}\right),$$

so *n* is the eigenvector of *M* with the largest eigenvalue. For the stored template this axis deviates by 0.26°
from *z*, and the CuA site rotated by 180° lies 1.37 Å RMSD from the CuB site: the site is only *pseudo*-symmetric,
which is why the axis is defined by a fit and not by any single pair of atoms.

**Why this frame.** Because the frame moves with the site, every coordinate in it is invariant to the rigid-body
motion of the protein. A residue that keeps its place in the pocket keeps its coordinates, whether the protein
rotates in the simulation box or two crystal structures are deposited in different orientations.

## 3. Cylindrical coordinates

Every point (*x*, *y*, *z*) in the SOLARIS frame is reported as

$$h = z, \qquad r = \sqrt{x^2 + y^2}, \qquad \varphi = \mathrm{atan2}(y, x) \bmod 360^\circ .$$

- *h*: height above (+) or below (−) the plane of the two coppers, towards the substrate entrance.
- *r*: distance from the pseudo-C2 axis.
- φ: azimuth around the axis; 0° points towards CuA and 180° towards CuB.

**φ near the axis.** φ is undefined on the axis, and close to it a small displacement changes it a lot: a
displacement δ perpendicular to the radius changes φ by about δ/*r* radians. At *r* = 2.5 Å, 0.5 Å changes φ by
about 11°; at *r* = 1 Å by about 29°. SOLARIS therefore sets φ to NaN when *r* < `r_min_phi` (2.5 Å by default):
those points have no φ in the tables, the summaries or the figures. The substrate OH and ring, which sit almost on
the axis, are described by *h*, *r* and distances instead.

## 4. Pocket landmarks

### 4.1 Roles and numbering

A landmark is a pocket residue defined by its *role* (conserved Glu, conserved Asn, carboxylate anchor,
gatekeeper, Ser/Ala next to the site, or any role defined in the configuration). Roles are written once per enzyme
in a reference numbering. A role that does not exist in an enzyme (for example the carboxylate anchor in JrTYR) is
set to `null` and skipped.

Simulation topologies are often renumbered. When the six His are named HD1–HD6 (as in MCPB.py topologies),
SOLARIS detects the offset between the reference numbering and the file from them; all six must give the same
offset. Otherwise `resid_offset` from the configuration is used. The residue types of the His, the conserved Glu
and the conserved Asn are checked, and every role must map to the same residue type in all the states of an
enzyme (section 10).

### 4.2 Three points per residue

| Point | Atoms averaged |
|---|---|
| CA | The Cα atom |
| SC | All heavy side-chain atoms (Cα for Gly) |
| TIP | The functional group: carboxylate for Glu/Asp, amide for Asn/Gln, guanidinium for Arg, NZ for Lys, the ring for His/Phe/Tyr/Trp, etc. (`TIP_ATOMS`; can be redefined per enzyme with `tip_atoms`) |

For each point SOLARIS reports *x*, *y*, *h*, *r*, φ and the distances to CuA and CuB. Each landmark has a
*plotted point* (`plot`: TIP for residues whose functional group matters, SC for steric residues such as the
gatekeeper); the figures and the SUMMARY use that point.

The orientation of the side chain is summarised by

$$dr = r_{\text{TIP}} - r_{\text{CA}}, \qquad dh = h_{\text{TIP}} - h_{\text{CA}} :$$

*dr* < 0 means that the side chain points towards the axis (into the pocket), *dh* > 0 that it points up, towards
the entrance.

## 5. Substrate descriptors (L-Tyr)

The substrate is L-tyrosine. Its atoms are mapped from the names in the topology to standard names (N, CA, C, CB,
CG, CD1, CD2, CE1, CE2, CZ, OH; GAFF2 names by default). If several copies are present, the one whose OH is
closest to the Cu–Cu midpoint in the first frame is used.

| Descriptor | Definition | Meaning |
|---|---|---|
| `L.OH.*`, `L.ring.*`, `L.CA.*` | Cylindrical coordinates of the OH, the ring centroid and the Cα | Position of the substrate in the pocket |
| `L.d_OH_CuA`, `L.d_OH_CuB` | Distances from the OH to each copper | Approach of the phenol to the site; their difference shows whether it is closer to one copper |
| `L.ring_tilt` | Angle between the ring normal and *z*, folded to 0–90° | 0° ring parallel to the Cu plane, 90° ring perpendicular to it |
| `L.ring_normal_vs_CuCu` | Angle between the ring normal and *x* (the Cu–Cu direction), 0–90° | Orientation of the ring relative to the two coppers |
| `L.angle_COH_vs_minus_z` | Angle between the CZ→OH bond and −*z* | 0°: the OH points straight down to the dicopper plane |
| `L.angle_CA_OH_vs_minus_z` | Angle between the Cα→OH vector and −*z* | Overall orientation of the molecule |
| `L.d_Cortho_Operox` | Shortest distance from an *ortho* carbon (CE1, CE2) to a peroxide oxygen | Geometry of the *ortho*-hydroxylation step, where a peroxide oxygen is transferred to the *ortho* carbon |
| `L.d_carboxylate_anchor` | Distance from the carboxylate carbon to the TIP of the anchor residue | Whether the substrate carboxylate is held by the anchor |
| `L.d_ring_gatekeeper` | Distance from the ring centroid to the SC of the gatekeeper | Steric contact with the gatekeeper |
| `L.d_ring_HisB2` and `L.HisB2_*` | See section 7 | π-stacking with HisB2 |

Ring normals are computed as the direction of least variance of the ring atoms: the right singular vector with the
smallest singular value of the centred coordinates. The π-stacking descriptors use the PLIP normal instead
(section 7).

## 6. The Cu2O2 core

| Descriptor | Definition | Meaning |
|---|---|---|
| `Cu.d_CuA_CuB` | Cu–Cu distance | Opening of the dicopper site |
| `Cu.max_d_Cu_NE2` | Largest of the six Cu–NE2 distances | Coordination check: a His that has left its copper |
| `Px.d_O_O` | O–O distance of the peroxide | Bound peroxide (about 1.4–1.5 Å) or cleaved |
| `Px.dihedral_CuA_O1_O2_CuB` | Absolute CuA–O1–O2–CuB dihedral | 180° for a planar side-on core; smaller values mean a butterfly distortion |
| `Px.Cu2O2_tilt` | Angle between the normal of the plane fitted to Cu, Cu, O, O and *z*, 0–90° | Orientation of the Cu2O2 core in the frame |

## 7. π-stacking of the substrate with HisB2

### 7.1 Geometry

With *c*<sub>T</sub> the centroid of the L-Tyr ring and *c*<sub>H</sub> the centroid of the HisB2 imidazole, the
ring normals are computed as in PLIP (Salentin et al., 2015): with the ring atoms listed in cyclic order
(CG, CD1, CE1, CZ, CE2, CD2 for L-Tyr; CG, ND1, CE1, NE2, CD2 for His), the normal is the cross product
(*a*<sub>2</sub> − *a*<sub>0</sub>) × (*a*<sub>0</sub> − *a*<sub>4</sub>), normalised. *n*<sub>H</sub> is oriented
towards the substrate side (*v* · *n*<sub>H</sub> ≥ 0 with *v* = *c*<sub>T</sub> − *c*<sub>H</sub>):

| Descriptor | Definition |
|---|---|
| `L.d_ring_HisB2` (*d*) | ‖*v*‖, distance between the ring centroids |
| `L.HisB2_h` (*h*<sub>⊥</sub>) | *v* · *n*<sub>H</sub>, height of the L-Tyr ring centre above the HisB2 plane |
| `L.HisB2_offset` | min(‖*v* − (*v* · *n*<sub>H</sub>) *n*<sub>H</sub>‖, ‖*v* − (*v* · *n*<sub>T</sub>) *n*<sub>T</sub>‖), the smaller sideways shift of one ring centre from the normal of the other ring |
| `L.HisB2_angle` (θ) | Angle between *n*<sub>T</sub> and *n*<sub>H</sub>, folded to 0–90° |

The first term of the offset is the arrow drawn in figure F4; the offset is that length or shorter. An offset of 0
means that one ring is exactly above the other (face-to-face); 1–2 Å is the parallel-displaced geometry, the most common and most stable
arrangement of stacked aromatic rings (Hunter and Sanders, 1990); beyond about 2 Å the rings no longer overlap
(the radius of an aromatic ring is about 1.4 Å).

### 7.2 Criterion

`L.HisB2_stack` is 1 when

$$d \le 5.5\ \text{Å}, \qquad \text{offset} \le 2.0\ \text{Å}, \qquad \theta \le 30^\circ\ (\text{parallel}) \ \text{or}\ \theta \ge 60^\circ\ (\text{T-shaped}),$$

and 0 otherwise. Normals, offset and thresholds are those of PLIP (Salentin et al., 2015). States close to a
threshold (an offset near 2 Å, an angle near 30°) can change class with the way the ring normal is computed;
check the values, not only the flag. The whole HisB2 imidazole is required; if atoms are missing, the
π-stacking descriptors and `L.d_ring_HisB2` are not computed and SUMMARY reports it.

Over a trajectory, the mean of `L.HisB2_stack` is the fraction of frames with stacking.

### 7.3 View of figure F4

F4 shows the rings looking down *n*<sub>H</sub>: origin at *c*<sub>H</sub>, first axis *u* = the HisB2 NE2 (the
nitrogen bound to CuB) projected on the HisB2 plane, second axis *n*<sub>H</sub> × *u*. All points are projected on
this plane, so the shift of the L-Tyr centre (the arrow) is a true length in the figure. For a trajectory, the regions holding 50 % and 90 % of
the L-Tyr ring centres are drawn, together with a representative frame: the one whose ring centre is closest to
the median centre.

## 8. Trajectories

Trajectories are read with MDAnalysis (topology + one or more trajectory files); multi-MODEL PDB files are also
read as trajectories. Frames are analysed every `stride` frames (10 by default), between `start` and `stop`.

**Periodic boundaries.** If a frame carries a usable box (all three box lengths > 40 Å), every needed atom is
placed at its minimum image relative to CuA before the fit. The box (*a*, *b*, *c*, α, β, γ) is converted to a
triclinic box matrix *M*; a vector *v* is reduced in fractional coordinates, *f* = *v* *M*<sup>−1</sup> −
round(*v* *M*<sup>−1</sup>), and the shortest of the 27 neighbouring images of *f* *M* is kept. A site split across
the box therefore gives the same result as an intact one.

**GaMD.** Gaussian accelerated MD (Miao et al., 2015) adds a boost potential, so the sampled frames are not a
Boltzmann ensemble. All trajectory statistics and figures in SOLARIS are **unweighted**: percentages are fractions
of frames, not populations. Equilibrium populations require reweighting with the boost energy of each frame, for
example by a second-order cumulant expansion (Miao et al., 2014), which SOLARIS does not do.

## 9. Statistics

**Per state** (`state_summary.csv`): number of values, mean, standard deviation (*n* − 1), 10th percentile, median
and 90th percentile of every descriptor.

**Angles φ** are circular quantities (359° and 1° are 2° apart), so their mean and spread are circular (Mardia and
Jupp, 2000). With the mean resultant length

$$\bar C = \tfrac{1}{n}\sum \cos\varphi_i, \quad \bar S = \tfrac{1}{n}\sum \sin\varphi_i, \quad R = \sqrt{\bar C^2 + \bar S^2},$$

the circular mean is atan2(*S̄*, *C̄*) and the circular standard deviation is √(−2 ln *R*), in degrees.

**Static structures within their trajectory** (`percentile_in_trajectory.csv`). For a static structure with value
*v* and the values *x<sub>j</sub>* of a trajectory of the same enzyme (at least 20 frames), the mid-rank percentile is

$$P = 100\left[\Pr(x_j < v) + \tfrac12 \Pr(x_j = v)\right].$$

A static structure is in the tail of its trajectory when *P* < 2.5 or *P* > 97.5 (`TAIL_PCT`). Angles φ and the
0/1 stacking flag are not ranked. This answers whether a minimum or a docking pose is a conformation that the
simulation actually visits.

**Landmark dispersion.** With RMS(*P*) = √(mean ‖*p* − *p̄*‖²) over the 3D positions of the plotted point,
`landmark_dispersion.csv` gives the RMS fluctuation of each landmark within each state. `family_dispersion.csv`
compares, per landmark, the spread within enzymes (their minima, or their GaMD frames) with the spread between
enzymes: *within* is pooled, √(Σ<sub>s</sub> *n*<sub>s</sub> RMS<sub>s</sub>² / Σ<sub>s</sub> *n*<sub>s</sub>),
i.e. every point around the mean of its own enzyme; *between* is the RMS of the enzyme means around the family
mean; the ratio is between / pooled within.

## 10. Quality checks

| Check | Threshold | Reported as |
|---|---|---|
| A His far from its copper | Cu–NE2 > 2.6 Å (structure), or in > 1 % of frames | Dicopper core checks |
| Distorted site | Core RMSD > 1.0 Å (structure), or in > 1 % of frames | Dicopper core checks |
| Wrong numbering | A role maps to different residue types in the states of one enzyme | Error (stops with `--strict`) |
| Wrong residue | His, conserved Glu or conserved Asn of the wrong type | The state is skipped |
| Missing atoms | Side-chain or TIP atoms absent; HisB2 imidazole incomplete; ligand atoms missing | Warnings |
| Missing ligand | No residue with the ligand name | Warning; no substrate descriptors |

The residue actually used for each role in each state, with the offset and how it was obtained, is written to
`residue_QC.csv` and to SUMMARY.

## 11. Figures

| Figure | What it shows | Encoding |
|---|---|---|
| F1_landmark_summary | *r*, *h* and φ of every landmark (plotted point) in the static structures, one sub-row per enzyme | Colour and marker = enzyme; large filled = global minimum, small = local minima, hollow = docking |
| F2_unrolled_cylinder | The cylinder cut open: φ against *h*, colour = *r*, static structures | Common *h* axis in all panels; the cut is placed in the widest empty sector of φ |
| F3_GaMD_&lt;enzyme&gt; | The landmarks of one enzyme along its trajectory (all, or those in `F3_landmarks`), with the static states on top | Top view and (*r*, *h*) view |
| F4_pistacking_HisB2 | The L-Tyr ring on HisB2, seen perpendicular to HisB2 (section 7.3) | Arrow = shift of the L-Tyr centre in the HisB2 plane; dotted circle = 2 Å offset limit |

**50 % and 90 % regions.** The points are binned in a 2D histogram, smoothed with a Gaussian kernel (σ = 1.4 bins),
and the smallest regions that hold 50 % and 90 % of the points are drawn (highest-density regions). With fewer than
30 points the points themselves are drawn.

## 12. Assumptions and limitations

- **One template for the family.** All sites are compared against one consensus core. Sites that deviate strongly
  from it (core RMSD > 1 Å) are flagged, but their coordinates remain less reliable.
- **The substrate is L-Tyr.** The substrate descriptors and the π-stacking analysis assume a phenol ring with OH,
  CZ and CE1/CE2. Other substrates need their atoms mapped onto these names.
- **π-stacking only with HisB2.** Other aromatic contacts (for example a Phe gatekeeper) are not evaluated.
- **Conventional thresholds.** 2.6 Å, 1.0 Å, 2.5 Å for φ and the π-stacking criteria are conventions; values close
  to a threshold should be reported as borderline.
- **Unweighted GaMD statistics** (section 8).
- **Assignments from the first frame.** The CuA/CuB labels and the substrate copy are chosen on the first frame of
  each trajectory.

## 13. References

- Gower, J. C. (1975). Generalized Procrustes analysis. *Psychometrika* 40, 33–51.
- Hunter, C. A. and Sanders, J. K. M. (1990). The nature of π–π interactions. *J. Am. Chem. Soc.* 112, 5525–5534.
- Kabsch, W. (1976). A solution for the best rotation to relate two sets of vectors. *Acta Cryst.* A32, 922–923.
- Kabsch, W. (1978). A discussion of the solution for the best rotation to relate two sets of vectors. *Acta Cryst.*
  A34, 827–828.
- Mardia, K. V. and Jupp, P. E. (2000). *Directional Statistics*. Wiley.
- Miao, Y., Sinko, W., Pierce, L., Bucher, D., Walker, R. C. and McCammon, J. A. (2014). Improved reweighting of
  accelerated molecular dynamics simulations for free energy calculation. *J. Chem. Theory Comput.* 10, 2677–2689.
- Miao, Y., Feher, V. A. and McCammon, J. A. (2015). Gaussian accelerated molecular dynamics: unconstrained
  enhanced sampling and free energy calculation. *J. Chem. Theory Comput.* 11, 3584–3595.
- Salentin, S., Schreiber, S., Haupt, V. J., Adasme, M. F. and Schroeder, M. (2015). PLIP: fully automated
  protein–ligand interaction profiler. *Nucleic Acids Res.* 43, W443–W447.
- Solomon, E. I., Heppner, D. E., Johnston, E. M., Ginsbach, J. W., Cirera, J., Qayyum, M., Kieber-Emmons, M. T.,
  Kjaergaard, C. H., Hadt, R. G. and Tian, L. (2014). Copper active sites in biology. *Chem. Rev.* 114, 3659–3853.
