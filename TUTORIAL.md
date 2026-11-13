# SOLARIS tutorial

This tutorial takes you from a set of tyrosinase structures and GaMD trajectories to the SOLARIS tables and
figures. The repository includes the structures of seven tyrosinases (GaMD free-energy minima and docking
receptors) so you can run every step before using your own files.

0. [Quick start](#0-quick-start)
1. [Install](#1-install)
2. [Prepare the input files](#2-prepare-the-input-files)
3. [Write the configuration](#3-write-the-configuration)
4. [Check residues and numbering](#4-check-residues-and-numbering)
5. [Analyse static structures](#5-analyse-static-structures)
6. [Analyse GaMD trajectories](#6-analyse-gamd-trajectories)
7. [Read the results](#7-read-the-results)
8. [Redraw figures](#8-redraw-figures)
9. [Troubleshooting](#9-troubleshooting)

---

## 0. Quick start

```bash
git clone https://github.com/tomascd123/SOLARIS.git
cd SOLARIS
pip install -r requirements.txt
python3 solaris.py examples/TYR_static.json --check
python3 solaris.py examples/TYR_static.json
```

The second command analyses the 27 example structures in about 10 s and writes `examples/output_static/`.
The same output is already included in the repository, so you can compare it with yours.

### Example data

```
examples/
├── TYR_static.json          configuration for the 7 enzymes (minima + docking)
├── TYR_GaMD.json            configuration template for the GaMD trajectories (not included)
├── structures/
│   ├── ScTYR/  ScTYR_global_renum.pdb  ScTYR_local1..3_renum.pdb  ScTYR_dock_prepared.pdb
│   ├── VvTYR/  ...
│   └── JrTYR/  ...
└── output_static/           reference output of TYR_static.json
```

| Enzyme | PDB template | Minima | Docking |
|---|---|---|---|
| ScTYR | 1WX2 | global + 3 local | receptor |
| VvTYR | 2P3X | global + 2 local | receptor |
| AbTYR | 2Y9X | global + 1 local | receptor |
| BmTYR | 3NQ1 | global + 2 local | receptor |
| HsTYR | 5M8L (model) | global + 1 local | receptor |
| VsTYR | 8BBQ | global + 2 local | receptor |
| JrTYR | 5CE9 | global + 2 local | complex with L-Tyr |

The minima are GaMD free-energy minima already renumbered like the docking structures
(`tools/renumber_like_docking.py`). The docking files of six enzymes are the prepared receptors (Cu and peroxide,
no ligand), so their L-Tyr descriptors are empty; JrTYR's docking file includes the docked L-Tyr.

## 1. Install

```bash
git clone https://github.com/tomascd123/SOLARIS.git
cd SOLARIS
pip install -r requirements.txt
```

SOLARIS needs Python ≥ 3.8 with numpy, pandas and matplotlib. MDAnalysis is only required to read trajectories
(`.parm7` + `.nc`/`.dcd`/`.xtc`); static structures are read directly from PDB files.

Test the installation:

```bash
python3 solaris.py -h
```

## 2. Prepare the input files

SOLARIS works with two kinds of input:

| Input | Format | Typical use |
|---|---|---|
| Static structures | PDB (one model) | GaMD free-energy minima, docking poses, crystal structures |
| Trajectories | topology + trajectory (`.parm7` + `.nc`) or multi-MODEL PDB | GaMD or conventional MD |

Every file must contain:

- **the two Cu atoms** (atom or residue named `CU`, `CU1`, `CU2`, or element Cu);
- **the six first-shell His** with their imidazole rings (CG, ND1, CD2, CE1, NE2);
- the residues you want to follow (conserved Glu and Asn, anchor, gatekeeper, …);
- optionally, **the peroxide** (default residue `OY1`, atoms `O1`, `O2`) and **the substrate** (default residue
  `TY1`; see [Ligand atom names](#ligand-atom-names)).

Water and ions are ignored, so you can pass solvated files as they are.

### Residue numbering

Each enzyme is described in **one reference numbering** (usually the crystal or docking numbering). Files that
use a different numbering are handled in two ways:

- **MCPB.py topologies.** If the six His are named `HD1`…`HD6`, SOLARIS finds them and computes the numbering
  offset automatically. This covers GaMD topologies and minima extracted from them.
- **Other files.** Give the offset explicitly with `"resid_offset"` in the state entry (default 0).

If your minima come from the GaMD topology and you prefer them to share the docking numbering, renumber them
once:

```bash
python3 tools/renumber_like_docking.py examples/structures/ScTYR/ScTYR_dock_prepared.pdb \
        my_minima/ScTYR_global.pdb my_minima/ScTYR_local1.pdb
```

This writes `ScTYR_global_renum.pdb` (and a `.txt` report with the residue mapping) next to each minimum. Read
the report: any line starting with `different:` means the two sequences do not match at that position.

## 3. Write the configuration

The configuration is one JSON file. Paths inside it are relative to the JSON file. A minimal example for one
enzyme with a global minimum, a docking pose and a GaMD trajectory:

```json
{
  "output_dir": "output_static",
  "landmarks": {
    "Conserved_Glu":      {"role": "Glu",        "plot": "TIP"},
    "Conserved_Asn":      {"role": "Asn",        "plot": "TIP"},
    "Carboxylate_anchor": {"role": "anchor",     "plot": "TIP"},
    "Gatekeeper":         {"role": "gatekeeper", "plot": "SC"}
  },
  "systems": {
    "ScTYR": {
      "roles": {"HisA1": 38, "HisA2": 54, "HisA3": 63, "HisB1": 190, "HisB2": 194, "HisB3": 216,
                "Glu": 182, "Asn": 191, "anchor": 55, "gatekeeper": 204},
      "structures": {
        "global":  "structures/ScTYR/ScTYR_global_renum.pdb",
        "docking": {"pdb": "structures/ScTYR/ScTYR_dock_prepared.pdb",
                    "peroxide": {"resname": "PER", "atoms": ["O1", "O2"]}},
        "GaMD":    {"top": "gamd/1wx2_centered.parm7", "traj": ["gamd/1wx2_prod_centered.nc"], "stride": 10}
      }
    }
  }
}
```

### Roles

`roles` maps each structural role to a residue number in the reference numbering of that enzyme.

| Role | Meaning |
|---|---|
| `HisA1`–`HisA3` | His coordinating CuA (required) |
| `HisB1`–`HisB3` | His coordinating CuB (required) |
| `Glu`, `Asn` | conserved Glu and Asn (checked to be GLU/ASN) |
| `anchor` | residue that binds the substrate carboxylate |
| `gatekeeper` | gatekeeper residue above CuA |
| any other name | any residue you want to follow, if a landmark uses it |

Write `"A:182"` to select a chain, and `"TODO"` (or leave it out) for a role an enzyme does not have; landmarks
that use it are skipped for that enzyme. Mark a landmark `"optional": true` to skip it silently.

The His order matters: HisA*k* and HisB*k* must be the pairs related by the pseudo-C2 symmetry of the site,
because the template is fitted atom by atom.

### Landmarks

Each landmark is followed at three points: `CA`, `SC` (centroid of the side-chain heavy atoms) and `TIP`
(functional group: carboxylate of Glu/Asp, amide of Asn/Gln, NZ of Lys, guanidinium of Arg, ring of aromatic
residues…). `plot` chooses which point is drawn in the figures; all three are always written to the tables.
Use `"tip": ["CD", "OE1", "OE2"]` in a landmark, or `"tip_atoms": {"anchor": "NZ"}` in an enzyme, to override the
TIP atoms.

### States

The keys of `structures` are state names. The name sets the category used in the figures:

| Name starts with | Category | Marker in the figures |
|---|---|---|
| `global` | global minimum | large, filled |
| `local` | local minimum | small, filled, linked to the global |
| `docking` | docking pose | hollow |
| `gamd`, `md`, `traj` | trajectory | density and subsampled frames |

A state is a trajectory if its entry has `traj` (or if the PDB has several MODELs). Trajectory entries accept
`start`, `stop` and `stride` (default 10).

### Ligand atom names

SOLARIS describes the substrate as L-Tyr. By default it looks for residue `TY1` with the GAFF2 atom names below;
change them under `"ligand"` globally, per enzyme or per state.

```json
"ligand": {"resname": "TY1", "atoms": {"N": "N1", "CA": "C1", "C": "C2", "CB": "C3", "CG": "C4",
           "CD1": "C5", "CD2": "C6", "CE1": "C7", "CE2": "C8", "CZ": "C9", "OH": "O2"}}
```

The ring atoms and OH are required; N, CA, C and CB are optional (descriptors that need them are skipped).
Standard PDB names (CA, CZ, OH…) are also recognised.

### Other settings

| Key | Default | Effect |
|---|---|---|
| `output_dir` | `SOLARIS_out` | output folder |
| `r_min_phi` | 2.5 | points closer than this to the axis (Å) are not drawn in φ |
| `r_color` | [0, 12] | colour scale of *r* in F2, fixed so runs are comparable |
| `points_per_trajectory` | 400 | frames drawn per trajectory (all frames are used in the tables) |
| `F3_landmarks` | all | landmarks drawn in F3, in this order; pocket labels plus `"L-Tyr OH"` and `"L-Tyr ring"` |
| `dpi` | 300 | figure resolution |

## 4. Check residues and numbering

Always start with `--check`. It opens every file, finds the residues and stops before any calculation:

```bash
python3 solaris.py examples/TYR_static.json --check
```

```
[OK] ScTYR/global: 1 frame(s), offset +0 (HD1-HD6); Conserved_Glu=GLU182, Conserved_Asn=ASN191, Carboxylate_anchor=ARG55, Gatekeeper=GLY204, A/S_adj=SER206; ligand yes; peroxide yes
[OK] ScTYR/docking: 1 frame(s), offset +0 (none (resid_offset 0)); Conserved_Glu=GLU182, Conserved_Asn=ASN191, ...; ligand no; peroxide yes
...
[OK] JrTYR/docking: 1 frame(s), offset +0 (none (resid_offset 0)); Conserved_Glu=GLU235, ...; ligand yes; peroxide yes

Check finished.
```

The minima take their offset from the HD1–HD6 residue names; the docking receptors have plain HIS residues and
use `resid_offset` 0. A missing file appears as `[SKIP] ... missing <path>` and is left out of the run.

For each state, check that:

- every landmark shows the residue type you expect (e.g. `Conserved_Glu=GLU…`);
- `ligand yes` / `peroxide yes` appear where they should (`ligand no` is expected for a receptor without ligand);
- there is no `[ERROR]` line. An error such as `role Glu is GLU/LYS depending on the state` means one file uses
  another numbering: fix its `resid_offset` or renumber it.

## 5. Analyse static structures

```bash
python3 solaris.py examples/TYR_static.json
```

Useful variants:

```bash
python3 solaris.py examples/TYR_static.json --only ScTYR VvTYR       # some enzymes
python3 solaris.py examples/TYR_static.json --states global docking  # some states
python3 solaris.py examples/TYR_static.json --no-docking --out minima_only
```

Missing files are skipped and listed at the end (and in `SUMMARY.md`); add `--strict` to stop instead.

## 6. Analyse GaMD trajectories

The trajectories are too large for the repository. Copy them to your machine and replace the `TODO/` paths in
`examples/TYR_GaMD.json` with their location (absolute paths or relative to the JSON file):

```json
"GaMD": {"top": "/data/gamd/ScTYR/1wx2_centered.parm7", "traj": ["/data/gamd/ScTYR/1wx2_prod_centered.nc"], "stride": 10}
```

Several trajectory files of one simulation can be listed in `traj`; they are read in order. Trajectories are long,
so run one enzyme at a time, each in its own folder:

```bash
python3 solaris.py examples/TYR_GaMD.json --only ScTYR --states GaMD --out SOLARIS_GaMD/ScTYR
```

or all enzymes in sequence:

```bash
bash tools/run_gamd_per_enzyme.sh            # all enzymes
bash tools/run_gamd_per_enzyme.sh ScTYR      # one enzyme
```

Notes:

- Periodic boundaries are handled: the site is rebuilt around CuA with the minimum-image convention, so the
  trajectory does not need to be imaged beforehand.
- The default `stride` is 10. Lower it for short trajectories; the percentile table needs at least 20 frames.
- Trajectory figures are **unweighted** distributions of the sampled frames.

To compare minima and docking with the trajectory of the same enzyme in one run (percentile table and F3 figure),
put the static states and the GaMD entry in the same configuration.

## 7. Read the results

```
examples/output_static/
├── SUMMARY.md
├── tables/
│   ├── descriptors_per_frame.csv
│   ├── landmarks_long.csv
│   ├── state_summary.csv
│   ├── percentile_in_trajectory.csv
│   ├── residue_QC.csv
│   └── pistack_geometry.csv
└── figures/
    ├── F1_landmark_summary.png
    ├── F2_unrolled_cylinder.png
    ├── F3_GaMD_<enzyme>.png
    └── F4_pistacking_HisB2.png / F4_pistacking_HisB2_GaMD.png
```

### Start with `SUMMARY.md`

Open `examples/output_static/SUMMARY.md` as an example of a complete report.

It lists, in this order:

1. **States not processed** and why.
2. **Dicopper core checks.** A His farther than 2.6 Å from its Cu, or a core RMSD to the template above 1.0 Å,
   points to lost coordination or a wrong His assignment. Check these before interpreting anything else.
3. **Warnings** (missing atoms, duplicated residue numbers…).
4. **Residues used** for every role in every state.
5. **Key values per state** (core RMSD, Cu–Cu, OH–CuA, OH–CuB, C*ortho*–O, ring tilt) and *h*/*r*/φ of every
   landmark.
6. **Static structures in the tails of their trajectory** (percentile < 2.5 or > 97.5).

### Coordinates

All positions are in the common frame (Å, degrees):

- ***h***: height above the dicopper core along the pseudo-C2 axis (+*h* towards the substrate entry);
- ***r***: distance from the axis (the substrate OH and ring sit at small *r*; pocket residues form the wall);
- **φ**: azimuth around the axis, 0° towards CuA and 180° towards CuB.

φ is undefined on the axis, so points with *r* < `r_min_phi` have no φ (NaN in the tables, left out of the φ panels).

### Tables

| Table | One row per | Use |
|---|---|---|
| `descriptors_per_frame.csv` | structure/frame | all descriptors; input for your own statistics or plots |
| `landmarks_long.csv` | point (CA/SC/TIP, L-Tyr OH/ring/CA) per frame | easy filtering and plotting |
| `state_summary.csv` | state × descriptor | mean, SD, p10, median, p90 (circular statistics for φ) |
| `percentile_in_trajectory.csv` | static state × descriptor | where each minimum/docking pose falls in its trajectory |
| `residue_QC.csv` | state × role | residue actually used, offset and how it was obtained |
| `pistack_geometry.csv` | point per frame | L-Tyr ring and HisB2 seen down the HisB2 normal (used by F4) |

Descriptor names start with the group: `Frame.` (fit), `Cu.` and `Px.` (core), `<landmark>.` (pocket residues)
and `L.` (substrate). For example `Conserved_Glu.TIP.h`, `Gatekeeper.dr` or `L.d_OH_CuA`.

### Figures

- **F1_landmark_summary**: *r*, *h* and φ of every landmark in the static structures, one row per landmark and
  one sub-row per enzyme. Every value of the family in one figure.
- **F2_unrolled_cylinder**: the cylinder cut open (φ vs *h*), colour = *r*, static structures. All panels share
  the *h* axis, so a conserved residue looks as tight as it is next to a variable one. The cut is placed
  automatically in the widest empty sector.
- **F3_GaMD_&lt;enzyme&gt;**: the landmarks of one enzyme along its trajectory (regions holding 50 % and 90 % of
  the frames, subsampled frames and median), with the static states on top. `F3_landmarks` chooses which
  landmarks are drawn; all by default. The fluctuation of each landmark is in `landmark_dispersion.csv`.
- **F4_pistacking_HisB2**: the L-Tyr ring seen perpendicular to HisB2, one panel per enzyme. The arrow is the
  shift of the L-Tyr ring centre in the HisB2 plane and the dotted circle the 2 Å offset limit.
  π-stacking needs d ≤ 5.5 Å, offset ≤ 2 Å and an angle between the ring planes ≤ 30° (parallel) or ≥ 60°
  (T-shaped); ring normals, offset and thresholds as in PLIP, computed by SOLARIS itself (the offset is the
  smaller of the two projections of one ring centre onto the other ring plane, so it can be shorter than the
  arrow). The whole HisB2 imidazole is needed; otherwise these descriptors are skipped with a warning. Descriptors: `L.HisB2_angle`,
  `L.HisB2_offset`, `L.HisB2_h` and `L.HisB2_stack` (1/0; its mean over a trajectory is the fraction of
  frames stacked). The GaMD version shows the regions holding 50/90 % of the L-Tyr ring centres and a
  representative frame (the one closest to the median); the percentages are unweighted.

Each enzyme keeps the same colour and marker in every figure.

## 8. Redraw figures

Figures can be redrawn from the saved tables without recomputing, e.g. to change the resolution or plot a
subset:

```bash
python3 solaris.py examples/TYR_static.json --figures-only --dpi 600
python3 solaris.py examples/TYR_static.json --figures-only --only ScTYR HsTYR
```

`--out` must point to the folder that holds `tables/`.

## 9. Troubleshooting

| Message | Cause and fix |
|---|---|
| `expected 2 Cu, found N` | the file has no Cu or extra Cu atoms; check atom names |
| `HisA1 = XXX123 is not a His (numbering?)` | wrong numbering for this file: set `resid_offset` or renumber it |
| `HD1-HD6 do not match a single numbering offset` | the His roles are not the six His named HD1–HD6 in this topology |
| `resid N appears 2 times ...; write 'chain:resid'` | two residues share the number: give the chain, e.g. `"A:404"` |
| `ligand TY1...: missing ...` | ligand atom names differ: set `ligand.atoms` |
| `role X is GLU/LYS depending on the state` | one state uses a different numbering |
| `trajectories need MDAnalysis` | `pip install MDAnalysis` |
| `--figures-only: ... not found` | `--out`/`output_dir` does not point to a previous run |
