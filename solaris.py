#!/usr/bin/env python3
"""
SOLARIS: Symmetry-Oriented Ligand Approach Reference In type-3 copper Sites.

Places every structure or trajectory frame of a type-3 copper protein in a common frame attached to the
dicopper site and reports pocket residues, substrate and Cu2O2 core in cylindrical coordinates (h, r, phi)
around the pseudo-C2 axis. See README.md for the method, the configuration file and the outputs.

    python3 solaris.py config.json --check
    python3 solaris.py config.json
    python3 solaris.py config.json --figures-only

Author: Tomás Castro (tomas.castro@utalca.cl)
"""

import argparse
import json
import math
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as patheffects
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon


# ============================================================================ constants
HIS_ROLES = ["HisA1", "HisA2", "HisA3", "HisB1", "HisB2", "HisB3"]
HIS_RING = ["CG", "ND1", "CD2", "CE1", "NE2"]

# Family consensus of CuA, CuB and the six His rings (Å) in the SOLARIS frame:
# origin at the Cu-Cu midpoint, z along the pseudo-C2 axis, x towards CuA, +z towards the conserved Glu/Asn.
TEMPLATE_KEYS = [("CuA", "CU"), ("CuB", "CU")] + [(r, a) for r in HIS_ROLES for a in HIS_RING]
TEMPLATE = np.array([
    [1.713, 0.000, 0.012], [-1.713, 0.000, -0.012],
    [4.305, -3.264, 0.085], [3.829, -3.235, 1.348], [3.650, -2.305, -0.602],
    [2.919, -2.287, 1.409], [2.774, -1.677, 0.244],
    [3.079, 3.314, 2.267], [4.098, 2.447, 2.500], [2.164, 2.670, 1.505],
    [3.794, 1.324, 1.886], [2.629, 1.396, 1.248],
    [1.549, 1.473, -3.906], [2.098, 2.514, -3.236], [1.540, 0.425, -3.036],
    [2.407, 2.102, -2.011], [2.063, 0.838, -1.818],
    [-4.512, 3.001, 0.086], [-3.638, 3.402, 1.052], [-3.949, 1.941, -0.531],
    [-2.577, 2.597, 1.007], [-2.720, 1.679, 0.059],
    [-4.173, -1.550, 3.008], [-3.400, -2.618, 2.717], [-3.925, -0.634, 2.044],
    [-2.698, -2.337, 1.617], [-2.963, -1.122, 1.165],
    [-2.364, -1.369, -3.979], [-2.132, -0.368, -4.008], [-2.265, -1.458, -2.721],
    [-1.874, 0.122, -2.806], [-1.935, -0.509, -1.977],
])
TEMPLATE_ROW = {k: i for i, k in enumerate(TEMPLATE_KEYS)}

BACKBONE = {"N", "CA", "C", "O", "OXT"}
TIP_ATOMS = {
    "GLU": ["CD", "OE1", "OE2"], "ASP": ["CG", "OD1", "OD2"],
    "ASN": ["CG", "OD1", "ND2"], "GLN": ["CD", "OE1", "NE2"],
    "ARG": ["CZ", "NH1", "NH2"], "LYS": ["NZ"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ2", "CZ3", "CH2"],
    "VAL": ["CB", "CG1", "CG2"], "LEU": ["CG", "CD1", "CD2"], "ILE": ["CB", "CG1", "CG2", "CD1"],
    "THR": ["CB", "OG1", "CG2"], "SER": ["OG"], "CYS": ["SG"], "MET": ["SD", "CE"],
    "ALA": ["CB"], "GLY": ["CA"], "PRO": ["CG", "CD"],
}
EXPECTED_RESNAME = {"Glu": {"GLU"}, "Asn": {"ASN"}}
DEFAULT_LANDMARKS = {
    "Conserved_Glu": {"role": "Glu", "plot": "TIP"},
    "Conserved_Asn": {"role": "Asn", "plot": "TIP"},
    "Carboxylate_anchor": {"role": "anchor", "plot": "TIP"},
    "Gatekeeper": {"role": "gatekeeper", "plot": "SC"},
    "A/S_adj": {"role": "SerAla", "plot": "SC"},
}
LIGAND_DEFAULT = {"resname": "TY1", "atoms": {
    "N": "N1", "CA": "C1", "C": "C2", "CB": "C3", "CG": "C4", "CD1": "C5", "CD2": "C6",
    "CE1": "C7", "CE2": "C8", "CZ": "C9", "OH": "O2"}}
PEROXIDE_DEFAULT = {"resname": "OY1", "atoms": ["O1", "O2"]}
LIG_RING = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
LIG_REQUIRED = set(LIG_RING) | {"OH"}
SOLVENT = {"WAT", "HOH", "TIP3", "TIP4", "T3P", "T4P", "SOL", "SPC", "NA", "NA+", "CL", "CL-",
           "K", "K+", "MG", "CA2", "ZN2", "IP", "IM", "IB"}
RESNAME_ALIASES = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
                   "CYX": "CYS", "CYM": "CYS", "GLH": "GLU", "ASH": "ASP", "LYN": "LYS", "ARN": "ARG"}
KEY_DESCRIPTORS = ["Cu.d_CuA_CuB", "Cu.max_d_Cu_NE2", "Px.d_O_O", "Px.Cu2O2_tilt",
                   "L.d_OH_CuA", "L.d_OH_CuB", "L.d_Cortho_Operox", "L.ring_tilt",
                   "L.angle_CA_OH_vs_minus_z", "L.d_ring_HisB2", "L.d_carboxylate_anchor",
                   "L.d_ring_gatekeeper", "L.OH.h", "L.ring.h"]
META_COLS = ["system", "state", "kind", "category", "frame"]
QC_COLS = ["system", "state", "role", "resid_ref", "resid_used", "resname", "offset", "offset_from"]

# figure style (Okabe-Ito colours per enzyme; marker shape is the redundant encoding)
TYR_COLORS = {"ScTYR": "#CC79A7", "VvTYR": "#E69F00", "AbTYR": "#56B4E9", "BmTYR": "#009E73",
              "HsTYR": "#0072B2", "VsTYR": "#D55E00", "JrTYR": "#882255"}
TYR_MARKERS = {"ScTYR": "o", "VvTYR": "s", "AbTYR": "^", "BmTYR": "D",
               "HsTYR": "p", "VsTYR": "h", "JrTYR": "*"}
EXTRA_COLORS = ["#7f7f7f", "#b07aa1", "#9c755f", "#bab0ac"]
EXTRA_MARKERS = ["v", "<", ">", "P", "X", "8", "d"]
CU_COLOR = "#b5651d"
TRAJ_BLUES = ["#c6dbef", "#6baed6"]
INK, GRID = "#222222", "#e8e8e8"
HIS_FC, HIS_EC = "#e2e2e2", "#a9a9a9"
STATE_STYLE = {"global": dict(edge="black", lw=1.8, fill="black"),
               "local": dict(edge="white", lw=1.8, fill="white"),
               "docking": dict(edge="0.35", lw=1.8, fill="#e66101")}
STATE_LABEL = {"global": "global min.", "local": "local min.", "docking": "docking", "other": "other"}
STATE_RANK = {"global": 0, "local": 1, "docking": 2, "other": 3}
RES3 = {"ALA": "Ala", "ARG": "Arg", "ASN": "Asn", "ASP": "Asp", "CYS": "Cys", "GLN": "Gln", "GLU": "Glu",
        "GLY": "Gly", "HIS": "His", "ILE": "Ile", "LEU": "Leu", "LYS": "Lys", "MET": "Met", "PHE": "Phe",
        "PRO": "Pro", "SER": "Ser", "THR": "Thr", "TRP": "Trp", "TYR": "Tyr", "VAL": "Val"}
R_COLOR = (0.0, 12.0)
SITE_HALF = 5.5  # half-width (Å) of the top view needed to show the Cu and the six His rings
N_POINTS = 400
TAIL_PCT = 2.5  # a static structure is "in the tail" of its trajectory below this percentile or above 100 - it


# ============================================================================ geometry
def kabsch(P, Q):
    """R, t minimising |R p + t - q| over the rows of P and Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)) or 1.0
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def angle_deg(u, v):
    c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def axis_tilt(u, axis):
    a = angle_deg(u, axis)
    return min(a, 180.0 - a)


def plane_normal(P):
    return np.linalg.svd(P - P.mean(0))[2][2]


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def cylindrical(p):
    return float(np.hypot(p[0], p[1])), float(np.degrees(np.arctan2(p[1], p[0])) % 360.0)


def box_matrix(box):
    a, b, c, al, be, ga = [float(v) for v in box[:6]]
    al, be, ga = np.radians([al, be, ga])
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz = np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return np.array([[a, 0, 0], [b * np.cos(ga), b * np.sin(ga), 0], [cx, cy, cz]])


_SHIFTS = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float)


def minimum_image(vecs, box):
    M = box_matrix(box)
    f = vecs @ np.linalg.inv(M)
    f -= np.round(f)
    v = f @ M
    cand = v[:, None, :] + (_SHIFTS @ M)[None, :, :]
    best = np.argmin((cand ** 2).sum(-1), axis=1)
    return cand[np.arange(len(v)), best]


def usable_box(box, need=40.0):
    if box is None:
        return False
    b = np.asarray(box, float)
    return b.size >= 6 and np.all(np.isfinite(b[:6])) and np.all(b[:3] > need)


# ============================================================================ input
def resolve_path(base, value):
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() else (base / p).resolve()


def std_resname(rn):
    rn = str(rn).upper()
    if re.fullmatch(r"H[DEIP]\d{1,2}", rn):
        return "HIS"
    return RESNAME_ALIASES.get(rn, rn)


def guess_element(name, resname):
    n = re.sub(r"[^A-Za-z]", "", str(name)).upper()
    r = re.sub(r"[^A-Za-z]", "", str(resname)).upper()
    for el in ("CU", "ZN", "FE", "MG", "MN", "NA", "CL"):
        if n.startswith(el) and (r.startswith(el) or n == el):
            return el
    return n[:1] if n else ""


class Topology:
    def __init__(self, names, resnames, resids, chains, elements):
        self.names = np.array([str(x).strip().upper() for x in names], dtype=object)
        self.resnames = np.array([str(x).strip().upper() for x in resnames], dtype=object)
        self.resids = np.array(resids, dtype=int)
        self.chains = np.array([str(x).strip() for x in chains], dtype=object)
        els = []
        for n, r, e in zip(self.names, self.resnames, elements):
            e = str(e).strip().upper() if e is not None else ""
            els.append(e if e and e not in ("X", "DUMMY") else guess_element(n, r))
        self.elements = np.array(els, dtype=object)
        self.n = len(self.names)
        # a residue is a contiguous block with the same (chain, resid, resname)
        self.residues = []
        prev = None
        for i in range(self.n):
            key = (self.chains[i], int(self.resids[i]), self.resnames[i])
            if key != prev:
                self.residues.append({"chain": key[0], "resid": key[1], "resname": key[2], "atoms": {}})
                prev = key
            self.residues[-1]["atoms"].setdefault(self.names[i], i)
        self.by_resid = {}
        for k, r in enumerate(self.residues):
            self.by_resid.setdefault(r["resid"], []).append(k)

    def find_residue(self, resid, chain=None):
        ks = [k for k in self.by_resid.get(resid, []) if chain is None or self.residues[k]["chain"] == chain]
        if len(ks) > 1:
            ks = [k for k in ks if self.residues[k]["resname"] not in SOLVENT] or ks
        return ks

    def residues_named(self, resname):
        return [r for r in self.residues if r["resname"] == str(resname).upper()]


class PDBSource:
    """Single- or multi-MODEL PDB file."""

    def __init__(self, path):
        names, resn, resi, chains, els, cur, models = [], [], [], [], [], [], []
        box = None
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("CRYST1"):
                    try:
                        box = np.array([float(line[a:b]) for a, b in
                                        ((6, 15), (15, 24), (24, 33), (33, 40), (40, 47), (47, 54))])
                    except ValueError:
                        box = None
                elif line.startswith(("ATOM", "HETATM")):
                    if line[16:17] not in (" ", "A", "1"):
                        continue
                    if not models:
                        names.append(line[12:16])
                        resn.append(line[17:21])
                        chains.append(line[21:22])
                        resi.append(int(line[22:26]))
                        els.append(line[76:78] if len(line) >= 78 else "")
                    cur.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                elif line.startswith("ENDMDL"):
                    if cur:
                        models.append(np.array(cur, float))
                    cur = []
        if cur:
            models.append(np.array(cur, float))
        if not models:
            raise ValueError(f"{path}: no ATOM/HETATM records")
        bad = [k for k, m in enumerate(models) if len(m) != len(models[0])]
        if bad:
            raise ValueError(f"{path}: MODEL {bad[:5]} do not have {len(models[0])} atoms")
        self.topo = Topology(names, resn, resi, chains, els)
        self.models, self.box, self.n_frames = models, box, len(models)

    def frames(self, start=0, stop=None, stride=1):
        stop = self.n_frames if stop is None else min(int(stop), self.n_frames)
        for i in range(int(start), stop, max(1, int(stride))):
            yield i, self.models[i], self.box


class MDASource:
    """Topology + trajectories read with MDAnalysis."""

    def __init__(self, top, trajs):
        try:
            import MDAnalysis as mda
        except ImportError as exc:
            raise RuntimeError("trajectories need MDAnalysis (pip install MDAnalysis)") from exc
        self.u = mda.Universe(str(top), *[str(t) for t in trajs])
        ag = self.u.atoms
        try:
            chains = list(ag.chainIDs)
        except Exception:
            chains = [""] * len(ag)
        try:
            els = list(ag.elements)
        except Exception:
            els = [""] * len(ag)
        self.topo = Topology(list(ag.names), list(ag.resnames), list(ag.resids), chains, els)
        self.n_frames = len(self.u.trajectory)

    def frames(self, start=0, stop=None, stride=1):
        stop = self.n_frames if stop is None else min(int(stop), self.n_frames)
        for ts in self.u.trajectory[int(start):stop:max(1, int(stride))]:
            box = None if ts.dimensions is None else np.asarray(ts.dimensions, float)
            yield int(ts.frame), ts.positions.astype(float), box


def open_source(base, entry, is_traj):
    if is_traj:
        top = resolve_path(base, entry.get("top") or entry.get("parm7") or entry.get("prmtop"))
        trajs = entry.get("traj") or entry.get("trajectory")
        trajs = [resolve_path(base, t) for t in (trajs if isinstance(trajs, list) else [trajs])]
        missing = [p for p in [top] + trajs if not p.exists()]
        if missing:
            raise FileNotFoundError(", ".join(map(str, missing)))
        return MDASource(top, trajs)
    p = resolve_path(base, entry["pdb"])
    if not p.exists():
        raise FileNotFoundError(str(p))
    return PDBSource(p)


# ============================================================================ roles -> atoms
def parse_spec(spec):
    """42, "42" or "A:42" -> (chain or None, resid)."""
    if isinstance(spec, (int, np.integer)):
        return None, int(spec)
    s = str(spec).strip()
    if ":" in s:
        c, r = s.split(":", 1)
        return (c.strip() or None), int(r)
    return None, int(s)


def is_undefined(spec):
    return spec is None or (isinstance(spec, str) and (not spec.strip() or "TODO" in spec.upper()))


def numbering_offset(topo, roles):
    """Offset between the reference numbering and this file, from residues named HD1..HD6 (None if absent)."""
    hits = {}
    for k, role in enumerate(HIS_ROLES, start=1):
        ks = [i for i, r in enumerate(topo.residues) if r["resname"] == f"HD{k}"]
        if len(ks) > 1:
            raise ValueError(f"{len(ks)} residues named HD{k}: resid {[topo.residues[i]['resid'] for i in ks]}")
        if ks:
            hits[role] = topo.residues[ks[0]]["resid"]
    if not hits:
        return None
    if len(hits) < 6:
        raise ValueError(f"only {', '.join(sorted(hits))} found as HD1-HD6")
    offs = {role: hits[role] - parse_spec(roles[role])[1] for role in HIS_ROLES}
    if len(set(offs.values())) != 1:
        det = ", ".join(f"{r}: HD{k + 1}={hits[r]} vs {parse_spec(roles[r])[1]} ({offs[r]:+d})"
                        for k, r in enumerate(HIS_ROLES))
        raise ValueError(f"HD1-HD6 do not match a single numbering offset: {det}")
    return offs["HisA1"]


class SiteResolver:
    """Maps the roles of one enzyme onto the atoms of one structure or trajectory."""

    def __init__(self, topo, roles, entry):
        self.topo, self.roles, self.warnings, self.qc = topo, roles, [], []
        off = numbering_offset(topo, roles)
        if off is None:
            self.offset = int(entry.get("resid_offset", 0))
            self.offset_from = "resid_offset" if "resid_offset" in entry else "none (resid_offset 0)"
        else:
            if "resid_offset" in entry and int(entry["resid_offset"]) != off:
                raise ValueError(f"resid_offset={entry['resid_offset']} in the config, but HD1-HD6 give {off:+d}")
            self.offset, self.offset_from = off, "HD1-HD6"

    def residue(self, role):
        spec = self.roles.get(role)
        if is_undefined(spec):
            return None
        chain, rid = parse_spec(spec)
        rid_s = rid + self.offset
        ks = self.topo.find_residue(rid_s, chain)
        if not ks:
            raise ValueError(f"{role}: residue {rid_s} not found (reference {rid}, offset {self.offset:+d})")
        if len(ks) > 1:
            want = {"HIS"} if role in HIS_ROLES else EXPECTED_RESNAME.get(role)
            ok = [k for k in ks if (std_resname(self.topo.residues[k]["resname"]) in want if want
                                    else std_resname(self.topo.residues[k]["resname"]) in TIP_ATOMS)]
            where = [(self.topo.residues[k]["chain"], self.topo.residues[k]["resname"]) for k in ks]
            if len(ok) != 1:
                raise ValueError(f"{role}: resid {rid_s} appears {len(ks)} times {where}; write 'chain:resid'")
            self.warnings.append(f"{role}: resid {rid_s} appears {len(ks)} times {where}; using "
                                 f"{self.topo.residues[ok[0]]['resname']}{rid_s}")
            ks = ok
        r = self.topo.residues[ks[0]]
        self.qc.append((role, rid, r["resid"], r["resname"]))
        return r

    def histidines(self):
        his = {}
        for role in HIS_ROLES:
            r = self.residue(role)
            if r is None:
                raise ValueError(f"role {role} is missing")
            if std_resname(r["resname"]) != "HIS":
                raise ValueError(f"{role} = {r['resname']}{r['resid']} is not a His (numbering?)")
            if "NE2" not in r["atoms"]:
                raise ValueError(f"{role} = {r['resname']}{r['resid']} has no NE2")
            his[role] = {a: r["atoms"][a] for a in HIS_RING if a in r["atoms"]}
        return his

    def coppers(self):
        t = self.topo
        cu = sorted({i for i in range(t.n) if t.elements[i] == "CU" or re.fullmatch(r"CU\d?", t.names[i])})
        if len(cu) != 2:
            raise ValueError(f"expected 2 Cu, found {len(cu)}")
        return cu

    def landmarks(self, landmark_defs, tip_override, system):
        out = {}
        for label, ld in landmark_defs.items():
            role = ld["role"]
            if is_undefined(self.roles.get(role)):
                if not ld.get("optional", False):
                    self.warnings.append(f"{label}: role '{role}' not defined for {system}; skipped")
                continue
            r = self.residue(role)
            std = std_resname(r["resname"])
            exp = EXPECTED_RESNAME.get(role)
            if exp and std not in exp:
                raise ValueError(f"{label} ({role}) = {r['resname']}{r['resid']}: expected {'/'.join(sorted(exp))} "
                                 f"(numbering? offset {self.offset:+d} from {self.offset_from})")
            heavy = {a: i for a, i in r["atoms"].items() if self.topo.elements[i] != "H"}
            if "CA" not in heavy:
                raise ValueError(f"{label}: {r['resname']}{r['resid']} has no CA")
            sc = [i for a, i in heavy.items() if a not in BACKBONE] or [heavy["CA"]]
            want = tip_override.get(role, ld.get("tip"))
            want = want.split() if isinstance(want, str) else (want or TIP_ATOMS.get(std) or [])
            tip = [heavy[a] for a in want if a in heavy] or sc
            miss = [a for a in want if a not in heavy]
            if miss:
                self.warnings.append(f"{label}: {r['resname']}{r['resid']} lacks {' '.join(miss)}; "
                                     f"TIP from the atoms present (or the side chain)")
            out[label] = {"role": role, "resname": r["resname"], "resid": r["resid"],
                          "CA": [heavy["CA"]], "SC": sc, "TIP": tip, "plot": ld.get("plot", "TIP")}
        return out

    def ligands(self, lcfg):
        amap = {**LIGAND_DEFAULT["atoms"], **lcfg.get("atoms", {})}
        res = self.topo.residues_named(lcfg["resname"])
        if not res:
            self.warnings.append(f"no residue named {lcfg['resname']}: no L-Tyr descriptors")
        if len(res) > 1:
            self.warnings.append(f"{len(res)} residues {lcfg['resname']}: using the one closest to the Cu")
        cand = []
        for r in res:
            at = r["atoms"]
            got = {}
            for a, n in amap.items():
                for nm in (str(n).upper(), a):
                    if nm in at:
                        got[a] = at[nm]
                        break
            miss = [f"{a}->{n}" for a, n in amap.items() if a not in got]
            if any(m.split("->")[0] in LIG_REQUIRED for m in miss):
                self.warnings.append(f"ligand {lcfg['resname']}{r['resid']}: missing {', '.join(miss)}; "
                                     f"no L-Tyr descriptors (check 'ligand.atoms')")
                continue
            if miss:
                self.warnings.append(f"ligand {lcfg['resname']}{r['resid']}: missing {', '.join(miss)}; "
                                     f"descriptors that need them are skipped")
            cand.append(got)
            if len(cand) == 1:
                self.qc.append(("ligand", None, r["resid"], str(lcfg["resname"]).upper()))
        return cand

    def peroxide(self, pcfg):
        res = self.topo.residues_named(pcfg["resname"])
        if len(res) > 1:
            self.warnings.append(f"{len(res)} residues {pcfg['resname']}: no peroxide descriptors")
        if len(res) != 1:
            return None
        at = res[0]["atoms"]
        if not all(a.upper() in at for a in pcfg["atoms"]):
            return None
        return [at[a.upper()] for a in pcfg["atoms"]][:2]


def merged(key, default, *cfgs):
    out = dict(default)
    for c in cfgs:
        out.update(c.get(key, {}))
    return out


def resolve_site(system, state, topo, sys_cfg, cfg, entry, landmark_defs):
    res = SiteResolver(topo, sys_cfg["roles"], entry)
    site = SimpleNamespace(offset=res.offset, offset_from=res.offset_from)
    site.his = res.histidines()
    site.cu_pair = res.coppers()
    site.landmarks = res.landmarks(landmark_defs, sys_cfg.get("tip_atoms", {}), system)
    site.ligands = res.ligands(merged("ligand", LIGAND_DEFAULT, cfg, sys_cfg, entry))
    site.peroxide = res.peroxide(merged("peroxide", PEROXIDE_DEFAULT, cfg, sys_cfg, entry))
    site.warnings = res.warnings
    site.qc = [(system, state, role, ref, used, rn, res.offset, res.offset_from) for role, ref, used, rn in res.qc]
    return site


# ============================================================================ per-frame descriptors
class FrameAnalyzer:
    def __init__(self, system, state, kind, category, site, r_min_phi=2.5):
        self.meta = {"system": system, "state": state, "kind": kind, "category": category}
        self.site = site
        self.r_min_phi = r_min_phi
        self.ready = False

    def _prepare(self, X):
        s = self.site
        ne2 = {r: X[s.his[r]["NE2"]] for r in HIS_ROLES}
        c0, c1 = X[s.cu_pair[0]], X[s.cu_pair[1]]
        dist = lambda c, roles: sum(np.linalg.norm(ne2[r] - c) for r in roles)
        direct = dist(c0, HIS_ROLES[:3]) + dist(c1, HIS_ROLES[3:])
        swap = dist(c1, HIS_ROLES[:3]) + dist(c0, HIS_ROLES[3:])
        cuA, cuB = (s.cu_pair[0], s.cu_pair[1]) if direct <= swap else (s.cu_pair[1], s.cu_pair[0])
        center = 0.5 * (X[cuA] + X[cuB])
        lig = min(s.ligands, key=lambda g: np.linalg.norm(X[g["OH"]] - center)) if s.ligands else None

        need = {cuA, cuB}
        for r in HIS_ROLES:
            need.update(s.his[r].values())
        for lm in s.landmarks.values():
            need.update(lm["CA"] + lm["SC"] + lm["TIP"])
        if lig:
            need.update(lig.values())
        if s.peroxide:
            need.update(s.peroxide)
        self.atoms = np.array(sorted(need), int)
        pos = {g: k for k, g in enumerate(self.atoms)}
        loc = lambda idx: np.array([pos[i] for i in idx], int)

        self.i_cuA, self.i_cuB = pos[cuA], pos[cuB]
        self.i_ne2 = {r: pos[s.his[r]["NE2"]] for r in HIS_ROLES}
        self.i_hisB2_ring = loc([s.his["HisB2"][a] for a in HIS_RING if a in s.his["HisB2"]])
        rows, trows = [self.i_cuA, self.i_cuB], [0, 1]
        for r in HIS_ROLES:
            for a in HIS_RING:
                if a in s.his[r]:
                    rows.append(pos[s.his[r][a]])
                    trows.append(TEMPLATE_ROW[(r, a)])
        self.i_core, self.t_core = np.array(rows, int), TEMPLATE[np.array(trows, int)]
        self.i_lm = {lab: {p: loc(lm[p]) for p in ("CA", "SC", "TIP")} for lab, lm in s.landmarks.items()}
        self.i_lig = {a: pos[i] for a, i in lig.items()} if lig else None
        self.i_px = loc(s.peroxide) if s.peroxide else None
        self.anchor = next((lab for lab, lm in s.landmarks.items() if lm["role"] == "anchor"), None)
        self.gatekeeper = next((lab for lab, lm in s.landmarks.items() if lm["role"] == "gatekeeper"), None)
        self.ready = True

    def run(self, source, start=0, stop=None, stride=1):
        rows, long_rows = [], []
        for fi, X, box in source.frames(start, stop, stride):
            if not self.ready:
                self._prepare(X)
            P = X[self.atoms].astype(float)
            if usable_box(box):
                ref = P[self.i_cuA].copy()
                P = ref + minimum_image(P - ref, box)
            row, lr = self.frame(P)
            tag = {**self.meta, "frame": int(fi)}
            rows.append({**row, **tag})
            long_rows += [{**d, **tag} for d in lr]
        return rows, long_rows

    def frame(self, P):
        R, t = kabsch(P[self.i_core], self.t_core)
        L = P @ R.T + t
        out = {"Frame.core_RMSD": float(np.sqrt(((L[self.i_core] - self.t_core) ** 2).sum(1).mean()))}
        cuA, cuB = L[self.i_cuA], L[self.i_cuB]
        out.update(self._core(L, cuA, cuB))
        lm_out, long_rows, centroids = self._landmarks(L, cuA, cuB)
        out.update(lm_out)
        if self.i_lig is not None:
            lig_out, lig_rows = self._ligand(L, cuA, cuB, centroids)
            out.update(lig_out)
            long_rows += lig_rows
        return out, long_rows

    def _core(self, L, cuA, cuB):
        out = {"Cu.d_CuA_CuB": float(np.linalg.norm(cuA - cuB)),
               "Cu.max_d_Cu_NE2": max(float(np.linalg.norm(L[self.i_ne2[r]] - (cuA if r.startswith("HisA") else cuB)))
                                      for r in HIS_ROLES)}
        if self.i_px is not None:
            o1, o2 = L[self.i_px[0]], L[self.i_px[1]]
            out["Px.d_O_O"] = float(np.linalg.norm(o1 - o2))
            out["Px.Cu2O2_tilt"] = axis_tilt(plane_normal(np.vstack([cuA, cuB, o1, o2])), np.array([0.0, 0.0, 1.0]))
            out["Px.dihedral_CuA_O1_O2_CuB"] = abs(dihedral(cuA, o1, o2, cuB))
        return out

    def _point_row(self, label, point, plot, resname, resid, c, cuA, cuB):
        r, phi = cylindrical(c)
        if r < self.r_min_phi:
            phi = float("nan")
        return {"landmark": label, "point": point, "plot": plot, "resname": resname, "resid": resid,
                "x": float(c[0]), "y": float(c[1]), "h": float(c[2]), "r": r, "phi": phi,
                "dCuA": float(np.linalg.norm(c - cuA)), "dCuB": float(np.linalg.norm(c - cuB))}

    def _landmarks(self, L, cuA, cuB):
        out, rows, centroids = {}, [], {}
        for lab, pts in self.i_lm.items():
            lm = self.site.landmarks[lab]
            c = {p: L[pts[p]].mean(0) for p in ("CA", "SC", "TIP")}
            for p, v in c.items():
                row = self._point_row(lab, p, p == lm["plot"], lm["resname"], lm["resid"], v, cuA, cuB)
                rows.append(row)
                for k in ("x", "y", "h", "r", "phi"):
                    out[f"{lab}.{p}.{k}"] = row[k]
            out[f"{lab}.TIP.dCuA"] = float(np.linalg.norm(c["TIP"] - cuA))
            out[f"{lab}.TIP.dCuB"] = float(np.linalg.norm(c["TIP"] - cuB))
            out[f"{lab}.dr"] = out[f"{lab}.TIP.r"] - out[f"{lab}.CA.r"]
            out[f"{lab}.dh"] = out[f"{lab}.TIP.h"] - out[f"{lab}.CA.h"]
            centroids[lab] = c
        return out, rows, centroids

    def _ligand(self, L, cuA, cuB, centroids):
        z = np.array([0.0, 0.0, 1.0])
        g = {a: L[k] for a, k in self.i_lig.items()}
        ring = np.vstack([g[a] for a in LIG_RING])
        rc, nrm = ring.mean(0), plane_normal(ring)
        out, rows = {}, []
        for nm, p in [("OH", g["OH"]), ("ring", rc)] + ([("CA", g["CA"])] if "CA" in g else []):
            row = self._point_row(f"L-Tyr {nm}", nm, True, "LIG", -1, p, cuA, cuB)
            rows.append(row)
            for k in ("x", "y", "h", "r", "phi"):
                out[f"L.{nm}.{k}"] = row[k]
        out["L.d_OH_CuA"] = float(np.linalg.norm(g["OH"] - cuA))
        out["L.d_OH_CuB"] = float(np.linalg.norm(g["OH"] - cuB))
        out["L.ring_tilt"] = axis_tilt(nrm, z)
        out["L.ring_normal_vs_CuCu"] = axis_tilt(nrm, np.array([1.0, 0.0, 0.0]))
        out["L.angle_COH_vs_minus_z"] = angle_deg(g["OH"] - g["CZ"], -z)
        if "CA" in g:
            out["L.angle_CA_OH_vs_minus_z"] = angle_deg(g["OH"] - g["CA"], -z)
        if self.i_px is not None:
            out["L.d_Cortho_Operox"] = float(min(np.linalg.norm(g[c] - L[k]) for c in ("CE1", "CE2") for k in self.i_px))
        out["L.d_ring_HisB2"] = float(np.linalg.norm(rc - L[self.i_hisB2_ring].mean(0)))
        if self.anchor and "C" in g:
            out["L.d_carboxylate_anchor"] = float(np.linalg.norm(g["C"] - centroids[self.anchor]["TIP"]))
        if self.gatekeeper:
            out["L.d_ring_gatekeeper"] = float(np.linalg.norm(rc - centroids[self.gatekeeper]["SC"]))
        return out, rows


# ============================================================================ statistics
def is_angle(col):
    return col.endswith(".phi")


def circ_stats(deg):
    a = np.radians(np.asarray(deg, float))
    C, S = np.cos(a).mean(), np.sin(a).mean()
    R = min(max(np.hypot(C, S), 1e-12), 1.0)
    return float(np.degrees(np.arctan2(S, C)) % 360.0), float(np.degrees(np.sqrt(-2.0 * np.log(R))))


def summarize(df, desc_cols):
    out = []
    for (sy, st), g in df.groupby(["system", "state"], sort=False):
        for c in desc_cols:
            v = g[c].dropna().to_numpy(float)
            if not len(v):
                continue
            row = {"system": sy, "state": st, "category": g["category"].iloc[0], "kind": g["kind"].iloc[0],
                   "descriptor": c, "n": len(v)}
            if is_angle(c):
                mu, sd = circ_stats(v)
                row.update(mean=mu, sd=sd if len(v) > 1 else np.nan, p10=np.nan, median=np.nan, p90=np.nan)
            else:
                row.update(mean=v.mean(), sd=v.std(ddof=1) if len(v) > 1 else np.nan,
                           p10=np.percentile(v, 10), median=np.median(v), p90=np.percentile(v, 90))
            out.append(row)
    return pd.DataFrame(out)


def key_columns(df, landmark_defs):
    cols = [c for c in KEY_DESCRIPTORS if c in df]
    for lab, ld in landmark_defs.items():
        p = ld.get("plot", "TIP")
        cols += [c for c in (f"{lab}.{p}.h", f"{lab}.{p}.r", f"{lab}.dr") if c in df]
    return cols


def percentile_in_trajectory(df, desc_cols, key_cols=()):
    """Percentile of each static structure within the trajectories of the same enzyme."""
    key_cols, out = set(key_cols), []
    for sy, g in df.groupby("system", sort=False):
        traj = g[g["kind"] == "trajectory"]
        if traj.empty:
            continue
        for st, gs in g[g["kind"] == "structure"].groupby("state", sort=False):
            for tst, gt in traj.groupby("state", sort=False):
                for c in desc_cols:
                    if is_angle(c) or c not in gs or gs[c].isna().all():
                        continue
                    ref = gt[c].dropna().to_numpy(float)
                    if len(ref) < 20:
                        continue
                    v = float(gs[c].iloc[0])
                    pct = 100.0 * (np.mean(ref < v) + 0.5 * np.mean(ref == v))
                    out.append({"system": sy, "state": st, "trajectory": tst, "descriptor": c, "value": v,
                                "percentile": pct, "key": c in key_cols,
                                "in_tail": bool(pct < TAIL_PCT or pct > 100 - TAIL_PCT)})
    return pd.DataFrame(out)


def core_checks(df):
    out = []
    for (sy, st), g in df.groupby(["system", "state"], sort=False):
        dn, rm = g["Cu.max_d_Cu_NE2"].to_numpy(float), g["Frame.core_RMSD"].to_numpy(float)
        if len(g) == 1:
            if dn[0] > 2.6:
                out.append(f"{sy}/{st}: a His is {dn[0]:.2f} Å from its Cu (> 2.6 Å): lost coordination or wrong His?")
            if rm[0] > 1.0:
                out.append(f"{sy}/{st}: core RMSD to the template {rm[0]:.2f} Å (> 1.0 Å): distorted site")
        else:
            f1, f2 = float(np.mean(dn > 2.6)), float(np.mean(rm > 1.0))
            if f1 > 0.01:
                out.append(f"{sy}/{st}: {100 * f1:.1f} % of frames with a His > 2.6 Å from its Cu")
            if f2 > 0.01:
                out.append(f"{sy}/{st}: {100 * f2:.1f} % of frames with core RMSD > 1.0 Å")
    return out


# ============================================================================ figure helpers
def _pick(system, systems, fixed, extra, default):
    if system in fixed:
        return fixed[system]
    others = [s for s in systems if s not in fixed]
    return extra[others.index(system) % len(extra)] if system in others else default


def marker_for(system, systems):
    return _pick(system, systems, TYR_MARKERS, EXTRA_MARKERS, "o")


def color_for(system, systems=()):
    return _pick(system, systems, TYR_COLORS, EXTRA_COLORS, "#7f7f7f")


def lighten(col, f):
    r, g, b = matplotlib.colors.to_rgb(col)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def darker(col, f=0.78):
    r, g, b = matplotlib.colors.to_rgb(col)
    return (r * f, g * f, b * f)


def halo():
    return [patheffects.withStroke(linewidth=2.2, foreground="white")]


def state_style(category):
    return STATE_STYLE.get(category, dict(edge="0.6", lw=1.2, fill="0.8"))


def res_label(resname, resid):
    rn = str(resname).upper()
    if rn in ("LIG", "TY1"):
        return "L-Tyr"
    return f"{RES3.get(std_resname(rn), rn.title())}{int(resid)}"


def display_name(label):
    return str(label).replace("_", " ")


def file_name(label):
    return re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")


def panels_for(long_df, landmark_order):
    present = list(dict.fromkeys(long_df["landmark"]))
    return [l for l in landmark_order if l in present] + [l for l in ("L-Tyr OH", "L-Tyr ring") if l in present]


def thin_frames(df, n=N_POINTS):
    """Regular subsample of about n frames per enzyme and state."""
    if df.empty or not n or n <= 0:
        return df
    keep = []
    for _, g in df.groupby(["system", "state"], sort=False):
        fr = np.unique(g["frame"].to_numpy())
        keep.append(g[g["frame"].isin(fr[::max(1, math.ceil(len(fr) / n))])])
    return pd.concat(keep)


def frame_step(df, n=N_POINTS):
    nf = df.groupby(["system", "state"])["frame"].nunique().max() if not df.empty else 0
    return max(1, math.ceil(nf / n)) if n and nf else 1


def auto_seam(phis):
    """Cut of the unrolled cylinder in the widest empty gap of phi."""
    p = np.sort(np.asarray(phis, float) % 360.0)
    if len(p) < 2:
        return 0.0
    gaps = np.diff(np.r_[p, p[0] + 360.0])
    k = int(np.argmax(gaps))
    return float((p[k] + gaps[k] / 2.0) % 360.0)


def unwrap(phi, seam):
    return (np.asarray(phi, float) - seam) % 360.0 + seam


def style_plain(ax):
    ax.grid(alpha=0.15)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def style_clean(ax):
    ax.tick_params(labelsize=7, colors=INK, length=2.5, width=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.6)
        ax.spines[s].set_color("#9a9a9a")


def his_ring_xy(role):
    return np.array([TEMPLATE[TEMPLATE_ROW[(role, a)]] for a in ("CG", "ND1", "CE1", "NE2", "CD2")])


def draw_circles(ax, radii=(5, 10)):
    for rr in radii:
        ax.add_patch(plt.Circle((0, 0), rr, fill=False, lw=0.7, color="#cfcfcf", zorder=0.5))


def draw_floor_top(ax, light=False):
    """Template His rings, Cu-NE2 bonds and Cu atoms seen from the top."""
    fc, ec = ("#f1f1f1", "#c4c4c4") if light else (HIS_FC, HIS_EC)
    xl, yl = ax.get_xlim(), ax.get_ylim()
    mx, my = 0.04 * (xl[1] - xl[0]), 0.04 * (yl[1] - yl[0])
    for k, roles in ((0, HIS_ROLES[:3]), (1, HIS_ROLES[3:])):
        x, y, _ = TEMPLATE[k]
        for role in roles:
            ne = TEMPLATE[TEMPLATE_ROW[(role, "NE2")]]
            ax.plot([x, ne[0]], [y, ne[1]], color=ec, lw=0.6, zorder=1)
    for role in HIS_ROLES:
        ax.add_patch(Polygon(his_ring_xy(role)[:, :2], closed=True, fc=fc, ec=ec, lw=0.6, zorder=1.2))
    for nm, k in (("CuA", 0), ("CuB", 1)):
        x, y, _ = TEMPLATE[k]
        ax.scatter([x], [y], s=55, color=CU_COLOR, edgecolors="white", linewidths=0.8, zorder=2)
        if xl[0] + mx < x < xl[1] - mx and yl[0] + my < y - 0.7 < yl[1] - my:
            ax.text(x, y - 0.7, nm, fontsize=6.5, color=darker(CU_COLOR), ha="center", va="top", zorder=2.1,
                    path_effects=halo())


def draw_reference_top(ax, lim):
    for rr in (5, 10):
        if rr < lim * 1.45:
            ax.add_patch(plt.Circle((0, 0), rr, fill=False, lw=0.7, color="#cfcfcf", zorder=1))
    for r in HIS_ROLES:
        x, y, _ = TEMPLATE[TEMPLATE_ROW[(r, "NE2")]]
        ax.plot(x, y, marker="o", ms=3, color="0.6", zorder=2)
    for nm, row in (("CuA", 0), ("CuB", 1)):
        x, y, _ = TEMPLATE[row]
        ax.plot(x, y, marker="o", ms=7, color=CU_COLOR, mec="black", mew=0.6, zorder=2)
        ax.text(x, y - 1.0, nm, fontsize=6.5, color=CU_COLOR, ha="center", va="top", zorder=2)


def draw_states(ax, sub, xcol, ycol, systems, size=1.0, links=True):
    """Global min. large filled, local min. small filled (linked to the global), docking hollow."""
    marks = []
    for sy, g in sub.groupby("system", sort=False):
        col, mk = color_for(sy, systems), marker_for(sy, systems)
        glob = g[g["category"] == "global"]
        if links and len(glob):
            gx, gy = float(glob[xcol].iloc[0]), float(glob[ycol].iloc[0])
            for _, r in g[g["category"] != "global"].iterrows():
                ax.plot([gx, r[xcol]], [gy, r[ycol]], color=col, lw=0.8, alpha=0.6, zorder=3, solid_capstyle="round")
        for _, r in g.iterrows():
            cat = r["category"]
            if cat == "global":
                kw = dict(s=100 * size, color=col, edgecolors="white", linewidths=0.9, zorder=5)
            elif cat == "docking":
                kw = dict(s=62 * size, facecolors="white", edgecolors=darker(col, 0.85), linewidths=1.5, zorder=4.5)
            else:
                kw = dict(s=34 * size, color=col, edgecolors="white", linewidths=0.6, zorder=4)
            ax.scatter([r[xcol]], [r[ycol]], marker=mk, **kw)
            marks.append((float(r[xcol]), float(r[ycol])))
    return marks


def draw_frames(ax, g, xcol, ycol, system, systems, size=7, alpha=0.35, jitter=0.0, y0=None, rng=None):
    xs = g[xcol].to_numpy(float)
    ys = np.full(len(g), y0, float) if y0 is not None else g[ycol].to_numpy(float)
    if jitter:
        ys = ys + (rng or np.random.default_rng(0)).uniform(-jitter, jitter, len(ys))
    ax.scatter(xs, ys, s=size, marker=marker_for(system, systems), color=color_for(system, systems),
               alpha=alpha, lw=0, zorder=3, rasterized=True)


def smooth2d(H, s):
    if s <= 0:
        return H
    r = int(3 * s) + 1
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / s) ** 2)
    k /= k.sum()
    H = np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"), 0, H)
    return np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"), 1, H)


def density_fill(ax, x, y, extent, nb=70, sigma=1.4, color=None, zorder=1):
    """Regions holding 50 % and 90 % of the frames."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    fills = TRAJ_BLUES if color is None else [lighten(color, 0.72), lighten(color, 0.35)]
    lines = ["#4292c6", "#08519c"] if color is None else [lighten(color, 0.2), darker(color, 0.75)]
    if len(x) >= 30:
        H, xe, ye = np.histogram2d(x, y, bins=nb, range=extent)
        H = smooth2d(H, sigma).T
        v = np.sort(H.ravel())[::-1]
        c = np.cumsum(v)
        if c[-1] <= 0:
            return
        c /= c[-1]
        l50, l90 = v[min(np.searchsorted(c, 0.5), len(v) - 1)], v[min(np.searchsorted(c, 0.9), len(v) - 1)]
        if l90 < l50 < H.max():
            xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
            ax.contourf(xc, yc, H, levels=[l90, l50, H.max() * 1.001], colors=fills, alpha=0.9, zorder=zorder)
            ax.contour(xc, yc, H, levels=[l90, l50], colors=lines, linewidths=0.6, zorder=zorder + 0.1)
            return
    ax.scatter(x, y, s=6, color=fills[1], alpha=0.5, lw=0, zorder=zorder + 1)


def state_handles():
    return [Line2D([0], [0], marker="o", ls="None", ms=8.5, mfc="#7f7f7f", mec="white", mew=0.8, label="global min."),
            Line2D([0], [0], marker="o", ls="-", lw=0.8, color="#9a9a9a", ms=5, mfc="#7f7f7f", mec="white",
                   label="local min."),
            Line2D([0], [0], marker="o", ls="None", ms=7, mfc="white", mec="#555555", mew=1.4, label="docking")]


def traj_handles(step, density=False, frames=True):
    h = [Line2D([0], [0], marker="o", ls="None", ms=3.5, mfc="#7f7f7f", mec="none", alpha=0.6,
                label=f"GaMD frame (1 every {step})" if step > 1 else "GaMD frame")] if frames else []
    h += [Line2D([0], [0], marker="o", ls="None", ms=8.5, mfc="#7f7f7f", mec="white", mew=0.8, label="median")]
    if density:
        h += [Patch(fc=lighten("#7f7f7f", 0.35), ec="none", label="50 % of frames"),
              Patch(fc=lighten("#7f7f7f", 0.72), ec="none", label="90 % of frames")]
    return h


def enzyme_handles(df, systems):
    out = []
    for sy in [s for s in systems if s in set(df["system"])]:
        g = df[df["system"] == sy]
        out.append(Line2D([0], [0], marker=marker_for(sy, systems), ls="None", ms=7.5, mfc=color_for(sy, systems),
                          mec="white", mew=0.6, label=f"{sy} ({res_label(g['resname'].iloc[0], g['resid'].iloc[0])})"))
    return out


def family_legend(fig, systems, cats, extra=()):
    h = [Line2D([0], [0], marker=marker_for(s, systems), ls="None", ms=8, mfc="white", mec="black", label=s)
         for s in systems]
    h += [Line2D([0], [0], marker="o", ls="None", ms=8, mfc="0.75", mec=state_style(c)["edge"], mew=2.0,
                 label=STATE_LABEL[c]) for c in ("global", "local", "docking", "other") if c in cats]
    fig.legend(handles=h + list(extra), loc="upper left", bbox_to_anchor=(0.93, 0.92), frameon=True, fontsize=8,
               title="TYR / state", title_fontsize=9)


# ============================================================================ figures
def fig_landmark_detail(long_df, landmark, systems, path, dpi, traj=False, n_points=N_POINTS):
    """One landmark seen down the pseudo-C2 axis, framed on the dicopper site and on the data."""
    pts = long_df[(long_df["landmark"] == landmark) & long_df["plot"]].dropna(subset=["x", "y"]).copy()
    if pts.empty:
        return
    here = [s for s in systems if s in set(pts["system"])]
    xy = pts[["x", "y"]].abs().to_numpy(float)
    lim = max(SITE_HALF, (float(np.nanpercentile(xy, 99)) if traj else float(np.nanmax(xy))) + 1.2)

    fig, ax = plt.subplots(figsize=(6.6, 7.4))
    fig.subplots_adjust(left=0.11, right=0.97, top=0.92, bottom=0.19)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    draw_circles(ax)
    draw_floor_top(ax)
    if traj:
        for sy in here:
            g = pts[pts["system"] == sy]
            col, mk = color_for(sy, systems), marker_for(sy, systems)
            density_fill(ax, g["x"], g["y"], [[-lim, lim], [-lim, lim]], nb=90, color=col, zorder=3)
            ax.scatter([g["x"].median()], [g["y"].median()], s=60, marker=mk, color=col, edgecolors="white",
                       linewidths=0.8, zorder=5)
    else:
        pts["_rank"] = pts["category"].map(STATE_RANK).fillna(3)
        draw_states(ax, pts.sort_values(["system", "_rank"]), "x", "y", systems)
    ax.set_aspect("equal")
    ax.set_xlabel("x (Å)  → CuA", fontsize=8, color=INK)
    ax.set_ylabel("y (Å)", fontsize=8, color=INK)
    ax.set_title(display_name(landmark) + (" · GaMD" if traj else ""), loc="left", fontsize=12, fontweight="bold",
                 color=INK)
    style_clean(ax)
    his = Patch(fc=HIS_FC, ec=HIS_EC, lw=0.6, label="His")
    second = (traj_handles(frame_step(pts, n_points), density=True, frames=False) if traj else state_handles()) + [his]
    enz = enzyme_handles(pts, systems)
    kw = dict(loc="lower center", fontsize=7.5, frameon=False, handletextpad=0.3, columnspacing=1.2)
    fig.legend(handles=enz, bbox_to_anchor=(0.54, 0.075), ncol=min(4, len(enz)), **kw)
    fig.legend(handles=second, bbox_to_anchor=(0.54, 0.01), ncol=len(second), **kw)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def fig_landmark_summary(long_df, landmark_order, path, dpi, r_min, seam=None, traj=False, n_points=N_POINTS):
    """r, h and phi of every landmark: one row per landmark, one sub-row per enzyme."""
    d = long_df[long_df["plot"] & (long_df["landmark"] != "L-Tyr CA")].copy()
    if d.empty:
        return
    systems = list(dict.fromkeys(d["system"]))
    if seam is None:
        seam = auto_seam(d.loc[d["r"] >= r_min, "phi"])
    d["phi_u"] = unwrap(d["phi"], seam)
    rows, ylab, groups, y = [], [], [], 0.0
    for lab in panels_for(d, landmark_order):
        sub, top = d[d["landmark"] == lab], y
        for sy in [s for s in systems if s in set(sub["system"])]:
            g = sub[sub["system"] == sy]
            rows.append((sy, y, g))
            rl = res_label(g["resname"].iloc[0], g["resid"].iloc[0])
            ylab.append((y, "" if rl == "L-Tyr" else rl))
            y += 1.0
        groups.append((lab, top, y - 1.0))
        y += 0.9
    H = 0.75 + 0.165 * y
    fig, axes = plt.subplots(1, 3, figsize=(11.5, H), sharey=True,
                             gridspec_kw=dict(width_ratios=[1, 1.1, 1.5], wspace=0.08))
    cols = [("r", "r (Å)", (0, None)), ("h", "h (Å)", (None, None)), ("phi_u", "φ (°)", (seam, seam + 360))]
    rng = np.random.default_rng(0)
    for ax, (col, xl, lim) in zip(axes, cols):
        for k, (_, top, bot) in enumerate(groups):
            if k % 2 == 0:
                ax.axhspan(top - 0.55, bot + 0.55, color="#f6f6f6", zorder=0, lw=0)
        for sy, yy, g in rows:
            if col == "phi_u":
                g = g[g["r"] >= r_min]
            if g.empty:
                continue
            c_, mk = color_for(sy, systems), marker_for(sy, systems)
            vals = g[col].to_numpy(float)
            if traj:
                lo, med, hi = np.nanpercentile(vals, [5, 50, 95])
                ax.plot([lo, hi], [yy, yy], color=darker(c_, 0.85), lw=1.1, zorder=4, solid_capstyle="round")
                draw_frames(ax, thin_frames(g, n_points), col, None, sy, systems, size=6, alpha=0.3,
                            jitter=0.3, y0=yy, rng=rng)
                ax.scatter([med], [yy], s=48, marker=mk, color=c_, edgecolors="white", linewidths=0.7, zorder=5)
                continue
            if len(vals) > 1:
                ax.plot([np.nanmin(vals), np.nanmax(vals)], [yy, yy], color=c_, lw=1.0, alpha=0.6, zorder=2,
                        solid_capstyle="round")
            for _, r in g.iterrows():
                if r["category"] == "global":
                    ax.scatter([r[col]], [yy], s=48, marker=mk, color=c_, edgecolors="white", linewidths=0.7, zorder=5)
                elif r["category"] == "docking":
                    ax.scatter([r[col]], [yy], s=34, marker=mk, facecolors="white", edgecolors=darker(c_, 0.85),
                               linewidths=1.2, zorder=4)
                else:
                    ax.scatter([r[col]], [yy], s=16, marker=mk, color=c_, edgecolors="white", linewidths=0.4, zorder=4)
        if col == "phi_u":
            for ang in (0.0, 180.0):
                ax.axvline(float(unwrap(ang, seam)), color=CU_COLOR, lw=0.7, alpha=0.6, zorder=1)
            ticks = np.arange(math.ceil(seam / 60.0) * 60, seam + 360.1, 60)
            ax.set_xticks(ticks)
            ax.set_xticklabels([{0: "0 (CuA)", 180: "180 (CuB)"}.get(int(round(t % 360)), f"{t % 360:.0f}")
                                for t in ticks])
        if col == "h":
            ax.axvline(0, color=CU_COLOR, lw=0.7, alpha=0.6, zorder=1)
        ax.set_xlim(*lim)
        ax.set_xlabel(xl, fontsize=8, color=INK)
        ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
        ax.xaxis.set_ticks_position("top")
        ax.xaxis.set_label_position("top")
        style_clean(ax)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(ax is axes[0])
        ax.tick_params(axis="y", length=0)
    axes[0].set_ylim(y - 0.4, -1.3)
    axes[0].set_yticks([t for t, _ in ylab])
    axes[0].set_yticklabels([s for _, s in ylab], fontsize=6.5, color=INK)
    for lab, top, bot in groups:
        axes[0].annotate(display_name(lab), xy=(0, (top + bot) / 2), xycoords=("axes fraction", "data"),
                         xytext=(-44, 0), textcoords="offset points", ha="right", va="center", fontsize=8,
                         fontweight="bold", color=INK)
    hs = [Line2D([0], [0], marker=marker_for(s, systems), ls="None", ms=6.5, mfc=color_for(s, systems), mec="white",
                 label=s) for s in systems]
    hs += (traj_handles(frame_step(d, n_points)) + [Line2D([0], [0], color="#555555", lw=1.1, label="p5–p95")]
           if traj else state_handles())
    fig.legend(handles=hs, loc="lower center", ncol=len(hs), bbox_to_anchor=(0.595, 0.0), fontsize=7,
               frameon=False, handletextpad=0.3, columnspacing=1.0)
    fig.suptitle("Landmark positions" + (" · GaMD" if traj else ""), fontsize=11, fontweight="bold",
                 color=INK, x=0.02, ha="left", y=0.995)
    fig.subplots_adjust(left=0.2, right=0.99, top=1 - 0.95 / H, bottom=0.45 / H)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def fig_unrolled_cylinder(long_df, landmark_order, path, dpi, r_min, traj=False, n_points=N_POINTS, r_color=R_COLOR,
                          seam=None):
    """phi vs h for every pocket landmark; marker = enzyme, colour = r."""
    d = long_df[long_df["plot"]]
    panels = [p for p in panels_for(d, landmark_order) if not p.startswith("L-Tyr")]
    if not panels:
        return None
    systems = list(dict.fromkeys(d["system"]))
    dd = d[d["landmark"].isin(panels)]
    if traj:
        dd = thin_frames(dd, n_points)
    if seam is None:
        seam = auto_seam(dd.loc[dd["r"] >= r_min, "phi"].dropna())
    norm = Normalize(*r_color)
    ncol = min(3, len(panels))
    nrow = math.ceil(len(panels) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.4 * ncol, 4.2 * nrow), squeeze=False)
    for ax, lab in zip(axes.flat, panels):
        for (sy, _), g in dd[dd["landmark"] == lab].groupby(["system", "state"], sort=False):
            if traj:
                g = g[g["r"] >= r_min]
                ax.scatter(unwrap(g["phi"], seam), g["h"], c=g["r"], cmap="viridis", norm=norm, s=16,
                           marker=marker_for(sy, systems), alpha=0.7, linewidths=0, zorder=3, rasterized=True)
                continue
            g = g[g["r"] >= r_min]
            if g.empty:
                continue
            sty = state_style(g["category"].iloc[0])
            ax.scatter(unwrap(g["phi"], seam), g["h"], c=g["r"], cmap="viridis", norm=norm, s=85,
                       marker=marker_for(sy, systems), edgecolors=sty["edge"], linewidths=sty["lw"], zorder=4)
        ax.axhline(0, lw=0.6, color="0.7")
        for nm, ang in (("CuA", 0.0), ("CuB", 180.0)):
            xv = float(unwrap(ang, seam))
            ax.axvline(xv, ls="--", lw=0.8, color=CU_COLOR, alpha=0.7)
            ax.text(xv + 3, 0.02, nm, transform=ax.get_xaxis_transform(), fontsize=6.5, color=CU_COLOR,
                    ha="left", va="bottom")
        ticks = np.arange(math.ceil(seam / 60.0) * 60, seam + 360.1, 60)
        ax.set_xlim(seam, seam + 360)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t % 360:.0f}" for t in ticks], fontsize=7)
        ax.set_xlabel("φ (°)", fontsize=8)
        ax.set_title(display_name(lab), fontsize=11, fontweight="bold")
        ax.set_ylabel("h (Å)", fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
        style_plain(ax)
    for ax in axes.flat[len(panels):]:
        ax.axis("off")
    fig.subplots_adjust(left=0.06, right=0.86, bottom=0.09, top=0.9, wspace=0.3, hspace=0.5)
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap="viridis"), cax=fig.add_axes([0.885, 0.25, 0.012, 0.5]))
    cb.set_label("r (Å)", fontsize=9)
    if traj:
        step = frame_step(d[d["landmark"].isin(panels)], n_points)
        extra = [Line2D([0], [0], marker="o", ls="None", ms=4, mfc="#7f7f7f", mec="none", alpha=0.7,
                        label=f"GaMD frame (1 every {step})" if step > 1 else "GaMD frame")]
        family_legend(fig, systems, set(), extra)
        fig.suptitle("Unrolled cylinder · GaMD", fontsize=11, fontweight="bold", y=0.98)
    else:
        family_legend(fig, systems, set(d["category"]))
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return seam


def fig_traj_enzyme(system, long_df, landmark_order, path, dpi, systems, n_points=N_POINTS):
    """All landmarks of one enzyme along its trajectory, with the static structures on top."""
    d = long_df[(long_df["system"] == system) & long_df["plot"]]
    traj = d[d["kind"] == "trajectory"]
    if traj.empty:
        return
    thin = thin_frames(traj, n_points)
    panels = panels_for(d, landmark_order)
    lim = max(6.0, float(np.nanpercentile(np.abs(d[["x", "y"]].to_numpy()), 99.5)) + 1.5)
    hmin, hmax = float(np.nanpercentile(d["h"], 0.5)) - 1, float(np.nanpercentile(d["h"], 99.5)) + 1
    rmax = float(np.nanpercentile(d["r"], 99.5)) + 1
    mk, col = marker_for(system, systems), color_for(system, systems)
    fig, axes = plt.subplots(2, len(panels), figsize=(3.1 * len(panels), 6.4), squeeze=False)
    for j, lab in enumerate(panels):
        t = traj[traj["landmark"] == lab]
        st = d[(d["landmark"] == lab) & (d["kind"] != "trajectory")]
        a0, a1 = axes[0, j], axes[1, j]
        draw_reference_top(a0, lim)
        density_fill(a0, t["x"], t["y"], [[-lim, lim], [-lim, lim]])
        draw_frames(a1, thin[thin["landmark"] == lab], "r", "h", system, systems, size=7, alpha=0.35)
        a1.scatter([t["r"].median()], [t["h"].median()], s=55, marker=mk, color=col, edgecolors="white",
                   linewidths=0.8, zorder=5)
        for (sy, _), g in st.groupby(["system", "state"], sort=False):
            fill = state_style(g["category"].iloc[0])["fill"]
            for ax, xx, yy in ((a0, g["x"], g["y"]), (a1, g["r"], g["h"])):
                ax.scatter(xx, yy, s=60, marker=marker_for(sy, systems), color=fill, edgecolors="black",
                           linewidths=0.9, zorder=4)
        a0.set_xlim(-lim, lim)
        a0.set_ylim(-lim, lim)
        a0.set_aspect("equal")
        a0.set_title(display_name(lab), fontsize=10, fontweight="bold")
        a0.set_xlabel("x (Å)", fontsize=8)
        a1.set_xlim(0, rmax)
        a1.set_ylim(hmin, hmax)
        a1.set_xlabel("r (Å)", fontsize=8)
        for ax in (a0, a1):
            ax.tick_params(labelsize=7)
            style_plain(ax)
    axes[0, 0].set_ylabel("y (Å)", fontsize=8)
    axes[1, 0].set_ylabel("h (Å)", fontsize=8)
    cats = set(d.loc[d["kind"] != "trajectory", "category"])
    step = frame_step(traj, n_points)
    h = [Patch(color=TRAJ_BLUES[1], label="GaMD 50 %"), Patch(color=TRAJ_BLUES[0], label="GaMD 90 %"),
         Line2D([0], [0], marker=mk, ls="None", ms=4, mfc=col, mec="none", alpha=0.6,
                label=f"GaMD frame (1 every {step})" if step > 1 else "GaMD frame"),
         Line2D([0], [0], marker=mk, ls="None", ms=8, mfc=col, mec="white", label="median")]
    h += [Line2D([0], [0], marker=mk, ls="None", ms=8, mfc=state_style(c)["fill"], mec="black", label=STATE_LABEL[c])
          for c in ("global", "local", "docking") if c in cats]
    fig.subplots_adjust(wspace=0.35, hspace=0.3, top=0.86)
    fig.legend(handles=h, loc="lower center", ncol=len(h), bbox_to_anchor=(0.5, 0.9), fontsize=8, frameon=False)
    fig.suptitle(f"{system} · GaMD (unweighted)", fontsize=10, y=0.99)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def make_figures(out, cfg, df, long_df, landmark_order, dpi, r_min_phi):
    figdir = out / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    n_pts = int(cfg.get("points_per_trajectory", N_POINTS))
    r_color = tuple(cfg.get("r_color", R_COLOR))
    # one cut of the phi axis for static and GaMD figures, so that they can be compared
    seam = None
    if not long_df.empty:
        pocket = long_df[long_df["plot"] & ~long_df["landmark"].str.startswith("L-Tyr") & (long_df["r"] >= r_min_phi)]
        if not pocket.empty:
            seam = auto_seam(pocket["phi"].dropna())
    for kind, sfx in (("structure", ""), ("trajectory", "_GaMD")):
        part = long_df[long_df["kind"] == kind] if not long_df.empty else long_df
        if part.empty:
            continue
        traj = kind == "trajectory"
        fig_unrolled_cylinder(part, landmark_order, figdir / f"F2_unrolled_cylinder{sfx}.png", dpi, r_min_phi,
                              traj=traj, n_points=n_pts, r_color=r_color, seam=seam)
        fig_landmark_summary(part, landmark_order, figdir / f"F1_landmark_summary{sfx}.png", dpi, r_min_phi, seam,
                             traj=traj, n_points=n_pts)
        lm_dir = figdir / f"F1_per_landmark{sfx}"
        lm_dir.mkdir(exist_ok=True)
        for old in lm_dir.glob("F1_*.png"):
            old.unlink()
        part_sys = list(dict.fromkeys(part["system"]))
        need = min(2, len(part_sys))  # a landmark present in a single enzyme cannot be compared
        detail = [lab for lab in panels_for(part[part["plot"]], landmark_order)
                  if not lab.startswith("L-Tyr") and part.loc[part["landmark"] == lab, "system"].nunique() >= need]
        for k, lab in enumerate(detail, start=1):
            fig_landmark_detail(part, lab, part_sys, lm_dir / f"F1_{k:02d}_{file_name(display_name(lab))}.png", dpi,
                                traj=traj, n_points=n_pts)
    systems = list(dict.fromkeys(df["system"]))
    for sy in systems:
        if (long_df.loc[long_df["system"] == sy, "kind"] == "trajectory").any():
            fig_traj_enzyme(sy, long_df, landmark_order, figdir / f"F3_GaMD_{sy}.png", dpi, systems, n_pts)
    return seam


# ============================================================================ report
def group_warnings(warnings, df):
    states = {sy: list(dict.fromkeys(g["state"])) for sy, g in df.groupby("system", sort=False)}
    grouped = {}
    for sy, st, msg in warnings:
        grouped.setdefault((sy, msg), []).append(st)
    out = []
    for (sy, msg), sts in grouped.items():
        sts = list(dict.fromkeys(sts))
        where = "all states" if sy in states and set(states[sy]) <= set(sts) else ", ".join(sts)
        out.append(f"{sy} ({where}): {msg}")
    return out


def fmt_cell(v):
    v = v.dropna()
    if v.empty:
        return "—"
    if len(v) == 1:
        x = float(v.iloc[0])
        return f"{x:.2f}" if np.isfinite(x) else "—"
    return f"{np.median(v):.2f} [{np.percentile(v, 10):.2f}–{np.percentile(v, 90):.2f}]"


def write_report(path, df, qc, skipped, warnings, seam, pct, landmark_defs):
    L = ["# SOLARIS summary\n",
         "- Frame: Kabsch fit of the core (2 Cu + 6 His rings) onto the family template; z = pseudo-C2 axis",
         f"- States: {df.groupby(['system', 'state']).ngroups}; structures/frames: {len(df)}"]
    if seam is not None:
        L.append(f"- Cut of the unrolled cylinder (F2): φ = {seam:.0f}°")
    sections = [("States not processed", skipped), ("Dicopper core checks", core_checks(df)),
                ("Warnings", group_warnings(warnings, df) if warnings else [])]
    for title, items in sections:
        if items:
            L += [f"\n## {title}\n"] + [f"- {s}" for s in items]
    L.append("\n## Numbering and residues used\n")
    q = pd.DataFrame(qc, columns=QC_COLS)
    for (sy, st), g in q.groupby(["system", "state"], sort=False):
        roles = ", ".join(f"{r}={rn}{u}" for r, u, rn in zip(g["role"], g["resid_used"], g["resname"]))
        L.append(f"- **{sy} / {st}** (offset {int(g['offset'].iloc[0]):+d}, {g['offset_from'].iloc[0]}): {roles}")
    cols = [("Frame.core_RMSD", "core RMSD"), ("Cu.d_CuA_CuB", "Cu–Cu"), ("L.d_OH_CuA", "OH–CuA"),
            ("L.d_OH_CuB", "OH–CuB"), ("L.d_Cortho_Operox", "C ortho–O"), ("L.ring_tilt", "ring tilt")]
    cols = [c for c in cols if c[0] in df]
    L += ["\n## Key values per state\n", "| Enzyme | State | n | " + " | ".join(c[1] for c in cols) + " |",
          "|---|---|---|" + "---|" * len(cols)]
    for (sy, st), g in df.groupby(["system", "state"], sort=False):
        L.append(f"| {sy} | {st} | {len(g)} | " + " | ".join(fmt_cell(g[c]) for c, _ in cols) + " |")
    L.append("\nTrajectories: median [p10–p90]. Distances in Å, angles in degrees.")
    L.append("\n## Landmarks (plotted point): h / r / φ\n")
    for (sy, st), g in df.groupby(["system", "state"], sort=False):
        parts = []
        for lab, ld in landmark_defs.items():
            pt = ld.get("plot", "TIP")
            if f"{lab}.{pt}.h" not in g or g[f"{lab}.{pt}.h"].isna().all():
                continue
            phi = g[f"{lab}.{pt}.phi"].dropna()
            phi_s = f"{circ_stats(phi)[0]:.0f}°" if len(phi) else "—"
            parts.append(f"{lab} {g[f'{lab}.{pt}.h'].median():.1f}/{g[f'{lab}.{pt}.r'].median():.1f}/{phi_s}")
        if parts:
            L.append(f"- {sy} / {st}: " + "; ".join(parts))
    if pct is not None and not pct.empty:
        tail = pct[pct["key"] & pct["in_tail"]]
        L.append(f"\n## Static structures in the tails of their trajectory (key descriptors, percentile < {TAIL_PCT:g} "
                 f"or > {100 - TAIL_PCT:g})\n")
        if tail.empty:
            L.append("- None: minima and docking poses fall within what the trajectory samples.")
        for (sy, st), g in tail.groupby(["system", "state"], sort=False):
            L.append(f"- {sy} / {st}: " + ", ".join(f"{d} (p{p:.0f})" for d, p in zip(g["descriptor"], g["percentile"])))
    Path(path).write_text("\n".join(L) + "\n", encoding="utf-8")


# ============================================================================ command line
def parse_args():
    ap = argparse.ArgumentParser(description="SOLARIS: cylindrical description of type-3 copper sites.")
    ap.add_argument("config", help="JSON configuration file")
    ap.add_argument("--check", action="store_true", help="only check residues and numbering")
    ap.add_argument("--only", nargs="+", help="enzymes to process")
    ap.add_argument("--states", nargs="+", help="states to process (e.g. global local1 docking)")
    ap.add_argument("--no-gamd", action="store_true", help="skip trajectories")
    ap.add_argument("--no-docking", action="store_true", help="skip docking poses")
    ap.add_argument("--strict", action="store_true", help="stop if a file is missing")
    ap.add_argument("--dpi", type=int)
    ap.add_argument("--out", help="output folder (overrides output_dir)")
    ap.add_argument("--figures-only", action="store_true", help="redraw the figures from the saved tables")
    return ap.parse_args()


def classify(state, entry):
    n = state.lower()
    cat = ("global" if n.startswith("glob") else "local" if n.startswith("loc") else
           "docking" if n.startswith("dock") else "GaMD" if n.startswith(("gamd", "md", "traj")) else "other")
    return cat, isinstance(entry, dict) and any(k in entry for k in ("traj", "trajectory"))


def wanted(args, state, cat, is_traj):
    return not ((args.states and state not in args.states) or ((is_traj or cat == "GaMD") and args.no_gamd)
                or (cat == "docking" and args.no_docking))


def collect_jobs(cfg, base, args, landmark_defs):
    systems = [s for s in cfg["systems"] if not s.startswith("_") and (not args.only or s in args.only)]
    jobs, skipped, qc, warns = [], [], [], []
    for system in systems:
        sys_cfg = cfg["systems"][system]
        for state, spec in sys_cfg.get("structures", {}).items():
            if state.startswith("_"):
                continue
            cat, is_traj = classify(state, spec)
            if not wanted(args, state, cat, is_traj):
                continue
            entry = spec if isinstance(spec, dict) else {"pdb": spec}
            try:
                src = open_source(base, entry, is_traj)
                site = resolve_site(system, state, src.topo, sys_cfg, cfg, entry, landmark_defs)
            except Exception as exc:
                if args.strict:
                    raise
                msg = f"{system}/{state}: {'missing ' if isinstance(exc, FileNotFoundError) else ''}{exc}"
                skipped.append(msg)
                print(f"[SKIP] {msg}")
                continue
            kind = "trajectory" if (is_traj or src.n_frames > 1) else "structure"
            if kind == "trajectory" and cat not in ("GaMD", "other"):
                cat = "GaMD"
            jobs.append((system, state, kind, cat, src, site, entry))
            qc += site.qc
            warns += [(system, state, w) for w in site.warnings]
            print(f"[OK] {system}/{state}: {src.n_frames} frame(s), offset {site.offset:+d} ({site.offset_from}); "
                  + ", ".join(f"{lab}={lm['resname']}{lm['resid']}" for lab, lm in site.landmarks.items())
                  + f"; ligand {'yes' if site.ligands else 'no'}; peroxide {'yes' if site.peroxide else 'no'}")
    return jobs, skipped, qc, warns


def role_inconsistencies(qc):
    """A role must map to the same residue type in every state of an enzyme."""
    bad = []
    for (sy, role), g in pd.DataFrame(qc, columns=QC_COLS).groupby(["system", "role"]):
        types = {std_resname(x) for x in g["resname"]}
        if role != "ligand" and len(types) > 1:
            bad.append((sy, "several", f"role {role} is {'/'.join(sorted(types))} depending on the state ("
                        + ", ".join(f"{s}={r}{u}" for s, r, u in zip(g["state"], g["resname"], g["resid_used"]))
                        + "): check the numbering"))
    return bad


def filter_tables(tables, args):
    out = []
    for t in tables:
        if not t.empty:
            m = pd.Series(True, index=t.index)
            if args.only:
                m &= t["system"].isin(args.only)
            if args.states:
                m &= t["state"].isin(args.states)
            if args.no_docking:
                m &= t["category"] != "docking"
            if args.no_gamd:
                m &= t["kind"] != "trajectory"
            t = t[m]
        out.append(t)
    return out


def load_tables(out):
    tdir = out / "tables"
    paths = [tdir / "descriptors_per_frame.csv", tdir / "landmarks_long.csv"]
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise SystemExit(f"--figures-only: {', '.join(missing)} not found (check --out or output_dir)")
    df, long_df = pd.read_csv(paths[0]), pd.read_csv(paths[1])
    if long_df["plot"].dtype != bool:
        long_df["plot"] = long_df["plot"].astype(str).str.strip().str.lower().isin(["true", "1"])
    return df, long_df


def main():
    args = parse_args()
    cfg_path = Path(args.config).resolve()
    base = cfg_path.parent
    with open(cfg_path, encoding="utf-8") as fh:
        cfg = json.load(fh)
    out = Path(args.out).resolve() if args.out else resolve_path(base, cfg.get("output_dir", "SOLARIS_out"))
    landmark_defs = {k: v for k, v in cfg.get("landmarks", DEFAULT_LANDMARKS).items() if not k.startswith("_")}
    landmark_order = list(landmark_defs)
    dpi = args.dpi or int(cfg.get("dpi", 300))
    r_min_phi = float(cfg.get("r_min_phi", 2.5))

    if args.figures_only:
        df, long_df = filter_tables(load_tables(out), args)
        make_figures(out, cfg, df, long_df, landmark_order, dpi, r_min_phi)
        print(f"Figures: {out / 'figures'}")
        return

    jobs, skipped, qc, warns = collect_jobs(cfg, base, args, landmark_defs)
    bad = role_inconsistencies(qc)
    for sy, _, b in bad:
        print(f"[ERROR] {sy}: {b}")
    if args.check:
        print("\nCheck finished." + (" INCONSISTENCIES FOUND." if bad else ""))
        return
    if bad and args.strict:
        raise SystemExit("Inconsistent residues between states (see above).")
    warns += bad
    if not jobs:
        raise SystemExit("Nothing to process.")

    rows, long_rows = [], []
    for system, state, kind, cat, src, site, entry in jobs:
        stride = int(entry.get("stride", 10 if kind == "trajectory" else 1))
        stop = entry.get("stop")
        print(f"[RUN] {system}/{state} ({kind}{f', stride {stride}' if kind == 'trajectory' else ''})")
        r, lr = FrameAnalyzer(system, state, kind, cat, site, r_min_phi).run(
            src, int(entry.get("start", 0)), None if stop is None else int(stop), stride)
        rows += r
        long_rows += lr

    df = pd.DataFrame(rows)
    desc_cols = [c for c in df.columns if c not in META_COLS]
    df = df[META_COLS + desc_cols]
    long_df = pd.DataFrame(long_rows)
    tdir = out / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    df.to_csv(tdir / "descriptors_per_frame.csv", index=False, float_format="%.4f")
    long_df.to_csv(tdir / "landmarks_long.csv", index=False, float_format="%.4f")
    pd.DataFrame(qc, columns=QC_COLS).to_csv(tdir / "residue_QC.csv", index=False)
    summarize(df, desc_cols).to_csv(tdir / "state_summary.csv", index=False, float_format="%.4f")
    pct = percentile_in_trajectory(df, desc_cols, key_columns(df, landmark_defs))
    if not pct.empty:
        pct.to_csv(tdir / "percentile_in_trajectory.csv", index=False, float_format="%.3f")

    seam = make_figures(out, cfg, df, long_df, landmark_order, dpi, r_min_phi)
    write_report(out / "SUMMARY.md", df, qc, skipped, warns, seam, pct if not pct.empty else None, landmark_defs)
    print(f"\nDone.\nTables:  {tdir}\nFigures: {out / 'figures'}\nSummary: {out / 'SUMMARY.md'}")
    if skipped:
        print(f"Skipped states: {len(skipped)} (see SUMMARY.md)")


if __name__ == "__main__":
    main()
