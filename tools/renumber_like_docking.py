#!/usr/bin/env python3
"""
Renumber minimised structures (sequential topology numbering, His HD1-HD6, CU1/CU2, OY1, TY1) so that they
share the residue numbering of the docking structure of the same enzyme.

Protein residues are matched by a global sequence alignment; Cu, peroxide and ligand take the numbers used in
the docking file (or the first free number if that number belongs to a protein residue). Residue names, atom
names and coordinates are not changed.

    python3 renumber_like_docking.py DOCKING.pdb MIN1.pdb [MIN2.pdb ...] [--suffix _renum] [--out DIR]
"""
import argparse
import sys
from pathlib import Path

AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
      "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
      "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
ALIAS = {"HIE": "HIS", "HID": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "CYX": "CYS",
         "CYM": "CYS", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "ARN": "ARG", "TYM": "TYR",
         **{f"HD{k}": "HIS" for k in range(1, 10)}}
CU = {"CU", "CU1", "CU2", "CUA", "CUB"}
PEROXIDE = {"PER", "OXY", "OY1", "PEO", "O2"}
LIGAND = {"TY1", "LIG", "TYR_L"}


def parse_res(line):
    """(chain, resseq, icode, resname), tolerant to shifted columns such as 'CU  A 400'."""
    tok = line[16:27].split()
    resname = tok[0] if tok else line[17:20].strip()
    chain = tok[1] if len(tok) == 3 else line[21] if len(tok) < 2 else " "
    last = tok[-1] if len(tok) > 1 else line[22:27].strip()
    icode = last[-1] if last and last[-1].isalpha() else ""
    try:
        return chain, int(last[:-1] if icode else last), icode, resname
    except ValueError:
        return line[21], int(line[22:26]), line[26].strip(), line[17:21].strip()


def residues(path):
    out, seen = [], set()
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            key = parse_res(line)
            if key not in seen:
                seen.add(key)
                out.append(key)
        elif line.startswith("ENDMDL"):
            break
    return out


def one_letter(resname):
    return AA.get(ALIAS.get(resname, resname))


def align(a, b, match=2, mismatch=-1, gap=-2):
    """Needleman-Wunsch; returns aligned index pairs (None = gap)."""
    n, m = len(a), len(b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        S[i][0] = i * gap
    for j in range(1, m + 1):
        S[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            S[i][j] = max(S[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch),
                          S[i - 1][j] + gap, S[i][j - 1] + gap)
    i, j, pairs = n, m, []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and S[i][j] == S[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch):
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif i > 0 and S[i][j] == S[i - 1][j] + gap:
            pairs.append((i - 1, None))
            i -= 1
        else:
            pairs.append((None, j - 1))
            j -= 1
    return pairs[::-1]


def build_map(dock_path, min_path):
    dock, mini = residues(dock_path), residues(min_path)
    dprot = [r for r in dock if one_letter(r[3])]
    mprot = [r for r in mini if one_letter(r[3])]
    chain = dprot[0][0] if dprot else " "
    lab = lambda r: f"{r[3]}{r[1]}{r[2]}"
    mapping, notes = {}, []
    for i, j in align([one_letter(r[3]) for r in mprot], [one_letter(r[3]) for r in dprot]):
        if i is None:
            notes.append(f"docking {lab(dprot[j])} has no match in the minimum")
        elif j is None:
            notes.append(f"minimum {lab(mprot[i])} has no match in the docking: left unchanged")
        else:
            if one_letter(mprot[i][3]) != one_letter(dprot[j][3]):
                notes.append(f"different: minimum {lab(mprot[i])} <-> docking {lab(dprot[j])}")
            mapping[mprot[i]] = (dprot[j][1], dprot[j][2])

    prot_nums = {v[0] for v in mapping.values()}
    next_free = [max(prot_nums | {r[1] for r in dock}, default=899) + 1]

    def take(want, what, r):
        if want is not None and want not in prot_nums:
            return (want, "")
        num = next_free[0]
        next_free[0] += 1
        why = f"docking uses {want} for a protein residue too" if want is not None else "absent from the docking"
        notes.append(f"{what} {r[3]}: {num} ({why})")
        return (num, "")

    dcu = [r for r in dock if r[3] in CU]
    dper = [r for r in dock if r[3] in PEROXIDE]
    dlig = [r for r in dock if r[3] in LIGAND and r[1] > 0]
    for k, r in enumerate([r for r in mini if r[3] in CU]):
        mapping[r] = take(dcu[k][1] if k < len(dcu) else None, "copper", r)
    for r in [r for r in mini if r[3] in PEROXIDE]:
        mapping[r] = take(dper[0][1] if dper else None, "peroxide", r)
    for r in [r for r in mini if r[3] in LIGAND]:
        mapping[r] = take(dlig[0][1] if dlig else None, "ligand", r)
    notes += [f"residue {lab(r)} has no rule: left unchanged" for r in mini if r not in mapping]
    return mapping, chain, notes, len(mprot), len(dprot)


def write(min_path, out_path, mapping, chain):
    lines = []
    for line in open(min_path, encoding="utf-8", errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            key = parse_res(line)
            if key in mapping:
                num, ic = mapping[key]
                line = line[:21] + chain + f"{num:4d}" + (ic or " ") + line[27:]
        elif line.startswith("TER"):
            line = "TER\n"
        lines.append(line)
    out_path.write_text("".join(lines), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description="Renumber minimised structures like the docking structure.")
    ap.add_argument("docking")
    ap.add_argument("minima", nargs="+")
    ap.add_argument("--suffix", default="_renum")
    ap.add_argument("--out", help="output folder (default: next to each minimum)")
    args = ap.parse_args()
    ok = True
    for mp in map(Path, args.minima):
        mapping, chain, notes, nm, nd = build_map(args.docking, mp)
        outdir = Path(args.out) if args.out else mp.parent
        outdir.mkdir(parents=True, exist_ok=True)
        out = outdir / f"{mp.stem}{args.suffix}.pdb"
        write(mp, out, mapping, chain)
        report = [f"{mp} -> {out}", f"docking: {args.docking} (chain '{chain}')",
                  f"protein residues: minimum {nm}, docking {nd}", "",
                  "Notes:" if notes else "No notes: the sequences match residue by residue.", *notes, "",
                  "Mapping (minimum -> docking):"]
        report += [f"  {k[3]:>4} {k[1]:>5}{k[2]} -> {v[0]}{v[1]}" for k, v in mapping.items()]
        (outdir / f"{mp.stem}{args.suffix}.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
        bad = [n for n in notes if n.startswith(("different", "minimum ", "residue"))]
        ok &= not bad
        print(f"[{'OK' if not bad else 'WARN'}] {mp.name} -> {out.name} ({nm} protein residues; {len(notes)} note(s))")
        for n in notes:
            print("      " + n)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
