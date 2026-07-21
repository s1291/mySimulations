#!/usr/bin/env python3
"""Reconcile a uniform AeroDyn v15 blade's CHORD to the CAD/v14 planform.

Motivation: the AeroDyn v15 blade's chord is ~6-9% narrower than both the CAD
geometry and the AeroDyn v14 table over the OUTER span (r > ~40 m); from r~8 m
to ~38 m all three agree to 1-2 cm (see tools/compare_chord_sources.py). If the
CAD geometry is the reference for CFD coupling, this tool rewrites ONLY the
BlChord column over the outer span to the CAD/v14 chord, leaving span positions,
twist, prebend, sweep, airfoil IDs and all polar files untouched.

Scope of the change (deliberately limited):
  * r < R_LO  (default 38 m): chord left exactly as v15. The inboard already
    matches; the root/cylinder transition (r < ~7 m) has an ill-defined "chord"
    in the CAD (near-circular sections), so it must NOT be overwritten.
  * R_LO..R_HI: smooth cosine ramp from v15 to the CAD/v14 reference, so there
    is no kink where the correction switches on.
  * r > R_HI  (default 40 m): full CAD/v14 reference chord.

Reference chord = CAD (NURBS-evaluated, data/cad_chord_vs_radius.csv) where the
CAD covers the radius, else the v14 table (data/ad14_blade.csv). Radius =
BlSpn + HubRad.

Usage:  python tools/rechord_to_cad.py <in_model_dir> <out_model_dir>
        (defaults: AeroDyn15_uniform  ->  AeroDyn15_uniform_cadchord)
"""
import os, sys, csv, shutil, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
HUBRAD = 1.51
R_LO, R_HI = 38.0, 40.0          # blend band [m]: v15 below R_LO, CAD/v14 above R_HI


def _read_csv(path, cols):
    rows = [r for r in csv.reader(open(path)) if r and not r[0][0].isalpha()]
    a = np.array([[float(x) for x in r] for r in rows])
    return tuple(a[:, c] for c in cols)


def reference_chord_fn():
    rC, cC = _read_csv(os.path.join(REPO, "data", "cad_chord_vs_radius.csv"), (0, 1))
    r14, c14 = _read_csv(os.path.join(REPO, "data", "ad14_blade.csv"), (0, 3))
    lo, hi = rC.min(), rC.max()

    def chord(r):
        r = np.atleast_1d(np.asarray(r, dtype=float))
        out = np.empty_like(r)
        in_cad = (r >= lo) & (r <= hi)
        out[in_cad] = np.interp(r[in_cad], rC, cC)
        out[~in_cad] = np.interp(r[~in_cad], r14, c14)
        return out
    return chord


def rechord(src_dir, dst_dir):
    shutil.rmtree(dst_dir, ignore_errors=True)
    shutil.copytree(src_dir, dst_dir)
    blade = next(f for f in os.listdir(dst_dir) if f.endswith("_blade_uniform.dat")
                 or (f.endswith(".dat") and "blade" in f.lower()))
    path = os.path.join(dst_dir, blade)
    raw = [l.rstrip("\r\n") for l in open(path, errors="replace")]
    hdr = next(i for i, l in enumerate(raw) if "BlSpn" in l and "BlChord" in l)
    chord_of = reference_chord_fn()
    out, changed, maxd = list(raw[:hdr + 2]), 0, 0.0
    for l in raw[hdr + 2:]:
        p = l.split()
        if len(p) < 7:
            out.append(l); continue
        try:
            vals = [float(x) for x in p[:6]]
        except ValueError:
            out.append(l); continue
        span = vals[0]
        r = span + HUBRAD
        v15_c = vals[5]
        ref_c = float(chord_of(r)[0])
        if r <= R_LO:
            blend = 0.0
        elif r >= R_HI:
            blend = 1.0
        else:
            blend = 0.5 * (1 - np.cos(np.pi * (r - R_LO) / (R_HI - R_LO)))
        new_c = (1 - blend) * v15_c + blend * ref_c
        if blend > 0:
            maxd = max(maxd, abs(new_c - v15_c)); changed += 1
        rest = p[6:]
        out.append(f"{vals[0]: .6E} {vals[1]: .6E} {vals[2]: .6E} {vals[3]: .6E} "
                   f"{vals[4]: .6E} {new_c: .6E} " + "  ".join(rest))
    if len(out) > 1:
        out[1] = ("Chord reconciled to CAD/v14 planform by rechord_to_cad.py "
                  "(only BlChord changed; aero polars & all else identical)")
    open(path, "w", newline="").write("\r\n".join(out) + "\r\n")
    return changed, maxd, path


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "AeroDyn15_uniform")
    dst = sys.argv[2] if len(sys.argv) > 2 else os.path.join(REPO, "AeroDyn15_uniform_cadchord")
    n, maxd, path = rechord(src, dst)
    print(f"rewrote BlChord on {n} nodes -> {os.path.relpath(path, REPO)}")
    print(f"max chord change vs v15: {maxd*100:.1f} cm")


if __name__ == "__main__":
    main()
