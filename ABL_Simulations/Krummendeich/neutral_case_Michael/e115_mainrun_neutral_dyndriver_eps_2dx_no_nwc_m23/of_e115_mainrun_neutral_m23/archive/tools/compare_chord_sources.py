#!/usr/bin/env python3
"""Compare chord vs span/radius across three representations of the same blade:
  - AeroDyn v14 blade table (RNodes = radius from rotor centre, Chord)
  - AeroDyn v15 blade table (BlSpn = distance from blade root; +HubRad = radius)
  - the CAD (IGES) blade-skin geometry (chord measured per spanwise section)

The key convention difference: v14 RNodes is measured from the rotor centre,
v15 BlSpn from the blade root. We add HubRad (1.51 m) to BlSpn so all three are
plotted against radius from the rotor centre and are directly comparable.

CAD chord is read from data/cad_chord_vs_radius.csv (extracted once from the
IGES blade-skin surface, so this script needs no CAD library or the 5 MB IGES).
The v14 table path is given on the command line; v15 is read from the repo.
"""
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
HUBRAD = 1.51   # m, ElastoDyn HubRad: BlSpn (from root) + HubRad = radius from centre


def read_v14(path):
    rows = [r for r in csv.reader(open(path)) if r and not r[0].startswith("RNodes")]
    a = np.array([[float(x) for x in r] for r in rows])
    return a[:, 0], a[:, 3], a[:, 1]          # RNodes(radius), Chord, AeroTwst


def read_v15_blade(path):
    raw = [l.rstrip("\r\n") for l in open(path, errors="replace")]
    hdr = next(i for i, l in enumerate(raw) if "BlSpn" in l and "BlChord" in l)
    rows = []
    for l in raw[hdr + 2:]:
        p = l.split()
        if len(p) < 7:
            continue
        try:
            rows.append([float(p[0]), float(p[4]), float(p[5])])  # BlSpn, BlTwist, BlChord
        except ValueError:
            continue
    a = np.array(rows)
    return a[:, 0] + HUBRAD, a[:, 2], a[:, 1]     # radius, chord, twist


def read_cad(path):
    a = np.array([[float(x) for x in r] for r in list(csv.reader(open(path)))[1:] if r])
    return a[:, 0], a[:, 1], a[:, 2]              # radius, chord, t/c


def main():
    v14_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "data", "ad14_blade.csv")
    r14, c14, t14 = read_v14(v14_path)
    r15, c15, t15 = read_v15_blade(os.path.join(REPO, "AeroDyn15", "AeroDyn15_blade.dat"))
    cad_csv = os.path.join(REPO, "data", "cad_chord_vs_radius.csv")
    rC, cC, tcC = read_cad(cad_csv)

    rr = np.linspace(6, 56, 60)
    a14, a15, aC = np.interp(rr, r14, c14), np.interp(rr, r15, c15), np.interp(rr, rC, cC)
    print("chord agreement on radius 6-56 m (mean |diff|, max |diff|):")
    print(f"  v14 vs v15 : {np.abs(a14-a15).mean()*100:.1f} cm, {np.abs(a14-a15).max()*100:.1f} cm")
    print(f"  v15 vs CAD : {np.abs(a15-aC).mean()*100:.1f} cm, {np.abs(a15-aC).max()*100:.1f} cm")
    print(f"  v14 vs CAD : {np.abs(a14-aC).mean()*100:.1f} cm, {np.abs(a14-aC).max()*100:.1f} cm")

    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    ax[0].plot(r14, c14, "^-", color="#1f77b4", ms=5, label=f"AeroDyn v14 ({len(r14)} nodes)")
    ax[0].plot(r15, c15, "o-", color="0.25", ms=4, label=f"AeroDyn v15 ({len(r15)} nodes)")
    ax[0].plot(rC, cC, "s--", color="#d62728", ms=4, mfc="none", label=f"CAD skin ({len(rC)} sections)")
    ax[0].set(xlabel="radius from rotor centre [m]", ylabel="chord [m]",
              title="Chord vs radius: AeroDyn v14 / v15 / CAD")
    ax[0].grid(alpha=0.3); ax[0].legend()

    ax[1].plot(rr, (a14 - a15) * 100, color="#1f77b4", label="v14 − v15")
    ax[1].plot(rr, (a15 - aC) * 100, color="#d62728", label="v15 − CAD")
    ax[1].plot(rr, (a14 - aC) * 100, color="0.4", label="v14 − CAD")
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set(xlabel="radius from rotor centre [m]", ylabel="chord difference [cm]",
              title="Pairwise chord difference (same blade, ~few cm)")
    ax[1].grid(alpha=0.3); ax[1].legend()
    fig.suptitle("Same blade in three representations: AeroDyn v14, AeroDyn v15, and the CAD",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    out = os.path.join(REPO, "figures", "chord_v14_v15_cad.png")
    fig.savefig(out, dpi=145, bbox_inches="tight")
    print("\nwrote", os.path.relpath(out, REPO))


if __name__ == "__main__":
    main()
