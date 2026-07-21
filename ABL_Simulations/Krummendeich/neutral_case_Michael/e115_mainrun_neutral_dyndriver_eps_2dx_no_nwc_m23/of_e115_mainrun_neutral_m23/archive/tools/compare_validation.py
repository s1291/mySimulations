#!/usr/bin/env python3
"""Compare the OpenFAST 3.5.5 results of the original (44-node) and the uniform
blended (57-node) AeroDyn blade, to validate that the re-discretization
preserves the aerodynamics. Reports rotor power, thrust and torque, and the
spanwise distributions, averaged over a steady window.
"""
import os, re, importlib.util
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/tmp/run_base"
UNI = "/tmp/run_uniform"
T0 = 150.0           # start of steady averaging window (s)
HUB, TIP = 1.51, 58.02

s = importlib.util.spec_from_file_location("u", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "uniformize_blade.py"))
u = importlib.util.module_from_spec(s); s.loader.exec_module(u)


def read_out(path):
    L = open(path, errors="replace").read().splitlines()
    hi = next(i for i, l in enumerate(L) if l.split()[:1] == ["Time"])
    names = L[hi].split()
    units = L[hi + 1].split()
    data = np.array([[float(x) for x in ln.split()] for ln in L[hi + 2:] if ln.strip()])
    return names, units, data


def avg(names, data, name):
    t = data[:, 0]
    return data[t >= T0, names.index(name)].mean()


def spanwise(names, data, ch, nb):
    return np.array([avg(names, data, f"AB1N{n:03d}{ch}") for n in range(1, nb + 1)])


def main():
    nb_b, db = u.read_blade(f"{BASE}/AeroDyn15/AeroDyn15_blade.dat")
    nb_u, du = u.read_blade(f"{UNI}/AeroDyn15/AeroDyn15_blade_uniform.dat")
    rb = (db[:, 0] + HUB) / TIP
    ru = (du[:, 0] + HUB) / TIP

    nB, uB, B = read_out(f"{BASE}/Turbine_driver.out")
    nU, uU, U = read_out(f"{UNI}/Turbine_driver.out")

    print(f"steady window: t >= {T0}s  (base {int((B[:,0]>=T0).sum())} samples, "
          f"uniform {int((U[:,0]>=T0).sum())} samples)\n")
    print(f"{'quantity':18}{'original':>14}{'uniform':>14}{'diff':>12}{'rel %':>9}")
    rows = [("RtAeroPwr", "kW", 1e-3), ("RtAeroFxh", "kN", 1e-3),
            ("RtAeroMxh", "kN·m", 1e-3), ("RtAeroCp", "-", 1.0),
            ("RtAeroCt", "-", 1.0)]
    summary = {}
    for ch, un, sc in rows:
        a, b = avg(nB, B, ch) * sc, avg(nU, U, ch) * sc
        rel = 100 * (b - a) / a if a else 0.0
        summary[ch] = (a, b, rel)
        print(f"{ch+' ['+un+']':18}{a:14.4f}{b:14.4f}{b-a:12.4f}{rel:9.3f}")

    chans = ["Fn", "Ft", "Alpha", "Cl", "Cd", "AxInd"]
    un = {"Fn": "N/m", "Ft": "N/m", "Alpha": "deg", "Cl": "-", "Cd": "-", "AxInd": "-"}
    m = (rb >= 0.1) & (rb <= 0.95)        # aerodynamically loaded region
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    print(f"\n{'spanwise':8}{'global max|d|':>15}{'loaded max|d|':>15}{'loaded RMS':>13}"
          "   (loaded = r/R 0.1..0.95)")
    for k, ch in enumerate(chans):
        vb = spanwise(nB, B, ch, len(db))
        vu = spanwise(nU, U, ch, len(du))
        a = ax[k // 3, k % 3]
        a.axvspan(0.1, 0.95, color="0.92", zorder=0)
        a.plot(rb, vb, "o-", color="0.25", ms=4, label="original (44)")
        a.plot(ru, vu, "s-", color="#d62728", ms=3, mfc="none", label="uniform (57)")
        d = np.interp(rb, ru, vu) - vb
        a.set(title=f"{ch} [{un[ch]}]", xlabel="r/R")
        a.grid(alpha=0.3); a.legend(fontsize=8)
        print(f"{ch:8}{np.abs(d).max():15.4g}{np.abs(d[m]).max():15.4g}"
              f"{np.sqrt((d[m]**2).mean()):13.4g}")
    p = summary
    txt = (f"power {p['RtAeroPwr'][0]:.0f}->{p['RtAeroPwr'][1]:.0f} kW ({p['RtAeroPwr'][2]:+.2f}%)    "
           f"thrust {p['RtAeroFxh'][0]:.0f}->{p['RtAeroFxh'][1]:.0f} kN ({p['RtAeroFxh'][2]:+.2f}%)    "
           f"torque {p['RtAeroMxh'][0]:.0f}->{p['RtAeroMxh'][1]:.0f} kN.m ({p['RtAeroMxh'][2]:+.2f}%)")
    fig.suptitle("OpenFAST 3.5.5: original 44-node vs uniform 57-node blended blade "
                 "(rigid, steady 11 m/s, 13.2 rpm)\n" + txt,
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    here = os.path.dirname(os.path.abspath(__file__))
    for out in (os.path.join(os.path.dirname(BASE), "validation_spanwise.png"),
                os.path.join(here, "..", "figures", "uniform_validation.png")):
        fig.savefig(out, dpi=150, bbox_inches="tight")
    print("\nshaded band = loaded region used for the in-band statistics")
    return summary


if __name__ == "__main__":
    main()
