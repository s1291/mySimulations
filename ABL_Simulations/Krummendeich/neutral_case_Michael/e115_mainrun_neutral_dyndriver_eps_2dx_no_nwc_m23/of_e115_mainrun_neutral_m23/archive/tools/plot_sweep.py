#!/usr/bin/env python3
"""Comprehensive validation figure: original 44-node vs uniform 57-node blended
blade across a wind-speed sweep, run in OpenFAST 3.5.5 (rigid, fixed 13.2 rpm,
fixed pitch). Top row: rotor power/thrust/torque vs wind speed for both models
and the % difference. Bottom row: spanwise normal force, angle of attack and
axial induction at three wind speeds spanning attached to deep-stall flow.
Reads /tmp/sweep/{base,uniform}_<V>/Turbine_driver.out."""
import os, importlib.util
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

VS = [5, 8, 11, 14, 17, 20, 25]
SPAN_VS = [8, 14, 20]          # representative: below rated, near rated, deep stall
T0 = 120.0
HUB, TIP = 1.51, 58.02

s = importlib.util.spec_from_file_location("c", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "compare_validation.py"))
c = importlib.util.module_from_spec(s); s.loader.exec_module(c)
u = c.u


def grids():
    _, db = u.read_blade("/tmp/sweep/base_11/AeroDyn15/AeroDyn15_blade.dat")
    _, dd = u.read_blade("/tmp/sweep/uniform_11/AeroDyn15/AeroDyn15_blade_uniform.dat")
    return (db[:, 0] + HUB) / TIP, (dd[:, 0] + HUB) / TIP, len(db), len(dd)


def main():
    rb, ru, nb_n, nu_n = grids()
    P = {"b": {}, "u": {}}
    data = {}
    for V in VS:
        nB, _, B = c.read_out(f"/tmp/sweep/base_{V}/Turbine_driver.out")
        nU, _, U = c.read_out(f"/tmp/sweep/uniform_{V}/Turbine_driver.out")
        data[V] = (nB, B, nU, U)
        for ch, sc in [("RtAeroPwr", 1e-3), ("RtAeroFxh", 1e-3), ("RtAeroMxh", 1e-3)]:
            P["b"].setdefault(ch, []).append(c.avg(nB, B, ch) * sc)
            P["u"].setdefault(ch, []).append(c.avg(nU, U, ch) * sc)

    fig = plt.figure(figsize=(17, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.33, wspace=0.27)
    V = np.array(VS)
    titles = {"RtAeroPwr": "Aerodynamic power [kW]", "RtAeroFxh": "Rotor thrust [kN]",
              "RtAeroMxh": "Aerodynamic torque [kN·m]"}
    for k, ch in enumerate(["RtAeroPwr", "RtAeroFxh", "RtAeroMxh"]):
        ax = fig.add_subplot(gs[0, k])
        b = np.array(P["b"][ch]); uu = np.array(P["u"][ch])
        ax.plot(V, b, "o-", color="0.25", label="original (44)")
        ax.plot(V, uu, "s--", color="#d62728", mfc="none", label="uniform (57)")
        ax.set_title(titles[ch]); ax.set_xlabel("wind speed [m/s]"); ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="upper left")
        axr = ax.twinx()
        rel = 100 * (uu - b) / np.where(np.abs(b) < 1e-9, np.nan, b)
        axr.plot(V, rel, "^:", color="#1f77b4", ms=5, alpha=0.7)
        axr.set_ylabel("difference [%]", color="#1f77b4", fontsize=9)
        axr.tick_params(axis="y", labelcolor="#1f77b4", labelsize=8)
        axr.axhline(0, color="#1f77b4", lw=0.5, alpha=0.4)

    chans = [("Fn", "Normal force Fn [N/m]"), ("Alpha", "Angle of attack [deg]"),
             ("AxInd", "Axial induction [-]")]
    colors = plt.cm.viridis(np.linspace(0.1, 0.8, len(SPAN_VS)))
    m = (rb >= 0.1) & (rb <= 0.95)
    print(f"Spanwise max|Δ| in loaded span (r/R 0.1-0.95), by wind speed:")
    for k, (ch, lab) in enumerate(chans):
        ax = fig.add_subplot(gs[1, k])
        ax.axvspan(0.1, 0.95, color="0.94", zorder=0)
        for col, Vv in zip(colors, SPAN_VS):
            nB, B, nU, U = data[Vv]
            vb = c.spanwise(nB, B, ch, nb_n); vu = c.spanwise(nU, U, ch, nu_n)
            ax.plot(rb, vb, "-", color=col, lw=1.3, label=f"{Vv} m/s orig")
            ax.plot(ru, vu, "x", color=col, ms=4, label=f"{Vv} m/s unif")
            d = np.abs(np.interp(rb, ru, vu) - vb)[m]
            if k == 0:
                print(f"  {ch:6} {Vv:2d} m/s: max {d.max():8.3g}, RMS {np.sqrt((d**2).mean()):8.3g}")
        ax.set_title(lab); ax.set_xlabel("r/R"); ax.grid(alpha=0.3)
        if k == 0:
            ax.legend(fontsize=7, ncol=3, loc="upper left")
    for k, (ch, lab) in enumerate(chans):
        if k == 0: continue
        nB, B, nU, U = data[14]

    fig.suptitle("OpenFAST 3.5.5 validation across wind speed: original 44-node vs "
                 "uniform 57-node blended blade\n(rigid, fixed 13.2 rpm, fixed pitch; "
                 "UA params recomputed by OpenFAST from each blended polar)",
                 fontsize=13, fontweight="bold")
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "..", "figures", "uniform_validation_sweep.png")
    fig.savefig(out, dpi=145, bbox_inches="tight")
    print("\nwrote", os.path.normpath(out))


if __name__ == "__main__":
    main()
