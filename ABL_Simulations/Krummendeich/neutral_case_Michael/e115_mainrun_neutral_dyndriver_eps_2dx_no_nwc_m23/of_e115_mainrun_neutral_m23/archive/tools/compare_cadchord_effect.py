#!/usr/bin/env python3
"""Effect of the CAD/v14 outboard chord reconciliation on the OpenFAST loads.
Compares the v15-chord uniform model (/tmp/sweep/uniform_<V>) against the
CAD-reconciled-chord model (/tmp/sweep/cadchord_<V>) across the wind sweep,
everything else identical (only BlChord for r>38 m differs). Reports rotor
power/thrust/torque and the spanwise normal force in the reconciled region."""
import os, importlib.util
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

VS = [5, 8, 11, 14, 17, 20, 25]
T0 = 120.0
HUB, TIP = 1.51, 58.02
s = importlib.util.spec_from_file_location("c", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "compare_validation.py"))
c = importlib.util.module_from_spec(s); s.loader.exec_module(c)
u = c.u


def trio(tag, V):
    n, _, D = c.read_out(f"/tmp/sweep/{tag}_{V}/Turbine_driver.out")
    return (c.avg(n, D, "RtAeroPwr")*1e-3, c.avg(n, D, "RtAeroFxh")*1e-3,
            c.avg(n, D, "RtAeroMxh")*1e-3)


def main():
    print(f"{'V':>4} | {'P_v15':>8} {'P_cad':>8} {'dP%':>6} | "
          f"{'T_v15':>7} {'T_cad':>7} {'dT%':>6} | {'Q_v15':>8} {'Q_cad':>8} {'dQ%':>6}")
    print("-"*86)
    P = {"v15": [], "cad": []}
    Tg = {"v15": [], "cad": []}
    rows = []
    for V in VS:
        p15, t15, q15 = trio("uniform", V)
        pc, tc, qc = trio("cadchord", V)
        P["v15"].append(p15); P["cad"].append(pc)
        Tg["v15"].append(t15); Tg["cad"].append(tc)
        dP = 100*(pc-p15)/p15 if p15 else float("nan")
        dT = 100*(tc-t15)/t15 if t15 else float("nan")
        dQ = 100*(qc-q15)/q15 if q15 else float("nan")
        rows.append((V, dP, dT, dQ))
        print(f"{V:>4} | {p15:8.1f} {pc:8.1f} {dP:6.2f} | {t15:7.1f} {tc:7.1f} {dT:6.2f} | "
              f"{q15:8.1f} {qc:8.1f} {dQ:6.2f}")
    print("-"*86)
    print("(+ = reconciled CAD-chord gives MORE than v15-chord; outboard chord is wider, so expect +)")

    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    V = np.array(VS)
    for k, (key, lab, un) in enumerate([(P, "Aero power", "kW"), (Tg, "Rotor thrust", "kN")][:2] +
                                       [(None, "Power difference", "%")]):
        if key is None:
            base = np.array(P["v15"]); rec = np.array(P["cad"])
            ax[k].plot(V, 100*(rec-base)/np.where(np.abs(base)<1e-9, np.nan, base),
                       "o-", color="#9467bd", label="power")
            tb = np.array(Tg["v15"]); tr = np.array(Tg["cad"])
            ax[k].plot(V, 100*(tr-tb)/tb, "s-", color="#2ca02c", label="thrust")
            ax[k].axhline(0, color="k", lw=0.5)
            ax[k].set(xlabel="wind speed [m/s]", ylabel="difference [%]",
                      title="CAD-chord vs v15-chord (% change)")
            ax[k].grid(alpha=0.3); ax[k].legend()
        else:
            ax[k].plot(V, key["v15"], "o-", color="0.3", label="v15 chord")
            ax[k].plot(V, key["cad"], "D--", color="#d62728", mfc="none", label="CAD chord")
            ax[k].set(xlabel="wind speed [m/s]", ylabel=f"{lab} [{un}]", title=lab)
            ax[k].grid(alpha=0.3); ax[k].legend()
    fig.suptitle("Effect of reconciling the outboard chord to CAD/v14 (OpenFAST 3.5.5, "
                 "rigid, fixed 13.2 rpm)", fontsize=12, fontweight="bold")
    fig.tight_layout()
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "..", "figures", "cadchord_load_effect.png")
    fig.savefig(out, dpi=145, bbox_inches="tight")
    print("\nwrote", os.path.relpath(out, os.path.dirname(here)))

    # spanwise normal force in the reconciled region at a few speeds
    _, du = u.read_blade("/tmp/sweep/uniform_11/AeroDyn15/AeroDyn15_blade_uniform.dat")
    rr = (du[:, 0] + HUB) / TIP
    reg = (rr >= 0.65) & (rr <= 0.97)
    print("\nspanwise Fn change in reconciled region (r/R 0.65-0.97):")
    for V in [8, 14, 20]:
        nb, _, B = c.read_out(f"/tmp/sweep/uniform_{V}/Turbine_driver.out")
        na, _, A = c.read_out(f"/tmp/sweep/cadchord_{V}/Turbine_driver.out")
        fb = c.spanwise(nb, B, "Fn", len(du)); fa = c.spanwise(na, A, "Fn", len(du))
        d = (fa - fb)[reg]
        print(f"  {V} m/s: mean {d.mean():+.1f} N/m, max {np.abs(d).max():.1f} N/m "
              f"({100*np.abs(d).max()/np.abs(fb).max():.1f}% of peak Fn)")


if __name__ == "__main__":
    main()
