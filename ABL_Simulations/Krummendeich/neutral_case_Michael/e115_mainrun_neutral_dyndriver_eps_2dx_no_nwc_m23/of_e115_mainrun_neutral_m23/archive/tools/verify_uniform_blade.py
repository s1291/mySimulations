"""INDEPENDENT verification of the committed AeroDyn15_uniform model.
Re-implements the expected blade resampling and polar blend from scratch (does
NOT import the generator), then checks the committed files against it. Writes a
PASS/FAIL report to /tmp/verify_report.txt."""
import os, re, glob
import numpy as np

import os
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORIG = os.path.join(REPO, "AeroDyn15")
UNI = os.path.join(REPO, "AeroDyn15_uniform")
HUBRAD = 1.51
R = []                      # report lines
def rep(tag, ok, msg): R.append((tag, ok, msg))

# ---------- parsers (independent) ----------
def read_blade(path):
    raw = [l.rstrip("\r\n") for l in open(path, errors="replace")]
    hdr = next(i for i, l in enumerate(raw) if "BlSpn" in l and "BlChord" in l)
    rows = []
    for l in raw[hdr+2:]:
        p = l.split()
        if len(p) < 7: continue
        try: vals = [float(x) for x in p[:6]] + [int(float(p[6]))]
        except ValueError: continue
        rows.append(vals)
    return np.array([r[:6] for r in rows]), np.array([r[6] for r in rows], int)

def parse_polar(path):
    raw = [l.rstrip("\r\n") for l in open(path, errors="replace")]
    scal = {}
    for l in raw:
        m = re.match(r"\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+([A-Za-z]\w*)\b", l)
        if m: scal.setdefault(m.group(2), float(m.group(1)))
    na = int(scal["NumAlf"])
    # find table: NumAlf rows of >=4 floats
    tbl, started = [], False
    for l in raw:
        p = l.split()
        if len(p) >= 4:
            try: v = [float(x) for x in p[:4]]
            except ValueError:
                if started: break
                continue
            # skip lines that are scalar+name (2 tokens) already excluded by >=4
            started = True; tbl.append(v)
            if len(tbl) == na: break
    return scal, np.array(tbl), raw

# ---------- load original ----------
ob_geo, ob_af = read_blade(os.path.join(ORIG, "AeroDyn15_blade.dat"))
span0 = ob_geo[:, 0]
n_orig_af = ob_af.max()
opol = {a: parse_polar(os.path.join(ORIG, "Airfoils", f"Polars{a:02d}.dat")) for a in range(1, n_orig_af+1)}

# ---------- load committed uniform ----------
ub_geo, ub_af = read_blade(os.path.join(UNI, "AeroDyn15_blade_uniform.dat"))
n = len(ub_geo)

# CHECK 1: node count & uniform spacing
su = ub_geo[:, 0]
d = np.diff(su)
rep("nodes==57", n == 57, f"n={n}")
rep("span starts 0", abs(su[0]) < 1e-9, f"su0={su[0]:.6f}")
rep("span ends at original tip", abs(su[-1]-span0[-1]) < 1e-6, f"su_end={su[-1]:.4f} vs {span0[-1]:.4f}")
rep("strictly increasing", np.all(d > 0), f"min d={d.min():.4f}")
# spacing uniform to within the file's 6-sig-fig print rounding (~1e-5 m); compare
# each node to the ideal uniform grid rather than diff-of-diffs.
ideal = np.linspace(su[0], su[-1], n)
rep("nodes on ideal uniform grid (<1e-4 m)", np.abs(su-ideal).max() < 1e-4,
    f"max dev from ideal {np.abs(su-ideal).max():.2e} m (file is 6 sig figs)")
rep("spacing >= 1 m floor", d.mean() >= 1.0, f"spacing={d.mean():.4f} m")

# CHECK 2: BlAFID = 1..57 (each node its own airfoil)
rep("BlAFID == 1..n", np.array_equal(ub_af, np.arange(1, n+1)), f"afid {ub_af.min()}..{ub_af.max()}")

# CHECK 3: geometry = resample of original continuous blade (within tight tol)
# chord(1),twist? columns: BlSpn0 BlCrvAC1 BlSwpAC2 BlCrvAng3 BlTwist4 BlChord5
gmax = {}
for col, nm in [(5, "chord"), (4, "twist"), (1, "prebendAC"), (2, "sweepAC"), (3, "curveAng")]:
    exp = np.interp(su, span0, ob_geo[:, col])
    gmax[nm] = np.abs(ub_geo[:, col] - exp).max()
rep("chord = interp(original)", gmax["chord"] < 1e-4, f"max dev {gmax['chord']:.2e} m")
rep("twist = interp(original)", gmax["twist"] < 1e-4, f"max dev {gmax['twist']:.2e} deg")
rep("prebend = interp(original)", gmax["prebendAC"] < 1e-4, f"max dev {gmax['prebendAC']:.2e} m")
rep("sweep = interp(original)", gmax["sweepAC"] < 1e-4, f"max dev {gmax['sweepAC']:.2e} m")
rep("no NaN/Inf in geometry", np.isfinite(ub_geo).all(), "")

# CHECK 4: polar files
upol_files = sorted(glob.glob(os.path.join(UNI, "Airfoils", "Polars*.dat")))
rep("57 polar files", len(upol_files) == 57, f"found {len(upol_files)}")

UA8 = {"alpha0","alpha1","alpha2","C_nalpha","Cn1","Cn2","Cd0","Cm0"}
alpha_ref = opol[2][1][:, 0]   # reference alpha grid from an airfoil station
# verify all original share one alpha grid
same_grid = all(np.array_equal(opol[a][1][:,0], alpha_ref) for a in range(2, n_orig_af+1))
rep("all original polars share alpha grid", same_grid, f"{len(alpha_ref)} pts")

# independent expected blend for each uniform node
afid0 = ob_af
max_tab_err = 0.0
ua_present_count = 0
bad_grid = 0
nonfinite = 0
numalf_bad = 0
worst_node = None
for k in range(n):
    s = su[k]
    scal, tab, raw = parse_polar(upol_files[k])
    # alpha grid intact
    if tab.shape[0] != 175: numalf_bad += 1
    if not np.array_equal(tab[:, 0], alpha_ref): bad_grid += 1
    if not np.isfinite(tab).all(): nonfinite += 1
    # UA lines should be ABSENT
    for kw in UA8:
        if kw in scal: ua_present_count += 1
    # expected blend: find bracketing original airfoils by span
    j = int(np.clip(np.searchsorted(span0, s) - 1, 0, len(span0)-2))
    a0, a1 = afid0[j], afid0[j+1]
    w = (s - span0[j])/(span0[j+1]-span0[j]) if span0[j+1] > span0[j] else 0.0
    if a0 == a1:
        exp_tab = opol[a0][1][:, 1:4]
    else:
        exp_tab = (1-w)*opol[a0][1][:, 1:4] + w*opol[a1][1][:, 1:4]
    err = np.abs(tab[:, 1:4] - exp_tab).max()
    if err > max_tab_err: max_tab_err, worst_node = err, k+1

rep("every polar has 175 alpha rows", numalf_bad == 0, f"{numalf_bad} bad")
rep("alpha grid intact in all polars", bad_grid == 0, f"{bad_grid} differ")
rep("no NaN/Inf in any polar table", nonfinite == 0, f"{nonfinite} bad")
rep("8 UA lines absent in all polars", ua_present_count == 0, f"{ua_present_count} present")
rep("blended table == independent linear blend", max_tab_err < 1e-5,
    f"max |diff| {max_tab_err:.2e} at node {worst_node}")

# CHECK 5: physical sanity of Cl/Cd over operating AoA (-10..20) for all nodes
clmax=cdmin=99; bad_phys=0
for f in upol_files:
    _, tab, _ = parse_polar(f)
    m = (tab[:,0]>=-10)&(tab[:,0]<=20)
    if tab[m,2].min() < 0: bad_phys += 1     # Cd negative in operating range
rep("Cd >= 0 in AoA[-10,20] all nodes", bad_phys == 0, f"{bad_phys} nodes with negative Cd")

# write report
ok_all = all(ok for _, ok, _ in R)
with open("/tmp/verify_report.txt", "w") as f:
    f.write("INDEPENDENT VERIFICATION OF AeroDyn15_uniform\n"+"="*60+"\n")
    for tag, ok, msg in R:
        f.write(f"[{'PASS' if ok else 'FAIL'}] {tag:42} {msg}\n")
    f.write("="*60+f"\nOVERALL: {'ALL PASS' if ok_all else 'SOME FAILED'}\n")
print("ALL PASS" if ok_all else "SOME FAILED")
