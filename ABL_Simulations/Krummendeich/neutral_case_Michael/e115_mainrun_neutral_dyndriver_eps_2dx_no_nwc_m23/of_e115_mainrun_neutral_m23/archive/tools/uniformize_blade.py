#!/usr/bin/env python3
"""Re-discretize the AeroDyn blade onto a UNIFORM spanwise grid while preserving
the aerodynamics, by blending the airfoil polars along the span.

The blade is a continuously blended-airfoil design sampled at 39 stations. A
uniform grid is a different sampling of the same continuous blade, provided we
also reconstruct the airfoil data at the new stations. At each uniform node we
therefore blend the two original polars that bracket it in span (weight = span
fraction): the Cl/Cd/Cm table row by row (all 39 files share one 175-point
alpha grid) and every Beddoes-Leishman UA scalar, since the case runs
AFAeroMod=2. The blend is exact at the 39 original stations (weight 0 or 1) and
linear between them, so the continuous airfoil-vs-span is preserved and only the
sampling becomes uniform.

Outputs (into <out>/): AeroDyn15.dat, AeroDyn15_blade_uniform.dat, and
Airfoils/Polars01..NN.dat (one blended table per uniform node).
"""
from __future__ import annotations
import os
import re
import shutil
import numpy as np

CASE = "/tmp/the_case"
NUM = re.compile(r"[+-]?\d+\.?\d*(?:[eE][+-]?\d+)?")

# The 8 Beddoes-Leishman UA parameters that OpenFAST derives FROM the static
# polar (alpha0 = zero-lift crossing, C_nalpha = smoothed Cn slope, alpha1/2 =
# separation-function f=0.7, Cn1/2 = Cn at alpha1/2, Cd0 = min Cd, Cm0 =
# Cm(alpha0)); see AirfoilInfo.f90:CalculateUACoeffs and AeroDyn theory_ua docs.
# Because they are functions of the polar, the correct way to blend airfoils is
# to blend the static table and then RE-DERIVE these from it -- which is exactly
# what WISDEM/WEIS do (Polar.unsteadyParams) and what OpenFAST itself does when
# the field is OMITTED from the file. We therefore drop these lines so OpenFAST
# recomputes them from each blended table (audited in <root>.UA.sum). The
# remaining UA fields (T_f0, b1, A1, St_sh, ... ) are airfoil-INDEPENDENT
# literature defaults already written as "Default" in the source files; they are
# left untouched. We never linearly interpolate UA scalars.
UA_FROM_POLAR = {"alpha0", "alpha1", "alpha2", "C_nalpha",
                 "Cn1", "Cn2", "Cd0", "Cm0"}
KEEP_SCALAR = {"NumAlf", "NumTabs", "NumCoords", "Ctrl"}  # counts/flags, do not blend


# ----------------------------- blade file IO ------------------------------- #
def read_blade(path):
    raw = [l.rstrip("\r\n") for l in open(path, encoding="utf-8", errors="replace")]
    hdr = next(i for i, l in enumerate(raw) if "BlSpn" in l and "BlChord" in l)
    rows = []
    for l in raw[hdr + 2:]:
        p = l.split()
        if len(p) < 7:
            continue
        try:
            v = [float(x) for x in p[:6]]
        except ValueError:
            continue
        extra = [float(x) for x in p[7:10]] if len(p) >= 10 else [0.0, 0.0, 0.0]
        rows.append(v + [float(int(float(p[6])))] + extra)
    return raw[:hdr + 2], np.array(rows)


def write_blade(path, head, data, note):
    lines = list(head)
    for i, l in enumerate(lines):
        if "NumBlNds" in l:
            lines[i] = re.sub(r"^\s*\d+", f"{len(data):>11d}", l, count=1)
    if len(lines) > 1:
        lines[1] = note
    for r in data:
        lines.append(f"{r[0]: .6E} {r[1]: .6E} {r[2]: .6E} {r[3]: .6E} "
                     f"{r[4]: .6E} {r[5]: .6E} {int(round(r[6])):6d} "
                     f"{r[7]:8.1f} {r[8]:8.1f} {r[9]:8.1f}")
    open(path, "w", newline="").write("\r\n".join(lines) + "\r\n")


# ----------------------------- polar IO/blend ------------------------------ #
def parse_polar(path):
    raw = [l.rstrip("\r\n") for l in open(path, encoding="utf-8", errors="replace")]
    scal = {}
    for i, l in enumerate(raw):
        m = re.match(r"\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+([A-Za-z]\w*)\b", l)
        if m:
            scal[m.group(2)] = (i, float(m.group(1)))
    ni = scal["NumAlf"][0]
    na = int(scal["NumAlf"][1])
    tstart, tbl = None, []
    for i in range(ni + 1, len(raw)):
        p = raw[i].split()
        if len(p) >= 4:
            try:
                vals = [float(x) for x in p[:4]]
            except ValueError:
                if tstart is not None:
                    break
                continue
            if tstart is None:
                tstart = i
            tbl.append(vals)
            if len(tbl) == na:
                break
    return {"lines": raw, "scalars": scal, "tstart": tstart, "table": np.array(tbl)}


def _sub_first_num(line, val):
    m = re.match(r"(\s*)(\S+)(.*)", line)
    return f"{m.group(1)}{val:.6f}{m.group(3)}"


def drop_ua_from_polar(lines):
    """Remove the 8 polar-derived UA lines so OpenFAST recomputes them from the
    static table. No-op for the cylinder polar (it carries no UA block)."""
    out = []
    for l in lines:
        m = re.match(r"\s*\S+\s+([A-Za-z]\w*)\b", l)
        if m and m.group(1) in UA_FROM_POLAR:
            continue
        out.append(l)
    return out


def blend_polar(A, B, w):
    """Blend two polars at span fraction w (0 -> A, 1 -> B).

    The static Cl/Cd/Cm table is linearly interpolated at matched angle of
    attack -- the standard NREL method (AirfoilPrep.Polar.blend / welib /
    WISDEM: coef = (1-w)*A + w*B). The 8 polar-derived UA scalars are NOT
    interpolated; their lines are dropped so OpenFAST re-derives them from this
    blended table (CalculateUACoeffs). All other lines pass through unchanged.
    """
    Ta, Tb = A["table"], B["table"]
    if not np.allclose(Ta[:, 0], Tb[:, 0]):
        raise ValueError("alpha grids differ; cannot blend row-wise")
    out = []
    for idx, line in enumerate(A["lines"]):
        m = re.match(r"\s*\S+\s+([A-Za-z]\w*)\b", line)
        if m and m.group(1) in UA_FROM_POLAR:
            continue                       # omit -> OpenFAST recomputes from blended polar
        if A["tstart"] <= idx < A["tstart"] + len(Ta):
            k = idx - A["tstart"]
            a = Ta[k, 0]
            cl, cd, cm = (1 - w) * Ta[k, 1:] + w * Tb[k, 1:]
            out.append(f"{a:14.6f} {cl:13.6f} {cd:13.6f} {cm:13.6f}")
        else:
            out.append(line)
    return out


# --------------------------- AeroDyn primary file -------------------------- #
def make_primary(src, dst, n):
    raw = [l.rstrip("\r\n") for l in open(src, encoding="utf-8", errors="replace")]
    out, i = [], 0
    n_orig = None
    while i < len(raw):
        l = raw[i]
        if "NumAFfiles" in l:
            n_orig = int(re.match(r"\s*(\d+)", l).group(1))
            out.append(re.sub(r"^\s*\d+", f"{n:>11d}", l, count=1))
            i += 1
        elif "AFNames" in l:
            out.append(f'"Airfoils/Polars01.dat"   AFNames            - '
                       f"Airfoil file names ({n} lines) (quoted strings)")
            for k in range(2, n + 1):
                out.append(f'"Airfoils/Polars{k:02d}.dat"')
            i += n_orig  # skip the original AFNames block
        elif "ADBlFile" in l:
            out.append(re.sub(r'"[^"]*"', '"AeroDyn15_blade_uniform.dat"', l, count=1))
            i += 1
        else:
            out.append(l)
            i += 1
    open(dst, "w", newline="").write("\r\n".join(out) + "\r\n")


# --------------------------------- build ----------------------------------- #
def build(case, out, n):
    ad = os.path.join(case, "AeroDyn15")
    head, data = read_blade(os.path.join(ad, "AeroDyn15_blade.dat"))
    span0, afid0 = data[:, 0], data[:, 6].astype(int)
    L = span0[-1]
    su = np.linspace(0.0, L, n)

    new = np.zeros((n, data.shape[1]))
    new[:, 0] = su
    for c in (1, 2, 3, 4, 5, 7, 8, 9):
        new[:, c] = np.interp(su, span0, data[:, c])
    new[:, 6] = np.arange(1, n + 1)

    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(os.path.join(out, "Airfoils"))
    write_blade(os.path.join(out, "AeroDyn15_blade_uniform.dat"), head, new,
                f"Uniform {L/(n-1):.4f} m spacing, {n} nodes, blended polars "
                f"(aerodynamics preserved) by uniformize_blade.py")

    pol = {a: parse_polar(os.path.join(ad, "Airfoils", f"Polars{a:02d}.dat"))
           for a in range(1, afid0.max() + 1)}
    weights = []
    for k, s in enumerate(su):
        j = int(np.clip(np.searchsorted(span0, s) - 1, 0, len(span0) - 2))
        a0, a1 = afid0[j], afid0[j + 1]
        w = (s - span0[j]) / (span0[j + 1] - span0[j]) if span0[j + 1] > span0[j] else 0.0
        lines = drop_ua_from_polar(pol[a0]["lines"]) if a0 == a1 \
            else blend_polar(pol[a0], pol[a1], w)
        open(os.path.join(out, "Airfoils", f"Polars{k+1:02d}.dat"),
             "w", newline="").write("\r\n".join(lines) + "\r\n")
        weights.append((s, a0, a1, w))
    make_primary(os.path.join(ad, "AeroDyn15.dat"),
                 os.path.join(out, "AeroDyn15.dat"), n)
    return head, data, new, pol, weights


def verify(case, data, new, pol):
    span0, afid0 = data[:, 0], data[:, 6].astype(int)
    # geometry: uniform polyline vs original polyline on a dense grid
    sd = np.linspace(0, span0[-1], 2000)
    gmax = 0.0
    for c, nm in [(5, "chord"), (4, "twist"), (1, "prebend")]:
        d = np.abs(np.interp(sd, new[:, 0], new[:, c]) - np.interp(sd, span0, data[:, c]))
        gmax = max(gmax, d.max())
    # static Cl/Cd/Cm blend is exact at original stations: blend at each
    # original span and compare the table to that station's own table.
    emax = 0.0
    for i in range(1, len(span0) - 1):
        s = span0[i]
        j = int(np.clip(np.searchsorted(span0, s) - 1, 0, len(span0) - 2))
        a0, a1 = afid0[j], afid0[j + 1]
        w = (s - span0[j]) / (span0[j + 1] - span0[j])
        lines = drop_ua_from_polar(pol[a0]["lines"]) if a0 == a1 \
            else blend_polar(pol[a0], pol[a1], w)
        tb = parse_polar_lines(lines)["table"]
        emax = max(emax, np.abs(tb - pol[afid0[i]]["table"]).max())
    return gmax, emax


def parse_polar_lines(lines):
    import tempfile
    f = tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False, newline="")
    f.write("\r\n".join(lines)); f.close()
    return parse_polar(f.name)


def main():
    import sys
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 57
    out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/uniform_ad"
    head, data, new, pol, weights = build(CASE, out, n)
    L = data[-1, 0]
    print(f"uniform grid: {n} nodes, spacing {L/(n-1):.4f} m (>= 1 m: {L/(n-1)>=1.0})")
    print(f"blended {n} polar files into {out}/Airfoils")
    gmax, emax = verify(CASE, data, new, pol)
    print(f"geometry max deviation (uniform vs original polyline): {gmax:.2e}")
    print(f"polar blend error at the 39 original stations:        {emax:.2e}")
    print("OK" if (gmax < 0.02 and emax < 1e-6) else "CHECK")


if __name__ == "__main__":
    main()
