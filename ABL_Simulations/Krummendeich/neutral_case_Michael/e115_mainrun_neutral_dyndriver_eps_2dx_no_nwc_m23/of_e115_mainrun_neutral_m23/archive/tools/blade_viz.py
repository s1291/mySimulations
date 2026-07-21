#!/usr/bin/env python3
"""Visualize an AeroDyn15 wind-turbine blade with each airfoil ID mapped to its
spanwise blade section.

It parses the model files directly so the figures stay in sync with the data:

  AeroDyn15/AeroDyn15_blade.dat    per-node span, chord, twist, prebend, sweep, BlAFID
  AeroDyn15/AeroDyn15.dat          ordered list of Polars##.dat airfoil files
  AeroDyn15/Airfoils/Polars##.dat  airfoil names -> blend family -> derived t/c

and writes three views into ``figures/``:

  blade_overview.png      static multi-panel matplotlib figure
  blade_interactive.html  Plotly 2D explorer (hover a section for its airfoil)
  blade_3d.html           Plotly 3D lofted blade colored by airfoil ID

NOTE: the polar files carry no shape coordinates (``NumCoords 0``). The 3D loft
therefore uses a *generic* NACA-type thickness form scaled to a relative
thickness (t/c) that is *derived from the airfoil naming convention* (see
``BASE_TC``). It is illustrative of the planform/twist/taper, not exact section
geometry.
"""
from __future__ import annotations

import argparse
import os
import re

import numpy as np

# --------------------------------------------------------------------------- #
# Airfoil family -> relative thickness t/c [%], inferred from the EC1xx naming
# convention used in this model (Cylinder = 100%, EC1-45 = 45%, ... EC1-16 = 16%).
# Override here if the convention differs for your airfoils.
# --------------------------------------------------------------------------- #
BASE_TC = {
    "Cylinder": 100.0,
    "EC145": 45.0,
    "EC135": 35.0,
    "EC128": 28.0,
    "EC122": 22.0,
    "EC116": 16.0,
}

PITCH_AXIS = 0.375  # chordwise fraction (LE=0, TE=1) used as the loft reference


def repo_root() -> str:
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
def parse_blade(path: str) -> dict:
    """Read the AeroDyn blade definition table into arrays (handles CRLF)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [ln.rstrip("\r\n") for ln in fh]
    hdr = next(i for i, ln in enumerate(lines) if "BlSpn" in ln and "BlChord" in ln)
    rows = []
    for ln in lines[hdr + 1:]:
        parts = ln.split()
        if len(parts) < 7:
            continue
        try:
            rows.append([float(p) for p in parts[:7]])
        except ValueError:
            continue  # skips the "(m)  (m) ..." units line
    arr = np.array(rows)
    return {
        "span": arr[:, 0],
        "prebend": arr[:, 1],   # BlCrvAC, out-of-plane
        "sweep": arr[:, 2],     # BlSwpAC, in-plane
        "twist": arr[:, 4],     # deg
        "chord": arr[:, 5],
        "afid": arr[:, 6].astype(int),
    }


def parse_airfoil_files(aerodyn_path: str) -> list[str]:
    """Ordered list of Polars##.dat paths (index 0 == AFID 1)."""
    names: list[str] = []
    with open(aerodyn_path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            m = re.search(r'"([^"]*Polars\d+\.dat)"', ln)
            if m:
                names.append(m.group(1).replace("\\", "/"))
    return names


def airfoil_label(path: str) -> str:
    """The human-readable airfoil name comment near the top of a polar file."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [ln.strip() for ln in fh]
    try:
        idx = next(i for i, ln in enumerate(lines) if "NumTabs" in ln)
    except StopIteration:
        idx = 0
    for ln in lines[idx + 1:]:
        if ln.startswith("!"):
            body = ln.lstrip("!").strip()
            if body and any(c.isalnum() for c in body):
                return body
    return os.path.basename(path)


def first_reynolds(path: str) -> float | None:
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if "Reynolds number" in ln:
                try:
                    return float(ln.split()[0])
                except (ValueError, IndexError):
                    return None
    return None


def _base_tc(token: str) -> float | None:
    return BASE_TC.get(token.strip().rstrip("Ff"))


def derive_tc(label: str) -> float | None:
    """Relative thickness t/c [%] from a blended-airfoil name.

    Blend names look like ``0.96769_Cylinder_EC145F_...`` (fraction toward the
    second family); pure names look like ``Cylinder__1_...`` or ``EC116__...``.
    """
    toks = [t for t in label.split("_") if t]
    if not toks:
        return None
    try:
        frac = float(toks[0])
    except ValueError:
        return _base_tc(toks[0])  # pure section
    a = _base_tc(toks[1]) if len(toks) > 1 else None
    b = _base_tc(toks[2]) if len(toks) > 2 else None
    if a is not None and b is not None:
        return (1.0 - frac) * a + frac * b
    return a


def load_model(blade_path: str, aerodyn_path: str, airfoils_dir: str) -> dict:
    blade = parse_blade(blade_path)
    files = parse_airfoil_files(aerodyn_path)
    base = os.path.dirname(aerodyn_path)
    labels, tcs, res = [], [], []
    for rel in files:
        p = os.path.join(base, rel)
        if not os.path.exists(p):  # fall back to the airfoils dir
            p = os.path.join(airfoils_dir, os.path.basename(rel))
        lbl = airfoil_label(p)
        labels.append(lbl)
        tcs.append(derive_tc(lbl))
        res.append(first_reynolds(p))
    blade["af_label"] = labels          # by AFID-1
    blade["af_tc"] = tcs
    blade["af_re"] = res
    # per-node convenience arrays
    blade["tc"] = np.array(
        [tcs[i - 1] if tcs[i - 1] is not None else np.nan for i in blade["afid"]]
    )
    blade["label"] = [labels[i - 1] for i in blade["afid"]]
    blade["re"] = [res[i - 1] for i in blade["afid"]]
    return blade


# --------------------------------------------------------------------------- #
# Generic geometry (illustrative only)
# --------------------------------------------------------------------------- #
def naca_profile(tc_pct: float, n: int = 60):
    """Closed symmetric NACA-type section of relative thickness ``tc_pct`` [%]."""
    t = max(tc_pct, 1.0) / 100.0
    beta = np.linspace(0.0, np.pi, n)
    x = (1.0 - np.cos(beta)) / 2.0
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                  + 0.2843 * x**3 - 0.1036 * x**4)
    xx = np.concatenate([x, x[::-1]])
    yy = np.concatenate([yt, -yt[::-1]])
    return xx, yy


def section_3d(px, py, chord, twist_deg, span, prebend, sweep, xp=PITCH_AXIS):
    th = np.radians(twist_deg)
    u = (px - xp) * chord
    v = py * chord
    xe = u * np.cos(th) - v * np.sin(th) + sweep
    yf = u * np.sin(th) + v * np.cos(th) + prebend
    z = np.full_like(xe, span)
    return xe, yf, z


# --------------------------------------------------------------------------- #
# Shared colour mapping (same Turbo scale across all three views)
# --------------------------------------------------------------------------- #
def color_setup(afid):
    import matplotlib as mpl

    cmin, cmax = int(afid.min()), int(afid.max())
    norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
    return mpl.colormaps["turbo"], norm, cmin, cmax


# --------------------------------------------------------------------------- #
# View 1: static multi-panel matplotlib
# --------------------------------------------------------------------------- #
def _annotate(ax, xs, ys, texts, dy=5, fs=6.5, color="0.15", rotation=0):
    """Place a small text label above each finite (x, y) point."""
    for x, y, t in zip(xs, ys, texts):
        if np.isfinite(y):
            ax.annotate(t, (x, y), textcoords="offset points", xytext=(0, dy),
                        ha="center", va="bottom", fontsize=fs, color=color,
                        rotation=rotation, clip_on=True)


def plot_static(m: dict, out_path: str) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable

    plt.rcParams.update({"font.size": 12})
    span, chord, afid = m["span"], m["chord"], m["afid"]
    cmap, norm, cmin, cmax = color_setup(afid)

    fig, axes = plt.subplots(
        5, 1, figsize=(26, 15), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.1, 1.1, 1.1, 1.1], "hspace": 0.13},
    )
    ax_pf, ax_id, ax_sp, ax_tc, ax_tw = axes

    # --- planform coloured by airfoil ID ---
    y_le = PITCH_AXIS * chord
    y_te = -(1.0 - PITCH_AXIS) * chord
    for i in range(len(span) - 1):
        xs = [span[i], span[i + 1], span[i + 1], span[i]]
        ys = [y_le[i], y_le[i + 1], y_te[i + 1], y_te[i]]
        ax_pf.fill(xs, ys, color=cmap(norm(afid[i])), edgecolor="none", zorder=1)
    ax_pf.plot(span, y_le, color="k", lw=1.0, zorder=3)
    ax_pf.plot(span, y_te, color="k", lw=1.0, zorder=3)
    ax_pf.axhline(0, color="0.4", lw=0.8, ls="--", zorder=2)
    # vertical guide at every airfoil-ID transition
    for i in range(1, len(afid)):
        if afid[i] != afid[i - 1]:
            ax_pf.axvline(span[i], color="0.5", lw=0.5, ls=":", zorder=2)
    # BlAFID number above every node
    for i in range(len(span)):
        ax_pf.text(span[i], y_le[i], f"{afid[i]}", va="bottom", ha="center",
                   fontsize=7, rotation=90, zorder=4)
    ax_pf.set_ylim(y_te.min() - 0.25, y_le.max() + 0.8)
    ax_pf.set_ylabel("chordwise extent  [m]\n(about pitch axis)")
    ax_pf.set_title("Blade planform coloured by airfoil ID  "
                    "(number above each node = BlAFID)")

    # --- airfoil ID vs span (step) ---
    ax_id.step(span, afid, where="post", color="0.25", lw=1.2)
    ax_id.scatter(span, afid, c=afid, cmap=cmap, norm=norm, s=28,
                  edgecolor="k", lw=0.4, zorder=3)
    _annotate(ax_id, span, afid, [str(a) for a in afid], dy=5, fs=8)
    ax_id.set_ylabel("airfoil ID\n(BlAFID)")
    ax_id.set_ylim(cmin - 3, cmax + 5)
    ax_id.grid(True, alpha=0.3)

    # --- distance between blade nodes (Δspan) ---
    # A gap spans two adjacent airfoils, so it is labelled with the pair
    # (inboard->outboard) and coloured neutrally, not by a single airfoil ID.
    mid = (span[:-1] + span[1:]) / 2.0
    dspan = np.diff(span)
    pair = [f"{afid[i]}→{afid[i + 1]}" for i in range(len(span) - 1)]
    ax_sp.vlines(mid, 0, dspan, color="0.7", lw=1.2, zorder=1)
    ax_sp.scatter(mid, dspan, color="#4c78a8", s=22,
                  edgecolor="k", lw=0.3, zorder=3)
    ax_sp.plot(span, np.zeros_like(span), marker="|", ls="none", color="k",
               ms=9, mew=1.0, zorder=2)  # node positions (rug)
    _annotate(ax_sp, mid, dspan, pair, dy=4, fs=7, color="#33475b",
              rotation=90)
    ax_sp.axhline(dspan.mean(), color="0.4", ls="--", lw=0.8,
                  label=f"mean {dspan.mean():.2f} m  (max {dspan.max():.2f} m)")
    ax_sp.set_ylabel("node spacing\nΔspan  [m]\n(between two airfoils)")
    ax_sp.set_ylim(0, dspan.max() * 1.5)
    ax_sp.legend(loc="upper left", fontsize=7, frameon=False)
    ax_sp.grid(True, alpha=0.3)

    # --- derived relative thickness ---
    ax_tc.plot(span, m["tc"], "-o", color="#1f77b4", ms=4)
    _annotate(ax_tc, span, m["tc"], [str(a) for a in afid], dy=5, fs=7,
              color="#16466b")
    ax_tc.set_ylabel("relative\nthickness t/c  [%]")
    ax_tc.margins(y=0.18)
    ax_tc.grid(True, alpha=0.3)

    # --- twist ---
    ax_tw.plot(span, m["twist"], "-o", color="#d62728", ms=4)
    _annotate(ax_tw, span, m["twist"], [str(a) for a in afid], dy=5, fs=7,
              color="#8b1a1a")
    ax_tw.set_ylabel("twist  [deg]")
    ax_tw.set_xlabel("blade span  BlSpn  [m]")
    ax_tw.margins(y=0.18)
    ax_tw.grid(True, alpha=0.3)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.01,
                        ticks=range(cmin, cmax + 1, 2))
    cbar.set_label("airfoil ID")

    fig.suptitle("AeroDyn15 blade: airfoil ID per blade section", y=0.995,
                 fontsize=16, fontweight="bold")
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path


# --------------------------------------------------------------------------- #
# View 2: interactive 2D Plotly explorer
# --------------------------------------------------------------------------- #
def _hover(m, i):
    re = m["re"][i]
    re_s = f"{re:g} M" if re is not None else "n/a"
    tc = m["tc"][i]
    tc_s = f"{tc:.1f}%" if np.isfinite(tc) else "n/a"
    span = m["span"]
    if i < len(span) - 1:
        gap = f"Δ to next node: {span[i + 1] - span[i]:.3f} m"
    else:
        gap = f"Δ from prev node: {span[i] - span[i - 1]:.3f} m"
    return (f"<b>node {i + 1}</b>  (span {span[i]:.2f} m)<br>"
            f"airfoil ID: <b>{m['afid'][i]}</b><br>"
            f"name: {m['label'][i]}<br>"
            f"t/c: {tc_s} &nbsp; Re: {re_s}<br>"
            f"chord: {m['chord'][i]:.3f} m &nbsp; twist: {m['twist'][i]:.2f}°<br>"
            f"{gap}")


def plot_interactive(m: dict, out_path: str) -> str:
    import matplotlib as mpl
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go

    span, chord, afid = m["span"], m["chord"], m["afid"]
    cmap, norm, cmin, cmax = color_setup(afid)

    def hexcol(v):
        return mpl.colors.to_hex(cmap(norm(v)))

    fig = make_subplots(
        rows=5, cols=1, shared_xaxes=True, vertical_spacing=0.035,
        row_heights=[0.38, 0.17, 0.17, 0.14, 0.14],
        subplot_titles=("Planform coloured by airfoil ID (hover a node)",
                        "Airfoil ID vs span",
                        "Distance between blade nodes (Δspan) [m]",
                        "Relative thickness t/c [%]", "Twist [deg]"),
    )

    y_le = PITCH_AXIS * chord
    y_te = -(1.0 - PITCH_AXIS) * chord
    for i in range(len(span) - 1):
        fig.add_trace(go.Scatter(
            x=[span[i], span[i + 1], span[i + 1], span[i], span[i]],
            y=[y_le[i], y_le[i + 1], y_te[i + 1], y_te[i], y_le[i]],
            fill="toself", mode="lines", line=dict(width=0),
            fillcolor=hexcol(afid[i]), hoverinfo="skip",
            showlegend=False), row=1, col=1)
    fig.add_trace(go.Scatter(x=span, y=y_le, mode="lines",
                             line=dict(color="black", width=1.2),
                             hoverinfo="skip", showlegend=False), row=1, col=1)
    fig.add_trace(go.Scatter(x=span, y=y_te, mode="lines",
                             line=dict(color="black", width=1.2),
                             hoverinfo="skip", showlegend=False), row=1, col=1)
    # hover markers + BlAFID labels carrying the per-section info and colorbar
    ids = [str(a) for a in afid]
    hov = [_hover(m, i) for i in range(len(span))]
    fig.add_trace(go.Scatter(
        x=span, y=np.zeros_like(span), mode="markers+text",
        marker=dict(size=9, color=afid, colorscale="Turbo", cmin=cmin, cmax=cmax,
                    line=dict(color="black", width=0.5),
                    colorbar=dict(title="airfoil ID", len=0.42, y=0.8)),
        text=ids, textposition="top center", textfont=dict(size=11),
        hovertext=hov, hoverinfo="text", showlegend=False), row=1, col=1)

    fig.add_trace(go.Scatter(x=span, y=afid, mode="lines+markers+text",
                             line=dict(color="#444", shape="hv"),
                             marker=dict(size=6, color=afid, colorscale="Turbo",
                                         cmin=cmin, cmax=cmax, showscale=False),
                             text=ids, textposition="top center",
                             textfont=dict(size=12),
                             hovertext=hov, hoverinfo="text",
                             showlegend=False), row=2, col=1)
    # node spacing: bar height = gap to next node, labelled with the adjacent
    # airfoil pair (a gap belongs to two airfoils, so it is coloured neutrally),
    # plus a rug of node positions
    mid = (span[:-1] + span[1:]) / 2.0
    dspan = np.diff(span)
    fig.add_trace(go.Bar(
        x=mid, y=dspan, width=np.minimum(dspan * 0.9, 0.9),
        marker=dict(color="#4c78a8"),
        text=[f"{afid[i]}→{afid[i + 1]}" for i in range(len(span) - 1)],
        textposition="outside", textfont=dict(size=10),
        customdata=np.stack([span[:-1] + 1, span[1:] + 1,
                             afid[:-1], afid[1:]], axis=-1),
        hovertemplate="nodes %{customdata[0]:.0f}->%{customdata[1]:.0f}<br>"
                      "airfoils %{customdata[2]:.0f} and %{customdata[3]:.0f}<br>"
                      "Δspan = %{y:.3f} m<extra></extra>",
        showlegend=False), row=3, col=1)
    fig.add_trace(go.Scatter(
        x=span, y=np.zeros_like(span), mode="markers",
        marker=dict(symbol="line-ns-open", size=8,
                    line=dict(color="black", width=1)),
        hoverinfo="skip", showlegend=False), row=3, col=1)
    fig.add_trace(go.Scatter(x=span, y=m["tc"], mode="lines+markers+text",
                             text=ids, textposition="top center",
                             textfont=dict(size=11),
                             line=dict(color="#1f77b4"), hoverinfo="skip",
                             showlegend=False), row=4, col=1)
    fig.add_trace(go.Scatter(x=span, y=m["twist"], mode="lines+markers+text",
                             text=ids, textposition="top center",
                             textfont=dict(size=11),
                             line=dict(color="#d62728"), hoverinfo="skip",
                             showlegend=False), row=5, col=1)

    fig.update_yaxes(title_text="chord extent [m]", scaleanchor="x",
                     scaleratio=1, row=1, col=1)
    fig.update_yaxes(title_text="BlAFID", range=[cmin - 3, cmax + 5], row=2, col=1)
    fig.update_yaxes(title_text="Δspan [m]", rangemode="tozero", row=3, col=1)
    fig.update_yaxes(title_text="t/c [%]", row=4, col=1)
    fig.update_yaxes(title_text="twist [°]", row=5, col=1)
    fig.update_xaxes(title_text="blade span BlSpn [m]", row=5, col=1)
    fig.update_layout(template="plotly_white", height=1200, autosize=True,
                      font=dict(size=15),
                      title="AeroDyn15 blade: airfoil ID per blade section",
                      margin=dict(t=80))
    fig.write_html(out_path, include_plotlyjs="cdn")
    return fig


# --------------------------------------------------------------------------- #
# View 3: 3D lofted blade
# --------------------------------------------------------------------------- #
def plot_3d(m: dict, out_path: str, n: int = 60) -> str:
    import plotly.graph_objects as go

    span, chord, afid, twist = m["span"], m["chord"], m["afid"], m["twist"]
    prebend, sweep, tc = m["prebend"], m["sweep"], m["tc"]
    _, _, cmin, cmax = color_setup(afid)

    xs, ys, zs, hover = [], [], [], []
    for i in range(len(span)):
        tci = tc[i] if np.isfinite(tc[i]) else 25.0
        px, py = naca_profile(tci, n=n)
        xe, yf, z = section_3d(px, py, chord[i], twist[i], span[i],
                               prebend[i], sweep[i])
        xs.append(xe); ys.append(yf); zs.append(z)
        hover.extend([_hover(m, i)] * len(xe))
    X = np.concatenate(xs); Y = np.concatenate(ys); Z = np.concatenate(zs)

    # One flat colour per spanwise strip (constant per section, no gradient):
    # each strip [i, i+1] takes its inboard node's airfoil ID, mirroring the
    # step-function used everywhere else.
    M = 2 * n
    I, J, K, face_id = [], [], [], []
    for s in range(len(span) - 1):
        for q in range(M):
            a = s * M + q
            b = s * M + (q + 1) % M
            c = (s + 1) * M + q
            d = (s + 1) * M + (q + 1) % M
            I += [a, a]; J += [b, d]; K += [d, c]
            face_id += [float(afid[s]), float(afid[s])]

    mesh = go.Mesh3d(
        x=X, y=Y, z=Z, i=I, j=J, k=K, intensity=face_id,
        colorscale="Turbo", cmin=cmin, cmax=cmax, intensitymode="cell",
        text=hover, hoverinfo="text", opacity=1.0, flatshading=True,
        colorbar=dict(title="airfoil ID"), name="blade",
    )
    # crisp section outlines at every node
    outlines = []
    for i in range(len(span)):
        tci = tc[i] if np.isfinite(tc[i]) else 25.0
        px, py = naca_profile(tci, n=n)
        xe, yf, z = section_3d(px, py, chord[i], twist[i], span[i],
                               prebend[i], sweep[i])
        outlines.append(go.Scatter3d(
            x=xe, y=yf, z=z, mode="lines",
            line=dict(color="rgba(0,0,0,0.35)", width=1.5),
            hoverinfo="skip", showlegend=False))

    fig = go.Figure(data=[mesh, *outlines])
    # Cross-section is exaggerated vs span so the planform/twist are visible
    # (a real 56 m blade is too slender to read at true proportions).
    fig.update_layout(
        title="AeroDyn15 blade: 3D loft coloured by airfoil ID "
              "(generic sections from derived t/c, illustrative; "
              "cross-section exaggerated)",
        scene=dict(
            xaxis_title="edgewise [m]", yaxis_title="flapwise [m]",
            zaxis_title="span [m]", aspectmode="manual",
            aspectratio=dict(x=0.5, y=0.32, z=2.4),
            camera=dict(eye=dict(x=1.9, y=1.5, z=0.5)),
        ),
        margin=dict(t=60, l=0, r=0, b=0), template="plotly_white",
        height=850, autosize=True, font=dict(size=15),
    )
    fig.write_html(out_path, include_plotlyjs="cdn")
    return fig


# --------------------------------------------------------------------------- #
def main() -> None:
    root = repo_root()
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["all", "static", "interactive", "3d"],
                    default="all")
    ap.add_argument("--blade",
                    default=os.path.join(root, "AeroDyn15", "AeroDyn15_blade.dat"))
    ap.add_argument("--aerodyn",
                    default=os.path.join(root, "AeroDyn15", "AeroDyn15.dat"))
    ap.add_argument("--airfoils",
                    default=os.path.join(root, "AeroDyn15", "Airfoils"))
    ap.add_argument("--outdir", default=os.path.join(root, "figures"))
    ap.add_argument("--snapshots", action="store_true",
                    help="also write PNG snapshots of the HTML views "
                         "(needs plotly image export / Chrome)")
    ap.add_argument("--scale", type=float, default=4.0,
                    help="pixel-density multiplier for the PNG snapshots "
                         "(default 4 => blade_interactive.png is 8800x4800). "
                         "Raise this for higher resolution; it keeps text "
                         "proportions (unlike inflating width/height).")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    m = load_model(args.blade, args.aerodyn, args.airfoils)
    print(f"Loaded {len(m['span'])} blade nodes, "
          f"{len(m['af_label'])} airfoils, "
          f"span 0 to {m['span'][-1]:.2f} m.")

    def snapshot(fig, path, **kw):
        if not args.snapshots:
            return
        try:
            fig.write_image(path, **kw)
            print("wrote", path)
        except Exception as exc:  # missing kaleido/Chrome, non-fatal
            print(f"  (skipped PNG snapshot {os.path.basename(path)}: {exc})")

    if args.mode in ("all", "static"):
        out = os.path.join(args.outdir, "blade_overview.png")
        plot_static(m, out)
        print("wrote", out)
    if args.mode in ("all", "interactive"):
        out = os.path.join(args.outdir, "blade_interactive.html")
        fig = plot_interactive(m, out)
        print("wrote", out)
        # Resolution is driven by `scale` (keeps text proportions); the fixed
        # logical width/height set how large the fonts read on the canvas.
        snapshot(fig, os.path.join(args.outdir, "blade_interactive.png"),
                 width=2200, height=1200, scale=args.scale)
    if args.mode in ("all", "3d"):
        out = os.path.join(args.outdir, "blade_3d.html")
        fig = plot_3d(m, out)
        print("wrote", out)
        snapshot(fig, os.path.join(args.outdir, "blade_3d.png"),
                 width=2000, height=1100, scale=args.scale)


if __name__ == "__main__":
    main()
