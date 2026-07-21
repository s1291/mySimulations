# WiValdi-115 blended-airfoil blade: uniform AeroDyn model for ALM

AeroDyn v15 aerodynamic model of the DLR **WiValdi-115** wind-turbine blade,
re-discretized onto a uniform spanwise grid and chord-matched to the CAD
geometry, for **OpenFAST coupled to CFD actuator-line (ALM)** simulations.

## Use this

```
"AeroDyn15/AeroDyn15.dat"   AeroFile
```

[`AeroDyn15/`](AeroDyn15/) is the production model: 57 nodes, uniform 1.009 m
spacing, 57 blended polars, CAD-matched chord. See
[`AeroDyn15/README.md`](AeroDyn15/README.md) for its properties and
[`archive/PROVENANCE.md`](archive/PROVENANCE.md) for how it was built and
validated (OpenFAST 3.5.5, 4-15 m/s).

> The polars omit the 8 polar-derived Beddoes-Leishman UA parameters on purpose:
> OpenFAST recomputes them from each blended table at init. That is the standard
> method, not a missing input.

## Repository layout

| Path | What |
|------|------|
| [`AeroDyn15/`](AeroDyn15/) | **Production blade** + 57 polars (point your `.fst` here). |
| [`InflowWind/`](InflowWind/) | InflowWind input for the OpenFAST + PALM coupling (`WindType=10`) and the PALM TCP connection file. |
| [`archive/`](archive/) | Everything used to build and validate the blade: superseded models, the generation and comparison tools, reference data, validation figures, and full provenance. Not needed to run. |

## Regenerate or revisualize the blade

The generation and visualization tools live in [`archive/tools/`](archive/tools/),
with reference data in [`archive/data/`](archive/data/):

```bash
# regenerate the uniform blade (57 nodes ~1 m; use a larger count for finer spacing)
python archive/tools/uniformize_blade.py 57 <out_dir>

# visualize the blade layout (static PNG + interactive HTML)
pip install -r archive/tools/requirements.txt
python archive/tools/blade_viz.py
```

The uniform model preserves the aerodynamics of the original case. To then match
the chord to the CAD geometry, see [`archive/PROVENANCE.md`](archive/PROVENANCE.md)
(Step 2). Visualization outputs and validation figures are in
[`archive/figures/`](archive/figures/).
