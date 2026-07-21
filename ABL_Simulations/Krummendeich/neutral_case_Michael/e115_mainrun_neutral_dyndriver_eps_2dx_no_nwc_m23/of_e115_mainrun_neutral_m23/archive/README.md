# archive/

This folder holds **superseded** material kept for traceability and provenance.
None of it is needed to run simulations with the current blade. The production
model lives at the top level in [`../AeroDyn15/`](../AeroDyn15/).

It is archived (not deleted) so the history of how the production blade was
derived, and the validation behind it, remains auditable.

## Why each item is here

| Item | What it is | Why it was superseded |
|------|------------|-----------------------|
| `AeroDyn15_original_44node/` | The blade exactly as originally delivered: 44 AeroDyn nodes, 39 distinct blended-airfoil polars, non-uniform spacing. | Replaced by a uniform-spacing re-discretization for CFD actuator-line coupling. Kept as the reference baseline every validation compares against. |
| `AeroDyn15_uniform_v15chord/` | Intermediate model: the original blade re-sampled onto a uniform 57-node grid with the polars re-blended, but keeping the original (v15) chord. | Its outboard chord (r > 40 m) is ~6-9% narrower than the CAD geometry. For ALM, the chord must match the physical blade, so the production model uses the CAD chord instead. This intermediate documents the chord change and the "aerodynamics-preserving" validation. |
| `tools/` | All scripts: the blade generator (`uniformize_blade.py`), the visualizer (`blade_viz.py`), and the comparison/validation tools (original-vs-uniform, v14/v15/CAD chord, sweep plots, the independent integrity verifier, the chord-reconciliation tool). | They were used to build and validate the production blade; `uniformize_blade.py` and `blade_viz.py` can still regenerate or visualize it (see the top-level README). |
| `data/` | Reference inputs: CAD chord-vs-radius (`cad_chord_vs_radius.csv`) and the AeroDyn v14 blade table (`ad14_blade.csv`). | Used by the generation and comparison tools; not read at simulation time. |
| `figures/` | Blade-layout views (overview, 3D, interactive) plus all comparison/validation figures (OpenFAST sweeps, v14/v15/CAD chord overlays, the 4-15 m/s verification, the chord-reconciliation and load-effect plots). | Illustrative; describe the blade and the old-vs-new comparison, not needed to run. |

See [`PROVENANCE.md`](PROVENANCE.md) for the full step-by-step history of how the
production blade was produced from `AeroDyn15_original_44node/`.
