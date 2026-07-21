# AeroDyn15/ production blade (WiValdi-115, uniform, CAD chord)

The current AeroDyn v15 aerodynamic model of the DLR WiValdi-115 blade, prepared
for **OpenFAST-CFD actuator-line (ALM) coupling**.

## Contents

| File | Description |
|------|-------------|
| `AeroDyn15.dat` | AeroDyn v15 primary input (points at the blade + 57 polars). |
| `AeroDyn15_blade.dat` | Blade table: 57 nodes, **uniform 1.009 m spacing**, CAD-matched chord. |
| `Airfoils/Polars01..57.dat` | One blended polar per node (Cl/Cd/Cm + Beddoes-Leishman UA block). |

## Key properties

- **Uniform spacing** (57 nodes, 1.009 m) - for clean ALM force projection.
- **Chord matches the CAD geometry** - the physical blade the CFD represents.
  In ALM the projected force scales with chord, so this is the geometrically
  faithful choice (independently corroborated by the AeroDyn v14 table).
- **Polars** are linear blends of the original 39 stations onto the uniform grid
  (exact at the originals). The 8 polar-derived UA parameters are intentionally
  **omitted**, so OpenFAST recomputes them from each blended table at init
  (audit them in the `<root>.UA.sum` file). The literature-default UA constants
  stay `Default`. Do not "fix" the missing UA lines - that is by design.

## Usage

Point the AeroDyn file in your OpenFAST `.fst` at:

```
"AeroDyn15/AeroDyn15.dat"   AeroFile
```

## Solver settings (configured for actuator-line CFD coupling)

`AeroDyn15.dat` is set up so AeroDyn computes only the **sectional airfoil
forces** and lets the CFD resolve the wake and tower:

| Setting | Value | Why |
|---------|-------|-----|
| `WakeMod` | **0** | No internal induction/wake - the CFD solver provides the induced velocity. Leaving BEMT/DBEMT on would double-count induction. |
| `TwrPotent` / `TwrAero` | **0 / False** | Tower influence is represented in the CFD domain; AeroDyn's tower model would double-count it. |
| `TipLoss` / `HubLoss` | **False** | Tip/hub-loss are BEMT induction corrections; disabled explicitly (already inactive when `WakeMod=0`). |
| `AFAeroMod` | 2 | Beddoes-Leishman unsteady airfoil aero - a *sectional* model, kept on; it does not conflict with the CFD wake. |

If you instead want a **standalone** (non-CFD) BEM run, set `WakeMod=2`,
`TwrPotent=1`, `TwrAero=True`, `TipLoss=True`, `HubLoss=True` - that is the
configuration used for the validation in `../archive/PROVENANCE.md`.

## Validation & provenance

Validated against the originally delivered blade in OpenFAST 3.5.5 over 4-15 m/s
(power/thrust within <0.6 %/0.2 % from the uniform re-discretization; the chord
change is a deliberate geometric correction). Full derivation and validation:
[`../archive/PROVENANCE.md`](../archive/PROVENANCE.md).
