# Provenance of the production blade

How `../AeroDyn15/` (the current production model) was derived from the originally
delivered blade, with the validation at each step. Everything referenced here is
preserved in this `archive/` folder.

## The blade

DLR **WiValdi-115** research turbine blade (~58 m tip radius), a continuously
*blended-airfoil* design: 39 of the original 44 AeroDyn nodes carry a distinct
airfoil, and each airfoil station is itself a named fractional blend of a small
set of base airfoils (e.g. `0.275_EC135F_EC128`).

## Step 1 - Uniform re-discretization (aerodynamics preserved)

**Why:** OpenFAST coupled to a CFD actuator-line solver needs uniform blade-node
spacing to avoid force-projection artifacts. But the blade cannot simply be
re-spaced, because the node positions *are* the airfoil-definition stations:
moving a node moves an airfoil.

**What:** `archive/tools/uniformize_blade.py` re-samples the
*continuous* blended blade onto a uniform 57-node grid (1.009 m spacing) and
rebuilds the airfoil data at each node:

- geometry (chord, twist, prebend, sweep) interpolated onto the uniform grid;
- one blended polar per node - the static Cl/Cd/Cm table linearly interpolated
  at matched angle of attack between the two bracketing original stations
  (`coef = (1-w)*lo + w*hi`), the standard NREL method (AirfoilPrep.py
  `Polar.blend`; identical in welib / openfast_toolbox / WISDEM). Exact at the
  39 original stations.
- the 8 polar-derived Beddoes-Leishman UA parameters (`alpha0`, `C_nalpha`,
  `alpha1/2`, `Cn1/2`, `Cd0`, `Cm0`) are *omitted*, so OpenFAST recomputes them
  from each blended table at init (`AirfoilInfo.f90:CalculateUACoeffs`) - what
  WISDEM/WEIS do. The airfoil-independent literature constants stay `Default`.

This produced `AeroDyn15_uniform_v15chord/`.

**Validation (OpenFAST 3.5.5, built from source; rigid, variable rotor speed,
4-15 m/s):** the uniform 57-node blade reproduces the original 44-node blade to
**|dP| <= 0.57 %, |dT| <= 0.20 %**, spanwise `Fn` within 1.5 % of peak. A grid
refinement at 10 m/s shows the small power deficit is pure discretization,
halving as the spacing halves (57: -0.47 %, 85: -0.20 %, 113: -0.10 %).
Independent re-derivation of the model (`archive/tools/verify_uniform_blade.py`) passed
all integrity checks. Figures: `figures/verify_4_15ms.png`,
`figures/verify_spanwise_4_15ms.png`, `figures/uniform_validation*.png`.

## Step 2 - Chord reconciled to the CAD geometry (for ALM)

**Why:** comparing chord across three independent representations of the blade
(AeroDyn v14 table, AeroDyn v15 model, and the CAD (IGES) blade-skin geometry)
showed that from r ~ 8-38 m all three agree to 1-2 cm, but **for r > 40 m the
v15 chord is ~6-9 % narrower than both the CAD and v14** (which agree with each
other to ~2.5 cm). The CAD chord was verified by evaluating the actual NURBS
blade-skin surface, so the deficit is real, not an extraction artifact.

In **Actuator Line Method** the sectional force projected into the CFD domain is
`F = 1/2 rho V_rel^2 c (Cl, Cd)` - it scales directly with chord `c`. For the
actuator lines to represent the physical blade that the CFD mesh/domain is built
around, `c` must be the **physical (CAD) chord**. Using the v15 chord would
inject ~6-9 % too little force over the outer span.

**What:** `archive/tools/rechord_to_cad.py` rewrites **only** the `BlChord` column for
r > 38 m to the CAD/v14 planform (smooth cosine ramp over 38-40 m), leaving
span, twist, prebend, sweep, airfoil IDs and **every polar file identical** to
Step 1. The inboard chord - including the cylinder/root transition, where the
CAD chord of a near-circular section is ill-defined - is left exactly as v15.

This produced the model now promoted to `../AeroDyn15/`. The **polars are
byte-identical** between the v15-chord and CAD-chord models; only the blade
chord differs. (Cl/Cd are dimensionless airfoil properties, valid regardless of
chord; the t/c distribution matches the CAD.)

**Load effect** (`figures/cadchord_load_effect.png`): the wider outboard chord
raises thrust at every wind speed and raises power/torque toward and above
rated - a deliberate geometric correction, not an error.

## Net result

`../AeroDyn15/` = uniform 57-node grid (1.009 m) + 57 blended polars (UA
recomputed by OpenFAST) + CAD-matched physical chord. Runs cleanly in OpenFAST
3.5.5; intended for OpenFAST-CFD actuator-line coupling of the WiValdi-115 blade.

## How to reproduce from scratch

```bash
# 1. uniform re-discretization (writes a model dir from the original case)
python archive/tools/uniformize_blade.py 57 <out_dir>
# 2. reconcile outboard chord to the CAD/v14 planform
python archive/tools/rechord_to_cad.py <out_dir> <out_dir_cadchord>
```
(`uniformize_blade.py` reads the original case; the CAD chord reference lives in
`data/cad_chord_vs_radius.csv`, extracted from the WiValdi-115 IGES.)
