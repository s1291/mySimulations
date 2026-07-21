# Simulation setup

## Precursor

- Neutral atmospheric boundary layer with zero surface heat flux.
- Domain extent in the x direction: 4.8 km.
- Domain extent in the y direction: 3.2 km.
- Domain extent in the z direction: 4.0 km.
- Grid spacing: 8 m in the x, y, and z directions.
- Target wind at 90 m: 8.0 m/s in the positive x direction.
- Surface roughness length: 0.1 m.
- Latitude: 53.83 degrees.
- Potential temperature is uniform below 500 m, increases by 2 K from 500 to 600 m, and has a stable gradient of 10 K/km above 600 m.
- Precursor duration: 18 hours.
- Developed inflow near hub height: about 8.01 m/s mean wind speed, 9.7 percent turbulence intensity, and 0.452 m/s friction velocity.
- Developed boundary-layer depth: about 620 m.

## Main run

- The main run uses the developed turbulent precursor inflow.
- Three grid refinement levels are used, with grid spacings of 8 m, 4 m, and 2 m.
- The wind turbine is located in the finest refinement level, where the grid spacing is 2 m.
- Total simulation time: 50 minutes.
- The first 20 minutes are discarded.
- Flow data from the final 30 minutes are saved every 1 second.
- Saved variables:
  - `u`: wind velocity in the x direction.
  - `v`: wind velocity in the y direction.
  - `w`: vertical wind velocity.
  - `p`: grid-resolved pressure perturbation. Only pressure differences are physically meaningful.
- The tower axis is located at x = 488 m and y = 176 m.
- The OpenFAST overhang is 5.0191 m upstream of the tower axis.
- The rotor-plane center is therefore located at x = 482.981 m, y = 176 m, and z = 90 m.
- The saved x range is 46 to 740 m. Relative to the rotor plane, this covers 437.0 m, or 3.47 rotor diameters, upstream and 257.0 m, or 2.04 rotor diameters, downstream.
- The saved y range is 2 to 350 m. Relative to the rotor centerline, this covers 174 m, or 1.38 rotor diameters, on each side.
- The saved z range is 2 to 160 m. The rotor extends from z = 27 to 153 m, so the saved region includes 25 m, or 0.20 rotor diameters, below the rotor and 7 m, or 0.06 rotor diameters, above it.
- The complete saved region is 694 m by 348 m by 158 m, equivalent to 5.51 by 2.76 by 1.25 rotor diameters in the x, y, and z directions.

## Wind turbine

- Turbine: NREL 5 MW.
- The turbine is modeled with the Actuator Line Method through PALM and OpenFAST coupling.
- Hub height: 90 m.
- Rotor diameter: 126 m.
- Rotor speed: fixed at 9.155 rpm.
- Collective pitch: 0 degrees.
- Precone: 0 degrees.
- Shaft tilt: 0 degrees.
- Mean yaw error: 0 degrees.
- Fixed platform, rigid blades, rigid tower, and rigid drivetrain.
- Generator and drivetrain speed changes are disabled.
- Steady airfoil aerodynamics are used.
- Rotor induction and wake development are resolved by PALM.
