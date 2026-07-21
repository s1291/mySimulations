# Comparison of four large-eddy simulation research codes and effects of model coefficient and inflow turbulence in actuator-line-based wind turbine modeling 

Luis A. Martínez-Tossas; Matthew J. Churchfield; Ali Emre Yilmaz © ; Hamid Sarlak; Perry L. Johnson; Jens N. Sørensen; Johan Meyers; Charles Meneveau

Check for updates
J. Renewable Sustainable Energy 10, 033301 (2018)
https://doi.org/10.1063/1.5004710

- CHORUS
![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-01.jpg?height=108&width=109&top_left_y=822&top_left_x=1539)

View Online
![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-01.jpg?height=127&width=97&top_left_y=822&top_left_x=1685)

## Articles You May Be Interested In

Collaborating teams compare wind turbine simulation codes
Scilight (May 2018)
Vortex induced vibrations of wind turbine blades: Influence of the tip geometry
Physics of Fluids (June 2020)
Effects of turbulent inflow time scales on wind turbine wake behavior and recovery
Physics of Fluids (September 2023)

# Special Topics Open for Submissions 

# Comparison of four large-eddy simulation research codes and effects of model coefficient and inflow turbulence in actuator-line-based wind turbine modeling 

Luis A. Martínez-Tossas, ${ }^{1,2}$ Matthew J. Churchfield, ${ }^{2}$ Ali Emre Yilmaz, ${ }^{3}$ Hamid Sarlak, ${ }^{4}$ Perry L. Johnson, ${ }^{1}$ Jens N. Sørensen, ${ }^{4}$ Johan Meyers, ${ }^{3}$ and Charles Meneveau ${ }^{1}$<br>${ }^{1}$ Department of Mechanical Engineering, Johns Hopkins University, Baltimore, Maryland 21218, USA<br>${ }^{2}$ National Renewable Energy Laboratory, Golden, Colorado 80401-3305, USA<br>${ }^{3}$ Department of Mechanical Engineering B3001, KU Leuven, Leuven, Belgium<br>${ }^{4}$ Fluid Mechanics Section, Department of Wind Energy, Technical University of Denmark, Lygby, Denmark

(Received 14 September 2017; accepted 23 February 2018; published online 16 May 2018)


#### Abstract

Large-eddy simulation (LES) of a wind turbine under uniform inflow is performed using an actuator line model (ALM). Predictions from four LES research codes from the wind energy community are compared. The implementation of the ALM in all codes is similar and quantities along the blades are shown to match closely for all codes. The value of the Smagorinsky coefficient in the subgrid-scale turbulence model is shown to have a negligible effect on the time-averaged loads along the blades. Conversely, the breakdown location of the wake is strongly dependent on the Smagorinsky coefficient in uniform laminar inflow. Simulations are also performed using uniform mean velocity inflow with added homogeneous isotropic turbulence from a public database. The time-averaged loads along the blade do not depend on the inflow turbulence. Moreover, and in contrast to the uniform inflow cases, the Smagorinsky coefficient has a negligible effect on the wake profiles. It is concluded that for LES of wind turbines and wind farms using ALM, careful implementation and extensive cross-verification among codes can result in highly reproducible predictions. Moreover, the characteristics of the inflow turbulence appear to be more important than the details of the subgrid-scale modeling employed in the wake, at least for LES of wind energy applications at the resolutions tested in this work. Published by AIP Publishing. https://doi.org/10.1063/1.5004710


## I. INTRODUCTION

In recent years, large-eddy simulation (LES) has become a prominent tool for numerical studies of wind turbine wakes and wind farms. ${ }^{1-3}$ One of the most accurate representations of a wind turbine in LES, apart from a fully resolved rotor, is the actuator line model (ALM). ${ }^{4-6}$ The wind energy computational research community has implemented ALM in a range of different LES frameworks (or numerical codes). To enhance the trustworthiness of LES-generated data using ALM and the different codes employed by different research groups, it is crucial to perform both a detailed cross-code comparison and sensitivity analysis for model parameters.

Recent studies have compared different numerical codes used in the community. The Blind Test campaign has compared many codes and numerical methods to experimental measurements. ${ }^{7-9}$ Many differences in the results from all codes were observed. The differences are caused not only by the numerical discretization method, but mostly because of different simulation parameters, such as boundary conditions, lift and drag coefficient tables, value of smoothing scale $\epsilon$ in the ALM, and nacelle and tower models. The work of Lignarolo et al. ${ }^{10}$ shows a comparison of different codes using the actuator disk with experimental measurements of flow over a porous disk. Good agreement was observed between the different codes in predicting the
velocity deficits in the wakes. One of the main differences between the codes was observed in predicting turbulence intensities depending on the subgrid-scale models used. Sarlak et al. ${ }^{11}$ compare two finite-volume codes and observe agreement in the near wake, but differences were observed depending on subgrid-scale modeling in some cases, and in the far wake. Furthermore, Martínez-Tossas et al. ${ }^{6}$ present comparisons of two codes (finite-volume and pseudo-spectral), with emphasis on how numerical discretization changes the breakdown of the wake. Sarlak et al. ${ }^{12}$ also present a comparison of a pseudo-spectral and finite-volume code with good agreement in the near wake, and differences in the breakdown location of the wake are observed in the case of uniform inflow.

In this study, we focus on the effects of numerical discretization and the Smagorinsky coefficient in the subgrid-scale model by running simulations of a single wind turbine under uniform laminar inflow using four different codes widely used in the wind energy research community. "Numerical discretization" is used to refer to the numerical discretization method used to solve the fundamental set of equations, which in this case consists of a range of pseudo-spectral, finite-difference, and finite-volume algorithms. The goal is to document the robustness of the ALM approach with respect to numerical discretization used in different codes, by running exactly the same case using different codes. Simulations are also run using one of the codes with turbulent inflow. The inflow used is homogeneous isotropic turbulence with zero mean shear. ${ }^{13}$ The effects of the Smagorinsky coefficient on quantities computed by the ALM and wake profiles are also studied in the case of turbulent inflow.

## II. BRIEF OVERVIEW OF THE FOUR CODES

Four different LES codes being used by the wind energy community are compared. This study focuses on possible effects of the different numerical methods and subgrid-scale turbulence model parameters. The codes compared employ a combination of finite-volume, finite-difference, and pseudo-spectral numerical discretization. Sections II A-IID provide the specifics on the numerical discretization for each framework. Two codes use a combination of pseudo-spectral and finite-difference numerical discretizations. The finite-difference resolution is chosen to be twice the spectral resolution. ${ }^{6}$

All codes solve the filtered Navier-Stokes equations

$$
\begin{gather*}
\nabla \cdot \widetilde{\boldsymbol{u}}=0,  \tag{1}\\
\frac{\partial \widetilde{\boldsymbol{u}}}{\partial t}+\widetilde{\boldsymbol{u}} \cdot \nabla \widetilde{\boldsymbol{u}}=-\nabla \widetilde{p}-\nabla \cdot \tau+\widetilde{\mathbf{f}}_{\epsilon}, \tag{2}
\end{gather*}
$$

where $\widetilde{\boldsymbol{u}}$ is the filtered resolved velocity vector, $\widetilde{p}$ is the filtered pressure divided by density, $\tau$ is the subgrid-scale stress tensor, and $\widetilde{\mathbf{f}}_{\epsilon}$ is the body force that represents the wind turbine. In all codes, the subgrid-stress tensor is modeled using an eddy-viscosity model, ${ }^{14}$ in which the deviatoric part of the subgrid-scale stress tensor is

$$
\begin{equation*}
\tau_{i j}=-2 \nu_{S G S} \widetilde{S}_{i j}, \quad \nu_{S G S}=\left(C_{s} \Delta\right)^{2}\left(2 \widetilde{S}_{i j} \widetilde{S}_{i j}\right)^{1 / 2}, \tag{3}
\end{equation*}
$$

where $C_{S}$ is the Smagorinsky coefficient and $\widetilde{S}_{i j}$ is the symmetric part of the resolved velocity gradient tensor. Some of the codes used have a dynamic Smagorinsky model implementation that uses the Germano identity to calculate $C_{S}$ as a function of space and time. ${ }^{15-18}$ The standard Smagorinsky model assumes a constant value of the Smagorinsky coefficient. The theoretically derived value for homogeneous isotropic turbulence with a spectral cutoff filter is $C_{S}=0.16 .^{14,19}$ To be consistent, the codes used were run using the same value for the coefficient ( $C_{S}=0.16$ ). This value ( $C_{S}=0.16$ ) is appropriate for homogeneous isotropic turbulence, but it is not ideal for simulations of wind turbine wakes. ${ }^{6}$ Section V addresses this issue in the context of the Lagrangian-scale-dependent model. ${ }^{18}$ Here, we chose $C_{s}=0.16$ nonetheless, since this value corresponds to a canonical reference case, thus allowing us to avoid having to choose another value more arbitrarily.

TABLE I. Codes used in this study with the description of the numerical discretization.
| Johns Hopkins University (JHU) | National Renewable Energy Laboratory (NREL) | KU Leuven | Technical University of Denmark (DTU) |
| :--- | :--- | :--- | :--- |
| LESGO | SOWFA | SP-Wind | EllipSys3D |
| Scheme $x$ and $y$ : pseudo-spectral $z$ : second-order finite difference | Scheme $x, y$, and $z$ : second-order finite volume | Scheme $x$ and $y$ : pseudo-spectral $z$ : fourth-order finite difference | Scheme $x, y$, and $z$ : second-order finite volume |


Sections II A-II D show a brief summary of the numerical methods in every code. Table I gives an outline of the code names and their numerical methods. All codes have a similar implementation of the ALM.

## A. LESGO-The Johns Hopkins University Code

LESGO is the pseudo-spectral LES code used by the Turbulence Research Group at Johns Hopkins University (JHU). ${ }^{20}$ The numerics are based on the early work of Moeng ${ }^{21}$ and Albertson. ${ }^{22}$ The code is pseudo-spectral, implying periodic boundary conditions in the streamwise and one of the spanwise directions. The other spanwise direction uses second-order, centered finite difference. The boundary conditions in this direction are zero shear stress with no penetration. Dealiasing for the nonlinear term is done using the $3 / 2$ rule. ${ }^{23}$ Time integration is done using the second-order Adams-Bashforth method. A fringe region (7.5\% of the domain) is used to smoothly transition the end of the domain to uniform inflow. ${ }^{24}$ The subgrid-scale models implemented are based on the standard Smagorinsky ${ }^{14}$ including several variants, such as dynamic, scale-dependent, and Lagrangian-averaged versions. ${ }^{16,18}$

## B. SOWFA-The National Renewable Energy Laboratory Code

Simulator for Wind Farm Applications (SOWFA) is an LES solver developed by the National Renewable Energy Laboratory (NREL) meant for wind farm simulations implemented under the OpenFOAM framework. ${ }^{25,26}$ It is a finite-volume code with second-order numerical discretization in both space and time. The boundary conditions are set to uniform inflow with zero normal pressure gradient at the inlet, zero normal gradient of velocity, and fixed pressure at the outlet. The lateral boundary conditions are set to zero gradient with no penetration. Time advancement is done using second-order backward differentiation. The subgrid-scale model used in the code for this study is of the standard Smagorinsky type with a fixed $C_{S}$ coefficient.

## C. SP-Wind-The KU Leuven Code

SP-Wind is a pseudo-spectral LES code from the Turbulent Flow Simulation and Optimization Group at KU Leuven. ${ }^{27,28}$ It too uses pseudo-spectral discretization in the streamwise and spanwise directions. The other spanwise direction uses fourth-order finite differencing. The boundary conditions are periodic in the pseudo-spectral directions and zero stress with no penetration in the finite-difference direction. A fringe region is used to drive the flow to uniform inflow. ${ }^{24}$ The subgrid-scale model used is the standard Smagorinsky model with a prescribed coefficient. ${ }^{14}$ Time integration is performed using a standard fourth-order Runge-Kutta method.

## D. EllipSys3D-The Technical University of Denmark Code

EllipSys3D is the finite-volume code from the Technical University of Denmark (DTU). ${ }^{29,30}$ EllipSys3D is a general-purpose finite volume solver on multiblock structured grids. In EllipSys3D, diffusive and convective terms are discretized using second-order central differencing schemes (CDSs) and a blend of CDS ( $10 \%$ ) and third-order QUICK scheme

TABLE II. Domain dimensions.
| $L_{x}: 24 D$ | $L_{y}: 6 D$ | $L_{z}: 6 D$ |
| :--- | :---: | :---: |
| Turbine location: $3 D$ | Inflow: $8 \mathrm{~m} / \mathrm{s}$ | $\rho=1.0 \mathrm{~kg} / \mathrm{m}^{3}$ |


( $90 \%$ ). Boundary conditions are set to symmetry on the walls with inflow and convective outflow boundary conditions. Temporal discretization is performed using a second-order backward Euler scheme, and the solution is marched in time using inner time stepping. Pressure checkerboarding is prevented by using Rhie-Chow interpolation on a collocated grid arrangement, and the pressure correction equation is solved using the Pressure-Implicit with Splitting of Operators (PISO) algorithm.

## III. DESCRIPTION OF SIMULATED FLOW

A standard test case is defined and simulated with the different codes. The turbine used is the NREL 5-MW reference turbine, which has a rotor diameter ( $D$ ) of $126 \mathrm{~m} .^{31}$ The dimensions of the computational box are shown in Table II, and a schematic of the domain is shown in Fig. 1. The rotational speed of the turbine is 9.155 rpm , which corresponds to a tip-speed ratio (TSR) of 7.55.

The value of $\epsilon$ is a parameter that establishes how the forces are smeared onto the grid. 4,32 The known numerical limit of the value is related to the grid size as $\epsilon \geq 2 \Delta x .^{4,5,32}$ Turbine quantities, such as lift and drag on the blades, are very dependent on $\epsilon$; for this reason, all cases presented have been run using a fixed value of $\epsilon=10 \mathrm{~m} .^{32}$ This is far from the optimum value recently found $\left(\epsilon / c \approx 0.25\right.$, where $c$ is the chord) ${ }^{33}$ but allows the use of a bigger domain with a uniform grid, which is needed in this case to study the wakes without local grid refinement. Running simulations using the optimal value of $\epsilon / c \approx 0.25$ requires very fine resolutions and local refinement near the rotor. ${ }^{34,35}$ This study focuses on the effect of numerical discretiza- tion and the Smagorinsky coefficient on the wake and its transition to turbulence in the context of LES using uniform grids and a practically affordable number of grid points. The grid resolution used in all codes is $\Delta / D=0.03125$ in the spectral directions and $\Delta / D=0.015625$ in the finite difference directions. The finite resolution in the spectral directions is twice the resolution in the finite difference directions. ${ }^{6}$ Initial tests were performed for different grid ratios, and a ratio of $\Delta_{\text {Spectral }} / \Delta_{\mathrm{FD}}=2$ gave converged Reynolds stresses in both directions.

The chord and twist angle as a function of blade radius are linearly interpolated using tabulated data. ${ }^{31}$ The actuator points in the tip of the blade are extrapolated from the tabulated data. The time stepping ensures that the tip of the blade does not go through more than one grid cell in a time step. ${ }^{5,36}$ The number of actuator points is $N=64$ for all codes, and no tip correction is used in any of the codes.

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-05.jpg?height=414&width=1055&top_left_y=2043&top_left_x=534)
FIG. 1. Schematic of the computational domain. Reproduced with permission from Martínez-Tossas et al., J. Phys.: Conf. Ser. 625, 012024 (2015). ${ }^{6}$ Copyright 2015 Author(s), Licensed under a Creative Commons Attribution 3.0 Unported License.

## IV. CODE COMPARISON RESULTS

Simulations using a constant Smagorinsky coefficient of $C_{S}=0.16$ have been performed to compare only the effects of numerical discretization for all codes. The simulations use uniform inflow without turbulence. The simulations were run for one flow-through time to give the wake time to develop, and then time averaging was done for $1-2$ flow-through time.

## A. Along the blade quantities

Time-averaged quantities along the blades are shown in Fig. 2. The figure compares the predictions of the various codes and plots blade element momentum (BEM) calculations as a reference. BEM is a well-established method that uses momentum theory together with lift and drag coefficient tables to predict aerodynamic loads along the blades. ${ }^{37}$ The BEM equations are solved iteratively, and the results without a tip loss correction are shown for reference.

Figure 2 shows that the codes agree very well in terms of the ALM implementation. The main differences come from the axial velocity [Fig. 2(b)]. These differences are small, and they do not change the angle of attack and loads significantly [Figs. 2(a), 2(c), and 2(d)]. The mean and standard deviation are computed at every actuator point based on the results from the four LES codes. The maximum standard deviation from all actuator points of the plotted quantities in the outer portion of the blade ( $r / R \geq 0.2$ ) is $0.33^{\circ}$ for the angle of attack ( $\alpha$ ) and 0.03 for lift in the nondimensional lift force units $\left(F_{L} / l D \rho U_{\infty}^{2}\right.$, where $D$ is the rotor diameter, $l$ is the width of each blade section, $\rho$ is the density, and $U_{\infty}$ is the inflow velocity). These results show that the implementation of the ALM in every code is similar and differences because of the numerical discretization method are negligible. It is important to note that the value of $\epsilon$ used in all codes is the same, and changing this value will change the prediction of quantities along the blades. ${ }^{32}$ The differences with BEM are also due to the chosen value of $\epsilon$. ${ }^{6,32}$ There

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-06.jpg?height=992&width=1288&top_left_y=1503&top_left_x=415)
FIG. 2. Angle of attack, nondimensional axial velocity, nondimensional lift, and drag along the blade for a case with $C_{S}=0.16$. The blade is presented for reference.

are discontinuities in the lift and drag radial distributions of the JHU, NREL, and KU Leuven codes. These discontinuities are from the tabulated data, which change abruptly from one airfoil type to the next. In reality, the blade sections change smoothly. This smoothness can be achieved by interpolating the lift and drag tables between nearby airfoil sections. The DTU code EllipSys3D does the interpolation, which is why the lift and drag curves from this code are smoother. This approach can lead to small differences in the axial velocity. In this case, the shear layer of the wake is sharper, leading to an earlier transition to turbulence, as will be shown in Sec. IV B. Even though the quantities along the blade are discrete, the loads implemented in the LES by the ALM are smoothened using a Gaussian kernel, which makes the force field in all codes smooth, even though the quantities computed along the blades may have discontinuities.

## B. Velocity distributions

The loads implemented by the ALM are similar in all of the codes. In uniform inflow conditions, the flow field close to the turbine is also similar for all codes. The near wake in all codes is the same as shown in Fig. 3. This region of the flow is governed by the inviscid equations, and turbulence is not triggered until a later stage in the wake. Differences are observed downstream once turbulence is triggered. For representative visualizations of instantaneous snapshots of vorticity resulting from such simulations, the reader is referred to Figs. 15 and 18 from Sarlak et al. ${ }^{12}$ In Fig. 3, the far wake differs among the various codes, although the differences are more evident in Fig. 4, where the $\left\langle u^{\prime} u^{\prime}\right\rangle$ Reynolds stress component more clearly shows when the flow becomes turbulent. We observe that in EllipSys3D, the transition to turbulence occurs slightly earlier. This is because of the differences in subtle ALM implementation at the rotor, which cause a sharper shear layer near the rotor tip, as shown in Fig. 5. There are always small differences in implementation of the ALM, which are difficult to track. The most common one is what happens when defining points near the tip of the blade. The last point for which there are tabulated airfoil data is not at the tip. So, there is no definite guideline on what to do in that case. Some codes extrapolate using the last two sections closest to the tip. Some others use interpolation based on zero chord at the tip. As commonly known, the numerical discretization can play a role in the location where breakdown to turbulent flow occurs. The finite-volume codes (EllipSys3D and SOWFA) use a collocated grid arrangement. The accuracy for the second-order finite difference can be improved significantly by using a staggered arrangement so as to enable discrete energy conservation. ${ }^{38}$ However, the finite-volume codes presented here do not take advantage of these features as of yet and lower accuracy numerical discretization (e.g., second-order finite-volume compared to pseudo-spectral) is thus expected to delay the transition to turbulence. The finite-difference schemes damp the higher wave numbers in the derivatives. As a result, turbulence that would be triggered by high wave numbers may not be triggered if the high wave

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-07.jpg?height=543&width=1056&top_left_y=1947&top_left_x=536)
FIG. 3. Mean streamwise velocity ( $u / U_{\infty}$ ) contours for simulations from all codes using $C_{S}=0.16$. The vertical black line denotes the location of the rotor.

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-08.jpg?height=544&width=1060&top_left_y=342&top_left_x=534)
FIG. 4. Streamwise component of the Reynolds stress tensor $\left(\left\langle u^{\prime} u^{\prime}\right\rangle / U_{\infty}^{2}\right)$ contours for simulations from all codes using $C_{S}=0.16$. The vertical black line denotes the location of the rotor.

number modes are damped. In pseudo-spectral numerical discretization, these higher wave numbers are present, facilitating the transition to turbulence. This is why the codes with lower-order numerical discretization take longer to transition. It is important to note that, as shown in Sec. V, this sensitivity to numerical discretization is observed only for uniform inflow, where there is no turbulence in the inflow condition.

Figure 5 shows profiles of the mean streamwise velocity and streamwise component of the Reynolds stress tensor at different distances downstream. The mean streamwise velocity is very similar for all codes in the near and far wake. The main differences can be observed at

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-08.jpg?height=1020&width=1269&top_left_y=1479&top_left_x=426)
FIG. 5. Streamwise mean velocity (top) and streamwise Reynolds stress component (bottom) at different distances downstream. Notice the change of scale in the bottom plot for $x / D>9$, where the wake has become turbulent.

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-09.jpg?height=417&width=1056&top_left_y=342&top_left_x=536)
FIG. 6. Mean streamwise velocity ( $u / U_{\infty}$ ) contours for different Smagorinsky coefficients in the JHU code.

intermediate downstream locations where the wakes become turbulent at various locations for different codes. The streamwise component of the Reynolds stress tensor has significant differences in the near wake. The differences are smaller in the far wake after the flow has become turbulent in all codes. Once turbulence is triggered in the wake, the Reynolds stress term becomes several orders of magnitude higher.

## V. EFFECT OF THE SUBGRID-SCALE MODEL PARAMETER

Here, we focus on comparing the effects of choosing different Smagorinsky coefficients for the turbulence model using only the JHU code. We apply three variants-with the first being the most general approach, namely, the scale-dependent Lagrangian averaged model. The other two cases correspond to the traditional Smagorinsky model with constant coefficients of $C_{s}=0.08$ and 0.16 .

We find that the time-averaged quantities along the blade do not depend on the Smagorinsky coefficient (consistent with Ref. 6). However, the wake profiles in the far wake are strongly dependent on the Smagorinsky coefficient. Figures 6 and 7 show that the location of transition to turbulence is strongly influenced by the Smagorinsky coefficient. As may be expected, a higher Smagorinsky coefficient delays the transition to turbulence. The scaledependent Lagrangian model computes values, which, on average, are close to $C_{S}=0.08$. When running a simulation with a fixed $C_{S}=0.08$, the wake profiles are closer to those computed by the Lagrangian-scale-dependent model. This strong dependence on the $C_{S}$ coefficient is present in cases with uniform inflow without turbulence. Note from Fig. 7 that the turbulence intensity in the far wake is larger for the $C_{s}=0.16$ case than for the cases with lower $C_{s}$. A possible reason is that $C_{s}=0.16$ increases the damping of the smallest resolved turbulent eddies, which in turn reduces the rate of cascading energy and leads to an accumulation of kinetic energy at the large scales that dominate the $\left\langle u^{\prime} u^{\prime}\right\rangle$ component of the Reynolds stress. This explanation is only a conjecture, and more research on this topic is needed to draw conclusions about the energy cascade in this flow.

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-09.jpg?height=416&width=1053&top_left_y=2074&top_left_x=536)
FIG. 7. Streamwise component of the Reynolds stress tensor ( $\left\langle u^{\prime} u^{\prime}\right\rangle / U_{\infty}^{2}$ ) contours for simulations with different Smagorinsky coefficients in the JHU code.

## VI. EFFECTS OF TURBULENT INFLOW

The first part of this study focused on differences in numerical discretization and the Smagorinsky coefficient, all for laminar, uniform inflow. It has been shown before that the ALM and its wake profiles depend weakly on the subgrid-scale turbulence model in a finitevolume code (EllipSys3D) with inflow turbulence. ${ }^{39}$ Here, we extend this study using the JHU pseudo-spectral code LESGO, with homogeneous isotropic turbulence as inflow. A fringe region was used to prescribe a desired turbulent inflow. ${ }^{24}$

## A. Use of public database for inflow prescription

In this work, we develop a new approach to adapt data from a public database of forced isotropic turbulence as inflow boundary conditions. The data were stored with high spatial and temporal resolutions from a direct numerical simulation at $R e_{\lambda} \approx 420$ in a $1024^{3}$ periodic box for a duration of five large-eddy turnover times. ${ }^{13,40}$ The box size for the isotropic simulation is four to five times the integral length scale of the turbulence. The velocity field from the database is filtered (using the available box filtering web service getBoxFiltered ${ }^{13}$ ) at a length scale that is two times the LES grid resolution. The transformation of database velocities to velocities to be used in LES is done as follows. The subscript (LES) denotes the quantities from the simulations. The subscript DB denotes the quantities from the database. A desired turbulence intensity ( $T I$ ) is prescribed. It must obey the following conditions:

$$
\begin{equation*}
T I=\frac{u_{\mathrm{LES}}^{\prime}}{U_{\mathrm{LES}}}=\frac{u_{\mathrm{DB}}^{\prime}}{U_{\mathrm{DB}}}, \tag{4}
\end{equation*}
$$

where $u_{\text {LES }}^{\prime}$ is the desired RMS of the LES, $U_{\text {LES }}$ is the mean desired inflow velocity in the LES, $T I$ is the desired turbulence intensity, and $u_{\mathrm{DB}}^{\prime}$ is the RMS of the database (note that $u^{\prime}$ here denotes the RMS values rather than instantaneous fluctuations). The database simulates isotropic turbulence with no mean velocity, and so, any arbitrary mean velocity can be selected. To obtain the appropriate timescales, we must determine the "database" mean velocity $U_{\mathrm{DB}}$ that will yield the desired turbulence intensity. Thus, data from the database are obtained by sweeping through the domain at a sweep velocity of

$$
\begin{equation*}
U_{\mathrm{DB}}=\frac{u_{\mathrm{DB}}^{\prime}}{T I} . \tag{5}
\end{equation*}
$$

A time-evolving field is then extracted from the database and used as the inflow condition. The domain sizes for the database and LES simulation are $L_{\mathrm{DB}}=2 \pi$ and $L_{\mathrm{LES}}=6 D$, where $D$ is the rotor diameter. The database velocities are rescaled according to

$$
\begin{gather*}
u_{\mathrm{LES}}\left(\mathbf{x}_{\mathrm{LES}}, t_{\mathrm{LES}}\right)=U_{\mathrm{LES}}+u_{\mathrm{DB}}\left(\mathbf{x}_{\mathrm{DB}}, t_{\mathrm{DB}}\right) \frac{U_{\mathrm{LES}}}{U_{\mathrm{DB}}}  \tag{6}\\
v_{\mathrm{LES}}\left(\mathbf{x}_{\mathrm{LES}}, t_{\mathrm{LES}}\right)=v_{\mathrm{DB}}\left(\mathbf{x}_{\mathrm{DB}}, t_{\mathrm{DB}}\right) \frac{U_{\mathrm{LES}}}{U_{\mathrm{DB}}}  \tag{7}\\
w_{\mathrm{LES}}\left(\mathbf{x}_{\mathrm{LES}}, t_{\mathrm{LES}}\right)=w_{\mathrm{DB}}\left(\mathbf{x}_{\mathrm{DB}}, t_{\mathrm{DB}}\right) \frac{U_{\mathrm{LES}}}{U_{\mathrm{DB}}} \tag{8}
\end{gather*}
$$

where for a desired location and time in the LES ( $\mathbf{x}_{\text {LES }}, t_{\text {LES }}$ ), we must use the data available at the following location and time in the database (recalling that any position outside $[0,2 \pi]$ can be obtained using the periodicity of the data):

$$
\begin{equation*}
x_{\mathrm{DB}}=x_{\mathrm{LES}} \frac{2 \pi}{6 D}-U_{\mathrm{DB}} t_{\mathrm{DB}}, \tag{9}
\end{equation*}
$$

$$
\begin{align*}
y_{\mathrm{DB}} & =y_{\mathrm{LES}} \frac{2 \pi}{6 D}  \tag{10}\\
z_{\mathrm{DB}} & =z_{\mathrm{LES}} \frac{2 \pi}{6 D} \tag{11}
\end{align*}
$$

and

$$
\begin{equation*}
t_{\mathrm{DB}}=t_{\mathrm{LES}} \times \frac{2 \pi}{6 D} \times \frac{U_{\mathrm{LES}}}{U_{\mathrm{DB}}} . \tag{12}
\end{equation*}
$$

The turbulence in the inflow is expected to decay as it travels downstream. Using Kolmogorov-type arguments, the length scale over which the turbulence decays can be estimated as

$$
\begin{equation*}
L_{\mathrm{decay}} \sim \frac{L_{\mathrm{int}}}{T I}, \tag{13}
\end{equation*}
$$

where $L_{\text {decay }}$ is the length scale over which turbulence decays, $L_{\text {int }}$ is the integral length scale, and $T I$ is the turbulence intensity. The turbine is placed at one $L_{\text {int }}$ from the inflow plane. This means that the turbulence is still very active by the time it reaches the rotor. In the case of 5\% $T I$, the turbulence has only experienced $5 \%$ of its characteristic decay because $L_{\mathrm{int}}$ is only $5 \%$ of $L_{\text {decay }}$.

## B. Effect of the subgrid-scale model parameter with turbulent inflow

The inflow turbulence and Smagorinsky coefficient have no effect on the time-averaged velocity field and loads computed by the ALM. This is shown in Fig. 8, where the lines completely overlap. This result is consistent with the expectation that quantities computed by the ALM are dominated by the mean velocity. The mean aerodynamic forces in Fig. 8 are the same as in the case of laminar inflow shown in Fig. 2. In the case of homogeneous isotropic turbulence as inflow, the mean velocity is the same as in the case of laminar inflow $\left(U_{\infty}\right)$. The mean aerodynamic forces computed by the ALM depend on the mean inflow velocity and are not affected by turbulent fluctuations.

The wake of a wind turbine under turbulent inflow is very different from that of laminar inflow. Figure 9 shows a volume rendering of instantaneous streamwise velocity for cases of laminar and turbulent inflow. The wake of the turbine under uniform inflow slowly evolves and eventually becomes turbulent far downstream. In the case of turbulent inflow, the wake becomes turbulent much faster and meandering is more noticeable.

The time-averaged streamwise velocity and Reynolds stress components are shown in Figs. 10 and 11. The contours show that the flow fields are very similar for all cases regardless of the Smagorinsky coefficient. In this case, the wake dynamics are governed mostly by the

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-11.jpg?height=464&width=1275&top_left_y=2035&top_left_x=426)
FIG. 8. Nondimensional lift and axial velocity along the blade for cases with different $C_{S}$ and turbulent inflow with a turbulence intensity of $5 \%$. The blade is presented for reference.

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-12.jpg?height=373&width=1249&top_left_y=351&top_left_x=428)
FIG. 9. Volume rendering of instantaneous streamwise velocity ( $u / U_{\infty}$ ) for a case with laminar inflow (left) and turbulent inflow with a turbulence intensity of 5\% (right).

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-12.jpg?height=418&width=1058&top_left_y=870&top_left_x=536)
FIG. 10. Mean streamwise velocity contours for simulations with an inflow turbulence intensity of 5\%.

inflow turbulence, and the Smagorinsky coefficient has a negligible effect. In the case of a very low turbulence intensity ( $1 \%$ ), the transition effect of $C_{S}$ is also negligibly small, as shown in Fig. 12. Even at such low turbulence levels, the inflow turbulence provides finite-amplitude per- turbations that trigger transition in the wake, as opposed to the uniform inflow cases, where the subgrid-scale model is the main factor determining where transition to turbulence occurs.

The dynamically computed Smagorinsky coefficient depends strongly on the inflow turbulence. Figure 13 shows the time-averaged Smagorinsky coefficient computed by the dynamic scale-dependent model for a case with uniform laminar inflow and a case with turbulent inflow with $5 \%$ turbulence intensity. In these figures, the dynamically computed value $C_{s}(x, y, z, t)$ is averaged in time and displayed as a function of $x$ and $z$. Differences are observed mostly in the thin shear layers, where the uniform inflow case takes a longer time to transition to turbulence. Stronger shear layers will produce higher $C_{s}$ values. The laminar inflow cases have a stronger shear layer, resulting in a higher $C_{s}$ coefficient. As shown, at the entrance, $C_{s}$ is close to zero also in the case of the Lagrangian model. The reason for this can be traced to several effects

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-12.jpg?height=416&width=1058&top_left_y=2074&top_left_x=534)
FIG. 11. Streamwise normal component of the Reynolds stress tensor $\left(\left\langle u^{\prime} u^{\prime}\right\rangle / U_{\infty}^{2}\right)$ for simulations with an inflow turbulence intensity of 5\%.

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-13.jpg?height=417&width=1056&top_left_y=342&top_left_x=536)
FIG. 12. Streamwise normal component of the Reynolds stress tensor $\left(\left\langle u^{\prime} u^{\prime}\right\rangle / U_{\infty}^{2}\right)$ for simulations with an inflow turbulence intensity of $1 \%$.

![](https://cdn.mathpix.com/cropped/fb8395c5-751e-4499-99c4-910340486ff1-13.jpg?height=324&width=1060&top_left_y=905&top_left_x=534)
FIG. 13. Smagorinsky coefficient for the Lagrangian-scale-dependent model for the case with laminar inflow (top) and 5\% turbulence intensity (bottom).

that in this case combine to produce a very low $C_{s}$ at the entrance. First, we use the code default $C_{s}=0$ (strictly speaking the numerator in the dynamic model) as the inflow boundary condition for the Lagrangian model. If there are turbulent motions at scales between the test filter and the grid-scale, previous experience shows that the dynamic coefficient increases signifi- cantly on timescales on the order of a grid-scale turnover time. However, in this case, the inflow has been filtered at scales twice the grid scale. Hence, there is essentially no turbulence that is sensed by the dynamic model between scales $\Delta$ and $2 \Delta$. In addition, the scale-dependent version of the dynamic model used here also uses test filtering at scale $4 \Delta$, which now does pick up some turbulence. The scale-dependent model then senses a finite value of the coefficient at scale $4 \Delta$ and a very small or zero value at scale $2 \Delta$. The model is based on extrapolating this trend down to the grid scale $\Delta$, thus resulting in the negligible values observed at the entrance. The scale-dependent model also responds at grid-turnover timescales. Estimating the grid-scale turnover timescale as $\left(L_{\text {int }} / u^{\prime}\right)\left(\Delta / L_{\text {int }}\right)^{2 / 3}$ yields that the distance traveled is given by $\left(L_{\text {int }} / T I\right)\left(\Delta / L_{\text {int }}\right)^{2 / 3}$. Hence, for $\mathrm{TI}=5 \%$ and $\Delta / D=0.0315$, we obtain that it takes at least a distance of about $2 D$ to build up turbulence and finite $C_{s}$ values.

## VII. CONCLUSIONS

Large-eddy simulations of a wind turbine under laminar uniform inflow using an actuator line model (ALM) were performed with four different large-eddy-simulation research codes. The parameters used with all the codes were matched as closely as possible. Excellent agreement was observed for the quantities along the blades and in the near-wake predicted by the ALM. It was concluded that numerical discretization did not noticeably affect the along-blade quantities predicted by ALM as long as implementation was carefully replicated among the codes. The Smagorinsky coefficient in the subgrid-scale turbulence model was shown to not affect the quantities along the blade computed by the ALM. In the case of uniform inflow, even though the implementation of ALM in every code is similar, the turbulent characteristics of the wake are different. The transition to turbulence depends strongly on the numerical discretization used and on the Smagorinsky coefficient used in the subgrid-scale model.

A series of simulations was performed with homogeneous isotropic turbulence from a database as inflow condition. A special methodology for rescaling a time-evolving, spatio-temporal data set was applied to match desired turbulence intensity and evolution timescales. The Smagorinsky coefficient had no noticeable effect on the wake with turbulent inflow. The wake is very different from the wake under uniform inflow; it becomes fully turbulent much faster, and strong meandering is observed in the case with inflow turbulence. This turbulence is not influenced by the Smagorinsky coefficient. In the case of flows with shear (e.g., atmospheric boundary layer), it has been shown that the Smagorinsky coefficient has a strong influence on the mean velocity and turbulent fluctuations. ${ }^{16-18}$ For this reason, for cases in which the inflow mean and turbulent statistics must be matched realistically to data, the use of more sophisticated subgridscale models (such as the Lagrangian-scale-dependent dynamic model) is recommended.

## ACKNOWLEDGMENTS

Simulations for the JHU code LESGO were performed using computational resources at the Maryland Advanced Research Computing Center and XSEDE. L.A.M.T. and C.M. thank the National Science Foundation for support (Grant Nos. 1230788 and 1243482, the WINDINSPIRE project). P.L.J. was supported by a National Science Foundation Graduate Research Fellowship Program under Grant No. DGE-1232825. The Alliance for Sustainable Energy, LLC (Alliance) is the manager and operator of the National Renewable Energy Laboratory (NREL). NREL is a national laboratory of the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy. This work was authored by the Alliance and supported by the U. S. Department of Energy under Contract No. DE-AC36-08GO28308. Funding was provided by the U.S. Department of Energy Office of Energy Efficiency and Renewable Energy Wind Energy Technologies Office. The views expressed in the article do not necessarily represent the views of the U.S. Department of Energy or the U.S. government. The U.S. government retains, and the publisher, by accepting the article for publication, acknowledges that the U.S. government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. government purposes.Simulations for the NREL code SOWFA were performed using NREL's Peregrine high-performance computing system. L.A.M.T. would also like to thank Sheri Anstedt of NREL for her help editing the manuscript.

[^0]${ }^{27}$ J. Meyers and P. Sagaut, Phys. Fluids 19, 095105 (2007).
${ }^{28}$ J. Meyers and C. Meneveau, "Large eddy simulations of large wind-turbine arrays in the atmospheric boundary layer," AIAA Paper No. 2010-827, 2010.
${ }^{29}$ J. Michelsen, "Basis3D-a platform for development of multiblock PDE solvers," Ph.D. thesis (Technical University of Denmark, Department of Fluid Mechanics, 1992); Technical Note AFM92-05.
${ }^{30}$ N. Sørensen, "General purpose flow solver applied to flow over hills," Ph.D. thesis (Technical University of Denmark, Risø National Laboratory for Sustainable Energy, 1995).
${ }^{31}$ J. Jonkman, S. Butterfield, W. Musial, and G. Scott, Technical Report No. NREL/TP-500-38060, National Renewable Energy Laboratory, Golden, CO, 2009.
${ }^{32}$ L. Martínez-Tossas, M. Churchfield, and S. Leonardi, Wind Energy 18, 1047 (2015).
${ }^{33}$ L. Martínez-Tossas, M. Churchfield, and C. Meneveau, Wind Energy 20, 1083 (2017).
${ }^{34}$ L. A. Martínez-Tossas, M. J. Churchfield, and C. Meneveau, J. Phys.: Conf. Ser. 753, 082014 (2016).
${ }^{35}$ M. J. Churchfield, S. Schreck, L. A. Martínez-Tossas, C. Meneveau, and P. R. Spalart, in 35th Wind Energy Symposium (2017), p. 1998.
${ }^{36}$ N. Troldborg, "Actuator line modeling of wind turbine wakes," Ph.D. thesis (Technical University of Denmark, 2009).
${ }^{37}$ M. O. L. Hansen, Aerodynamics of Wind Turbines, 2nd ed. (Routledge, London; Sterling, VA, 2007).
${ }^{38}$ P. Moin and R. Verzicco, "Vortical structures and wall turbulence," Eur. J. Mech. B 55, 242 (2016).
${ }^{39}$ H. Sarlak, C. Meneveau, and J. Sørensen, Renewable Energy 77, 386 (2015).
${ }^{40}$ E. Perlman, R. Burns, Y. Li, and, C. Meneveau, in Proceedings of the 2007 ACM/IEEE Conference on Supercomputing (ACM, 2007), p. 23.


[^0]:    ${ }^{1}$ M. Churchfield, S. Lee, P. Moriarty, L. Martínez Tossas, S. Leonardi, G. Vijayakumar, and J. Brasseur, "A large-Eddy simulation of wind-plant aerodynamics," AIAA Paper No. 2012-0537.
    ${ }^{2}$ Y. Wu and F. Porté-Agel, Boundary-Layer Meteorol. 146, 181 (2013).
    ${ }^{3}$ M. Calaf, C. Meneveau, and J. Meyers, Phys. Fluids 22, 015110 (2010).
    ${ }^{4}$ J. Sørensen and W. Shen, J. Fluids Eng. 124, 393 (2002).
    ${ }^{5}$ N. Troldborg, J. Sørensen, and R. Mikkelsen, Wind Energy 13, 86 (2010).
    ${ }^{6}$ L. Martínez-Tossas, M. Churchfield, and C. Meneveau, J. Phys.: Conf. Ser. 625, 012024 (2015).
    ${ }^{7}$ P. Krogstad and P. Eriksen, Renewable Energy 50, 325 (2013).
    ${ }^{8}$ F. Pierella, P. Krogstad, and L. Sætran, Renewable Energy 70, 62 (2014).
    ${ }^{9}$ P. Krogstad, L. Sætran, and M. Adaramola, J. Fluids Struct. 52, 65 (2015).
    ${ }^{10}$ L. Lignarolo, D. Mehta, R. Stevens, A. Yilmaz, G. van Kuik, S. Andersen, C. Meneveau, C. Ferreira, D. Ragni, J. Meyers et al., Renewable Energy 94, 510 (2016).
    ${ }^{11}$ H. Sarlak, F. Pierella, R. Mikkelsen, and J. N. Sørensen, J. Phys.: Conf. Ser. 524, 012145 (2014).
    ${ }^{12}$ H. Sarlak, T. Nishino, L. Martínez-Tossas, C. Meneveau, and J. N. Sørensen, Renewable Energy 93, 340 (2016).
    ${ }^{13}$ Y. Li, E. Perlman, M. Wan, Y. Yang, C. Meneveau, R. Burns, S. Chen, A. Szalay, and G. Eyink, J. Turbul. 9, N31 (2008).
    ${ }^{14}$ S. Pope, Turbulent Flows (Cambridge University Press, 2001).
    ${ }^{15}$ M. Germano, J. Fluid Mech. 238, 325 (1992).
    ${ }^{16}$ C. Meneveau, T. S. Lund, and W. Cabot, J. Fluid Mech. 319, 353 (1996).
    ${ }^{17}$ F. Porté-Agel, C. Meneveau, and M. Parlange, J. Fluid Mech. 415, 261 (2000).
    ${ }^{18}$ E. Bou-Zeid, C. Meneveau, and M. Parlange, Phys. Fluids 17, 025105 (2005).
    ${ }^{19}$ J. Meyers, C. Meneveau, and B. Geurts, Phys. Fluids 22, 125106 (2010).
    ${ }^{20}$ J. Graham and C. Meneveau, Phys. Fluids 24, 125105 (2012).
    ${ }^{21}$ C. Moeng, J. Atmos. Sci. 41, 2052 (1984).
    ${ }^{22}$ J. Albertson and M. Parlange, Water Resources Res. 35, 2121 (1999).
    ${ }^{23}$ C. Canuto, M. Hussaini, A. Quarteroni, and T. Zang, Fundamentals in single domains (2006).
    ${ }^{24}$ R. Stevens, J. Graham, and C. Meneveau, Renewable Energy 68, 46 (2014).
    ${ }^{25}$ M. Churchfield and S. Lee, SOWFA/NWTC Information Portal, https://nwtc.nrel.gov/SOWFA.
    ${ }^{26}$ H. Weller, G. Tabor, H. Jasak, and C. Fureby, Comput. Phys. 12, 620 (1998).

