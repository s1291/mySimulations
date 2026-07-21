# A vortex-based tip/smearing correction for the actuator line 

Alexander R. Meyer Forsting, Georg Raimund Pirrung, and Néstor Ramos-García<br>DTU Wind Energy, Technical University of Denmark, Frederiksborgvej 399, 4000 Roskilde, Denmark<br>Correspondence: Alexander R. Meyer Forsting (alrf@dtu.dk)

Received: 21 December 2018 - Discussion started: 21 January 2019
Revised: 13 May 2019 - Accepted: 16 May 2019 - Published: 28 June 2019


#### Abstract

The actuator line (AL) was intended as a lifting line (LL) technique for computational fluid dynamics (CFD) applications. In this paper we prove - theoretically and practically - that smearing the forces of the actuator line in the flow domain forms a viscous core in the bound and shed vorticity of the line. By combining a near-wake representation of the trailed vorticity with a viscous vortex core model, the missing induction from the smeared velocity is recovered. This novel dynamic smearing correction is verified for basic wing test cases and rotor simulations of a multimegawatt turbine. The latter cover the entire operational wind speed range as well as yaw, strong turbulence and pitch step cases. The correction is validated with lifting line simulations with and without viscous core, which are representative of an actuator line without and with smearing correction, respectively. The dynamic smearing correction makes the actuator line effectively act as a lifting line, as it was originally intended.


## 1 Introduction

The actuator line (AL) technique developed by Sørensen and Shen (2002) is a lifting line (LL) representation of the wind turbine rotor suitable for computational fluid dynamics (CFD) simulations. It captures transient physical features like shed and trailed vorticity (including root/tip vortices), without the computational cost associated with resolving the full rotor geometry. Thus, the AL model enables large-eddy simulations (LES) of large wind farms in realistic, turbulent atmospheric boundary layers (Vollmer et al., 2017).

However, in contrast to LL vortex formulations, the blade forces are dispersed in the flow domain - most commonly in form of a Gaussian projection - to avoid numerical instabilities. A length scale - also referred to as smearing coefficient - controls this force redistribution, the lower limit of which is linked to the grid size by numerical stability requirements (Troldborg et al., 2009). Mikkelsen (2003) observed a large sensitivity of the blade velocities to this length scale, which consequently also propagated to the blade forces. Especially in regions along the blade exhibiting stark load changes, such as around the root and tip, forces are substantially overpredicted, meaning that this effect is exacerbated by non-
tapered and low-aspect-ratio blades. As actuator disc formulations suffer from similar issues towards the blade tip, Glauert (1935) type tip corrections are also frequently applied to ALs (Shen et al., 2005). However, these correct discs for missing discrete blades and should therefore be unnecessary - strictly even invalid - for ALs. Shives and Crawford (2013) and Jha et al. (2014) achieved a reduction in the force over-prediction by varying the originally fixed smearing factor with respect to the blade chord. However, their methods cannot decouple the blade forces from the smearing length scale: a smeared force distribution in the flow domain unavoidably leads to lower induction at the blade - increasing lift and drag - compared with an actual LL with a concentrated, spatially singular force.

### 1.1 The vortex smearing hypothesis

Shives and Crawford (2013) noticed the similarity between the velocities induced across an actuator line and those predicted by a viscous vortex core model. These models include the limiting effect of viscous shear forces on the induced velocities around vortex cores. A similar comparison of the swirl velocities about an infinite vortex line is shown in

Fig. 1 - here with a Lamb-Oseen vortex core model (Lamb, 1932; Oseen, 1911). Without viscosity (inviscid) the velocities approach infinity towards the vortex centre. The startling agreement between the Lamb-Oseen and AL velocities was first demonstrated by Dag et al. (2017). The Gaussian body force smearing in the AL technique thus produces similar swirl velocities to a viscous vortex. Ignoring viscous effects, the AL should, in principle, induce the same velocities as a LL-equivalent to the inviscid solution. The missing induced velocity in the AL model (shaded area in Fig. 1) can be approximated following Dag et al. (2017) as follows:

$$
\begin{align*}
\Delta v_{\theta}(r) & =\overbrace{\frac{\Gamma}{2 \pi r}}^{\text {inviscid }}-\overbrace{\frac{\Gamma}{2 \pi r}\left[1-\exp \left(-r^{2} / \epsilon^{2}\right)\right]}^{\text {viscous core }} \\
& =\frac{\Gamma}{2 \pi r} \exp \left(-r^{2} / \epsilon^{2}\right) \tag{1}
\end{align*}
$$

where $\Gamma$ represents the vortex line's circulation, $r$ is the distance from the vortex core and $\epsilon$ is the length scale used in the force smearing. This formulation can be split into an inviscid and viscous/smearing contribution:

$$
\begin{align*}
& \Delta v_{\theta}(r)=\overbrace{v_{\theta}(r)}^{\text {inviscid }} \underbrace{f_{\epsilon}\left(r_{\epsilon}\right)}_{\text {smearing }} \text { with } v_{\theta}(r)=\frac{\Gamma}{2 \pi r}, \\
& f_{\epsilon}\left(r_{\epsilon}\right)=\exp \left(-r_{\epsilon}^{2}\right), \quad r_{\epsilon}=\frac{r}{\epsilon} . \tag{2}
\end{align*}
$$

If this viscous behaviour of the force smearing in AL simulations was limited to the bound vortex representing the blade, it would not influence the blade forces as long as the blade was straight. However, Dag et al. (2017) argued that the trailing vortices (in the wake) exhibit the same viscous core, as they originate from the bound vortex. Hence, the wake of an AL induces lower velocities at the blade than a LL. The missing velocity can be estimated from the viscous core equivalence and thus correct the velocities at the blade. This mostly impacts blade forces by changing the angle of attack at the blade sections.

### 1.2 Contributions of this paper

Dag et al. (2017) corrected AL simulations of a rectangular wing and two rotors with different aspect ratios by recuperating the missing induced velocity introduced by the viscous core. For all of their simulations they were able to show the beneficial effect of the correction on the blade load distribution - represented by more physical behaviour, especially towards the tip and root. However, their implementation of the correction did not fully couple the flow field with the blade forces and the induction correction.

The major contributions of this paper are as follows:

- The development of a tuning-free, dynamic and numerically robust smearing correction, which is fully coupled to the AL model.

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-02.jpg?height=770&width=606&top_left_y=225&top_left_x=1182)
Figure 1. Distribution of the tangential velocity component in a plane orthogonal to an infinite vortex line (along $x$ ) obtained from either an inviscid or viscous (Lamb-Oseen) theoretical vortex and an actuator line (AL) CFD simulation.

- A theoretical proof of the force smearing - vortex core equivalence.
- Proof of the vortex core inheritance in trailed vorticity.
- The confirmation of the missing velocity assumption by comparing LL simulations with/without viscous core and AL results with/without correction.

The test cases include constantly and elliptically loaded wings as well as rotor simulations of a multimegawatt turbine covering the entire operational wind speed range. As the AL model is especially attractive for wind farm simulations, the focus here is on coarsely resolved ALs. The correct dynamic behaviour of the new correction is verified through yawed inflow and pitch step simulations.

## 2 Proof of force smearing - vortex core equivalence

The equivalence between the velocity field induced by an AL and a viscous vortex can be derived directly from the incompressible Navier-Stokes equations. This proof follows the approach by Forsythe et al. (2015) that successfully connected an AL's vorticity field to its force projection. Starting by taking the curl of the incompressible momentum equation, the vorticity transport equation is obtained ( $\boldsymbol{\omega}=\nabla \times \boldsymbol{u}$ ):

$$
\begin{equation*}
\underbrace{\frac{\partial \boldsymbol{\omega}}{\partial t}}_{\substack{=0 \\ \text { steady }}}+(\boldsymbol{u} \nabla) \boldsymbol{\omega}=\underbrace{(\boldsymbol{\omega} \nabla) \boldsymbol{u}}_{\substack{=0 \\ 2 D}}+\underbrace{v \nabla^{2} \boldsymbol{\omega}}_{\substack{=0 \\ \text { inviscid }}}+\nabla \frac{\boldsymbol{f}}{\rho}, \tag{3}
\end{equation*}
$$

where $v$ is the viscosity, $\rho$ is density and $\boldsymbol{f}$ represents the body forces from the AL. Away from the root and tip, the flow around a high-aspect-ratio blade is nearly twodimensional, as the span-wise flow is negligible. Viscous effects are disregarded in light of the large Reynolds numbers encountered. Furthermore the relationship between the body force and flow field becomes quasi-steady assuming the flow is attached. Cancelling the respective terms and noting that in two-dimensional flow ( $y-z$ plane) $\boldsymbol{\omega}=\omega_{x} \hat{\boldsymbol{e}}_{x}$ and $\nabla=\left(0, \frac{\partial}{\partial y}, \frac{\partial}{\partial z}\right)$ :
$(\boldsymbol{u} \nabla) \omega_{x} \hat{\boldsymbol{e}}_{x}=\nabla \frac{\boldsymbol{f}}{\rho}$.

Assuming the drag to be negligible, the force - in the form of lift - exerted by the AL on the flow in terms of its circulation $\Gamma$ becomes

$$
\begin{equation*}
\boldsymbol{f}=-\boldsymbol{f}_{\text {aero }} g(r) \tag{5}
\end{equation*}
$$

$\boldsymbol{f}_{\text {aero }}=L=\rho \boldsymbol{u} \times \Gamma \hat{\boldsymbol{e}}_{x} \quad g(r)=\frac{1}{\pi \epsilon^{2}} \exp \left(-r^{2} / \epsilon^{2}\right)$.

Here $g$ represents a two-dimensional Gaussian force projection with $r$ indicating the distance from the AL. Inserting these expression into Eq. (4) and exploiting standard matrix transformation and mass conservation ${ }^{1}$

$$
\begin{equation*}
(\boldsymbol{u} \nabla) \omega_{x} \hat{\boldsymbol{e}}_{x}=\nabla\left(\Gamma g \hat{\boldsymbol{e}}_{x} \times \boldsymbol{u}\right)=(\boldsymbol{u} \nabla) \Gamma g \hat{\boldsymbol{e}}_{x} \tag{7}
\end{equation*}
$$

$(\boldsymbol{u} \nabla) \omega_{x}=(\boldsymbol{u} \nabla) \Gamma g$.

Due to mass conservation the $\boldsymbol{u} \nabla$ term can be inverted; this gives a direct relationship between the force projection and vorticity:
$\omega_{x}=\Gamma g=\frac{\Gamma}{\pi \epsilon^{2}} \exp \left(-r^{2} / \epsilon^{2}\right)$.

As the body force is axially symmetric, the vorticity only induces tangential velocities

$$
\begin{align*}
\omega_{x}(r) & =\frac{1}{r}(\frac{\partial r u_{\theta}}{\partial r}-\underbrace{\frac{\partial u_{r}}{\partial \theta}}_{\substack{=0 \\
\text { axisymmetry }}}) \\
& \Rightarrow u_{\theta}=\frac{1}{r} \int_{0}^{r} r \omega_{x}(r) \mathrm{d} r . \tag{10}
\end{align*}
$$

Inserting Eq. (9) and integrating gives the swirl velocity induced by a smeared body force
$u_{\theta}=\frac{\Gamma}{2 \pi r}\left[1-\exp \left(-r^{2} / \epsilon^{2}\right)\right]$.

[^0]This expression equals that of the Lamb-Oseen vortex, only with the viscous core radius replaced by the smearing coefficient ${ }^{2}$. This marks the theoretical confirmation of the observations by Dag et al. (2017), which additionally indicates that a viscous core behaviour with an AL requires inviscid, two-dimensional and locally steady flow conditions.

## 3 Numerical methodology

### 3.1 Actuator line simulations

The discretized incompressible Navier-Stokes equations are solved using DTU's CFD code EllipSys3D (Sørensen, 1995; Michelsen, 1994a, b). The flow is iteratively solved at each time instant by the SIMPLE algorithm (Patanker and Spalding, 1972). Depending on the turbulence model either the third-order accurate QUICK (Leonard, 1979) scheme (RANS) or a fourth-order CDS scheme (LES) discretizes the convective terms. As the flow variables are located at the cell centres, a modified Rhie and Chow (Réthoré and Sørensen, 2012) algorithm avoids pressure-velocity decoupling. Further details regarding the numerical techniques are given in Meyer Forsting et al. (2017). For all comparisons with the LL code, the RANS equations are solved using the $k-\omega$ shearstress transport turbulence closure of Menter (1993). Only the turbulent inflow cases in Sect. 5.2.4 are computed using the DES technique of Strelets (2001). The AL model was implemented by Mikkelsen (2003) in EllipSys3D. We employ a version utilizing three-dimensional Gaussian force projection, which follows the original formulation of Sørensen and Shen (2002). As the AL model is especially attractive for wind farm simulations, the focus here is on coarsely resolved ALs, with either 9 or 19 sections ( $N_{s}$ ) along the blade. They are uniformly spaced and discretize the blade starting from the root at 1.5 m to the tip at 63 m . The smearing length scale is connected to the number of sections, such that $\epsilon=2 R /\left(N_{s}+1\right)-R$ defining the rotor radius - which ensures that the forces in the domain change smoothly between sections (Nathan, 2018). The tower and nacelle are not modelled.

The numerical domain for the rotor simulations is discretized in a verified, standard manner (Meyer Forsting et al., 2017; Troldborg et al., 2009). It consists of a box with $25 R$ side length that contains an inner box with a uniformly spaced refined mesh of $3.2 R$ edge length at its centre surrounding the rotor (see Fig. 2). To capture the velocity gradients around the AL correctly the mesh spacing is $\Delta x= R / 40$. This is twice the recommended minimum (Troldborg et al., 2009); however, it delivers more accurate angle of attack estimates at the section centres (Shives and Crawford, 2013). In total, 256 cells discretize the flow domain along each dimension, resulting in $16.8 \times 10^{6}$ degrees of freedom. All variables, except pressure and its correction, which ne-

[^1]![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-04.jpg?height=783&width=836&top_left_y=225&top_left_x=164)
Figure 2. Numerical box domain with a structured mesh and uniform spacing around the rotor at its centre. Only every eighth grid point is shown.

cessitate special treatment (Sørensen, 1995), obey symmetry conditions on the lateral boundaries, whereas at the inflow and outflow faces they follow Dirichlet and Neumann conditions, respectively. The wing test cases follow the same approach only with $N_{s}=32$ and an inner box edge length of $3 b$, where $b$ is the wing's half-span. This results in 80 cells along each dimension and $5.1 \times 10^{5}$ degrees of freedom. To ensure that the blade tip remains inside a single cell during one time step $\Delta t<\Delta x /(\Omega R)$ - with mesh spacing $\Delta x$ and rotational speed $\Omega$. Without rotation, the term $\Omega R$ is replaced by the advection speed of the wake. The kinematic viscosity and air density are kept constant at $1.789 \times 10^{-5} \mathrm{~kg} \mathrm{~m}^{-1} \mathrm{~s}^{-1}$ and $1.225 \mathrm{~kg} \mathrm{~m}^{-3}$, respectively. Simulations are stopped when the thrust residual reaches $1 \times 10^{-5}$.

The sensitivity of the rotor thrust to the domain size, time step and grid size is explored in Fig. 3. The length of the domain edges is doubled to $50 R$, the time step is halved with respect to a set-up obeying the method described above. A simulation of the NREL 5-MW at $8 \mathrm{~ms}^{-1}$ with either 40 or 60 grid cells along the rotor depending on the smearing length scale acts as a reference. With $\epsilon=R / 10$ and $R / 20$, this represents 4 and 3 times the recommended resolution, respectively (Troldborg et al., 2009). Although non-zero, the sensitivity of the results is acceptable in code comparison and should impact AL simulations with and without correction similarly.

### 3.2 Free-wake lifting line rotor simulations

The in-house solver MIRAS has been employed to perform the free-wake lifting line simulations. MIRAS is a multi-

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-04.jpg?height=417&width=608&top_left_y=225&top_left_x=1180)
Figure 3. Thrust sensitivity of the NREL 5-MW AL simulations at a $8 \mathrm{~ms}^{-1}$ wind speed with respect to grid size, doubling the domain size and halving the time step at two smearing length scales ( $N_{s}= \{9,19\}, T_{\text {ref }}=\{4.20,4.06\} \times 10^{5} \mathrm{~N}$ ).

fidelity computational vortex model, which is mainly used for predicting the aerodynamic behaviour of wind turbines and their wakes. It has been developed at DTU over the last decade and been extensively validated for small to large size wind turbine rotors by Ramos-García et al. (2017, 2014a, b).

The free-wake vortex method essentially models the wake of a wind turbine using a bundle of infinitely thin vortex filaments. To avoid numerical singularities, a viscous core must be introduced, which represents a more physical distribution of the velocities induced by each vortex filament, desingularizing the Biot-Savart law near the centre of the filament. The velocity induced by each one of the elements is obtained directly by evaluating the Biot-Savart law, and by summing the velocity induced by all filaments, the total wake induction is obtained as follows:
$\boldsymbol{u}\left(\boldsymbol{x}_{i}\right)=\sum_{j=1}^{N} K_{i j} \frac{\gamma_{j}}{4 \pi} \frac{\boldsymbol{t}_{j} \times \boldsymbol{r}_{i j}}{r_{i j}^{3}}$,
where $K_{i j}=\frac{r_{i j}^{2}}{\left(\varepsilon_{j}^{2 z}+r_{i j}^{2 z}\right)^{1 / z}}$,
$N$ is total number of filaments that form the wake, $\boldsymbol{r}_{i j}= \boldsymbol{x}_{i}-\boldsymbol{y}_{j}$ is the distance vector from the vortex element $\boldsymbol{y}_{j}$ to the evaluation point $\boldsymbol{x}_{i}, \gamma_{j}$ is the circulation of the filament, $\boldsymbol{t}_{j}$ is the unit orientation vector of the $j$ th filament and $r_{i j}=\left|\boldsymbol{r}_{i j}\right|, \varepsilon_{j}$ is the vortex core radius of the filament, and $z$ defines the cut-off velocity profile where the LambOseen model (Lamb, 1932; Oseen, 1911), $z=2$, has been employed.

A viscous core model is applied to emulate the effect of viscosity by changing the vortex core radius as a function of time (Leishman et al., 2002):
$\varepsilon_{i}(t)=\sqrt{4 \alpha_{v} \delta_{v} v t_{i}}+\varepsilon_{0}$,
where $\alpha_{v}$ is a constant set to 1.25643 (Ananthan and Leishman, 2004), $v$ is the kinematic viscosity and $t_{i}$ is the time elapsed since the generation of the $i$ th filament. In order to
represent the diffusive timescales, the viscous core radius is set to change with the vortex age by adding a turbulence eddy viscosity, $\delta_{v}$, first proposed by Squire (1965), and in this work set to $1 \times 10^{-3}$. To avoid the singular behaviour of newly released vortex elements, an initial core radius, $\varepsilon_{0}$, is introduced. In accordance with Ramos-García et al. (2017), who found that a small core radius is necessary to have flow convergence, a core radius of $0.1 \%$ of the local chord at the release station is used.

For the sake of the present study, two different approaches to compute the angle of attack have been followed.

- Inviscid (LL), where the non-regularized Biot-Savart law is used to compute the induction from the wake filaments at the quarter-chord location. This is the standard method used in a lifting line solvers.
- Viscous (LL+core), where the regularized Biot-Savart law is used to compute the induction at the quarterchord location. A viscous core with a radius equal to the actuator line smearing factor is used for a direct comparison of the methods.

This enables a double validation of the models. On the one hand the corrected AL simulations can be validated against the LL calculations, and on the other hand the raw AL model, without tip correction, can be compared against the LL+core simulations which include the smearing effect in the freewake model.

## 4 Tip/smearing correction for the actuator line

Applying the velocity correction methodology introduced in Sect. 1.1 in three-dimensional space yields a velocity correction vector. The viscous core behaviour of the AL bound vorticity - proven in Sect. 2 to originate from the force smearing - is inherited by the trailing vortices, as will be demonstrated in Sect. 5.1.1. Therefore, the induction from the trailed vorticity at the blade is lower than without force smearing. Figure 4 shows the path of trailed vorticity shed from in-between two sections of a blade with a strength of $\Delta \Gamma=\Gamma_{s}-\Gamma_{s+1}$ with
$\Gamma_{s}=\frac{1}{2} \sqrt{v_{s}^{2}+w_{s}^{2}} C_{\mathrm{L}}(\alpha) c$.

Here $s$ defines the blade section index, $C_{\mathrm{L}}$ is the sectional lift coefficient, $c$ is the section chord and $\alpha$ is the angle of attack, which depends on the inflow angle in combination with blade pitch and twist at the section. The missing induction from this single trailed vortex at a point $C$ is obtained by integrating along the vortex line
$\boldsymbol{u}^{*}=\int_{0}^{\infty} f_{\epsilon} \delta \widetilde{\boldsymbol{u}} \mathrm{d} l$,

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-05.jpg?height=442&width=608&top_left_y=221&top_left_x=1180)
Figure 4. Trailed vorticity path. The blade rotates in the $x-y$ plane and $z$ points downstream. The vortex element $\delta l$ with a strength of $\Gamma_{s}-\Gamma_{s+1}$ is shed at $r$ and transported downstream by the local velocity. The distance from the shedding location $r$ to a point $C$ along the blade is $h\left(h=r-C_{x}\right)$, where $\delta l$ induces tangential and axial velocities.

where $l$ is the vortex following coordinate. Here $\delta \widetilde{\boldsymbol{u}}$ is the velocity induced by an infinitesimal element $\boldsymbol{\delta} \boldsymbol{l}$ of the vortex line at point $C$, which is given by the Biot-Savart law:
$\delta \widetilde{\boldsymbol{u}}=\frac{\Delta \Gamma}{4 \pi} \frac{\delta \boldsymbol{l} \times \boldsymbol{x}}{|\boldsymbol{x}|^{3}}$,
where $\boldsymbol{x}$ is the vector pointing from the element towards $C$. The smearing factor for this vortex element becomes
$f_{\epsilon}=\exp \left(-\frac{\left(\boldsymbol{x} \hat{\boldsymbol{e}}_{\perp}\right)^{2}}{\epsilon^{2}}\right)$.

The viscous core only acts in the plane orthogonal to the vortex element $\boldsymbol{\delta} \boldsymbol{l}$; hence, $\hat{\boldsymbol{e}}_{\perp}$ projects $\boldsymbol{x}$ onto this plane. This is different to using the distance, $|\boldsymbol{x}|$, as Dag et al. (2017) proposed, which violates the two-dimensional nature of the viscous core.

The total missing induction at a blade section $s$ is obtained by summing the contribution from all trailed vortices. The number of trailed vortices $N_{v}$ is directly related to the number of blade sections $N_{v}=N_{s}+1$. Discretizing the vortices in time, the missing induction at a certain blade section becomes
$\boldsymbol{u}_{s}^{*}=\sum_{v}^{N_{v}} \sum_{n}^{N_{t}} f_{s, v}^{n} \boldsymbol{\Delta} \widetilde{\boldsymbol{u}}_{s, v}^{n}$.

Here $v$ denotes the trailed vortex index, $n$ is the time index and $N_{t}$ is the number of time steps. Note that $n=1$ is the most recently shed vortex element. As a tip/smearing correction should remain computationally inexpensive, numerically solving the Biot-Savart law in Eq. (16) to obtain $\boldsymbol{\Delta} \widetilde{\boldsymbol{u}}_{s, v}^{n}$ is unfeasible. This would necessitate $N_{t} N_{v} N_{s}$ or $N_{t}\left(N_{s}+1\right) N_{s}$ evaluations. An accurate, yet fast, alternative to solving the Biot-Savart law directly is the near-wake model (NWM) for trailed vorticity by Pirrung et al. $(2016,2017 b)$, which also includes downwind convection. It performs well
for dynamic flow cases and exhibits great numerical stability, as it was originally developed to enhance the aerodynamic accuracy of blade element momentum (BEM) models. Its formulation is based on a lifting line representation of the blade's trailed vorticity as depicted in Fig. 4 and approximates the induced velocities from a single trailed vortex line by two indicial functions. The velocity induced by a vortex element is given in the NWM as
$\boldsymbol{\Delta} \widetilde{\boldsymbol{u}}_{s, v}^{n}=\left(X_{s, v}^{n}+Y_{s, v}^{n}\right)\left[\begin{array}{c}0 \\ \sin \left(\phi^{n}\right) \\ -\cos \left(\phi^{n}\right)\end{array}\right]$,
with $\phi$ representing the helix angle of the vorticity shed in the CFD domain (see Fig. 4). The indicial functions take the following form:

$$
\begin{align*}
\left\{X_{s, v}^{n}, Y_{s, v}^{n}\right\} & =a_{\{X, Y\}} \frac{r_{v}}{4 \pi h_{s}\left|h_{s}\right|} \Delta \Gamma_{v}^{n} \phi_{s, v}^{*^{n}} \\
& {\left[1-\exp \left(-b_{\{X, Y\}} \frac{\Delta \beta^{*^{n}}}{\phi_{s, v}^{*^{n}}}\right)\right] } \\
& \exp \left(-b_{\{X, Y\}} \sum_{i}^{n-1} \frac{\Delta \beta^{*^{i}}}{\phi_{s, v}^{*^{i}}}\right) . \tag{20}
\end{align*}
$$

The definitions of $a_{\{X, Y\}}, b_{\{X, Y\}}, \beta^{*}$ and $\phi^{*}$ are those of Pirrung et al. $(2016,2017 \mathrm{~b})$. The indicial functions allow for the solution to be time-advanced by a mere multiplication, considerably reducing the model evaluations to $N_{v} N_{s}+N_{v}$. In the original formulation this removes the need for bookkeeping; however, as the smearing factor also changes with the position of the vortex element, all previously shed elements are advanced individually in this specific implementation. This is only an experimental feature for testing the smearing correction and should be simple to remove in a future, more practical implementation.

Following the lifting line formulation of the NWM shown in Fig. 4 (refer to Appendix A for a detailed mathematical description) the perpendicular distance from the vortex element to $C$ becomes

$$
\begin{align*}
& \boldsymbol{x}_{\perp}=\frac{\boldsymbol{\delta} \boldsymbol{l}}{|\boldsymbol{\delta} \boldsymbol{l}|} \times \boldsymbol{x}= \\
& r \cos \phi\left(\begin{array}{c}
\tan \phi(\beta \cos \beta-\sin \beta) \\
-\tan \phi(-1+h / r+\cos \beta+\beta \sin \beta) \\
-1+(1-h / r) \cos \beta
\end{array}\right) . \tag{21}
\end{align*}
$$

Thus, the smearing factor becomes

$$
\begin{equation*}
f_{s, v}^{n}=\exp \left(-\frac{\left|\boldsymbol{x}_{\perp}\left(r_{v}, \beta^{n}, h_{s}, \phi^{n}\right)\right|^{2}}{\epsilon^{2}}\right) . \tag{22}
\end{equation*}
$$

When discretizing in time, $\beta^{n}$ is taken to be at the mid-point of the vortex element.

Finally the missing velocities computed in Eq. (18) correct the original velocities from the CFD simulations

$$
\begin{equation*}
\boldsymbol{u}_{s}=\boldsymbol{u}_{s}^{\mathrm{CFD}}+\boldsymbol{u}_{s}^{*} . \tag{23}
\end{equation*}
$$

Therefore, the correction influences the blade forces through the angle of attack and the velocity magnitude, although it is the former that dominates. It also changes the circulation at each blade section through Eq. (14) and, thus, the shed vorticity and its induction. Hence determining the correction velocity is an iterative procedure. The correction algorithm is executed after the flow field is solved and takes the following form:

1. Interpolate the velocity vector $\boldsymbol{u}^{\mathrm{CFD}}$ at the section centres from the CFD flow field.
2. Compute the helix angles $\phi$, where $\phi_{v}= -\tan ^{-1}\left(\frac{w_{v-1}^{\mathrm{CFD}}+w_{v}^{\mathrm{CFD}}}{v_{v-1}^{\mathrm{CFD}}+v_{v}^{\mathrm{CFD}}}\right)$ and $\phi_{\left\{1, N_{v}\right\}}=-\tan ^{-1}\left(\frac{w_{\left\{1, N_{s}\right\}}^{\mathrm{CFD}}}{v_{\left\{1, N_{s}\right\}}^{\mathrm{CFD}}}\right)$.
3. Combine the CFD velocities with the respective correction from the previous time step $n-1$, such that $\boldsymbol{u}_{n}=\boldsymbol{u}_{n}^{\mathrm{CFD}}+\boldsymbol{u}_{n-1}^{*}$.
4. Compute the smearing factor $f_{\epsilon}$ for all time steps, sections and elements (Eq. 22).
5. Determine the angle of attack and velocity magnitude from $\boldsymbol{u}_{n}$ to determine $\Gamma_{s}$ (Eq. 14) ${ }^{3}$,
6. Compute the velocities from the newly released vortex element $\boldsymbol{\Delta} \widetilde{\boldsymbol{u}}_{n}$ (Eq. 19).
7. At the first iteration of each time step, advance the previous elements in time.
8. Compute the velocity correction at the current time step $\boldsymbol{u}_{n}^{*}$ (Eq. 18).
9. Update the velocity at the sections with some form of relaxation $\boldsymbol{u}_{n}=\boldsymbol{u}_{n}^{\mathrm{CFD}}+\boldsymbol{u}_{n}^{*}$.
10. Repeat steps 5-9 until convergence is reached.

We use the technique by Pirrung et al. (2017a) established for the NWM to accelerate and ensure its convergence. Furthermore, the activation of the correction is delayed until the starting vorticity of the rotor has been transported at least one blade length away from the rotor plane. This enhances its numerical stability, as induction has already built up at the blades by its time of activation.

## 5 Results

### 5.1 Basic wing test cases

To verify the smearing hypothesis (Sect. 1.1) and the novel smearing correction (Sect. 4) two basic wing flow cases with known theoretical solutions are modelled using CFD. Either

[^2]![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-07.jpg?height=250&width=704&top_left_y=223&top_left_x=230)
Figure 5. Definition of the wing test cases with either a rectangular or elliptic planform. Vortices are trailed in-between sections and the actuator line forces are computed and exerted at the sections' centres.

Table 1. Input parameters common to both rectangular and elliptic wing simulations.
| $w_{\infty}\left(\mathrm{m} \mathrm{s}^{-1}\right)$ | $2 b(\mathrm{~m})$ | $r_{v=1}(\mathrm{~m})$ | $\Omega\left(\mathrm{rad} \mathrm{s}^{-1}\right)$ |
| :--- | ---: | ---: | ---: |
| 10 | 10 | 0.5 | 0 |


a rectangular or an elliptic wing is represented by an AL as shown in Fig. 5, where the coordinate system is unchanged from the definition in Fig. 4. The AL is discretized in uniformly spaced sections in-line with the underlying flow grid, and the smearing parameter is twice the section width, which ensures a continuous force distribution along the wing. The common simulation parameters are given in Table 1. Unless specifically stated the sectional lift coefficient $C_{\mathrm{L}}=1$ and drag is zero along the wing, independent of the angle of attack. The chord of the rectangular wing is set to 1 m and the elliptical chord distribution is
$c(x)=c_{0} \sqrt{1-\left(\frac{x-\left(b+r_{v=1}\right)}{b}\right)^{2}}$,
with the root chord $c_{0}=4 \mathrm{~m}$. All simulations are performed within the same computational domain, defined in Sect. 3.1.

The theoretical predictions of the velocity field are achieved by representing the vortex system of Fig. 5 with vortex filaments. The velocity induced by a filament with a viscous vortex core at an arbitrary point $C$ is

$$
\begin{align*}
\boldsymbol{u} & =f_{\epsilon}\left(x_{\perp}\right) \frac{\Gamma}{4 \pi} \frac{\left(x_{1}+x_{2}\right)\left(\boldsymbol{x}_{\mathbf{1}} \times \boldsymbol{x}_{\mathbf{2}}\right)}{x_{1} x_{2}+\boldsymbol{x}_{\mathbf{1}} \cdot \boldsymbol{x}_{\mathbf{2}}}  \tag{25}\\
x_{\perp} & =\frac{\left|\boldsymbol{x}_{\mathbf{1}} \times \boldsymbol{x}_{\mathbf{2}}\right|}{\left|x_{2}-x_{1}\right|} \quad \text { and } \quad x_{i}=\left|\boldsymbol{x}_{i}\right| \tag{26}
\end{align*}
$$

where $\boldsymbol{x}_{1}$ points from the start of the filament to $C$ and $\boldsymbol{x}_{2}$ from its end. For a definition of $f_{\epsilon}\left(x_{\perp}\right)$ refer to Eq. (17); without viscous core $f_{\epsilon}\left(x_{\perp}\right)=1$. The contribution from different segments is summed to give the overall velocity field. As the wing is lightly loaded, we assume that all vortex segments remain in the $x-z$ plane.

### 5.1.1 Trailed vorticity smearing

The vortex smearing hypothesis assumes the trailed vorticity inheriting the smeared velocity field from the bound vor-

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-07.jpg?height=523&width=838&top_left_y=223&top_left_x=1065)
Figure 6. Velocities induced perpendicular to a rectangular wing predicted by an actuator line (AL) without correction and by vortex segments with a viscous core (Vortex) with different smearing parameters. Velocities are shown along lines cutting the bound and trailed vortices at right angles and $y=0$. Only half of the horseshoe vortex is depicted as $x^{\prime}=x-\left(b+r_{v=1}\right)$.

tex, which was confirmed theoretically by Martínez-Tossas and Meneveau (2019) for straight wings. We test this further by simulating a rectangular wing without any correction. All vorticity is shed from the wing tips, creating the well-known horseshoe vortex. Hence, for the hypothesis to be valid, the velocity distribution in the plane orthogonal to the trailed vortices should be identical to that of the bound vortex .

Figure 6 compares the velocities induced by a rectangular wing predicted by an AL and three vortex segments (one bound, two trailed) for five different smearing parameters. Only half of the wing is presented, due to symmetry. Velocity distributions are shown for lines cutting the vortex segments at right angles for $y=0$. Clearly the velocity smearing is identical between trailed and bound vorticity, confirming the smearing hypothesis. Slight differences are linked to the numerical discretization of the Gaussian force projection (Shives and Crawford, 2013) and numerical diffusion.

### 5.1.2 Smearing correction verification

As mentioned in Sect. 4, the new smearing correction uses a lifting line representation of the trailed vorticity. Thus, the prediction of the velocity correction with our model or vortex segments should be identical. To simultaneously verify its numerical implementation, our model only receives the sampled velocities from the flow domain to compute the circulation at the sections. Furthermore, the body forces are not applied inside the domain to avoid influencing the trailed vortex paths, and the correction velocities are not added to the CFD velocities to keep the circulation unchanged. This holds the trailed vortices in the $x-z$ plane, simplifying the representation of the wake with vortex segments. The segments' circulation is exactly the same as in the smearing correction to avoid any numerical effects influencing the comparison. Fig-

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-08.jpg?height=403&width=847&top_left_y=217&top_left_x=157)
Figure 7. Analytical (An.) and corresponding model prediction of the velocity correction for varying force smearing ( $\epsilon=2 b /\left(N_{v}-1\right)$ ) along a rectangular and elliptic wing, where $x^{\prime}=x-\left(b+r_{v=1}\right)$.

ure 7 compares the velocity correction predicted by the analytical vortex segments and our model at each section along a rectangular and elliptic wing for different smearing factors (i.e. $\epsilon=2 b /\left(N_{v}-1\right)$ ). With decreasing force smearing the velocity correction concentrates towards the tips, as the induced velocity gradients increase. Therefore, even at higher resolution the smearing correction remains significant, although more localized. Generally the model slightly overpredicts the missing induction at the wing, becoming more prominent with increasing resolution. At $N_{v}=64$ the difference reaches a maximum of $6.7 \%$ (rectangular) and $1.7 \%$ (elliptic) with respect to the inflow velocity. The average error does not breach $0.5 \%$ in any case. The velocity jump towards the tip sections of the elliptical wing is related to the equidistant discretization of the wing (Pirrung et al., 2014).

### 5.1.3 Coupled AL-smearing correction verification

The coupling between velocity correction and the flow domain is verified by comparing the corrected downwash at an elliptical wing to the theoretical expectation. The downwash should be constant along the wing and is given by

$$
\begin{equation*}
v_{\mathrm{th}}=-\frac{\Gamma_{0}}{4 b}=-\frac{w_{\infty} c_{0} C_{\mathrm{L}}}{8 b}, \tag{27}
\end{equation*}
$$

where $\Gamma_{0}$ is the circulation at the wing root. Similar to Shives and Crawford (2013) the $C_{\mathrm{L}}$ was not fixed, but instead followed the theoretical lift curve slope for thin airfoils $C_{\mathrm{L}}=2 \pi$. For the wing to operate at a constant lift coefficient $C_{\mathrm{L}}=1$, its angle of attack needed to include the effect of the induced velocities:

$$
\begin{equation*}
\alpha=\alpha_{\mathrm{eff}}+\alpha_{i}=\frac{C_{\mathrm{L}}}{2 \pi}+\tan ^{-1}\left(\frac{c_{0} C_{\mathrm{L}}}{8 b}\right) . \tag{28}
\end{equation*}
$$

This represents a more rigorous test of the coupled system than prescribing the loading along the wing, as only the correct downwash leads to the desired, constant sectional lift coefficients.

Figure 8 shows the downwash predicted by AL simulations with different smearing parameters and active correction. The CFD components of the velocities are shown

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-08.jpg?height=525&width=608&top_left_y=221&top_left_x=1180)
Figure 8. Downwash at an elliptical wing predicted by AL simulations with different smearing factors and smearing correction. The CFD components of the velocities are shown (dashed) as well as the total downwash incorporating the correction (solid). The theoretical value acts as reference.

Table 2. Input parameters for the NREL 5-MW simulations.
| $V_{\infty}\left(\mathrm{m} \mathrm{s}^{-1}\right)$ | $\Omega(\mathrm{rpm})$ | Pitch $\left(^{\circ}\right)$ |
| :--- | ---: | ---: |
| 4 | 4.6 | 0.00 |
| 6 | 6.9 | 0.00 |
| 8 | 9.2 | 0.00 |
| 14 | 12.1 | 2.59 |
| 25 | 12.1 | 23.09 |


( $v^{\mathrm{CFD}}$ ) separately to emphasize the contribution of the correction to arrive at the correct, constant downwash of $1 \mathrm{~ms}^{-1}$. Clearly without the correction, the induced velocities are a function of the smearing factor and only arrive at the theoretically expected value for $N_{v}=32$. Including the correction greatly reduces the dependence of the downwash on the force smearing. The insufficient correction towards the tips feeds back to the equidistant discretization of the AL (Pirrung et al., 2014), which is linked to the uniform spacing of the underlying flow grid.

### 5.2 Rotor simulations - NREL 5-MW

The validity of the smearing hypothesis and its correction in rotor applications is demonstrated with simulations of the NREL 5-MW turbine (Jonkman et al., 2009) using actuator line (AL) and lifting line (LL) models. The input parameters for these simulations are given in Table 2.

### 5.2.1 Uniform inflow

Figure 9 compares the AL results with and without the novel smearing correction to the LL with and without viscous core. At this wind speed of $8 \mathrm{~m} \mathrm{~s}^{-1}$ the thrust coefficient is highest $\left(C_{\mathrm{T}}=0.84\right)-$ and hence induction is highest - thus lending itself as a strong verification case. Clearly, there is an equiva-

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-09.jpg?height=384&width=847&top_left_y=217&top_left_x=157)
Figure 9. Normal and tangential forces on the NREL 5-MW blades at $8 \mathrm{~ms}^{-1}$ predicted by AL simulations with/without smearing correction and LL with/without viscous core ( $N_{s}=19, \epsilon=0.1 R$ ).

lence between the original AL and the LL with a viscous core and the corrected AL with the LL. Therefore, the smearing correction makes the AL effectively act as a LL, as originally intended by Sørensen and Shen (2002). The impact of the viscous core is most prominent toward the blade root and tip. The sudden drop in the forces predicted by the AL/LL+core for the tip section of the blade - located at $r / R=0.97-$ is not triggered by any aerodynamic tip effects, but relates to a pronounced reduction in chord. Figure 10 shows the diminishing effect of the correction on the angle of attack towards the root and tip. As the correction velocities are negligible with respect to the rotational velocity - they impact the velocity magnitude by less than $0.1 \%$ in the lifting region of the blade - it is ultimately the change in the angle of attack that explains the observations in the force distributions. While not greatly affecting the magnitude of the forces in the mid-section of the blade, the viscous core does introduce greater fluctuations in the force distribution. Hence, the missing induction introduced by the viscous core reduces the coupling between neighbouring blade sections. The smearing correction also recovers this behaviour of the LL. Surpassing rated wind speed, forces increase inboard until cutout. Thus, just before cut-out at $25 \mathrm{~ms}^{-1}$ loading reaches a maximum towards the root, causing an equally pronounced influence of the smearing correction in this region, as demonstrated in Fig. 11. Again, the equivalence of the AL and LL implementations is remarkable. This high wind speed case also demonstrates our correction is not only a tip correction.

The comparison of AL and LL is summarized in Fig. 12 in the form of local thrust and power distributions at different wind speeds. Note that for the wind speeds below rated $\left(<11.4 \mathrm{~m} \mathrm{~s}^{-1}\right)$, the coefficients are identical. The results are only presented for simulations with $N_{s}=9$ for visibility, but compare equally well at higher resolution. As mentioned earlier, the smearing correction predominantly acts towards the tip and root. An additional overview of all results is given in Table 3. Here the total rotor thrust $T$ and power $P$ predicted by the corrected actuator line (AL*) and the lifting line (LL) are listed as well as the influence of adding the viscous core relative to AL* and LL, respectively. The AL and LL solu-

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-09.jpg?height=500&width=608&top_left_y=223&top_left_x=1178)
Figure 10. Angle of attack with/without smearing correction on the NREL 5-MW blades at $8 \mathrm{~m} \mathrm{~s}^{-1}\left(N_{s}=19, \epsilon=0.1 R\right)$.

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-09.jpg?height=380&width=842&top_left_y=875&top_left_x=1063)
Figure 11. Normal and tangential forces on the NREL 5-MW blades at $25 \mathrm{~ms}^{-1}$ predicted by AL simulations with/without smearing correction and LL with/without viscous core ( $N_{s}=9, \epsilon= 0.2 R$ ).

tions are not directly compared to avoid including any mean bias in the comparison. The impact of the correction on the AL forces is nearly identical to removing the viscous core in the LL simulations at any wind speed, which further supports our correction methodology. In light of the large errors incurred without any correction, unsurprisingly, some form of tip correction is usually applied in AL simulations.

### 5.2.2 Yawed inflow

As the smearing correction does not include yaw effects the wake is assumed to advect normal to the rotor plane we tested its influence at yaw angles $\chi$ of 15,30 and $45^{\circ}$ at $8 \mathrm{~m} \mathrm{~s}^{-1}$. Again the LL with and without viscous core acted as a reference. The time steps remained the same as in uniform inflow. Here only the results for the most extreme case at $45^{\circ}$ yaw are shown, as the differences are most severe in this case. Figure 13 presents the normal and tangential force variation during one rotation, averaged over three distinct regions of the blade, at a wind speed of $8 \mathrm{~ms}^{-1}$. Whilst the agreement is best towards the blade tip, the force variation with azimuthal position is similar between AL and LL simulations across all sections. AL results are shifted downwards with respect to the LL predictions at the inner sections, hint-

Table 3. An overview of the simulation inputs and results for the NREL 5-MW in uniform inflow. Results are grouped by blade/grid resolution. For the actuator line (AL) and lifting line (LL) the simulation time step $\Delta t$, the total thrust $T$ and power $P$ as well as the relative change in these quantities caused by the correction/removing the viscous core are listed. Note that AL* represents the corrected AL results, and the change is expressed relative to the rotor thrust and power.
| $N_{s}$ | $\epsilon / R$ | $V_{\infty}\left[\mathrm{m} \mathrm{s}^{-1}\right]$ | $\Delta t \times 10^{-2}[\mathrm{~s}]$ |  | $T \times 10^{5}[\mathrm{~N}]$ |  | $\Delta T[\%]$ |  | $P \times 10^{6}[\mathrm{~W}]$ |  | $\Delta P$ [\%] |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  |  | AL | LL | AL* | LL | $\Delta \mathrm{AL}$ | $\Delta \mathrm{LL}$ | AL* | LL | $\Delta$ AL | $\Delta \mathrm{LL}$ |
| 9 | 0.2 | 4 | 15.90 | 18.12 | 1.01 | 1.02 | 3.50 | 3.21 | 0.26 | 0.26 | 9.19 | 8.74 |
|  |  | 6 | 13.82 | 12.08 | 2.27 | 2.30 | 3.78 | 3.20 | 0.89 | 0.88 | 8.91 | 8.72 |
|  |  | 8 | 10.36 | 9.06 | 4.08 | 4.09 | 2.84 | 3.20 | 2.13 | 2.08 | 7.51 | 8.71 |
|  |  | 14 | 7.87 | 6.89 | 4.64 | 4.77 | 4.00 | 3.19 | 5.39 | 5.54 | 5.64 | 4.90 |
|  |  | 25 | 7.87 | 6.89 | 2.84 | 2.99 | 2.77 | 0.82 | 5.40 | 5.68 | 3.08 | 1.65 |
| 19 | 0.1 | 4 | 10.37 | 18.12 | 0.98 | 0.99 | 2.80 | 2.13 | 0.25 | 0.25 | 7.48 | 5.81 |
|  |  | 6 | 6.90 | 12.08 | 2.21 | 2.22 | 2.79 | 2.12 | 0.83 | 0.85 | 7.45 | 5.79 |
|  |  | 8 | 5.17 | 9.06 | 3.92 | 3.95 | 2.81 | 2.13 | 1.98 | 2.02 | 7.49 | 5.81 |
|  |  | 14 | 3.93 | 6.89 | 4.52 | 4.61 | 3.09 | 5.98 | 5.22 | 5.35 | 4.71 | 7.50 |


![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-10.jpg?height=474&width=839&top_left_y=994&top_left_x=161)
Figure 12. Local thrust and power coefficients along the NREL 5MW blades at different wind speeds predicted by AL simulations with/without smearing correction and LL with/without core ( $N_{s}= 9, \epsilon=0.2 R$ ).

ing at the AL experiencing higher induction in this region. However, for the verification of the smearing correction this shift is irrelevant, instead its impact on the AL forces needs to be assessed relative to the difference between LL with and without core. In this respect the smearing correction behaves correctly, increasing forces in a similar fashion as a LL without core in the mid-section of the blade and reducing them towards the root and tip.

### 5.2.3 Pitch step

The pitch step is defined as

$$
\begin{equation*}
\psi=\psi_{0}+\frac{\Delta \psi}{2}\left[1+\tanh \left(k\left(t-t_{0}\right)\right)\right], \tag{29}
\end{equation*}
$$

with $\psi_{0}$ defining the pitch angle before the step, $t_{0}$ representing the time instant of the step and $\Delta \psi$ denoting the pitch change. Here an extremely violent step is chosen - determined by $k$ - to encourage an equally pronounced blade

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-10.jpg?height=727&width=836&top_left_y=997&top_left_x=1065)
Figure 13. Normal ( $\mathbf{a , b , c}$ ) and tangential (d, e,f) force variation during one rotation as a function of azimuthal position of the NREL 5-MW blades at $8 \mathrm{~ms}^{-1}$ and $45^{\circ}$ yaw. The forces are averaged over three sections (inner: $\mathbf{a , d}$; middle: $\mathbf{b , ~} \mathbf{e}$; outer: $\mathbf{c , f}$ ) of the blade and are predicted by AL simulations with/without smearing correction and LL with/without viscous core. The blade is facing upwind at $\chi=0^{\circ}$ and is pointing vertically up at $\chi=90^{\circ}\left(N_{s}=9, \epsilon=0.2 R\right)$.

force response and test the numerical stability of our correction. The parameters governing this comparison are given in Table 4, which realize a pitch step of $\pm 2^{\circ}$ in 0.14 s ( $10 \%$ to $90 \%$ pitch). To capture the swift change in pitch, the time step is adjusted in both AL and LL simulations to $3.94 \times 10^{-2} \mathrm{~s}$. The blade force response is normalized as

$$
\begin{equation*}
\hat{F}(t)=\frac{F(t)-F_{0}}{F_{\infty}-F_{0}}, \tag{30}
\end{equation*}
$$

Table 4. Inputs defining the pitch step.
| $V_{\infty}\left[\mathrm{m} \mathrm{s}^{-1}\right]$ | $\Omega[\mathrm{rpm}]$ | $\psi_{0}\left[{ }^{\circ}\right]$ | $\Delta \psi\left[{ }^{\circ}\right]$ | $k$ |
| :--- | ---: | ---: | ---: | ---: |
| 14 | 12.1 | 2.59 | $\pm 2$ | 16 |


![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-11.jpg?height=361&width=834&top_left_y=477&top_left_x=166)
Figure 14. Normalized tangential force response at two blade sections (middle and tip) of the NREL 5-MW following a pitch step of $+2^{\circ}$ in 0.14 s at $14 \mathrm{~m} \mathrm{~s}^{-1}$ predicted by AL simulations with/without smearing correction and LL with/without viscous core ( $N_{s}=9$, $\epsilon=0.2 R$ ).

with $F_{0}$ and $F_{\infty}$ denoting the steady-state values before and after the pitch step, respectively.

Here only the tangential force response after a $+2^{\circ}$ step for the mid and tip blade sections are shown in Fig. 14, as they capture the main features of the response. As the definition here is positive pitch to feather, the force decreases along the blade for positive pitch changes. The AL simulations exhibit a faster response with a kink at 0.14 s , coinciding with the pitch change reaching $99 \%$ of the step. Therefore, the AL seems to capture the pitch rate lift. The LL does not show this feature so - as in yaw - the correct behaviour of the smearing correction on the AL force response should be assessed relative to the influence of removing the viscous core in the LL model. Overall, the correction has limited effect on the dynamic response, which is also confirmed by the LL simulations, but the correction essentially acts on the forces in the same fashion as removing the core in the LL. In the midsection it reduces the forces by a maximum of $1 \%$ during the first 2 s , dropping to $0.5 \%$ afterwards. At the tip section it intensifies the response by a maximum of $1 \%$, diminishing to $0.3 \%$ at 4 s .

### 5.2.4 Turbulent inflow

Highly turbulent inflow should challenge the numerical stability of the new smearing correction by introducing strong and abrupt changes in the angle of attack. Comparing simulations with and without inflow turbulence should also reveal whether turbulence alters the nature of the correction. Figure 15 shows the impact of the smearing correction on the time-averaged normal and tangential blade forces at an $8 \mathrm{~m} \mathrm{~s}^{-1}$ mean wind speed for AL simulations with uniform and turbulent inflow. At a turbulence intensity (TI) of $15 \%$,

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-11.jpg?height=378&width=838&top_left_y=221&top_left_x=1065)
Figure 15. Time-averaged normal and tangential forces on the NREL 5-MW blades at an $8 \mathrm{~ms}^{-1}$ mean wind speed and changing inflow turbulence predicted by the AL model with and without smearing correction ( $N_{s}=35, \epsilon=R / 8$ ).

![](https://cdn.mathpix.com/cropped/0083d7ae-da3b-4a17-a3e2-0169bfae444f-11.jpg?height=361&width=836&top_left_y=828&top_left_x=1067)
Figure 16. Variation of normal and tangential forces on the NREL 5-MW blades at an $8 \mathrm{~ms}^{-1}$ mean wind speed and turbulence intensity of $15 \%$ predicted by the AL model with and without smearing correction ( $N_{s}=35, \epsilon=R / 8$ ).

the forces are unsurprisingly slightly larger ( $\approx 20 \mathrm{Nm}^{-1}$ ) than in uniform inflow. However, the change in forces introduced by the correction is nearly identical ( $<2 \mathrm{Nm}^{-1}$ ). When comparing the standard deviation of the forces with and without correction in Fig. 16 the smoothing and dampening effect of the smearing correction on the forces in highly turbulent inflow also becomes apparent. Madsen et al. (2018) observed a corresponding reduction of the force variations on the whole rotor blade when comparing near-wake model results against BEM results for the NM 80 rotor in turbulent inflow. This illustrates that the smearing correction leads to the same dynamic coupling between neighbouring blade sections as a lifting line model.

## 6 Conclusions

The actuator line was intended as a lifting line technique for CFD applications. In this paper we prove - theoretically and practically - that smearing the forces of the actuator line in the flow domain leads to smeared velocity fields. For the typical Gaussian force projection, the widely known LambOseen (Lamb, 1932; Oseen, 1911) viscous core appears in both bound and trailed vorticity. This core reduces the velocities approaching the vortex centre compared with the inviscid solution of the lifting line. Thus, the trailed vorticity of
an actuator line induces lower velocities at the blade owing to the force projection. We recover this missing induction by combining a near-wake model of the trailed vorticity with the Lamb-Oseen viscous core model and coupling it with the actuator line model. Basic wing test cases with theoretical solutions verify the correction, as it recovers nearly all induction independent of the severity of the force smearing. Furthermore, rotor simulations show the applicability and strength of the correction over the entire operational wind speed range as well as in yaw, strong turbulence and undergoing pitch steps. Here the correction is validated with lifting line simulations with and without viscous core, which are representative of an actuator line with and without smearing correction, respectively. The agreement between the respective actuator line and lifting line results is remarkable.

The current implementation of the smearing correction relies on heavy bookkeeping. In future versions the latter will be removed without jeopardizing stability or accuracy, making it suitable for wind farm simulations in realistic atmospheric flows. Potentially, the correction might also enable accurate rotor simulations at lower discretization.

Code and data availability. All data and parts of the code covering the smearing correction are available upon request. Commercial and research licenses for EllipSys3D can be purchased from DTU.

## Appendix A: Fixed-wake equations

The equations governing the fixed-wake approach underlying the smearing correction (see Fig. 4) are summarized for completeness.

Direction vector from the vortex element to a control point along the blade:
$\boldsymbol{x}=\left(\begin{array}{c}-r \cos \beta+r-h \\ r \sin \beta \\ -v_{h} \beta / \Omega\end{array}\right)$.

Definition of the vortex element:

$$
\boldsymbol{\delta} \boldsymbol{l}=\delta l \cos \phi\left(\begin{array}{c}-\sin \beta  \tag{A2}\\ -\cos \beta \\ v_{h} /(\Omega r)\end{array}\right)
$$

with $\delta l=\frac{r \delta \beta}{\cos \phi}$. Note that

$$
\begin{equation*}
\cos \phi=\frac{\Omega r}{\sqrt{(\Omega r)^{2}+v_{h}^{2}}}=\frac{1}{\sqrt{1+\left(\frac{v_{h}}{\Omega r}\right)^{2}}} \tag{A3}
\end{equation*}
$$

with $\tan \phi=\frac{v_{h}}{\Omega r}$.

Incremental velocity induced by the vortex element

$$
\begin{align*}
& \delta \widetilde{\boldsymbol{u}}=\frac{\Delta \Gamma}{4 \pi} \frac{\delta \boldsymbol{l} \times \boldsymbol{x}}{|\boldsymbol{x}|^{3}}= \\
& \quad \operatorname{Ar}\left(\begin{array}{c}
\left(v_{h} / \Omega r\right)(\beta \cos \beta-\sin \beta) \\
-\left(v_{h} / \Omega r\right)\left(-1+h_{r}+\cos \beta+\beta \sin \beta\right) \\
{\left[-1+\left(1-h_{r}\right) \cos \beta\right]}
\end{array}\right) \tag{A4}
\end{align*}
$$

with
$A=\frac{\Delta \Gamma r \delta \beta}{4 \pi}$

$$
\begin{aligned}
& \underbrace{\left[r^{2}\left(1+\left(1-h_{r}\right)^{2}-2\left(1-h_{r}\right) \cos \beta+\left(\frac{v_{h} \beta}{\Omega r}\right)^{2}\right)\right]^{-\frac{3}{2}}}_{|\boldsymbol{x}|^{2}} \\
& h_{r}=\frac{h}{r}
\end{aligned}
$$

Author contributions. ARMF was responsible for writing the paper, except for Sect. 3.2 which was authored by NRG. The original numerical implementation of the near-wake model by GRP was adapted by ARMF to incorporate the smearing correction with the help of GRP. NRG performed all simulations with the lifting line, and ARMF performed all simulations with the actuator line. All authors decided on the flow cases to simulate and commented on the paper.

Competing interests. The authors declare that they have no conflict of interest.

Acknowledgements. We would like to acknowledge DTU Wind Energy's internal project "Virtual Atmosphere" for partially funding this research. Furthermore many thanks to senior researcher Mac Gaunaa for his insights on vortex aerodynamics and senior researcher Niels Troldborg for his input on actuator line simulations/modelling (both from DTU Wind Energy). Thanks also to Ang Li for his help regarding the near-wake model.

Review statement. This paper was edited by Alessandro Bianchini and reviewed by two anonymous referees.

## References

Ananthan, S. and Leishman, J. G.: Role of Filament Strain in the Free-Vortex Modeling of Rotor Wakes, J. Am. Helicopter Soc., 9, 176-191, 2004.
Dag, K., Sørensen, J., Sørensen, N., and Shen, W.: Combined pseudo-spectral/actuator line model for wind turbine applications, PhD thesis, DTU Wind Energy, Denmark, 2017.
Forsythe, J. R., Lynch, E., Polsky, S., and Spalart, P.: Coupled Flight Simulator and CFD Calculations of Ship Airwake using Kestrel, in: 53rd AIAA Aerospace Sciences Meeting, AIAA, Kissimmee, Florida, USA, https://doi.org/10.2514/6.2015-0556, 2015.
Glauert, H.: Airplane Propellers, Springer Berlin Heidelberg, Berlin, Heidelberg, 169-360, https://doi.org/10.1007/978-3-642-91487-4_3, 1935.
Jha, P. K., Churchfield, M. J., Moriarty, P. J., and Schmitz, S.: Guidelines for Volume Force Distributions Within Actuator Line Modeling of Wind Turbines on Large-Eddy Simulation-Type Grids, J. Sol. Energ.-T. ASME, 136, 031003, https://doi.org/10.1115/1.4026252, 2014.
Jonkman, J., Butterfield, S., Musial, W., and Scott, G.: Definition of a 5-MW reference wind turbine for offshore system development, Tech. rep., NREL/TP-500-38060, National Renewable Energy Laboratory (NREL), Colorado, USA, 2009.
Lamb, H.: Hydrodynamics, C.U.P, 6th Edn., Cambridge University Press, Cambridge, 1932.
Leishman, J. G., Bhagwat, M. J., and Bagai, A.: Free-Vortex Filament Methods for the Analysis of Helicopter Rotor Wakes, J. Aircraft, 39, 759-775, 2002.
Leonard, B.: A stable and accurate convective modelling procedure based on quadratic upstream interpolation, Comput. Method. Appl. M., 19, 59-98, 1979.

Madsen, H. A., Sørensen, N. N., Bak, C., Troldborg, N., and Pirrung, G.: Measured aerodynamic forces on a full scale 2 MW turbine in comparison with EllipSys3D and HAWC2 simulations, J. Phys. Conf. Ser., 1037, 022011, https://doi.org/10.1088/17426596/1037/2/022011, 2018.
Martínez-Tossas, L. A. and Meneveau, C.: Filtered lifting line theory and application to the actuator line model, J. Fluid Mech., 863, 269-292, https://doi.org/10.1017/jfm.2018.994, 2019.
Menter, F. R.: Zonal two equation $k-\omega$ turbulence models for aerodynamic flows, in: 23rd Fluid Dynamics, Plasmadynamics, and Lasers Conference, Fluid Dynamics and Co-located Conferences, Orlando,FL, https://doi.org/10.2514/6.1993-2906, 1993.
Meyer Forsting, A., Troldborg, N., Bechmann, A., and Réthoré, P.-E.: Modelling Wind Turbine Inflow: The Induction Zone, PhD thesis, DTU Wind Energy, Denmark, https://doi.org/10.11581/DTU:00000022, 2017.
Michelsen, J.: Basis3D - a platform for development of multiblock PDE solvers, Tech. rep., Dept. of Fluid Mechanics, Technical University of Denmark, DTU, 1994a.
Michelsen, J.: Block structured multigrid solution of 2D and 3D elliptic PDE's, Tech. rep., Dept. of Fluid Mechanics, Technical University of Denmark, DTU, 1994b.
Mikkelsen, R.: Actuator Disc Methods Applied to Wind Turbines, PhD thesis, Technical University of Denmark, 2003.
Nathan, J.: Application of Actuator Surface Concept in LES Simulations of the Near Wake of Wind Turbines, PhD thesis, École de Technologie Supérieure, Montreal, Canada, 2018.
Oseen, C.: Über Wirbelbewegung in einer reibenden Flüssigkeit, Arkiv för matematik, astronomi och fysik, Ark. Mat. Astron. Fys., 7, 14-21, 1911.
Patanker, S. and Spalding, D.: A calculation procedure for heat, mass and momentum transfer in three-dimensional parabolic flows, Int. J. Heat Mass Tran., 15, 59-98, 1972.
Pirrung, G. R., Hansen, M. H., and Madsen, H. A.: Improvement of a near wake model for trailing vorticity, J. Phys. Conf. Ser., 555, 012083, https://doi.org/10.1088/1742-6596/555/1/012083, 2014.
Pirrung, G., Madsen, H. A., Kim, T., and Heinz, J.: A coupled near and far wake model for wind turbine aerodynamics, Wind Energy, 19, 2053-2069, https://doi.org/10.1002/we.1969, 2016.
Pirrung, G., Riziotis, V., Madsen, H., Hansen, M., and Kim, T.: Comparison of a coupled near- and far-wake model with a free-wake vortex code, Wind Energ. Sci., 2, 15-33, https://doi.org/10.5194/wes-2-15-2017, 2017a.
Pirrung, G. R., Madsen, H. A., and Schreck, S.: Trailed vorticity modeling for aeroelastic wind turbine simulations in standstill, Wind Energ. Sci., 2, 521-532, https://doi.org/10.5194/wes-2-521-2017, 2017b.
Ramos-García, N., Shen, W. Z., and Sørensen, J. N.: Validation of a three-dimensional viscous-inviscid interactive solver for wind turbine rotors, Renew. Energ., 70, 78-92, https://doi.org/10.1016/j.renene.2014.04.001, 2014a.
Ramos-García, N., Shen, W. Z., and Sørensen, J. N.: Three-dimensional viscous-inviscid coupling method for wind turbine computations, Wind Energy, 19, 67-93, https://doi.org/10.1002/we.1821, 2014b.
Ramos-García, N., M. Mølholm, J. N. S., and Walther, J. H.: Hybrid vortex simulations of wind turbines using a three-dimensional viscous-inviscid panel method, Wind Energy, 20, 1187-1889, 2017.

Réthoré, P.-E. and Sørensen, N.: A discrete force allocation algorithm for modelling wind turbines in computational fluid dynamics, Wind Energy, 15, 915-926, https://doi.org/10.1002/we.525, 2012.

Shen, W., Sørensen, J., and Mikkelsen, R.: Tip loss correction for actuator/Navier Stokes computations, J. Sol. Energ.-T. ASME, 59, 209-213, 2005.
Shives, M. and Crawford, C.: Mesh and load distribution requirements for actuator line CFD simulations, Wind Energy, 16, 657669, https://doi.org/10.1002/we.1546, 2013.
Sørensen, J. N. and Shen, W. Z.: Numerical modelling of wind turbine wakes, J. Fluid. Eng.-T. ASME, 124, 393-399, https://doi.org/10.1115/1.1471361, 2002.

Sørensen, N.: General purpose flow solver applied to flow over hills, PhD thesis, Risø National Laboratory, 1995.
Squire, H. B.: The growth of a vortex in turbulent flow, Aeronaut. Quart., 16, 302-306, 1965.
Strelets, M.: Detached eddy simulation of massively separated flows, in: 39th AIAA Aerospace Sciences Meeting and Exhibit, AIAA Paper 2001-0879, Reno, NV, 2001.
Troldborg, N., Sørensen, J., and Mikkelsen, R.: Actuator Line Modeling of Wind Turbine Wakes, PhD thesis, Technical University of Denmark, 2009.
Vollmer, L., Steinfeld, G., and Kühn, M.: Transient LES of an offshore wind turbine, Wind Energ. Sci., 2, 603-614, https://doi.org/10.5194/wes-2-603-2017, 2017.


[^0]:    $$
    \begin{equation*}
    { }^{1} \nabla \times\left(\hat{\boldsymbol{e}}_{x} \times \boldsymbol{u}\right)=(\boldsymbol{u} \nabla) \hat{\boldsymbol{e}}_{x}-\left(\hat{\boldsymbol{e}}_{x} \nabla\right) \boldsymbol{u}+\hat{\boldsymbol{e}}_{x}(\nabla \boldsymbol{u})-\boldsymbol{u}\left(\nabla \hat{\boldsymbol{e}}_{x}\right)= \tag{11}
    \end{equation*}
    $$

    $(\boldsymbol{u} \nabla) \hat{\boldsymbol{e}}_{x}+0+0+0$

[^1]:    ${ }^{2}$ Note that in the $x-y$ plane the circulation would be $-\Gamma$.

[^2]:    ${ }^{3}$ Strictly, the influence of the shed vorticity on the velocity at the AL should be removed as remarked by Martínez-Tossas and Meneveau (2019), however its influence is negligible.

