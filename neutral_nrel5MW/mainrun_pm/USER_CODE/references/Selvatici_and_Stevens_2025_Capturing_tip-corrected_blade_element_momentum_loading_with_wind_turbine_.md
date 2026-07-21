# Capturing tip-corrected blade element momentum loading with wind turbine models 

Davide Selvatici *, Richard J.A.M. Stevens<br>Physics of Fluids Group, Max Planck Center Twente for Complex Fluid Dynamics, J.M. Burgers Center for Fluid Dynamics, University of Twente, P.O. Box 217, Enschede, 7500 AE, The Netherlands

## ARTICLE INFO

## Keywords:

Actuator line method
Actuator disk model
Blade element momentum
Large eddy simulation


#### Abstract

Wind turbine models integrated into flow solvers often overestimate blade loading in the tip region. This overestimation arises from overestimating the sampled velocity and angle of attack near the blade tip. We use the ratio of tip-corrected to non-tip-corrected axial and tangential velocity components from Blade Element Momentum (BEM) theory to consistently correct the sampled velocity along the blade. This correction accurately captures both the angle of attack and tip loading and it does not require any inter-processor communication, a desirable feature in modern computing systems like GPUs. The method has been validated against a vortex-based smearing correction consistent with lifting line theory and BEM with tip correction, using both NREL-5MW and DTU-10MW turbines and adopting different inflow conditions and grid discretizations. The correction is demonstrated for the Actuator Line Method. It requires only basic turbine data, such as tip speed ratio and airfoil characteristics, making it adaptable to other turbine models, such as the actuator disk with rotation.


## 1. Introduction

Understanding the interaction between wind turbines and the surrounding flow is critical for their design, operation, and optimization. Computational fluid dynamics is an effective tool for studying such complex aerodynamics under controlled conditions. However, bladeresolved simulations remain challenging due to the need for precise blade geometry details - often proprietary and not publicly available - and the immense computational requirements. Therefore, wind turbine models representing blades through distributed forces offer a more efficient and practical approach.

Three commonly used wind turbine models are the actuator disk, actuator disk with rotation, and actuator line (ALM) models. Each model differs in its specific details, offering a unique balance between computational efficiency and simplicity versus accuracy [1-3]. These wind turbine models are frequently employed in Reynolds-Averaged Navier-Stokes (RANS) and large-eddy simulation (LES), with the choice of method depending on the specific objectives of the study. The actuator disk model is the simplest and most computationally efficient, requiring only the specification of the turbine thrust coefficient, with [4] proposing resolution corrections for this model. In contrast, the ALM and actuator disk with rotation use tabulated lift and drag coefficients to calculate instantaneous turbine blade loading, enabling a more precise representation of wake dynamics. The ALM effectively
captures unsteady and turbulent structures within wind turbine wakes, while the actuator disk model with rotation, although being more computationally efficient and suitable for coarser meshes, does not capture the blade motion [1,2].

Accurate force calculations require accurate determination of the local velocity near the blade tips. Given the complexity that this problem entails, the most effective method to accurately capture the blade tip loading in ALM is a subject of ongoing debate. The ALM is thought to inherently generate tip vortices, which reduce loads near the blade tips [5,6]. High-resolution simulations have indeed shown that ALM can accurately capture tip loadings when blade geometry is taken into account [7-9]. However, achieving this level of accuracy requires extremely fine grid resolution - up to 900 nodes per diameter at the blade tip - rendering it impractical for many applications. On coarser, more practical meshes, ALM tends to overestimate blade loading near the tips, underscoring the need for correction in many practical scenarios [10,11].

To address the overestimation of blade loading near the tip at practical grid resolutions, various models have been developed to enhance accuracy in the tip region [10,12]. Both ALM and the actuator disk with rotation are based on the fundamental principles of blade element momentum (BEM) theory, which utilizes similar aerodynamic force calculations. Therefore, Shen et al. [12] proposed a practical solution:

[^0]directly incorporating BEM tip corrections into ALM. While this approach improves ALM results on coarser meshes, it does not correct for mismatches in angle of attack and induction factors, which are essential for consistency with BEM theory. Other approaches involve correcting the velocity sampled from the flow field. While the velocity is typically sampled at the actuator point, alternative methods include using a weighted average in a volume near the actuator point [9,13], employing Lagrangian-averaged velocity sampling [14], and other sampling methods [15,16].

Dağ [11] and Dağ and Sørensen [17] observed that the bound vortex created by ALM exhibits a Gaussian vorticity distribution similar to a Lamb-Oseen vortex model. This observation led to the development of a vortex-based smearing correction method that approximates the velocity in ALM to that induced by the singular vorticity distribution predicted by Prandtl lifting line theory [18,19]. This approach was theoretically confirmed by Martínez-Tossas and Meneveau [20]. To recover the missing induction, the Biot-Savart law [18] would ideally be used. However, it is very computationally expensive, making wake models a more efficient alternative for correcting vortex smearing [19]. This approach has shown excellent agreement with results from lifting line theory. The accuracy of vortex smearing corrections has been further validated by an iterative method based on filtered lifting line theory [20] and a non-iterative method based on linearized lifting line theory [21]. Although being recognized for its accuracy, vortex-based smearing corrections could come with considerable complexity and computational demands, thus possibly limiting their wide adoption, especially when computational efficiency is relevant.

In the present work, we present a novel way to use BEM to selfconsistently correct the sampled velocity at the blade location, relying solely on locally available information. This approach does not require communication between different processors, enhancing computational efficiency especially in modern systems like GPUs. It offers a simple method for correcting velocity sampling in ALM, effectively addressing tip-loading distribution by recovering the missing induction and accounting for angle of attack differences using only the turbine's tip speed ratio and blade airfoil data as inputs. Additionally, this approach can be easily adapted to other models, such as the actuator disk with rotation, that utilize similar information. Our results demonstrate that this method accurately estimates thrust and power for the NREL-5MW and DTU-10MW turbines across a wide range of operating conditions, underscoring its potential in applications where computational efficiency and simplicity are crucial [16].

The remainder of this manuscript is structured as follows: Section 2 discusses tip corrections within BEM theory, and Section 3 reviews the classic ALM implementation. Section 4 introduces the proposed correction method. Section 5 details the LES used to validate the proposed method against BEM with tip correction and the vortex-based smearing correction by [18,19] for the NREL-5MW (Section 6) and DTU-10MW (Section 7) wind turbine rotors. Section 8 presents some concluding remarks regarding the computational efficiency of the method. Finally, Section 9 summarizes the findings of this study and highlights the broader impact of the proposed method on wind turbine simulation practices.

## 2. Tip-corrections in Blade Element Momentum theory

Based on the equivalence between the blade element and axial momentum theory, Blade Element Momentum theory is an analytic approach that has been developed to evaluate the induction velocity and loads at the rotor [6,22]. Experimental investigations and bladeresolved simulations have shown that neglecting tip corrections leads to inaccuracies in simulating the pressure equalization observed at the blade tip [23-26], hence including tip corrections has become standard practice in BEM to address these discrepancies. The main feature that BEM with tip correction introduces is the ability to model the deflection of flow streamlines in the rotor wake, enabling accurate reconstruction

![](https://cdn.mathpix.com/cropped/ec93945f-e863-4bc4-af5e-0180bd348132-2.jpg?height=431&width=683&top_left_y=185&top_left_x=1158)
Fig. 1. Sketch of the velocity components, angles, and forces as defined in BEM theory [22]. Flow direction is in positive x-direction.

of rotor induction. This is an improvement over standard BEM, which assumes a rigid wake and tends to overestimate induction [6,22]. This effect was modeled by Prandtl and Glauert using potential flow theory, who analytically derived the tip-loss function:

$$
\begin{equation*}
F(r)=\frac{2}{\pi} \arccos \left(\exp \left(-\frac{B(R-r)}{2 r \sin \phi}\right)\right), \tag{1}
\end{equation*}
$$

where $B$ is number of blades, $R$ is the blades radius, $r$ is the radial coordinate, and $\phi$ is the wind angle, defined as

$$
\begin{equation*}
\phi=\arctan \left(\frac{U}{V}\right), \tag{2}
\end{equation*}
$$

with $U$ and $V$ the axial and tangential components of velocity at the blade (see Fig. 1):

$$
\begin{equation*}
U=U_{\infty}(1-a) \quad \text { and } \quad V=\Omega r\left(1+a^{\prime}\right) \tag{3}
\end{equation*}
$$

where $a$ and $a^{\prime}$ are the axial and tangential induction factors, and $\Omega$ is the angular velocity of the rotor. Other tip correction functions exist to take into account the blades prebend or sweep angles [see for example 27]. For a broader view on tip correction functions and BEM theory we refer the reader to the books by [6,22].

When the tip-loss function is adopted, the induction factors become dependent on the blade location and the resulting axial velocity is reduced [6,22,28]:

$$
\begin{equation*}
a(r)=\frac{1}{\frac{4 F(r) \sin ^{2} \phi}{\sigma(r) C_{n}}+1} \text { and } a^{\prime}(r)=\frac{1}{\frac{4 F(r) \sin \phi \cos \phi}{\sigma(r) C_{t}}-1}, \tag{4}
\end{equation*}
$$

in which $\sigma(r)=B c(r) /(2 \pi r)$ is the radially dependent solidity factor, $c(r)$ is the chord length, and $C_{n}(r)$ and $C_{t}(r)$ are the normal and tangential load coefficients, which are defined as

$$
\begin{equation*}
C_{n}(r)=C_{L}(r) \cos \phi+C_{D}(r) \sin \phi \quad \text { and } \quad C_{t}(r)=C_{L}(r) \sin \phi-C_{D}(r) \cos \phi \tag{5}
\end{equation*}
$$

where $C_{L}(r)$ and $C_{D}(r)$ are the radially dependent lift and drag coefficients. The increase of induction is reflected in the tendency of $F(r)$ to 0 towards the tip. Given that the normal component of the velocity ( $C_{n}$ ) at the blade tip significantly exceeds the tangential component ( $C_{t}$ ), the relative increase in tangential velocity is smaller than the decrease in axial velocity. This causes the load reduction towards the blade tip.

The aerodynamic coefficients $C_{L}$ and $C_{D}$ are obtained from look-up tables that document the coefficients as a function of the angle of attack $\alpha$, which is determined as (see Fig. 1)

$$
\begin{equation*}
\alpha(r)=\phi(r)-(\gamma(r)+\theta(r)), \tag{6}
\end{equation*}
$$

where $\gamma(r)$ and $\theta(r)$ are the radially dependent twist and pitch angles of the blade.

The normal and tangential forces at the blade are obtained from the relative incoming velocity $U_{\text {rel }}=\sqrt{U^{2}+V^{2}}$ as:

$$
\begin{equation*}
F_{n}(r)=\frac{1}{2} \rho c U_{\mathrm{rel}}^{2} C_{n}(r) \quad \text { and } \quad F_{t}(r)=\frac{1}{2} \rho c U_{\mathrm{rel}}^{2} C_{t}(r) . \tag{7}
\end{equation*}
$$

![](https://cdn.mathpix.com/cropped/ec93945f-e863-4bc4-af5e-0180bd348132-3.jpg?height=479&width=1194&top_left_y=183&top_left_x=434)
Fig. 2. Comparison of axial (a) and tangential (b) velocity components for BEM with and without tip correction adopting the NREL-5MW rotor.

The rotor thrust and power coefficients follow from integration over the blades:

$$
\begin{equation*}
C_{T}=\frac{\int_{r_{\mathrm{hub}}}^{R} F_{n}(r) B d r}{\frac{1}{2} \rho U_{\infty}^{2} \pi R^{2}} \quad \text { and } \quad C_{P}=\frac{\int_{r_{\mathrm{hub}}}^{R} F_{t}(r) r B \Omega d r}{\frac{1}{2} \rho U_{\infty}^{3} \pi R^{2}}, \tag{8}
\end{equation*}
$$

where $r_{\text {hub }}$ is the hub radius.
We use the BEM implementation by Ning [29], which guarantees convergence employing a one-dimensional residual function. Furthermore, for inductions greater than 0.4 , we include the thrust-induction correction by Buhl [30]. To illustrate the effect of the tip correction, we consider the NREL-5MW rotor in a uniform inflow of $8 \mathrm{~m} / \mathrm{s}$, rotating at an angular velocity of 9.1552 rpm , at a tip-speed ratio (TSR) of 7.55 [31]. Fig. 2 reports the radial distribution of velocity components along the blade, showing that tip correction effects are most evident in the axial velocity. Without tip correction, the normalized velocity increases towards the tip. With tip correction, the axial velocity decreases towards the tip, reflecting the impact of three-dimensional effects such as tip vortices that deform and decelerate the flow in the tip-blade region [22,23]. Fig. 2b illustrates that the tip correction has a limited impact on the tangential velocity $V$, which is primarily determined by the blade's rotational speed. The power and thrust coefficients with tip correction are $C_{P}=0.4898$ and $C_{T}=0.7892$. These values are $6.7 \%$ and $2.5 \%$ lower than those obtained without tip correction, where $C_{P}$ is 0.5249 and $C_{T}$ is 0.8096 .

## 3. Classic actuator line method

In the classic ALM the blades are represented by lines discretized by actuator points. By interpolation of the velocity field, a velocity vector is assigned to each blade point, and a smearing function distributes the computed forces to the grid points of the flow solver. Typically, a Gaussian projection

$$
\begin{equation*}
\eta(\boldsymbol{d})=\frac{1}{\varepsilon^{3} \pi^{3 / 2}} \exp \left(-\frac{|\boldsymbol{d}|^{2}}{\varepsilon^{2}}\right), \tag{9}
\end{equation*}
$$

centered around the blade point, with $\boldsymbol{d}$ the distance from the blade point and $\varepsilon=2 \Delta_{\text {grid }}$ the smearing radius, is employed [32-34]. For a detailed description of the results obtained with this formulation we refer the reader to Section 6.1 and to the works by $[5,6,35]$. The sampled velocity plays a crucial role in the determination of the blade loads, as it is used to compute both the angle of attack and the loads themselves. Moreover, since the loadings are smeared to the flow solver, a non-linear interaction between the loads and the sampled velocity is generated. The consequence is that when a relatively coarse grid is employed, the smearing width extends to a distance much larger than the chord at the blade tip (about 1.42 m for the NREL-5MW rotor, and 1.14 m for the DTU-10MW rotor), causing a more rigid wake and therefore a missing induction at the rotor that resembles the one of BEM without tip correction (see also Section 2). In Section 6.1, we provide a detailed comparison between the sampled velocity of classic ALM and the induction predicted by BEM without tip correction, finding good agreement between the two (see Fig. 4).

## 4. Correcting the sampled velocity

The proposed correction aims at ensuring that the ALM results capture the ones obtained using BEM with tip correction, which are often used as benchmark for ALM [see for example 17,25]. This approach also ensures consistency with Lifting Line theory (see Section 6.2). The algorithm leverages BEM to compute the ratio between tip-corrected and non-tip-corrected velocity components, which only requires the turbine's TSR and blade airfoil data. Indeed, by incorporating in the wind angle definition Eq. (2) both the TSR definition, i.e.

$$
\begin{equation*}
\mathrm{TSR}=\frac{\Omega R}{U_{\infty}}, \tag{10}
\end{equation*}
$$

and Eq. (3), we obtain

$$
\begin{equation*}
\phi=\arctan \left(\frac{U_{\infty}}{\Omega r} \frac{(1-a)}{\left(1+a^{\prime}\right)}\right)=\arctan \left(\frac{R / r}{\operatorname{TSR}} \frac{(1-a)}{\left(1+a^{\prime}\right)}\right), \tag{11}
\end{equation*}
$$

which is used as first step for the BEM iterative algorithm. Then the two velocity ratios, $U_{\text {tip }} / U_{\text {notip }}$ for the axial velocity component and $V_{\text {tip }} / V_{\text {notip }}$ for the tangential one, can be found from the obtained induction factors (as shown in Fig. 3) without the need of specifying the inflow velocity and with guaranteed convergence [29]. In simulations with uniform inflows, the TSR is constant, making the correction constant over time and therefore not requiring any additional computational cost. For cases with variable TSR, as shown by Diaz et al. [36], the TSR can be estimated iteratively, but the computational overhead of the correction is still negligible. In particular, from axial momentum theory, by using the definition of axial induction factor $a=\left(U_{\infty}-U_{d}\right) / U_{\infty}$, with $U_{d}$ being the disk-averaged velocity and the thrust coefficient $C_{T}=4 a(1-a)$ [22, equation 9.27 therein], one can derive that

$$
\begin{equation*}
U_{\mathrm{ref}}=\frac{2\left\langle U_{\mathrm{d}}\right\rangle}{1+\sqrt{1-C_{T}}} . \tag{12}
\end{equation*}
$$

In Eq. (12), $C_{T}$ is function of $U_{\text {ref }}$ as per turbine definition (see Jonkman et al. [31]); thus, Eq. (12) can be solved iteratively to find the tip-speed ratio TSR $=\Omega R / U_{\text {ref }}$, in which $\Omega=\Omega\left(U_{\text {ref }}\right)$ is the angular velocity. This enables the current method to be used even when the free-stream velocity is not directly accessible, for example, within a wind farm. Fig. 3 also shows how the velocity ratios are incorporated in the classic ALM. The full algorithm consists of the following steps:

1. BEM is used to compute $U_{\text {tip }} / U_{\text {notip }}$ and $V_{\text {tip }} / V_{\text {notip }}$ using the TSR as the only input;
2. Compute the tip-corrected velocity components: $U_{\mathrm{c}}=U \times U_{\text {tip }} / U_{\text {notip }}$ and $V_{\mathrm{c}}=V \times V_{\text {tip }} / V_{\text {notip }}$;
3. Determine the relative velocity to the blade point: $U_{\text {rel }}= \sqrt{U_{\mathrm{c}}^{2}+V_{\mathrm{c}}^{2}} ;$
4. Determine the wind angle and angle of attack at the blade point using: $\phi=\arctan \left(U_{\mathrm{c}} / V_{\mathrm{c}}\right)$ and $\alpha=\phi-(\gamma+\theta)$;

![](https://cdn.mathpix.com/cropped/ec93945f-e863-4bc4-af5e-0180bd348132-4.jpg?height=647&width=1367&top_left_y=185&top_left_x=350)
Fig. 3. Illustration of the algorithm and its incorporation into the traditional ALM framework. Light-blue boxes illustrate the iterative loop to calculate the velocity ratios. Green boxes indicate the classic ALM steps. The velocity ratios that capture the tip corrections are incorporated in step 2.

![](https://cdn.mathpix.com/cropped/ec93945f-e863-4bc4-af5e-0180bd348132-4.jpg?height=994&width=1369&top_left_y=983&top_left_x=348)
Fig. 4. Comparison of the radial distributions of (a) normal and (b) tangential loads, (c) axial velocity component, and (d) angle of attack between classic ALM, the present study ALM (utilizing 64 grid points per diameter) and BEM, both with and without tip correction.

5. By utilizing $\alpha$ and $\phi$, the aerodynamic coefficients $C_{L}, C_{D}$, and $C_{n}, C_{t}$ can be calculated Eq. (5), which are used to calculate the loads Eq. (7).

For a more detailed discussion about the computational efficiency of the method, we refer to Section 8.

## 5. Large Eddy simulation

To numerically study the wind turbine we use an in-house pseudospectral LES method. We select LES because using ALM in this framework is still challenging and LES guarantees to resolve small-scale structures such as tip vortices, which are known to have an impact
on the tip-blade loading. The governing equations are the filtered continuity and momentum conservation equations [37-39]:

$$
\begin{align*}
\nabla \cdot \tilde{\boldsymbol{u}} & =0 \\
\frac{\partial \tilde{\boldsymbol{u}}}{\partial t}+\tilde{\boldsymbol{u}} \cdot \nabla \tilde{\boldsymbol{u}} & =\boldsymbol{f}_{\mathrm{ALM}}-\nabla \tilde{p}-\nabla \cdot \boldsymbol{\tau} \tag{13}
\end{align*}
$$

where $\boldsymbol{u}$ is the velocity vector, $\boldsymbol{f}_{\text {ALM }}$ the forces on the blades calculated by the ALM, $p$ is the kinematic pressure, and $\tau$ is the sub-grid stresses (SGS) tensor. The tilde ${ }^{\sim}$ indicates that a filtered velocity field is considered. The anisotropic minimum dissipation model is used for the SGS deviatoric stresses [39,40]. The trace of the SGS stress tensor is absorbed into the filtered modified pressure $\tilde{p}^{*}=\tilde{p} / \rho_{0}-p_{\infty} / \rho_{0}+ \operatorname{Tr}(\boldsymbol{\tau}) / 3$, where $\tilde{p}$ is the kinematic pressure and $\rho_{0}$ is the air density

![](https://cdn.mathpix.com/cropped/ec93945f-e863-4bc4-af5e-0180bd348132-5.jpg?height=469&width=1191&top_left_y=180&top_left_x=437)
Fig. 5. Rotor power (a) and thrust (b) for different inflow velocities computed with the present study ALM, classic ALM (utilizing 128 grid points per diameter), BEM with and without tip correction, the Effective Velocity Model by Muscari et al. [16].

Table 1
Parameters adopted for the simulations with varying inflow velocities. From top to bottom: inflow velocity, angular velocity, tip-speed ratio, and pitch angle.
| $\mathrm{U}_{\infty}[\mathrm{m} / \mathrm{s}]$ | 4.0 | 5.5 | 7.0 | 8.5 | 10.0 | 11.5 | 13.0 | 14.5 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| $\Omega$ [rpm] | 7.190 | 7.707 | 8.465 | 9.727 | 11.444 | 12.073 | 12.073 | 12.073 |
| TSR | 11.86 | 9.25 | 7.98 | 7.55 | 7.55 | 6.92 | 6.13 | 5.49 |
| $\theta\left[^{\circ}\right]$ | 0 | 0 | 0 | 0 | 0 | 1.283 | 6.560 | 9.529 |


here $1.225 \mathrm{~kg} / \mathrm{m}^{3}$. We use uniform inflow and free-slip boundary conditions with zero vertical velocity at the bottom and top of the domain. The time integration is performed using a third-order accurate Adams-Bashforth scheme. Spatial derivatives in the vertical direction are calculated using a second-order central finite difference scheme. The pseudo-spectral method is applied in the horizontal directions, resulting in periodic boundary conditions in the streamwise and lateral directions.

For the baseline simulations used to validate the method we adopt the NREL-5MW reference turbine rotor [31], which has a rotor diameter $D$ of 126 m . It is simulated considering a uniform inflow with a velocity of $8 \mathrm{~m} / \mathrm{s}$. At this inflow velocity, the rotational speed of the rotor is set to 9.1552 rpm , corresponding to a TSR of 7.55 . The number of blade points is kept constant at 60 , which ensures convergence of the obtained results across the adopted mesh spacings since the criterion of $\Delta_{\text {blade }} / \Delta_{\text {grid }}<3, \Delta_{\text {blade }}$ being the spacing between actuator points, is always satisfied [7,35]. When different simulation setups have been used for validation, these are outlined in the specific sections.

The computational domain is discretized with $n_{x}, n_{y}$, and $n_{z}$ points in the streamwise, spanwise, and vertical directions. A uniform grid is used in every direction, with a corresponding grid spacing of $\Delta_{\text {grid }}= L_{x} / n_{x}=L_{y} / n_{y}=L_{z} / n_{z}$, where $L_{x}, L_{y}, L_{z}$ are the dimensions of the computational domain. The computational grid is staggered in the vertical direction. The simulation domain is $L_{x} \times L_{y} \times L_{z}=8 D \times 6 D \times 6 D$. The turbine rotor hub is placed at three diameters after the inlet and in the middle of the $\mathrm{y}-\mathrm{z}$ plane. This computational setup is based on previous ALM studies [8,17,34,35,41,42]. We ensure that the statistical stationary state is reached and that the blade tip moves less than one grid point per time step. Blades are set perpendicular to the incoming flow, and the chord, twist, and aerodynamic coefficient data are interpolated in the radial direction to smooth out discontinuities in the turbine characteristics [6].

## 6. Validation using the NREL-5MW turbine rotor

### 6.1. Validation against Blade Element Momentum theory

In Fig. 4 we compare the present study ALM with the classic ALM and BEM, with and without tip correction. For this analysis, a grid discretization of 64 grid points per diameter, corresponding to $\Delta_{\text {grid }}=$

2 m , is employed for ALM. Panels a and b of Fig. 4 illustrate that the classic ALM consistently overestimates the blade loadings of BEM with tip correction, and this is due to the increasing trend of axial velocity (panel c) towards the tip, which causes the angle of attack to increase in the same region. Furthermore, the results of classic ALM resemble the ones of BEM without tip correction, a point highlighted in Section 3. In contrast, when the proposed correction is applied to ALM, the results approach the desired distribution of BEM with tip correction (green lines in Fig. 4) thanks to the reduction of axial velocity near the tip, leading also to an improved distribution of angle of attack. This trend finds confirmation in experimental observations [43] and bladeresolved simulations [44]. There are small differences between BEM and ALM within the mid-blade region, originating from the different assumptions in BEM and ALM. Specifically, BEM assumes independence among blade points, whereas in ALM blade points influence each other through the flow solver [35]. Since the ratio of tip-corrected to non-tip-corrected velocities approaches to 1 in the mid-blade region. This difference does not affect the effectiveness of the proposed method, which is instead focused on the tip region only. To study the effectiveness of the proposed method under varying inflow conditions, we perform eight simulations with inflow velocity ranging from 4 to $14.5 \mathrm{~m} / \mathrm{s}$. We use 128 grid points per diameter, or $\Delta_{\text {grid }}=1 \mathrm{~m}$. The chosen velocities correspond to a range of Reynolds number between $R e_{D}=3.4 \times 10^{7}$ and $1.2 \times 10^{8}$. The TSR, angular velocity and pitch angle $\theta$ are reported in Table 1 and are obtained from the turbine definition Jonkman et al. [31]. Fig. 5 shows the comparison among the present study, BEM with tip correction, classic ALM, BEM without tip correction and the Effective Velocity Model (EVM) by [16], based on sampling the velocity on a line upstream of the rotor. Here, we report the aerodynamic thrust, i.e. exempt from the rotor's weight as in the report by [31], and the power has been multiplied by the constant electrical generator efficiency, $94.4 \%$. The agreement between the proposed method and BEM with tip correction is excellent for all inflow velocities. Moreover, the difference between the present study ALM and the classic ALM increases at higher velocities due to the cubic dependence of the power on the sampled velocity, highlighting that a correction of tip-blade loading is even more important when the inflow velocity is increased. The figure also demonstrates the robustness of the proposed algorithm across a wide range of flow conditions, underlying the possibility to employ the proposed correction in various flow scenarios. The present method significantly improves on both classic ALM and EVM at all velocities, especially in rotor power, where the maximum difference between ALM and the present study ALM is reduced to $4.6 \%$ from $10.2 \%$ at $11.5 \mathrm{~m} / \mathrm{s}$. The rotor thrust is also improved, although at this discretization, the classic ALM deviates less from BEM already.

### 6.2. Validation against vortex-based smearing correction

To further validate our results, we integrated the vortex-based smearing correction method proposed by Forsting et al. [18,19] into

![](https://cdn.mathpix.com/cropped/ec93945f-e863-4bc4-af5e-0180bd348132-6.jpg?height=1007&width=1361&top_left_y=180&top_left_x=353)
Fig. 6. Comparison of the radial distributions of (a-b) normal and tangential loads, (c) axial velocity component, and (d) angle of attack obtained using present study ALM and the vortex-based smearing correction method (utilizing 64 grid points per diameter), against BEM with tip correction.

our computational framework. The integration implemented in our code guarantees that the different ALM versions can be compared accurately, isolating the effect of the methodology. To analyze the effect of the vortex-based smearing correction method on the axial velocity we use 64 grid points per diameter, or $\Delta_{\text {grid }}=2 \mathrm{~m}$. As illustrated in Fig. 6c, the axial velocity distribution at the blade closely matches both the one obtained with the present method and the one with BEM with tip correction. In the mid-blade region the tangential loading and axial velocity retrieve the difference with respect to BEM discussed in the previous section. Moreover, Figs. 6a and 6b show that the normal and tangential loads calculated by both methods are nearly identical along the entire blade. Since the method by [18] is designed to be consistent with Lifting Line theory, we can conclude that also the proposed method is consistent with the latter theory, and demonstrating the reliability of BEM as a reference for the results.

Finally, it is worth remarking that the method by [18] was developed in a RANS framework obtaining similar results, hence given the consistency between the two methods, we expect that the effectiveness of the proposed algorithm would be unchanged when a different flow solver is adopted, as it mainly depends on the rotor induction, which in this case is less related to the turbulence characteristics of the inflow.

### 6.3. Validation at different grid discretizations

The sampled velocity within ALM is known to be influenced by numerical discretization, particularly near the blade tip. This dependency arises from the use of a variable smearing width equal to $\varepsilon= 2 \Delta_{\text {grid }}$ to ensure that the integrated force remains conserved as the grid resolution varies. This effect is well established in the literature (see, e.g. Dağ and Sørensen [17]). Fig. 7 shows the radial distributions of the corrected loading, axial velocity and angle of attack when the grid resolution is varied from 32 to 128 grid points per diameter. As the figure illustrates, the present study ALM approaches the results of BEM with tip correction as the grid spacing ( $\Delta_{\text {grid }}$ ) decreases. Instead, the results of classic ALM overestimate the blade loading of BEM with tip

Table 2
Convergence of the thrust and power coefficients from the present study ALM and the vortex-based smearing correction method by [18] with mesh refinement.
| Simulation | $C_{P}$ | Error \% | $C_{T}$ | Error \% |
| :--- | :--- | :--- | :--- | :--- |
| BEM with tip correction | 0.4898 |  | 0.7892 |  |
| Present study, $D / \Delta_{\text {grid }}=32$ | 0.5443 | 11.1 | 0.8181 | 3.7 |
| $D / \Delta_{\text {grid }}=64$ | 0.5241 | 7.0 | 0.8050 | 2.0 |
| $D / \Delta_{\text {grid }}=128$ | 0.5122 | 4.6 | 0.7949 | 0.7 |
| Vortex-smearing correction, |  |  |  |  |
| $D / \Delta_{\text {grid }}=32$ | 0.5164 | 5.4 | 0.8015 | 1.6 |
| $D / \Delta_{\text {grid }}=64$ | 0.5148 | 5.1 | 0.8010 | 1.5 |
| Classic ALM, $D / \Delta_{\text {grid }}=32$ | 0.5767 | 17.8 | 0.8337 | 5.6 |
| $D / \Delta_{\text {grid }}=64$ | 0.5554 | 13.4 | 0.8209 | 4.0 |
| $D / \Delta_{\text {grid }}=128$ | 0.5400 | 10.2 | 0.8106 | 2.7 |


correction even at the finest grid adopted here (see [35]).
Table 2 presents the relative difference in the power ( $C_{P}$ ) and thrust $\left(C_{T}\right)$ coefficients obtained with the proposed correction, BEM theory with tip correction and the vortex-based smearing correction method by [18]. At the finest grid spacing, the difference between the present study ALM and BEM reduces to $4.6 \%$ for the power coefficient ( $C_{P}$ ) and $0.7 \%$ for the thrust coefficient ( $C_{T}$ ), while the difference is higher in the case of classic ALM ( $10.2 \%$ and $2.7 \%$, respectively) as also mentioned earlier. These aerodynamic coefficients do improve when the correction is used as the tip region has a larger swept area and relative velocity, underscoring the importance to correctly model the tip-blades loading. It is however important to recognize that perfect convergence is not anticipated due to the previously discussed inherent differences between BEM and ALM [see also 35].

## 7. Validation using the DTU- 10 MW turbine rotor

To further validate the versatility of the proposed method, we adopt a different rotor geometry, employing the DTU-10MW reference

![](https://cdn.mathpix.com/cropped/ec93945f-e863-4bc4-af5e-0180bd348132-7.jpg?height=994&width=1367&top_left_y=180&top_left_x=350)
Fig. 7. Comparison of the radial distribution of the (a) normal loads, (b) tangential loads, (c) axial velocity, and (d) the axial velocity ratio obtained using the present study ALM (using 128, 64 and 32 grid points per diameter) against BEM with tip correction.

![](https://cdn.mathpix.com/cropped/ec93945f-e863-4bc4-af5e-0180bd348132-7.jpg?height=554&width=1373&top_left_y=1294&top_left_x=346)
Fig. 8. Comparison of the radial distributions of (a) axial velocity and (b) angle of attack between the present study ALM (utilizing 90.6 grid points per diameter) against classic ALM and BEM, both with and without tip correction, using the DTU- 10 MW wind turbine rotor.

rotor [45]. Here, we use 90.6 grid points per diameter (or $\Delta_{\text {grid }}=2$ m ) and the same numerical domain used in the previous sections. The rotor diameter is 178.3 m , and the inflow velocity is set to $9 \mathrm{~m} / \mathrm{s}$, which corresponds to a Reynolds number $\operatorname{Re}_{D}=1.1 \times 10^{8}$. The angular velocity of the rotor is 7.229 rpm , the corresponding TSR is 7.67 and the number of actuator points per blade used is again set to 60 for both ALM and BEM.

Fig. 8a shows the resulting axial velocity distribution of the present study ALM and classic ALM with respect to the one computed with BEM. While the distribution of axial velocity in the case of BEM without tip correction tends to increase towards the tip and classic ALM follows the same trend (as expected from Section 6.1), the present study ALM agrees well with the distribution of BEM with tip correction in the tip region. This demonstrates that the present method is suitable for different turbine geometries and sizes. The difference between ALM and BEM distributions in the mid-blade is independent on the correction (see Sections 6.1 and 6.2). Fig. 8b shows that the angle of attack
distribution is by consequence in very good agreement with BEM with tip correction as well.

## 8. Computational considerations

A key benefit of our approach is that it is fully parallelizable due to the independence of actuator points within the BEM framework. Additionally, efficient BEM algorithms are readily available, and its computational overhead is negligible compared to the flow solver calculations. This guarantees that the method is highly efficient, scalable, and energy-effective for simulations on modern heterogeneous computing platforms, especially those utilizing accelerators like GPUs. Vortex-based smearing corrections, instead, require computing induction corrections in three-dimensional space for each blade point, making their cost increase with the number of actuator points; they also necessitate additional communication between parallel processes. While non-iterative approaches and optimized algorithms based
on wake parameterizations exist, these limitations are not entirely resolved.

We note that using local data, however, leads to some trade-off in accuracy compared to methods that use inter-processor communication and non-local physical information. Our approach to use BEM to correct the sampled velocity is more effective than tip-correction methods like Shen et al. [46], which is focused only on matching loading conditions a posteriori. The simplicity of our approach also allows it to be adopted in RANS or URANS settings where computational simplicity, efficiency, and robustness are essential. Indeed, when different flow solvers are used-as shown by Nathan et al. [47] (using DES and LES), Meyer Forsting et al. [18] (using RANS), and Stevens et al. [48] (using LES)consistently overestimate loads compared to BEM with tip correction. Remarkably, Nathan et al. [47] show that axial velocity obtained using DES and LES agree well and increase towards the tip, consistent with our analysis in Section 6.1.

## 9. Conclusions

This work aims to address the issue of tip-loading overestimation commonly observed in simulations with wind turbine models, focusing specifically on the Actuator Line Method (ALM). Initial approaches, such as the one by Shen et al. [46] used Prandtl tip corrections to correct loading at the blades, but do not recover the correct angle of attack and induction at the blade. We show that the loading overestimation results from an overestimation of the sampled velocity in the tip region. The key novelty of the present approach is to use BEM to selfconsistently correct the sampled velocity at the blade location, relying only on local information. This approach eliminates communication between processors, enhancing computational efficiency in modern computing systems like GPUs.

The proposed algorithm relies only on basic turbine informationspecifically, the tip speed ratio (TSR), airfoil blade data, and the blade sampled velocity. Since the TSR can be estimated from local flow properties, this approach can be readily applied across various flow scenarios, including wind farm simulations. Its consistency with the widely used BEM theory enables seamless integration with tools built on this extensively utilized framework [see for example $6,24,27,49]$. The proposed method has undergone extensive validation against BEM with tip correction and a vortex-smearing correction technique by Meyer Forsting et al. [18], both consistent with Lifting Line theory, across diverse inflow velocities, grid discretizations, and using two turbine rotors, namely the NREL-5MW and DTU-10MW. It has demonstrated reliable accuracy in estimating turbine power and thrust, making it well-suited for applications that require a simple yet dependable turbine modeling approach.

## CRediT authorship contribution statement

Davide Selvatici: Writing - review \& editing, Writing - original draft, Validation, Software, Methodology, Investigation, Formal analysis, Data curation, Conceptualization. Richard J.A.M. Stevens: Writing - review \& editing, Supervision, Software, Project administration, Funding acquisition, Conceptualization.

## Declaration of competing interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Acknowledgments

This project has received funding from the European Research Council Horizon Europe program (Grant No. 101124815). This work was partly carried out on the Dutch national e-infrastructure with the support of SURF Cooperative. We acknowledge the EuroHPC Joint Undertaking for awarding the project EHPC-REG-2023R03-178 access to the EuroHPC supercomputer Discoverer, hosted by Sofia Tech Park (Bulgaria).

## References

[1] R.J.A.M. Stevens, C. Meneveau, Flow structure and turbulence in wind farms, Annu. Rev. Fluid Mech. 49 (1) (2017) 311-339.
[2] F. Porté-Agel, M. Bastankhah, S. Shamsoddin, Wind-turbine and wind-farm flows: A review, Bound.-Layer Meteorol. 174 (1) (2020) 1-59.
[3] Z. Li, X. Liu, X. Yang, Review of turbine parameterization models for large-eddy simulation of wind turbine wakes, Energies 15 (18) (2022).
[4] C.R. Shapiro, D.F. Gayme, C. Meneveau, Filtered actuator disks: Theory and application to wind turbine models in large eddy simulation, Wind Energy 22 (10) (2019) 1414-1420.
[5] L.A. Martínez-Tossas, M.J. Churchfield, S. Leonardi, Large eddy simulations of the flow past wind turbines: actuator line and disk modeling, Wind Energy 18 (6) (2015) 1047-1060.
[6] J.N. Sørensen, General Momentum Theory for Horizontal Axis Wind Turbines, Springer, 2016.
[7] P.K. Jha, M.J. Churchfield, P.J. Moriarty, S. Schmitz, Guidelines for volume force distributions within actuator line modeling of wind turbines on large-eddy simulation-type grids, J. Sol. Energy Eng. 136 (3) (2014).
[8] P.K. Jha, S. Schmitz, Actuator curve embedding - an advanced actuator line model, J. Fluid Mech. 834 (2017).
[9] M.J. Churchfield, S. Schreck, L.A. Martínez-Tossas, C. Meneveau, P.R. Spalart, An Advanced Actuator Line Method for Wind Energy Applications and Beyond, AIAA SciTech Forum, 2017.
[10] J.N. Sørensen, K.O. Dag, N. Ramos-García, A refined tip correction based on decambering, Wind Energy 19 (5) (2016) 787-802.
[11] K.O. Dağ, Combined pseudo-spectral / actuator line model for wind turbine applications, 2017, DTU Wind Energy PhD Vol. 67.
[12] W.Z. Shen, R. Mikkelsen, J.N. Sørensen, C. Bak, Tip loss corrections for wind turbine computations, Wind Energy 8 (4) (2005) 457-475.
[13] X. Yang, F. Sotiropoulos, R.J. Conzemius, J.N. Wachtler, M.B. Strong, Largeeddy simulation of turbulent flow past wind turbines/farms: the Virtual Wind Simulator (VWiS), Wind Energy 18 (12) (2015) 2025-2045.
[14] S. Xie, An actuator-line model with Lagrangian-averaged velocity sampling and piecewise projection for wind turbine simulations, Wind Energy 24 (10) (2021) 1095-1106.
[15] A.G. Sanvito, A. Firpo, P. Schito, V. Dossena, A. Zasso, G. Persico, A novel vortex-based velocity sampling method for the actuator-line modeling of floating offshore wind turbines in windmill state, Renew. Energy 231 (2024) 120927.
[16] C. Muscari, P. Schito, A. Viré, A. Zasso, J.-W. van Wingerden, The effective velocity model: An improved approach to velocity sampling in actuator line models, Wind Energy 27 (5) (2024) 447-462.
[17] K.O. Dağ, J.N. Sørensen, A new tip correction for actuator line computations, Wind Energy 23 (2) (2019) 148-160.
[18] A.R.M. Forsting, G.R. Pirrung, N. Ramos-García, A vortex-based tip/smearing correction for the actuator line, Wind Energy Sci. 4 (2) (2019) 369-383.
[19] A.R.M. Forsting, G.R. Pirrung, N. Ramos-García, Brief communication: A fast vortex-based smearing correction for the actuator line, Wind Energy Sci. 5 (1) (2020) 349-353.
[20] L.A. Martínez-Tossas, C. Meneveau, Filtered lifting line theory and application to the actuator line model, J. Fluid Mech. 863 (2019) 269-292.
[21] V.G. Kleine, A. Hanifi, D.S. Henningson, Non-iterative vortex-based smearing correction for the actuator line method, J. Fluid Mech. 961 (2023).
[22] E. Branlard, Wind Turbine Aerodynamics and Vorticity-Based Methods, Springer, 2017.
[23] H. Glauert, Airplane Propellers, Springer, 1935.
[24] H.A. Madsen, C. Bak, M. Døssing, R. Mikkelsen, S. Øye, Validation and modification of the blade element momentum theory based on comparisons with actuator disc simulations, Wind Energy 13 (4) (2010) 373-389.
[25] G. Bangga, Comparison of blade element method and CFD simulations of a 10 MW wind turbine, Fluids 3 (4) (2018).
[26] G. Bangga, T. Lutz, Aerodynamic modeling of wind turbine loads exposed to turbulent inflow and validation with experimental data, Energy 223 (2021).
[27] A. Li, G.R. Pirrung, M. Gaunaa, H.A. Madsen, S.G. Horcas, A computationally efficient engineering aerodynamic model for swept wind turbine blades, Wind Energy Sci. 7 (1) (2022) 129-160.
[28] T. Burton, D. Sharpe, N. Jenkins, E. Bossanyi, Wind Energy Handbook, Wiley, 2001.
[29] S.A. Ning, A simple solution method for the blade element momentum equations with guaranteed convergence, Wind Energy 17 (2013) 1327-1345.
[30] M.L. Buhl, A new empirical relationship between thrust coefficient and induction factor for the turbulent windmill state, Technical Report, National Renewable Energy Laboratory, 2005.
[31] J. Jonkman, S. Butterfield, W. Musial, G. Scott, Definition of a 5-MW reference wind turbine for offshore system development, Technical Report, National Renewable Energy Laboratory (NREL), 2009.
[32] J.N. Sørensen, W.Z. Shen, Numerical modeling of wind turbine wakes, J. Fluids Eng. 124 (2) (2002) 393-399.
[33] N. Troldborg, Actuator Line Modeling of Wind Turbine Wakes (Ph.D. thesis), Technical University of Denmark, 2008.
[34] L.A. Martínez-Tossas, M.J. Churchfield, A.E. Yilmaz, H. Sarlak, P.L. Johnson, J.N. Sørensen, J. Meyers, C. Meneveau, Comparison of four large-eddy simulation research codes and effects of model coefficient and inflow turbulence in actuator-line-based wind turbine modeling, Renew. Sustain. Energy Rev. 10 (3) (2018) 033301.
[35] L. Liu, L. Franceschini, D.F. Oliveira, F.C.C. Galeazzo, B.S. Carmo, R.J.A.M. Stevens, Evaluating the accuracy of the actuator line model against blade element momentum theory in uniform inflow, Wind Energy 25 (6) (2022) 1046-1059.
[36] G.P.N. Diaz, A.D. Otero, H. Asmuth, J.N. Sørensen, S. Ivanell, Actuator line model using simplified force calculation methods, Wind Energy Sci. 8 (3) (2023) 363-382.
[37] R.J.A.M. Stevens, J. Graham, C. Meneveau, A concurrent precursor inflow method for large eddy simulations and applications to finite length wind farms, Renew. Energy 68 (2014) 46-50.
[38] R.J.A.M. Stevens, M. Wilczek, C. Meneveau, Large-eddy simulation study of the logarithmic law for second- and higher-order moments in turbulent wall-bounded flow, J. Fluid Mech. 757 (2014) 888-907.
[39] S.N. Gadde, A. Stieren, R.J.A.M. Stevens, Large-eddy simulations of stratified atmospheric boundary layers: Comparison of different subgrid models, Bound.-Layer Meteorol. 178 (3) (2021) 363-382.
[40] W. Rozema, R.W.C.P. Verstappen, A.E.P. Veldman, J.C. Kok, Low-dissipation simulation methods and models for turbulent subsonic flow, Arch. Comput. Methods Eng. 27 (1) (2018) 299-330.
[41] M.H.A. Madsen, F. Zahle, N.N.S. rensen, J.R.R.A. Martins, Multipoint high-fidelity CFD-based aerodynamic shape optimization of a 10 MW wind turbine, Wind Energy Sci. 4 (2) (2019) 163-192.
[42] Z. Sun, W.J. Zhu, W.Z. Shen, W. Zhong, J. Cao, Q. Tao, Aerodynamic analysis of coning effects on the DTU 10 MW wind turbine rotor, Energies 13 (21) (2020) 5753.
[43] S.I. Green, Wing Tip Vortices. Fluid Mechanics and Its Applications, Springer, 1995.
[44] P.F. Melani, O.S. Mohamed, S. Cioni, F. Balduzzi, A. Bianchini, An insight into the capability of the actuator line method to resolve tip vortices, Wind Energy Sci. 9 (3) (2024) 601-622.
[45] C. Bak, F. Zahle, R. Bitsche, T. Kim, A. Yde, L.C. Henriksen, A. Natarajan, M. Hansen, Description of the DTU 10 MW reference wind turbine, Technical Report, DTU Wind Energy, 2013.
[46] W.Z. Shen, J.N. Sørensen, R. Mikkelsen, Tip loss correction for actuator/NavierStokes computations, J. Sol. Energy Eng. 127 (2) (2005) 209-213.
[47] J. Nathan, A.R.M. Forsting, N. Troldborg, C. Masson, Comparison of OpenFOAM and EllipSys3D actuator line methods with (NEW) MEXICO results, J. Phys.: Conf. Series 854 (2017).
[48] R.J.A.M. Stevens, L.A. Martínez-Tossas, C. Meneveau, Comparison of wind farm large eddy simulations using actuator disk and actuator line models with wind tunnel experiments, Renew. Energy 116 (2018) 470-478.
[49] C. Crawford, Re-examining the precepts of the blade element momentum theory for coning rotors, Wind Energy 9 (5) (2006) 457-478.


[^0]:    * Corresponding author.

    E-mail address: d.selvatici@utwente.nl (D. Selvatici).

