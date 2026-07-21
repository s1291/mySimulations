# Brief communication: A fast vortex-based smearing correction for the actuator line 

Alexander R. Meyer Forsting, Georg R. Pirrung, and Néstor Ramos-García<br>DTU Wind Energy, Technical University of Denmark, Frederiksborgvej 399, 4000 Roskilde, Denmark<br>Correspondence: Alexander R. Meyer Forsting (alrf@dtu.dk)

Received: 18 September 2019 - Discussion started: 10 October 2019
Revised: 27 January 2020 - Accepted: 20 February 2020 - Published: 23 March 2020


#### Abstract

The actuator line is a lifting line representation of aerodynamic surfaces in computational fluid dynamics applications but with non-singular forces, which reduces the self-induced velocities at the line. The vortex-based correction by Meyer Forsting et al. (2019a) recovers this missing induction and thus the intended lifting line behaviour of the actuator line. However, its computational cost exceeds that of existing tip corrections and quickly grows with blade discretization. Here we present different methods for reducing its computational cost to the level of existing corrections without jeopardizing the stability or accuracy of the original method. The cost is reduced by at least $98 \%$, whereas the power is maximally affected by $0.8 \%$ with respect to the original formulation. This accelerated smearing correction remains a dynamic correction by modelling the variation in trailed vorticity over time. The correction is openly available (Meyer Forsting et al., 2019b).


## 1 Introduction

The actuator line (AL) Sørensen and Shen (2002) is a lifting line (LL) representation of aerodynamic surfaces in Eulerian computational fluid dynamics (CFD) applications. It allows simulating the interaction between the atmosphere and wind farms, as it captures all the important flow features of fully resolved rotors, at a fraction of the computational cost. However transferring a LL into the CFD domain requires dispersing the concentrated blade forces of the LL over a certain region - most commonly in the form of a Gaussian projection to avoid causing numerical instabilities. This force smearing leads to the formation of a viscous core in the released vorticity, which subsequently reduces the induced velocity at the blade (Dag, 2017; Meyer Forsting et al., 2019a; MartínezTossas and Meneveau, 2019). Lower induction implies larger angles of attack and thus increased blade forces. Especially in regions presenting large load changes, as around the root and tip of the blade, does the AL thus overestimate the forces.

Meyer Forsting et al. (2019a) - following the approach proposed by Dag (2017) - presented a correction to the AL that combines the fast and dynamic near-wake model by Pirrung et al. (2016, 2017a, b) with a viscous core model (Lamb, 1932; Oseen, 1911) to recover the missing induction. With
the correction, the AL truly functions as a LL, which was verified over the entire operational wind speed range of modern turbines as well as in yaw and for dynamic pitch steps (Meyer Forsting et al., 2019a). The numerical stability of the correction was not challenged by any of those flow cases not even by extreme inflow turbulence.

The only disadvantage of the new smearing correction is its computational cost. Though it is incorrect to apply conventional tip corrections to ALs - they correct actuator discs for missing discrete blades - their low cost makes them attractive. In this paper we present different methods that reduce the computational cost of the new correction to that of existing corrections without jeopardizing the stability or accuracy of the method.

## 2 Methods for increasing speed

Computing the missing induction requires re-evaluating the velocity contribution from each previously released vortex element at each time step. The velocity contribution from a single trailed vortex at some point along the blade is obtained by integrating along the vortex line

$$
\begin{equation*}
\boldsymbol{u}^{\star}=\int_{0}^{\infty} f_{\epsilon} \delta \tilde{\boldsymbol{u}} \mathrm{d} l \tag{1}
\end{equation*}
$$

Here $\boldsymbol{\delta} \tilde{\boldsymbol{u}}$ is the velocity induced by an infinitesimal element $\delta l$ of a vortex line and $f_{\epsilon}$ represents the smearing factor, originating from the presence of a viscous core in the released vorticity. Integrating over the vortex length is equivalent to integrating over time, as at each time step an element is released. Originally, the near-wake model by Pirrung et al. (2016, 2017a, b) provides directly the integrated velocities $\tilde{\boldsymbol{u}} .^{1}$ It was only broken into elements, as $f_{\epsilon}$ is a function of the perpendicular distance from the vortex to the blade element, which varies in time. As the distance changes at each time step, the velocity contribution from each vortex element also needs to be updated each time step. Hence the more vortex lines, the costlier the correction becomes.

### 2.1 Reduce wake length (orig. $\beta_{\text {max }}=\pi / 2$ )

In the work verifying the smearing correction by Meyer Forsting et al. (2019a), the integration along the vortex lines was performed until $\beta_{\text {max }}=2 \pi$, where $\beta$ defines the rotation angle, to ensure most induction is captured. However, the near-wake model is devised to provide only the induction from the vortex lines until $\beta=\pi / 2$. Considering that the vortex core effect is only active in the near-wake, $\beta_{\text {max }}$ could equally be set to $\pi / 2$, thus reducing the number of vortex elements significantly.

### 2.2 Reduce inner loops (cut loops)

The computational cost of vortex methods grows with the square of the blade elements, which could lead to escalating costs with increasing discretization. Usually, the induction of each vortex line on each blade section needs to be determined. Yet the limited size of the viscous core allows shortcutting this procedure by considering only the blade sections closest to the vortex line. The velocity missing in AL simulations in two dimensions is given by

$$
\begin{equation*}
v^{\star}(r, \epsilon)=\tilde{v} \overbrace{\exp \left(-r^{2} / \epsilon^{2}\right)}^{f_{\epsilon}}, \tag{2}
\end{equation*}
$$

with $r$ representing the distance from the vortex core and $\epsilon$ the force smearing length scale. To determine the size of the vortex core $r_{\text {max }}$, the ratio between cut and fully resolved vortex core is computed:

$$
\begin{equation*}
I=\frac{\int_{0}^{r_{\max }} v^{\star}(r, \epsilon) \mathrm{d} r}{\int_{0}^{\infty} v^{\star}(r, \epsilon) \mathrm{d} r} \tag{3}
\end{equation*}
$$

[^0]Different ratios were tested; however $I=0.99$ - corresponding to $r_{\text {max }}=1.83 \epsilon$-provides a beneficial balance between accuracy and speed.

### 2.3 Constant smearing factor, $f_{\epsilon}$ (fixed $x_{\perp}$ )

A more radical approach than just reducing the wake length, as described in Sect. 2.1, is fixing the perpendicular distance between the vortex and blade element and thus the smearing factor. In the three-dimensional formulation the smearing factor is given by (Meyer Forsting et al., 2019a)

$$
\begin{equation*}
f_{\epsilon}=\exp \left(-\frac{\left|\boldsymbol{x}_{\perp}(r, \beta, h, \phi)\right|^{2}}{\epsilon^{2}}\right) \tag{4}
\end{equation*}
$$

with the perpendicular distance

$$
\boldsymbol{x}_{\perp}=r \cos \phi\left(\begin{array}{c}
\tan \phi(\beta \cos \beta-\sin \beta)  \tag{5}\\
-\tan \phi(-1+h / r+\cos \beta+\beta \sin \beta) \\
-1+(1-h / r) \cos \beta
\end{array}\right) .
$$

The greatest simplification is achieved by setting $\beta=0$, such that $\boldsymbol{x}_{\perp}$ becomes the distance between vortex trailing point and blade section, which is a geometric constant for rigid blades.

$$
\begin{equation*}
\left|x_{\perp}(\beta=0)\right|=h \tag{6}
\end{equation*}
$$

The smearing factor no longer needs to be updated for all vortex elements at each time step, and the velocity correction in Eq. (1) simply becomes

$$
\begin{equation*}
\boldsymbol{u}^{\star}=f_{\epsilon} \tilde{\boldsymbol{u}}, \tag{7}
\end{equation*}
$$

where $\tilde{\boldsymbol{u}}$ is directly determined by the near-wake model. Thus it is very computationally efficient and does not require saving and integrating the induced velocities from discretized vortex arcs. At each time step and for each blade section, the influence from each previously trailed vortex arc can be simply updated by multiplying with an exponential decay factor and adding the influence of the newly trailed element (Pirrung et al., 2016). In this paper this method is run in conjunction with the cutting loops approach.

## 3 Results

This section compares the influence of the different speedup methods presented in Sect. 2 with the original results of Meyer Forsting et al. (2019a). All results are obtained with exactly the same computational set-up as presented in Sect. 3 of the same paper. The AL models the NREL 5 MW (Jonkerman et al., 2009) under uniform inflow. For the inflow wind speed specific turbine parameters refer to Table 2 in Meyer Forsting et al. (2019a).

Figure 1 compares the force distributions for the NREL 5 MW at two wind speeds obtained with the speed-up methods presented in Sect. 2 to those obtained with the original

![](https://cdn.mathpix.com/cropped/14c31a1d-b470-4013-af94-0b8adf07fa88-3.jpg?height=514&width=1200&top_left_y=217&top_left_x=432)
Figure 1. Normal and tangential forces on the NREL 5 MW blades at 8 and $25 \mathrm{~ms}^{-1}$ predicted by AL simulations (blades discretized by 19 sections) with smearing correction and different computational speed-up methods. The reference is the original formulation (orig.) with $\beta_{\text {max }}=2 \pi$.

![](https://cdn.mathpix.com/cropped/14c31a1d-b470-4013-af94-0b8adf07fa88-3.jpg?height=676&width=1202&top_left_y=920&top_left_x=430)
Figure 2. Difference in normal and tangential forces over the NREL 5MW blades at 8,14 and $25 \mathrm{~m} \mathrm{~s}^{-1}$ predicted by AL simulations (blades discretized by 19 sections) with different smearing correction speed-up methods with respect to simulations without speed-up.

model. The influence of reducing the wake length is only shown for a wind speed of $25 \mathrm{~m} \mathrm{~s}^{-1}$, but it is similar at lower wind speeds. With increasing wind speed, the peak force moves clearly from the tip to the root whilst the smearing correction ensures the smooth behaviour towards the blade ends. From pure visual inspection there is no change in the forces when applying any of the speed-up options. To highlight their impact, only the change in the force distributions with respect to the unmodified model is shown in Fig. 2 here additionally the results for a wind speed of $14 \mathrm{~ms}^{-1}$ are presented. Reducing the wake length has a negligible effect on the forces as does reducing the inner loops, except close to the root. Fixing the smearing factor additionally to cutting the loops has the largest influence. However even at $25 \mathrm{~m} \mathrm{~s}^{-1}$ the deviation does not exceed $17 \mathrm{Nm}^{-1}$. With respect to the local force the difference remains below $1 \%$.

A full result overview - the impact of the speed-up methods on thrust and power as well as their influence on the computational cost per blade - is given in Table 1. Results are shown for rotors discretized by 9 and 19 blade sections. Firstly, the greatest change in thrust or power across all methods occurs when fixing the smearing constant, yet never by more than $0.8 \%$ and only at the highest wind speed. The positive influence of cutting the inner loops on performance grows with increasing resolution. However, the largest reduction in the computational cost comes from limiting the wake length and ultimately fixing the smearing factor. With the latter approach the longest of all smearing correction iterations lasted $8 \times 10^{-4} \mathrm{~s}$.

Table 1. An overview of the influence of the computational speed-up methods on thrust, power and computational cost per blade for two different blade discretizations -9 and 19 blade sections. Only for the original model are the nominal values shown, otherwise the relative change to the original is given in percent.
| $N_{s}$ |  | $V_{\infty}\left(\mathrm{m} \mathrm{s}^{-1}\right)$ | Orig. $\beta_{\text {max }}=2 \pi$ | Orig. $\beta_{\text {max }}=\pi / 2$ (\%) | Cut loops (\%) | Fixed $x_{\perp}$ (\%) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | Thrust | 8 | 406 kN | - | $1.96 \times 10^{-2}$ | $-3.83 \times 10^{-2}$ |
|  |  | 14 | 466 kN | - | $1.35 \times 10^{-2}$ | $-1.79 \times 10^{-1}$ |
|  |  | 25 | 286 kN | $4.43 \times 10^{-3}$ | $-2.06 \times 10^{-2}$ | $-5.65 \times 10^{-1}$ |
|  | Power | 8 | 2.11 MW | - | $4.37 \times 10^{-2}$ | $-1.37 \times 10^{-1}$ |
|  |  | 14 | 5.43 MW | - | $1.88 \times 10^{-2}$ | $-2.75 \times 10^{-1}$ |
|  |  | 25 | 5.47 MW | $6.78 \times 10^{-3}$ | $-2.28 \times 10^{-2}$ | $-7.95 \times 10^{-1}$ |
|  | Cost | 8 | $8.68 \times 10^{-3} \mathrm{~s}$ | - | -43.5 | -99.3 |
|  |  | 14 | $9.91 \times 10^{-3} \mathrm{~s}$ | - | -44.9 | -98.3 |
|  |  | 25 | $1.04 \times 10^{-2} \mathrm{~s}$ | -83.4 | -45.9 | -99.4 |
| 19 | Thrust | 8 | 394 kN | - | $9.13 \times 10^{-2}$ | $2.81 \times 10^{-2}$ |
|  |  | 14 | 456 kN | - | $-7.40 \times 10^{-4}$ | $-9.81 \times 10^{-2}$ |
|  |  | 25 | 277 kN | $8.73 \times 10^{-5}$ | $-3.30 \times 10^{-2}$ | $-4.75 \times 10^{-1}$ |
|  | Power | 8 | 2.00 MW | - | $1.71 \times 10^{-1}$ | $1.43 \times 10^{-2}$ |
|  |  | 14 | 5.28 MW | - | $-1.80 \times 10^{-3}$ | $-1.56 \times 10^{-1}$ |
|  |  | 25 | 5.29 MW | $1.32 \times 10^{-4}$ | $-4.36 \times 10^{-2}$ | $-6.83 \times 10^{-1}$ |
|  | Cost | 8 | $4.35 \times 10^{-2} \mathrm{~s}$ | - | -63.3 | -99.0 |
|  |  | 14 | $4.45 \times 10^{-2} \mathrm{~s}$ | - | -63.7 | -98.2 |
|  |  | 25 | $4.37 \times 10^{-2} \mathrm{~s}$ | -82.9 | -62.6 | -99.0 |


## 4 Conclusions

The smearing correction by Meyer Forsting et al. (2019a) recovered the lifting line behaviour of the actuator line, however at a larger computational cost than existing actuator disc tip corrections. This paper presents different methods for reducing the cost of the smearing correction to those levels. The number of wake elements manifests itself as the key cost driver. Reducing the wake length therefore significantly reduces the computational cost without negatively impacting the blade forces. The greatest speed-up comes from utilizing the near-wake model to avoid recomputing the contributions from each element at each time step, leading to a fall in the cost of at least $98 \%$. This is accompanied by changes in thrust and power of maximally $0.8 \%$ and $0.7 \%$, respectively. Still, with respect to the great gain in performance this is acceptable and lies well within CFD simulation uncertainty. Furthermore, the new, faster method avoids any form of bookkeeping, greatly simplifying the implementation of the smearing correction. It also remains a dynamic correction that takes into account how the trailed vorticity changes over time and moves away from the blades. This faster and simpler version of the smearing correction is openly available (Meyer Forsting et al., 2019b).

Code availability. All data are available on request. Commercial and research licences for EllipSys3D can be purchased from DTU. The source code of the fast smearing correction is openly available (Meyer Forsting et al., 2019b).

Competing interests. The authors declare that they have no conflict of interest.

Special issue statement. This article is part of the special issue "Wind Energy Science Conference 2019". It is a result of the Wind Energy Science Conference 2019, Cork, Ireland, 17-20 June 2019.

Acknowledgements. We would like to acknowledge DTU Wind Energy's internal project Virtual Atmosphere for partially funding this research.

Financial support. This research has been supported by the DTU Wind Energy (project Virtual Atmosphere).

Review statement. This paper was edited by Rebecca Barthelmie and reviewed by David Wood and Claudio Balzani.

## References

Dag, K.: Combined pseudo-spectral/actuator line model for wind turbine applications, PhD thesis, DTU Wind Energy, Denmark, 2017.

Jonkman, J., Butterfield, S., Musial, W., and Scott, G.: Definition of a 5-MW reference wind turbine for offshore system development, Tech. rep., NREL/TP-500-38060, National Renewable Energy Laboratory (NREL), Colorado, USA, 2009.
Lamb, H.: Hydrodynamics, C.U.P., 6th Edn., Cambridge University Press, Cambridge, 1932.
Martínez-Tossas, L. A. and Meneveau, C.: Filtered lifting line theory and application to the actuator line model, J. Fluid Mech., 863, 269-292, https://doi.org/10.1017/jfm.2018.994, 2019.
Meyer Forsting, A. R., Pirrung, G. R., and Ramos-García, N.: A vortex-based tip/smearing correction for the actuator line, Wind Energ. Sci., 4, 369-383, https://doi.org/10.5194/wes-4-369-2019, 2019a.
Meyer Forsting, A. R., Pirrung, G. R., and Ramos-García, N.: Actuator-Line-Smearing-Correction, DTU Data, Technical University of Denmark, https://doi.org/10.11583/DTU.9752285.v1, 2019b.

Oseen, C.: Über Wirbelbewegung in einer reibenden Flüssigkeit, Arkiv för matematik, astronomi och fysik, Ark. Mat. Astron. Fys., 7, 14-21, 1911.
Pirrung, G., Madsen, H. A., Kim, T., and Heinz, J.: A coupled near and far wake model for wind turbine aerodynamics, Wind Energy, 19, 2053-2069, https://doi.org/10.1002/we.1969, 2016.
Pirrung, G., Riziotis, V., Madsen, H., Hansen, M., and Kim, T.: Comparison of a coupled near- and far-wake model with a free-wake vortex code, Wind Energ. Sci., 2, 15-33, https://doi.org/10.5194/wes-2-15-2017, 2017a.
Pirrung, G. R., Madsen, H. A., and Schreck, S.: Trailed vorticity modeling for aeroelastic wind turbine simulations in standstill, Wind Energ. Sci., 2, 521-532, https://doi.org/10.5194/wes-2-521-2017, 2017b.
Sørensen, J. N. and Shen, W. Z.: Numerical modelling of wind turbine wakes, J. Fluid. Eng.-T. ASME, 124, 393-399, https://doi.org/10.1115/1.1471361, 2002.


[^0]:    ${ }^{1}$ Note that the integration only covers the near-wake region, from 0 to $\pi / 2$.

