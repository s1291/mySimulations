# Generalised grid requirements minimizing the actuator line angle-of-attack error 

Forsting, AR Meyer; Troldborg, Niels

## Published in:

Journal of Physics: Conference Series

Link to article, DOI:
10.1088/1742-6596/1618/5/052001

## Publication date:

2020

## Document Version

Publisher's PDF, also known as Version of record

Link back to DTU Orbit

## Citation (APA):

Forsting, AR. M., \& Troldborg, N. (2020). Generalised grid requirements minimizing the actuator line angle-ofattack error. Journal of Physics: Conference Series, 1618(5), Article 052001. https://doi.org/10.1088/17426596/1618/5/052001

## General rights

Copyright and moral rights for the publications made accessible in the public portal are retained by the authors and/or other copyright owners and it is a condition of accessing publications that users recognise and abide by the legal requirements associated with these rights.

- Users may download and print one copy of any publication from the public portal for the purpose of private study or research.
- You may not further distribute the material or use it for any profit-making activity or commercial gain
- You may freely distribute the URL identifying the publication in the public portal


# Generalised grid requirements minimizing the actuator line angle-ofattack error 

To cite this article: AR Meyer Forsting and N Troldborg 2020 J. Phys.: Conf. Ser. 1618052001

View the article online for updates and enhancements.
![](https://cdn.mathpix.com/cropped/f88616bd-9680-438f-a8a6-d8213c97b2ad-02.jpg?height=497&width=657&top_left_y=2256&top_left_x=127)

## IOP ebooks ${ }^{\text {™ }}$

Bringing together innovative digital publishing with leading authors from the global scientific community.

Start exploring the collection-download the first chapter of every title for free.

# Generalised grid requirements minimizing the actuator line angle-of-attack error 

AR Meyer Forsting ${ }^{1}$, N Troldborg ${ }^{1}$<br>${ }^{1}$ DTU Wind Energy, Technical University of Denmark, Ris $\varnothing$ Campus, DK-4000 Roskilde, Denmark<br>E-mail: alrf@dtu.dk


#### Abstract

The actuator line (AL) is a lifting line (LL) representation of aerodynamic surfaces in computational fluid dynamics (CFD) applications. The AL blade forces are computed from 2D airfoil polars and the CFD velocity vector extracted at the line position, as the self-induction at the very centre of the bound vortex should, following vortex theory, be nil. Yet, this is not the case in CFD, which leads to errors in the angle-of-attack computation. We derive an expression for the error in the lift force from vortex considerations and show it to be a function of chord, the smearing length scale used in distributing the AL forces over the numerical domain and the number of grid cells per smearing length scale. Thereby demonstrating that the required number of grid cells - contrary to current belief - needs to grow faster than the inverse of the smearing length scale refinement to maintain the error level. We additionally show that the error can be large for the commonly used ratio of 2 grid cells per length scale, especially if the latter is relatively small with respect to the rotor radius. Ultimately, the recommendation is to always run with the largest smearing length scale possible for the specific application in conjunction with a smearing correction, as this minimizes the error in the blade forces whilst reducing the computational resources required.


## 1. Introduction

The actuator line (AL) [1] is a lifting line (LL) representation of aerodynamic surfaces in Eulerian computational fluid dynamics (CFD) applications. It allows simulating the aerodynamic behaviour of wind turbines by capturing all major flow features of fully resolved rotors, at a fraction of the computational cost. However transferring a LL into the CFD domain requires dispersing the concentrated blade forces of the LL over a certain region - most commonly in form of a Gaussian projection - to avoid causing numerical instabilities. Smearing the forces creates a Lamb-Oseen [2,3] viscous core [4,5,6], which reduces the velocity gradients of the vortex forming about the AL and thus allows representing the velocity field without instabilities in the CFD domain. The AL blade forces are computed from 2D airfoil polars and the CFD velocity vector extracted at the line position, as the self-induction at the very centre of the bound vortex should, following vortex theory, be nil. However, Shives and Crawford [7] found that in CFD simulations it is in fact always non-zero - its magnitude strongly depending on grid resolution. There are some methods trying to avoid relying solely on the AL centre velocity (cf. [8]) by extracting velocities at several locations around the AL core, however those rely on steady-state assumptions which are violated in realistic atmospheric flows. Without sufficient resolution, the large velocity gradient at the vortex centre cannot be accurately captured, leading to non-zero
velocities and hence errors in the angles-of-attack. When accepting an angle-of-attack below $0.5^{\circ}$, Shives and Crawford find the necessary ratio between the force smearing length scale, $\epsilon$, and the grid size, $\Delta x$, to exceed 4 for their particular case. Martìnez-Tossas et al. [9] similarly investigated the interdependency of grid resolution and smearing length scale but unlike Shives and Crawford for rotors instead of wings. Although the focus of their study lies on the wake deficit, they clearly indicate an increase in power with grid refinement at constant $\epsilon$ and a drop in power across all grid sizes when reducing $\epsilon$. Interestingly, the power output shows no convergence with $\epsilon$ but simply keeps falling with it. The smearing correction by Meyer Forsting et al. [ $5,10,11$ ] corrects for the missing induction caused by the formation of the viscous core in the trailed vorticity for ALs, yet even with this correction, there remained a tendency for the power (and thrust) to fall with $\epsilon[12]$. On the left in figure 1 the reduction in the tangential forces for the NREL 5MW [13] at $8 \mathrm{~m} / \mathrm{s}$ is reproduced from this investigation where $n_{\epsilon}=R / \epsilon(R=63 \mathrm{~m})$. Additionally, the ratio between length scale and grid size, $\epsilon / \Delta x$, is given for reference. From smallest $n_{\epsilon}$ to the largest, the power drops by roughly $6 \%$ and the thrust by $3 \%$. The right plot shows the corresponding difference in the angle-of-attack, $\alpha$, with respect to $n_{\epsilon}=16$. The change in $\alpha$ lies below the level deemed acceptable by Shives and Crawford, however in the context of power and load evaluation for rotors it is clearly not. Furthermore there remains large uncertainty as to which combination of $n_{\epsilon}$ and $\epsilon / \Delta x$ will yield accurate rotor loads. This

![](https://cdn.mathpix.com/cropped/f88616bd-9680-438f-a8a6-d8213c97b2ad-04.jpg?height=696&width=1578&top_left_y=1290&top_left_x=239)
Figure 1. The influence of smearing length scale, $\epsilon$, on the NREL 5 MW tangential blade forces (left) and the corresponding difference in the angle-of-attack, $\alpha$, (right) with respect to $n_{\epsilon}=16$.

work expands the investigation of Shives and Crawford [7], giving a more general description of the grid requirements with respect to the Gaussian force smearing length scale. For this purpose the self-induction errors of an AL representing an infinite wing at different circulation strengths and inflow velocities are related to vortex properties. A general formulation for estimating the error in the rotor loads with respect to $n_{\epsilon}$ and $\epsilon / \Delta x$ is established.

## 2. Methodology

The finite-volume solver, EllipSys3D, discretises the Navier-Stokes equations over a blockstructured domain [14]. The turbulence is modelled by a Reynolds-averaged Navier-Stokes formulation with a Menter $k-\omega$ shear-stress transport closure. The convective terms are discretised by the QUICK [15] scheme. The body forces of the AL are either directly applied in
the momentum equations or applied as pressure jumps (PJ) [16] in a Rhie-Chow like approach [17].

The blade length of the NREL 5 MW rotor acted as reference (yet the absolute value is of no significance) and ranges $n_{\epsilon}=[8,64]$ and $\epsilon / \Delta x=[2,16]$ were explored. The force smearing was 2D. The inflow velocity, $V_{0}$, was set to $10 \mathrm{~m} / \mathrm{s}$ and the air density, $\rho_{0}$, to $1 \mathrm{~kg} / \mathrm{m}^{3}$ unless otherwise stated. With a dynamic viscosity of $10^{-5} \mathrm{~kg} / \mathrm{m} / \mathrm{s}$ the rotor-based Reynolds number ${ }^{1}$ is above the recommended $10^{5}$ [18] and the flow is essentially inviscid. The AL force was applied perpendicular to the inflow direction. The reference circulation, $\Gamma_{0}$, about the AL was $10 \mathrm{~m}^{2} / \mathrm{s}$. Note that the lift force on the AL is directly related to the circulation i.e. $L=\rho_{0} \Gamma V_{0}$. Velocities are extracted at the AL position and the steady-state computations are converged to $10^{-12}$. The AL was generally placed on a mesh node, but the effect of placing it at the cell-centre (CC) was also tested.

The numerical domain follows a box topology with a uniformly spaced refined mesh region in the centre with side lengths $6 \epsilon$. The lateral boundary conditions are periodic, thus creating a quasi 2D domain. The other boundaries off the main flow direction are of the symmetry type, whereas the inflow and outflow faces obey Dirichlet and Neumann conditions, respectively. As the domain boundaries seemed to influence the results significantly, they were placed minimally 10 km from the AL. The box domain was thus $20 \mathrm{~km} \times 20 \mathrm{~km} \times 0.3 \mathrm{~km}$ or equivalently $2540 \epsilon_{\max } \times 2540 \epsilon_{\max } \times 38 \epsilon_{\max }$.

## 3. Results and Discussion

### 3.1. Velocity profiles

In figure 2 the streamwise variation in the vertical velocity through the AL core for $n_{\epsilon}=8$ (left) and 16 (right) and its dependency on grid resolution/force application is shown. In the upper row the absolute velocity is compared to the theoretically correct Lamb-Oseen velocity profile (note that for $x / \epsilon \geq 0, v \leq 0$ ) and in the lower the velocity difference with respect to the Lamb-Oseen solution is shown. The behaviour with respect to grid resolution is similar whether $n_{\epsilon}=8$ or 16 , it is first and foremost the magnitude that changes. Yet the velocity difference furthermore reveals that at $n_{\epsilon}=16$ minor instabilities form downstream of the AL for $\epsilon / \Delta x=4$. Whereas the difference in velocities with active PJ is nearly symmetric about the AL, it is asymmetric without it. Most importantly for the evaluation of the angle-of-attack the velocity at the AL is always non-zero and negative. At the same grid resolution the velocity is smaller with active PJ. These trends remain for increasing $n_{\epsilon}$. Increasingly negative velocities at the AL also signify a continues reduction in the angle-of-attack, thereby demonstrating the underlying mechanism leading to reducing blade forces with growing $n_{\epsilon}$.

## 3.2. $A L$ non-zero self-induction error

This section establishes a relationship between grid resolution and smearing length scale. Figure 3 demonstrates the variation in the velocity error at the AL position with grid resolution and smearing length scale and some additional parameters. Velocities are normalised by circulation as this determines its magnitude, as demonstrated by cases with $\Gamma=10 \Gamma_{0}$ falling on the same lines. Clearly the larger $n_{\epsilon}$ - or equivalently the smaller $\epsilon$ and more acute the AL forces - the larger the error. With the Rhie-Chow like PJ algorithm active, the error can be significantly reduced, as already reported by Shives and Crawford [7]. The error is insensitive to the location of the AL, as nearly no change is seen in the error when the AL is moved from the cell node to the cell centre (CC). Some additional sensitivities are explored in figure 4. As these results are obtained for a smaller domain ( $200 \mathrm{~m} \times 200 \mathrm{~m} \times 300 \mathrm{~m}$ ) the error is not converging towards zero,

[^0]![](https://cdn.mathpix.com/cropped/f88616bd-9680-438f-a8a6-d8213c97b2ad-06.jpg?height=1457&width=1593&top_left_y=392&top_left_x=230)
Figure 2. Streamwise variation in vertical velocity through AL core for $n_{\epsilon}=8$ (left) and 16 (right) and its dependency on grid resolution/force application. In the upper row the absolute velocity is compared to the Lamb-Oseen velocity profile (note that for $x / \epsilon \geq 0, v \leq 0$ ) and in the lower the velocity difference with respect to the Lamb-Oseen solution is given.

however the trends are equivalent to those obtained for the larger domain. It reinforces that the magnitude of the error is determined by $\Gamma$ alone and that lower order schemes lead to larger errors, as expected.

To obtain a more general estimate of the velocity error an adequate normalisation is necessary. As the error varies with circulation and smearing length scale, the drop in vorticity over one grid cell could be used as an error measure. The drop in vorticity as function of the radial coordinate $r$ is given by:

$$
\begin{equation*}
\Delta \omega=\frac{\Gamma}{\pi \epsilon^{2}}\left[1-\exp \left(-r^{2} / \epsilon^{2}\right)\right] \tag{1}
\end{equation*}
$$

![](https://cdn.mathpix.com/cropped/f88616bd-9680-438f-a8a6-d8213c97b2ad-07.jpg?height=668&width=780&top_left_y=431&top_left_x=242)
Figure 3. Convergence of the self-induction error at the AL position with grid resolution at different smearing length scales.

![](https://cdn.mathpix.com/cropped/f88616bd-9680-438f-a8a6-d8213c97b2ad-07.jpg?height=675&width=771&top_left_y=424&top_left_x=1055)
Figure 4. Sensitivity of self-induction error. UDS: upwind differencing scheme; CDS: central scheme. Note these results were obtained for a smaller numerical domain and thus have a constant offset.

Integrating from the AL centre over one grid cell $\Delta x$

$$
\begin{equation*}
\int_{0}^{\Delta x} \Delta \omega / \Gamma d r=\frac{2 \Delta x-\sqrt{\pi} \epsilon \operatorname{erf}(\Delta x / \epsilon)}{2 \pi \epsilon^{2}}=K \tag{2}
\end{equation*}
$$

where $\operatorname{erf}()$ is the error function. Multiplying both $\Delta v / \Gamma$ and $K$ by $\epsilon$ yields two dimensionless quantities. The previous results are presented by those measures in figure 5, which shows that they now collapse for different smearing length scales and exhibit a linear relationship using log scales. With PJ the error is generally lower and the error convergence is about first order. Those normalised measures can be rewritten in terms of more tangible quantities. The velocity error can be reformulated in terms of an angle-of-attack error assuming small angles and $\Gamma$ replaced by lift:

$$
\begin{equation*}
\frac{\Delta v \epsilon}{\Gamma}=\frac{\Delta \alpha V_{\infty} \epsilon}{\Gamma}=\frac{2 \Delta \alpha}{C_{l}} \frac{\epsilon}{c} \tag{3}
\end{equation*}
$$

Here $C_{l}$ is the sectional lift coefficient and $c$ the chord. The latter can be replaced by the normalised chord i.e. $\tilde{c}=c / R$ and introducing $n_{\epsilon}$ :

$$
\begin{equation*}
\frac{\Delta v \epsilon}{\Gamma}=\frac{2}{C_{l}} \frac{\Delta \alpha}{\tilde{c} n_{\epsilon}} \tag{4}
\end{equation*}
$$

The integral measure in equation 2 can also be simplified by introducing an approximation to the error function:

$$
\begin{equation*}
\operatorname{erf}(x) \approx \frac{2}{\sqrt{\pi}}\left[x-\frac{x^{3}}{3}+\mathcal{O}\left(x^{5}\right)\right] \tag{5}
\end{equation*}
$$

This approximation is sufficiently accurate for the commonly used values of $\epsilon / \Delta x$ exceeding 2 .

![](https://cdn.mathpix.com/cropped/f88616bd-9680-438f-a8a6-d8213c97b2ad-08.jpg?height=674&width=787&top_left_y=434&top_left_x=237)
Figure 5. Normalised self-induction error.

![](https://cdn.mathpix.com/cropped/f88616bd-9680-438f-a8a6-d8213c97b2ad-08.jpg?height=666&width=780&top_left_y=438&top_left_x=1037)
Figure 6. Correlation between grid refinement and angle-of-attack error, and trend lines fitted to data with and without active PJ force allocation.

Introducing this simplification into equation 2 and multiplying by $\epsilon$, we obtain:

$$
\begin{align*}
K \epsilon & =\frac{1}{2 \pi \epsilon}\left[\frac{2 \epsilon}{3}\left(\frac{\Delta x}{\epsilon}\right)^{3}\right]  \tag{6}\\
& =\frac{1}{3 \pi}\left(\frac{\Delta x}{\epsilon}\right)^{3} \tag{7}
\end{align*}
$$

From the previously established relationship between both quantities:

$$
\begin{equation*}
\log _{10}\left[\frac{\Delta \alpha}{\tilde{c} n_{\epsilon} C_{l}}\right] \propto \log _{10}\left[(\Delta x / \epsilon)^{3}\right] \tag{8}
\end{equation*}
$$

A linear fit to this relationship takes the form:

$$
\begin{equation*}
\log _{10}\left[\frac{\Delta \alpha}{\tilde{c} n_{\epsilon} C_{l}}\right]=a_{1} \log _{10}\left[(\Delta x / \epsilon)^{3}\right]+a_{2} \tag{9}
\end{equation*}
$$

,which allows to establish an equation for estimating the angle-of-attack error:

$$
\begin{equation*}
\Delta \alpha=a_{3} \tilde{c} n_{\epsilon} C_{l}(\epsilon / \Delta x)^{-3 a_{1}} \tag{10}
\end{equation*}
$$

,where $a_{\text {. }}$ represent constants. Figure 6 shows the relationship between those alternative measures and the corresponding fits. The error in the lift force furthermore can be established from the expression for $\Delta \alpha$ in equation 10 :

$$
\begin{equation*}
\tilde{L}=\frac{\Delta L}{L}=\frac{\Delta C_{l}}{C_{l}} \approx \frac{2 \pi \Delta \alpha}{C_{l}}=a_{4} \tilde{c} n_{\epsilon}(\epsilon / \Delta x)^{-3 a_{1}} \tag{11}
\end{equation*}
$$

This leads to an interesting conclusion that the grid resolution does not simply grow proportionally with $n_{\epsilon}$, in fact the grid resolution needs to increase more quickly to keep $\tilde{L}$ constant:

$$
\begin{align*}
n_{\epsilon_{2}}\left(\epsilon_{2} / \Delta x_{2}\right)^{-3 a_{1}} & =n_{\epsilon_{1}}\left(\epsilon_{1} / \Delta x_{1}\right)^{-3 a_{1}}  \tag{12}\\
\Delta x_{2} & =\Delta x_{1}\left(\epsilon_{2} / \epsilon_{1}\right)^{\frac{1+3 a_{1}}{3 a_{1}}} \approx \Delta x_{1}\left(\epsilon_{2} / \epsilon_{1}\right)^{4 / 3} \tag{13}
\end{align*}
$$

Thus halving the smearing length scale i.e. $\epsilon_{2}=0.5 \epsilon_{1}$ requires $\Delta x_{2}=0.4 \Delta x_{1}$.
For completeness all relevant constants are summarised in table 1. While those constants might vary with the specific numerical method employed, they are easily established by extracting the self-induction error at different $\epsilon / \Delta x$ ratios and following the procedure laid out in this section.

Table 1. Induction error fit constants.
|  | $a_{1}$ | $a_{2}$ | $a_{4}$ |
| ---: | :--- | :--- | :--- |
| w/o PJ | $9.700 \mathrm{e}-1$ | -1.259 | $3.465 \mathrm{e}-1$ |
| PJ | $8.669 \mathrm{e}-1$ | -2.047 | $5.637 \mathrm{e}-2$ |


### 3.3. Estimating the error in blade forces for rotors

The error in the lift is a function of the normalised chord size. In figure 7 the normalised chord distributions are shown for five very different rotors, from model scale (MEXICO, $R=2.25 \mathrm{~m}$ ) to futuristic designs with radii beyond 100 m (AVATAR, $R=103 \mathrm{~m}$ ). Only chord lengths for aerofoils with thickness ratios below $40 \%$ are deemed relevant, as even thicker aerofoils are less sensitive to angle-of-attack changes and it is the outer part of the blades that are more heavily loaded. The dot indicates the rotor-averaged chord length. Even though the turbines vary extremely in size, their normalised chord ratios are similar - their averages lying in-between $3.5 \%$ to $5.8 \%$. Following the expression in equation 11 contours of the relative lift error as function of $\epsilon / \Delta x$ and normalised chord length are given in figure 8 either with inactive PJ (left) or active PJ algorithm (right). As reference the rotor-averaged chord lengths of the five turbines are marked as well. Multiplying the contour value by $n_{\epsilon}$ yields the final error. Taking the AVATAR average chord for instance, $n_{\epsilon}=10$ and the commonly used ratio for $\epsilon / \Delta x=2$, the error is $1.61 \%(\mathrm{w} / \mathrm{o} \mathrm{PJ})$ and $0.32 \%(\mathrm{PJ})$, respectively, and for the MEXICO case $2.66 \%$ and $0.54 \%$. It is easily seen that increasing $n_{\epsilon}$, whilst beneficial for resolving the near-wake rotor wake and tip vortices, renders the blade forces inaccurate if $\epsilon / \Delta x$ stays constant.

## 4. Conclusion

The well-known over-prediction of forces by actuator lines (AL), especially towards the ends of lifting-surfaces, has been successfully addressed by recent publications and attributed to missing vorticity in the AL wake.

Yet another, commonly overlooked issue of the AL method lies in the assumption of negligible self-induction at its very centre, at which velocities are extracted to determine the angle-ofattack. We show that, though small in magnitude, the self-induction is sufficiently large to have an detrimental impact on the estimated power of multi-MW turbines through small changes in the angle-of-attack. By correlating the self-induction error with the vorticity drop across a

![](https://cdn.mathpix.com/cropped/f88616bd-9680-438f-a8a6-d8213c97b2ad-10.jpg?height=680&width=759&top_left_y=424&top_left_x=646)
Figure 7. Normalised chord distribution for five different rotors, in regions where the relative aerofoil thickness lies below $40 \%$. Dots indicate the rotor-average.

![](https://cdn.mathpix.com/cropped/f88616bd-9680-438f-a8a6-d8213c97b2ad-10.jpg?height=732&width=1593&top_left_y=1254&top_left_x=242)
Figure 8. Contours of relative lift error as function of grid cells per smearing length scale and normalised chord length, without PJ (left) and with active PJ algorithm (right). Lines represent the rotor-averaged chord lengths from figure 7.

grid cell - as given by the AL force distribution - we establish that the error in the lift force is a function of chord, the smearing length scale used in distributing the AL forces over the numerical domain and the number of grid cells per smearing length scale:

$$
\frac{\Delta L}{L} \propto \frac{c}{\epsilon}\left(\frac{\epsilon}{\Delta x}\right)^{-3}
$$

Thus the error - contrary to common belief - does not remain constant as long as the ratio of smearing length scale to grid size equally remains constant, demonstrating that the frequently cited ratio of 2 grid cells per length scale is not a sufficient condition for minimizing the AL blade force error.

Ultimately, the recommendation is to always run with the largest smearing length scale possible for the specific application in conjunction with a smearing correction [10], as this minimizes the error in the blade forces whilst reducing the computational resources required and applying the blade forces inside the CFD domain in the form of pressure jumps, when using an incompressible solver with collocated grids.

## References

[1] Sørensen J N and Shen W Z 2002 Journal of Fluids Engineering 124 393-399
[2] Lamb H 1932 Hydrodynamics (C.U.P)
[3] Oseen C 1911 Über Wirbelbewegung in einer reibenden Flüssigkeit Arkiv för matematik, astronomi och fysik (Almqvist \& Wiksells)
[4] Dağ K O and Sørensen J N 2020 Wind Energy 23 148-160
[5] Meyer Forsting A R, Pirrung G R and Ramos-García N 2019 Wind Energy Science 4 369-383
[6] Martínez-Tossas L A and Meneveau C 2019 Journal of Fluid Mechanics 863 269-292
[7] Shives M and Crawford C 2013 Wind Energy 16 657-669 ISSN 1099-1824
[8] Forsythe J R, Lynch E, Polsky S and Spalart P 2015 53rd AIAA Aerospace Sciences Meeting
[9] Martínez-Tossas L A, Churchfield M J and Leonardi S 2015 Wind Energy 18 1047-1060
[10] Meyer Forsting A R, Pirrung G R and Ramos-García N 2019 Actuator-line-smearing-correction
[11] Meyer Forsting A, Pirrung G and Ramos García N 2020 Wind Energy Science ISSN 2366-7443
[12] Meyer Forsting A, Pirrung G and Ramos García N 2019 Journal of Physics: Conference Series 1256012020
[13] Jonkerman J, Butterfield S, Musial W and Scott G 2009 Definition of a $5-\mathrm{mw}$ reference wind turbine for offshore system development Tech. rep. NREL
[14] Sørensen N 1995 General purpose flow solver applied to flow over hills Ph.D. thesis Ris $\varnothing$ National Laboratory
[15] Leonard B 1979 Computer Methods in Applied Mechanics and Engineering 19 59-98
[16] Troldborg N, Sørensen N, Réthoré P E and van der Laan M 2015 Computers \& Fluids 119 197-203
[17] Rhie C M and Chow W L 1983 AIAA Journal
[18] Troldborg N, Sørensen J and Mikkelsen R 2009 Actuator Line Modeling of Wind Turbine Wakes Ph.D. thesis Technical University of Denmark


[^0]:    ${ }^{1}$ An artificial reference quantity in this 2D setup, but allows comparing it to the usual 3D applications of the AL.

