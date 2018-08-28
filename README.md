# turbpy
## Offline land-surface turbulence schemes
Package designed to allow simulation of surface turbulent fluxes outside a land-model context.

![Example of the stability functions available in the turbpy package.](https://github.com/klapo/turbpy/blob/master/TF.OfflineTurb.Idealized.jpg)

a) The stability corrections for the bulk aerodynamic methods

b) the conductance for all methods

c) the resulting sensible heat fluxes for each method assuming a surface temperature of 265K, wind speed of 1 m/s,
surface roughness of 0.005, and varying the air temperature from 253K (Ri = -1) to 293K (Ri = 2).

See the jupyter notebook [TF.turbpyDemonstration](./TF.turbpyDemonstration.ipynb) for an example of how the turbpy package can be used. The image above is generated in this notebook.

## Available conductance schemes

# turbpy: simulating sensible heat fluxes

### Bulk methods
**1) Standard/Anderson/Choudhury and Monteith (1988)/Jordan (1991)**

**2) Louis (unmodified)**

**2a) Louis (modified)**

Adjusts the b' variable to 12 from 4.7. 

**2b) Louis (capped)**

Implemented the CROCUS capping. CROCUS allows multiple methods for capping conductance, I chose the one with a publication to cite. In this case, $Ri_b \le 0.026$.

**3) Mahrt**

### Monin-Obukhov methods

This publication describes both options 1 and 2.  
Beljaars, A. C. M., & Holtslag, A. A. M. (1991). Flux Parameterization over Land Surfaces for Atmospheric Models. Journal of Applied Meteorology. https://doi.org/10.1175/1520-0450(1991)030


**1)	Holtslag and de Bruin**

A Monin-Obukhov based stability parameterization. This was the original stability correction used in the study, but its behavior was hidden by the capping. I've modified the code to accept various cappings (described below), allowing me to easily access this stability parameterization along with others.

Note, the original citation for this method is incorrect. It was only described in Launiainen and vihma (1990). The stability correction does not appear in the originally cited paper.

This method has a critical Richardson number ($Ri_c = 1.43$)

**2)	Beljaars and Holtslag**

An adjustmebnt to the Holtslag and de Bruin method. Also a Monin-Obukhov similarity theory parameterization for the stability corrections. Makes an important assumptions that $\Psi_H \neq \Psi_M$, which required some further code refactoring.

The expression for $\Psi_H$ (stability correction for heat) comes from BH91, equation 32.

$$
\Psi_H = -(1 + \frac{2}{3}a\zeta)^(\frac{3}{2}) - b(\zeta - \frac{c}{d}) e^{-d\zeta} - \frac{bc}{d} + 1 $$

to solve for the Obukhov length iteratively, I need to find $\frac{d\Psi_h}{d\zeta}$. This is the expression I've used.

$$
\frac{d\Psi_h}{d\zeta} = (b e^{-d\zeta}(d\zeta - 1 - c) - a(\frac{2a\zeta}{3} + 1)^{\frac{1}{2}}
$$

**3) Modified MO scheme **

Jiménez, P. A., Dudhia, J., González-Rouco, J. F., Navarro, J., Montávez, J. P., & García-Bustamante, E. (2012). A Revised Scheme for the WRF Surface Layer Formulation. Monthly Weather Review, 140(3), 898–918. https://doi.org/10.1175/MWR-D-11-00056.1

#### Log-linear methods

**5)	Log-linear/Webb/Paulson** 

These are the fifth (also called the classical theory) sixth (Webb, 1970), and seventh (Paulson 1970) names for this method I have found. The Anderson/standard method is assumed to approximate these MO stability functions. I explore this assumption a bit at the end of the notebook.

$$
\Psi_H = \Psi_M = -\alpha\zeta
$$

where alpha is assumed to be a constant with a value of 5. This leads to a critical Richardson number of ~0.2. We know this is non-physical, and likely plays a role in the substantial model ad hocery. I implement this method with $\alpha$ being a tunable parameter.

**5a) Webb as implemented in CLMv4.5**

Requires that the horizontal wind speed remain greater than one ($V \ge 1$) to "stop runaway cooling". 

For $0\le \zeta \le 1 $
$$\phi_h = 1 + \alpha \zeta $$

For $\zeta > 1$
$$\phi_h = \alpha + \zeta $$

The relationship between $\phi$ and $\Psi$ is

$$
\Psi = \int_{0}^{\zeta} \frac{1 - \phi(x)}{x}dx 
$$

So if we integrate to find $\Psi$ for the two cases

For $0\le \zeta \le 1 $
$$\Psi_h = -\alpha \zeta $$

For $\zeta > 1$
$$\Psi_h = (1 - \alpha)\log{\zeta} + \zeta $$

The citation that justifies these equations does not contain an expression for the stability correction or $\Psi$ (from which we can determine the stability correction). I am unsure what assumptions go into these expressions as the expressions do not converge at a value of $\zeta = 1$. This causes unstable behavior near that threshold value and leads to big jumps in the simulated $Q_h$.

**5b) Webb as implemented in NoahMP v1.1 (used in WRFv3.4)**

NoahMP has four options for parameterizing turbulence. Here, we will only examine the Webb method. Other methods include substantial capping/limiting behavior (e.g., the original Chen 1997 flux parameterization). For instance, the NoahMP includes a Paulson scheme, with an identical stability correction to the Webb scheme, but with adjustments to $u*$, $z_0$, and additional capping of certain values like $\zeta$. The Webb scheme is the most straighforward.

One of the options also has identical stability corrections but calls them the Paulson scheme.

$$
\Psi_H = \Psi_M = -\alpha\zeta
$$

with the limit of $\zeta <= 1$

### Capping
**I have catalogued four capping behaviors used in Monin-Obukhov schemes.**

Each capping method is implemented with its intended corresponding stability correction ($L$ and $\zeta$ capping with the Webb correction and the windless exchange coefficient with the Holtslag and de Bruin correction).

**1) $\zeta$ capping**

Implemented in CLMv4.5,  $\zeta$ is capped at a value of 2, whereas in NoahMP-Webb it is capped at 1. This was also implemented in conjunction with the Webb stability correction.

**2) Windless exchange coefficient**

If the effective conductance is below a set amount (2, but note this is not a physically meaningful value) the coductance is set to a constant. This was implemented with the Holtslag and de Bruin stability correction.


NoahMP does this, but in a round about fashion (from line 3291 in module_sf_noahmplsm.F of NoahMPv1.1).

**3) Minimum wind speed**

CLM invokes a minimum wind speed of 1m/s. Note that this does not actually keep the sensible heat fluxes from reaching zero. NoahMP has a minimum $u*$ of 0.1m/s.

** Capping in bulk schemes **

**1) Louis (CROCUS; Martin and Lejuene 1998; Lafarsse et al., 2017)**

Applies a critical Ri threshold above which the conductance is held as a constant. Justifies this using the same general line of reasoning as Mahrt 1987, but with no data. My personal bet was this justification was determined a posteriori when trying to stop runaway cooling and to justify the constant conductance.



