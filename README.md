# turbpy
## Offline land-surface turbulence schemes
Package designed to allow simulation of surface turbulent fluxes outside of a land-model context. See the jupyter notebook [TF.turbpyDemonstration](./TF.turbpyDemonstration.ipynb) for an example of how the turbpy package can be used. The example plots are generated in this notebook.

![Example of the stability functions available in the turbpy package.](https://github.com/klapo/turbpy/blob/master/turbpy.idealized.jpg)

## Available conductance schemes

### Bulk methods
**1) Standard/Anderson/Choudhury and Monteith (1988)/Jordan (1991)**

A bulk aerodynamic approximation to the Monin-Obukhov based stability correction of Webb (1970). These expressions are identical except for using the Obukhov length or the bulk Richardson number (below). Sets $Q_h$ to zero above a critical Richardson number. This number does not exist in nature. This method is likely the most common method implemented within land models.

![Demonstration of the Webb/log-linear methods being approximately equal.](https://github.com/klapo/turbpy/blob/master/turbpy.idealized_loglinear.jpg)

**2) Louis (unmodified)**

$$
b' = 4.7
$$

**2a) Louis (modified)**

Adjusts the b' variable to 12 from 4.7. Some very early parameter testing found that varying the b' value had almost not discernible impact on the mean bias. Based on the plots below, this seems likely to be the case again.

**2b) Louis (capped)**

Implemented the CROCUS capping where C_h is set to a constant value for Ri greater than some threshold, $Ri_{st}$ such that

$$
C_h(Ri_{st} = 0.026)
$$

**3) Mahrt**

Attempts to represent the role of spatial variability in the bulk Richardson number through an empirical parameterization.

## Monin-Obukhov methods

This publication describes both options 1 and 2.  
Beljaars, A. C. M., & Holtslag, A. A. M. (1991). Flux Parameterization over Land Surfaces for Atmospheric Models. Journal of Applied Meteorology. https://doi.org/10.1175/1520-0450(1991)030


**1) Holtslag and de Bruin**

A Monin-Obukhov based stability parameterization. We implement the method as described in Launiainen and vihma (1990)

This method has a critical Richardson number ($Ri_c = 1.43$)

**2)	Beljaars and Holtslag**

An adjustmebnt to the Holtslag and de Bruin method. Includes the assumptions that $\Psi_H \neq \Psi_M$.

The expression for $\Psi_H$ (stability correction for heat) comes from BH91, equation 32.

$$
\Psi_H = -(1 + \frac{2}{3}a\zeta)^(\frac{3}{2}) - b(\zeta - \frac{c}{d}) e^{-d\zeta} - \frac{bc}{d} + 1 $$

to solve for the Obukhov length iteratively, I need to find $\frac{d\Psi_h}{d\zeta}$. This is the expression I've used.

$$
\frac{d\Psi_h}{d\zeta} = (b e^{-d\zeta}(d\zeta - 1 - c) - a(\frac{2a\zeta}{3} + 1)^{\frac{1}{2}}
$$

**4) Cheng and Brutsaert 2005 **



**5) Log-linear/Webb/Paulson** 

The log-linear method appears under a variety of names, but with identical form. The Anderson/standard method is assumed to approximate these MO stability functions. I explore this assumption a bit at the end of the notebook.

$$
\Psi_H = \Psi_M = -\alpha\zeta
$$

where alpha is assumed to be a constant with a value of 5. This leads to a critical Richardson number of ~0.2. We know this is non-physical. We allow $\alpha$ to be a tunable parameter.

In addition, I try to mimic CLMv4.5 and NoahMPv1.1 as closely as possible as described below.

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

The citation that justifies these equations does not contain an expression for the stability correction or $\Psi$ (from which we can determine the stability correction). I am unsure what assumptions go into these expressions as the expressions do not converge at a value of $\zeta = 1$. This causes unstable behavior near that threshold value and leads to big jumps in the simulated $Q_h$. We do not show this method in the example figure, as the instability is likely a result of something we do not know, rather than the behavior actually in the model.

**5b) Webb as implemented in NoahMP v1.1 (used in WRFv3.4)**

We implement the Webb method as found in NoahMPv1.1 since the model is relatively well-documented with descriptive comments and citations. A technical note, while not perfectly matching the code, is also available and useful for understanding some of the implementation.

NoahMPv1.1 has four options for parameterizing turbulence. Here, I will only examine the Webb method. Other methods include substantial capping/limiting behavior (e.g., the original Chen 1997 flux parameterization). For instance, the NoahMP includes a Paulson scheme, with an identical stability correction to the Webb scheme, but with adjustments to $u*$, $z_0$, and additional capping of certain values like $\zeta$. The Webb scheme is the most straighforward.

One of the options also has identical stability corrections but calls them the Paulson scheme.

$$
\Psi_H = \Psi_M = -\alpha\zeta
$$

with the limit of $\zeta <= 1$

### Capping
**I have catalogued four capping behaviors used in Monin-Obukhov schemes.**

**1) Obukhov lenth capping**

L is not allowed to increase beyond some critical value during very stable conditions.

**2) $\zeta$ capping**

Implemented in CLMv4.5 (possibly other model versions as well). $\zeta$ is capped at a value of 2, whereas in NoahMP-Webb it is capped at 1. This was also implemented in conjunction with the Webb stability correction. Again, to reiterate, the Webb model has a critical Richardson number above which is shuts off turbulence, which is non-physical.

**3) Windless exchange coefficient**

If the effective conductance is below a set amount (2, but note this is not a physically meaningful value) the coductance is set to a constant. This was implemented with the Holtslag and de Bruin stability correction, but completely masks the behavior of that method.

Each capping method is implemented with its intended corresponding stability correction ($L$ and $\zeta$ capping with the Webb correction and the windless exchange coefficient with the Holtslag and de Bruin correction).

NoahMP does this, but in a round about fashion (from line 3291 in module_sf_noahmplsm.F of NoahMPv1.1).

> RAHC = MAX(1.,1./(CH*UR))

The translation: the aerodynamic resistance for heat = maximum of (1, conductance for heat * horizontal wind speed).

Let's take an example value for $C_h$ from below, $C_h = 0.002$ with a wind speed of $U = 1.5$ that gives a value of 0.003. So NoahMPv1.1 effectively only uses conductance values greater than 1.

**4) Minimum wind speed**

CLM invokes a minimum wind speed of 1m/s. Note that this does not actually keep the sensible heat fluxes from reaching zero. NoahMP has a minimum $u*$ of 0.1m/s.

**5) Louis (CROCUS; Martin and Lejuene 1998; Lafarsse et al., 2017)**

Applies a critical Ri threshold above which the conductance is held as a constant. Justifies this using the same general line of reasoning as Mahrt 1987, but with no data. My personal bet was this justification was determined a posteriori when trying to stop runaway cooling and to justify the constant conductance.



