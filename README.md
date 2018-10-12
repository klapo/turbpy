# turbpy
## Offline land-surface turbulence schemes
Package designed to allow simulation of surface turbulent fluxes outside of a land-model context. See the jupyter notebook [TF.turbpyDemonstration](./TF.turbpyDemonstration.ipynb) for an example of how the turbpy package can be used. The example plots are generated in this notebook.

![Example of the stability functions available in the turbpy package.](https://github.com/klapo/turbpy/blob/master/turbpy.idealized.jpg)

# Available conductance schemes

See the demonstration notebook for details (above).

## Bulk methods
**1) Standard/Anderson/Choudhury and Monteith (1988)/Jordan (1991)**

A bulk aerodynamic approximation to the Monin-Obukhov based stability correction of Webb (1970). These expressions are identical except for using the Obukhov length or the bulk Richardson number (below). Sets $Q_h$ to zero above a critical Richardson number. This number does not exist in nature. This method is likely the most common method implemented within land models.

![Demonstration of the Webb/log-linear methods being approximately equal.](https://github.com/klapo/turbpy/blob/master/turbpy.idealized_loglinear.jpg)

**2) Louis (unmodified)**

Includes a tunable parameter, b', and the option for capping the conductance.

**3) Mahrt**

Attempts to represent the role of spatial variability in the bulk Richardson number through an empirical parameterization.

## Monin-Obukhov methods

**1) Holtslag and de Bruin**

**2)	Beljaars and Holtslag**

**3) Cheng and Brutsaert 2005**

**4) Log-linear/Webb/Paulson** 

The log-linear method appears under a variety of names, but with identical form. The Anderson/standard method is assumed to approximate these MO stability functions. I explore this assumption a bit at the end of the notebook.

### Capping

If the effective conductance is below a set amount (2, but note this is not a physically meaningful value and could be adjusted) the coductance is set to a constant. This can be implemented with any of the Monin-Obukhov schemes.
