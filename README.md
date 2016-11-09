# turbpy
## Offline land-surface turbulence schemes 
Package designed to allow simulation of surface turbulent fluxes outside a land-model context.

![Example of the stability functions available in the turbpy package.](https://github.com/klapo/turbpy/blob/master/TF.OfflineTurb.Idealized.pdf)

a) The stability corrections for the bulk aerodynamic methods

b) the conductance for all methods

c) the resulting sensible heat fluxes for each method assuming a surface temperature of 265K, wind speed of 1 m/s, 
surface roughness of 0.005, and varying the air temperature from 253K (Ri = -1) to 293K (Ri = 2).

See the ipython notebook [TF.turbpyTest](./TF.turbpyTest.ipynb) for an example of how the turbpy package can be used. The image above is generated in [TF.turbpyDemonstration](./TF.turbpyDemonstration.ipynb).

## References
Many of these equations and definitions are derived from Chris Bretherton's notes on boundary layer meteorology, from Andreas 2001, and from SNTHERM revision 4.

Andreas, E. L. (2001), Parameterizing Scalar Transfer over Snow and Ice:
A Review, J. Hydrometeorol., 3, 417–432.

Chris Bretheron's notes on boundary layer meteorology: atmos.washington.edu/~breth/classes/AS547
