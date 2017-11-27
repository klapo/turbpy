# turbpy
## Offline land-surface turbulence schemes
Package designed to allow simulation of surface turbulent fluxes outside a land-model context.

![Example of the stability functions available in the turbpy package.](https://github.com/klapo/turbpy/blob/master/TF.OfflineTurb.Idealized.png)

a) The stability corrections for the bulk aerodynamic methods

b) the conductance for all methods

c) the resulting sensible heat fluxes for each method assuming a surface temperature of 265K, wind speed of 1 m/s,
surface roughness of 0.005, and varying the air temperature from 253K (Ri = -1) to 293K (Ri = 2).

See the ipython notebook [TF.turbpyTest](./TF.turbpyTest.ipynb) for an example of how the turbpy package can be used. The image above is generated in [TF.turbpyDemonstration](./TF.turbpyDemonstration.ipynb).

## Available conductance schemes
Conductance schemes are broken into two categories.

### 1) Bulk aerodynamic
- Anderson: Assumes turbulence ceases at a critical bulk Richardson number [Anderson 1976]. This scheme can be found under a variety of names including Choudhury and Monteith [1988] or Jordan [1991], but each scheme is functionally identical. We refer to this scheme as the Anderson or "standard" method.  
Anderson, E. A. (1976). A point energy and mass balance model of a snow cover. NOAA Technical Report NWS 19. https://doi.org/10.1016/S0074-6142(99)80039-4  
Choudhury, B. J., & Monteith, J. L. (1988). A four-layer model for the heat budget of homogeneous land surfaces. Quarterly Journal of the Royal Meteorological Society, 114, 373–398. https://doi.org/10.1002/qj.49711448006  
Jordan, R. (1991). A One-Dimensional Temperature Model for a Snow Cover: Technical Documentation for SNTHERM.  

- Louis: Bulk aerodynamic approximation to the Monin-Obukhov similarity theory [Louis 1979].  
Louis, J. F. (1979). A parametric model of vertical eddy fluxes in the atmosphere. Boundary-Layer Meteorology, 17(2), 187–202. https://doi.org/10.1007/BF00117978

- Mahrt: Parameterizes spatial variability in the stability [Mahrt 1987].  
Mahrt, L. (1987). Grid-Averaged Surface Fluxes. Monthly Weather Review. https://doi.org/10.1175/1520-0493(1987)115<1550:GASF>2.0.CO;2  

### 2) Monin-Obukhov Similarity theory
We implement the numerical approximation from SNTHERM revision 4, as described in Andreas 2001 and Jordan and Andreas 2004.

Andreas, E. L. (2001), Parameterizing Scalar Transfer over Snow and Ice: A Review, J. Hydrometeorol., 3, 417–432.  
Andreas, E. L., Jordan, R. E., & Makshtas, A. P. (2004). Simulations of Snow, Ice, and Near-Surface Atmospheric Processes on Ice Station Weddell. Journal of Hydrometeorology, 5(4), 611–624. https://doi.org/10.1175/1525-7541(2004)005<0611:SOSIAN>2.0.CO;2  
Chris Bretheron's [notes on boundary layer meteorology](atmos.washington.edu/~breth/classes/AS547)
