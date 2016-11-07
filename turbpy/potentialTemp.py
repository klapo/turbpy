import numpy as np

from turbpy.bulkRichardson import bulkRichardson
import turbpy.multiConst as mc


def potentialTemp(mHeight, airTemp, correctionMethod, airVaporPress=np.nan):
    '''
    Converts absolute temperature to potential temperature using different
    specified assumptions.

    INPUTS:
    mHeight                 measurement height above the surface (m)
    airTemp                 air temperature (K)
    airVaporPress           vapor pressure of air (Pa)
    '''

# ------------------------------------------------------------------------------
# Sub-functions (potential temperature assumptions)
# ------------------------------------------------------------------------------
    def marht2004():
        '''
        Gives a LOCAL potential temperature. Note that this is different than
        the actual potential temperature.

        Mahrt, L., and D. Vickers (2004), Bulk formulation of the surface heat flux,
            Boundary-Layer Meteorol., 110(3), 357–379,
            doi:10.1023/B:BOUN.0000007244.42320.1e.
        '''
        theta = airTemp + .01 * mHeight
        return(theta)

# Adding a non-local potential temperature will be an involved process.
# Using the loca potential temperature for the time being
    # def bolton1980():
    #     # Estimate pressure here
    #     theta = airTemp * (1e5 / p) ** (mc.R_da / mc.Cp_air)
    #
    # ' Stability Error Message '
    # def stabErrMess():
    #     raise ValueError(
    #         'Unrecognized choice for calculating potential temperature: '
    #         + correctionMethod
    #         + '\n'
    #         + 'Valid options: '
    #         + stabilityCase.keys
    #         )
    #
    # def pressIdealGas():
    #     Tv/(1+(1/0.622-1)*rv)
    #     p = rho * R_d * T_v
