import numpy as np

import multiConst as mc
from aStability import aStability
from surfFluxCalc import moninObukhov


def turbFluxes(
        airTemp,                            # air temperature at some height above the surface (K)
        airPres,                            # air pressure of the air above the vegetation canopy (Pa)
        airVaporPress,                      # vapor pressure of the air above the vegetation canopy (Pa)
        windspd,                            # wind speed above the canopy (m s-1)
        sfcTemp,                            # ground temperature (K)
        sfcVaporPress,                      # Vapor pressure at the surface (Pa)
        snowDepth,                          # depth of snow (m)
        mHeight,                            # height of observations (m)
        groundSnowFraction=1,               # Fraction of ground covered by snow (-)
        ixDerivMethod=False,                 # choice of method used to compute derivative (analytical or numerical)
        ixStability='mahrtExponential',     # method for calculating stability
        ixStabParam=mc.stabParams,          # Stability params from mc library
        z0Ground=.005                       # Surface roughness length for log-layer (m)
        ):
    '''
    turbFluxes.py

    Offline turbulent fluxes for non-canopy surfaces. Python companion code for
    the SUMMA model. https://github.com/NCAR/summa

    INPUT:
        airTemp (K), airPres (Pa), airVaporPress (Pa), windspd (m s-1), sfcTemp (K),
        soilRelHumidity [0-1], snowDepth (m)
            - scalars of turbulence boundary conditions
        mHeight (m)
            - height of atmospheric observations. All observations assumed to come
             from same height at the moment.
        ixDerivMethod (analytical, numerical, or False)
            - model control on how derivatives are computeDerivative
        ixStability (see multiConst.py for options)
            - stability scheme method (string)
        ixStabParam (see multiConst.py for default values)
            - Parameter values for the stability scheme.
        z0Ground (m)
            - surface roughness length

    OUTPUT:
        groundConductanceSH - ground conductance for sensible heat (m s-1)
        groundConductanceLH - ground conductance for latent heat (m s-1)
        senHeatGround - sensible heat flux from ground surface (W m-2)
        latHeatGround - latent heat flux from ground surface (W m-2)
        turbFluxGround - net turbulent heat fluxes at the ground surface (W m-2)
        dTurbFluxGround_dTGround - derivative in net turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
    '''
# --------------------------------------------------------------------------------------------------------------------
    if np.isnan(airTemp) or np.isnan(airPres) or np.isnan(airVaporPress) or np.isnan(sfcTemp):
        raise ValueError('Input data includes nan value.')
    # check that measurement height above the ground surface is above the roughness length
    if mHeight < snowDepth + z0Ground:
        raise ValueError('measurement height < snow depth + roughness length')
    # define height above the snow surface
    heightAboveGround = mHeight - snowDepth

    ########
    # Local variables
    volHeatCapacityAir = mc.iden_air * mc.Cp_air              # volumetric heat capacity of air (J m-3)
    latentHeatConstant = mc.iden_air * mc.w_ratio / airPres   # latent heat constant for (kg m-3 Pa-1)

    #######
    # Latent Heat - Vaporization or sublimation.
    # NOTE: The physics implied here may be wrong. Point to bring up with Martyn/Bart
    if groundSnowFraction > 0. and sfcTemp < mc.Tfreeze:
        latHeatSubVapGround = mc.LH_sub  # sublimation from snow
    # case when the ground is snow-free
    # evaporation of water in the soil pores, this occurs even if frozen because of super-cooled water
    elif groundSnowFraction == 0.:
        latHeatSubVapGround = mc.LH_vap
    else:
        latHeatSubVapGround = mc.LH_sub

    ########
    # compute resistances
    if ixDerivMethod:
        derivDesired = ('analytical' in ixDerivMethod
                        or 'numerical' in ixDerivMethod)
    else:
        derivDesired = False

    stabOut = aStability(
        derivDesired,               # flag for derivatives
        ixStability,                # choice of stability function
        ixStabParam,                # Parameters dictionary
        heightAboveGround,          # measurement height (m)
        airTemp,                    # air temperature @ mHeight (K)
        airVaporPress,              # vapor pressure @ mHeight (Pa)
        sfcTemp,                    # ground temperature (K)
        sfcVaporPress,              # vapor pressure @ surface (Pa)
        windspd,                    # wind speed at @ mHeight (m s-1)
        z0Ground,                   # surface roughness length (m)
        )

    # Unpack conductances
    (RiBulkGround, stabilityCorrectionParameters,
        stabilityCorrectionDerivatives, conductanceSensible,
        conductanceLatent) = stabOut

    ########
    # compute sensible and latent heat fluxes (positive downwards)
    # Turbulent fluxes using bulk aerodynamic stability corrections.
    senHeatGround = -volHeatCapacityAir * conductanceSensible * \
        (sfcTemp - airTemp)
    latHeatGround = -latHeatSubVapGround * latentHeatConstant * \
        conductanceLatent * (sfcVaporPress - airVaporPress)
    if 'moninObukhov' in ixStability:
        # Turbulent fluxes for Monin-Obukhov similarity theory
        senHeatGround, latHeatGround = moninObukhov(airTemp, airVaporPress,
                                                    sfcTemp, sfcVaporPress, stabilityCorrectionParameters,
                                                    senHeatGround, latHeatGround,
                                                    conductanceSensible, conductanceLatent)

    # compute derivatives
    # if ixDerivMethod == 'analytical' and not ixStability == 'moninObukhov':
    # compute derivatives for the ground fluxes w.r.t. ground temperature
    # d(ground sensible heat flux)/d(ground temp)
    # dSenHeatGround_dTGround = (-volHeatCapacityAir *
    #     dGroundCondSH_dGroundTemp) *\
    #     (sfcTemp - airTemp) + \
    #     (-volHeatCapacityAir * groundConductanceSH)
    # d(ground latent heat flux)/d(ground temp)
    # dLatHeatGround_dTGround = (-latHeatSubVapGround * latentHeatConstant *
    #     dGroundCondLH_dGroundTemp) * \
    #     (satVP_GroundTemp * soilRelHumidity - airVaporPress) + \
    #     (-latHeatSubVapGround * latentHeatConstant * groundConductanceLH) *\
    #     dSVPGround_dGroundTemp * soilRelHumidity

    ########
    # net turbulent flux at the ground surface (W m-2)
    turbFluxGround = senHeatGround + latHeatGround

    # compute derivatives
    # if ixDerivMethod == 'analytical':
    #     # (energy derivatives)
    #     # derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
    #     dTurbFluxGround_dTGround = dSenHeatGround_dTGround + \
    #         dLatHeatGround_dTGround
    # else:
    #     dTurbFluxGround_dTGround = -9999

    return (
        conductanceSensible,  # ground conductance for sensible heat (m s-1)
        conductanceLatent,  # ground conductance for latent heat (m s-1)
        senHeatGround,  # sensible heat flux from ground surface (W m-2)
        latHeatGround,  # latent heat flux from ground surface (W m-2)
        turbFluxGround,  # net turbulent heat fluxes at the ground surface (W m-2)
        stabilityCorrectionParameters,  # Stability correction (0-1)
        stabilityCorrectionDerivatives,
        np.nan,  # derivative in net turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
        )
###########################################################################################################
