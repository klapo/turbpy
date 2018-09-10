import numpy as np

import turbpy.multiConst as mc
from turbpy.aStability import aStability
from turbpy.parameter_methods import get_params
from turbpy import surfFluxCalc


def turbFluxes(airTemp,  # air temperature at some height above the surface (K)
               airPres,  # air pressure of the air above the vegetation canopy (Pa)
               airVaporPress,  # vapor pressure of the air above the vegetation canopy (Pa)
               windspd,  # wind speed above the canopy (m s-1)
               sfcTemp,  # ground temperature (K)
               sfcVaporPress,  # Vapor pressure at the surface (Pa)
               snowDepth,  # depth of snow (m)
               mHeight,  # height of observations (m)
               groundSnowFraction=0,  # Fraction of surface covered by snow (-). Used in latent heat calculation.
               param_dict=None,  # dictionary of parameters and parameterizations
               z0Ground=.005,  # Surface roughness length for log-layer (m)
               ):
    '''
    turbFluxes.py

    Offline turbulent fluxes for non-canopy, relatively smooth surfaces (e.g., snow)

    INPUT:
        airTemp (K), airPres (Pa), airVaporPress (Pa), windspd (m s-1), sfcTemp (K),
        soilRelHumidity [0-1], snowDepth (m)
            - scalars of turbulence boundary conditions
        mHeight (m)
            - height of atmospheric observations. All observations assumed to come
             from same height at the moment.
        param_dict
            - Dictionary containing the choices for simulating turbulence.
        z0Ground (m)
            - surface roughness length

    OUTPUT:
        groundConductanceSH - ground conductance for sensible heat (m s-1)
        groundConductanceLH - ground conductance for latent heat (m s-1)
        senHeatGround - sensible heat flux from ground surface (W m-2)
        latHeatGround - latent heat flux from ground surface (W m-2)
        turbFluxGround - net turbulent heat fluxes at the ground surface (W m-2)
    '''

    # Check for the forcing data
    if np.isnan(airTemp) or np.isnan(airPres)\
       or np.isnan(airVaporPress) or np.isnan(sfcTemp):
        raise ValueError('Input data includes nan value.')
    # Check that measurement height above the ground surface is above
    # the roughness length.
    if mHeight < snowDepth + z0Ground:
        raise ValueError('measurement height < snow depth + roughness length')
    # Define height above the snow surface.
    heightAboveGround = mHeight - snowDepth

    # Get the parameters for running turbpy.
    param_dict = get_params(param_dict)

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
    # evaporation of water in the soil pores, this occurs even if frozen
    # because of super-cooled water
    elif groundSnowFraction == 0.:
        latHeatSubVapGround = mc.LH_vap
    else:
        latHeatSubVapGround = mc.LH_sub

    ########
    # compute resistances
    stabOut = aStability(
        param_dict,                # parameter defining dictionary
        heightAboveGround,         # measurement height (m)
        airTemp,                   # air temperature @ mHeight (K)
        airVaporPress,             # vapor pressure @ mHeight (Pa)
        sfcTemp,                   # ground temperature (K)
        sfcVaporPress,             # vapor pressure @ surface (Pa)
        windspd,                   # wind speed at @ mHeight (m s-1)
        z0Ground,                  # surface roughness length (m)
    )

    # Unpack conductances
    (RiBulkGround, stabilityCorrectionParameters,
     conductanceSensible, conductanceLatent) = stabOut

    ########
    # compute sensible and latent heat fluxes (positive downwards)
    # Turbulent fluxes using bulk aerodynamic stability corrections.
    if param_dict['stability_method'] == 'monin_obukhov':
        # MO uses a windless exchange coefficient for stable conditions
        # Pass in the "wind function" to calculate fluxes.
        senHeatGround = volHeatCapacityAir * conductanceSensible
        latHeatGround = latHeatSubVapGround * latentHeatConstant * \
            conductanceLatent

        # Determine capping behavior
        cap = param_dict['monin_obukhov']['capping']

        # Turbulent fluxes for Monin-Obukhov similarity theory
        (senHeatGround, latHeatGround,
         conductanceSensible, conductanceLatent) = \
            surfFluxCalc.moninObukhov(airTemp, airVaporPress,
                                      sfcTemp, sfcVaporPress, windspd,
                                      stabilityCorrectionParameters,
                                      senHeatGround, latHeatGround,
                                      conductanceSensible, conductanceLatent,
                                      volHeatCapacityAir, latHeatSubVapGround,
                                      latentHeatConstant, cap)
    else:
        senHeatGround = (-volHeatCapacityAir * conductanceSensible *
                         (sfcTemp - airTemp))
        latHeatGround = (-latHeatSubVapGround * latentHeatConstant *
                         conductanceLatent * (sfcVaporPress - airVaporPress))

    return (
        conductanceSensible,  # ground conductance for sensible heat (m s-1)
        conductanceLatent,  # ground conductance for latent heat (m s-1)
        senHeatGround,  # sensible heat flux from ground surface (W m-2)
        latHeatGround,  # latent heat flux from ground surface (W m-2)
        stabilityCorrectionParameters,  # Stability correction (0-1)
        param_dict  # Return the parameterizations used in this run.
    )
