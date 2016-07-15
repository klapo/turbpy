import numpy as np
from .multiConst import *
from .conversionTools import *

def turbFluxes(
        # input: atmosphere boundary conditions
        airTemp,                        # air temperature at some height above the surface (K)
        airPres,                        # air pressure of the air above the vegetation canopy (Pa)
        VPair,                          # vapor pressure of the air above the vegetation canopy (Pa)
        windspd,                        # wind speed above the canopy (m s-1)
        # input: surface boundary conditions
        groundTemp,                     # ground temperature (K)
        # input: diagnostic variables
        soilRelHumidity,                # relative humidity in the soil pores [0-1]
        snowDepth,                      # depth of snow (m)
        #
        mHeight,                        # height of observations (m)
        # input: model control
        ixDerivMethod=None,             # choice of method used to compute derivative (analytical or numerical)
        ixStability='mahrtExponential', # method for calculating stability
        ixStabParam=None,               # Variable holding stability scheme parameters (unitless)
        z0Ground=.005                   # Surface roughness length for log-layer (m)
        ):
    '''
    turbFluxes.py

    Offline turbulent fluxes for non-canopy surfaces. Python companion code for
    the SUMMA model. https://github.com/NCAR/summa

    INPUT:
        airTemp (K), airPres (Pa), VPair (Pa), windspd (m s-1), groundTemp (K),
        soilRelHumidity [0-1], snowDepth (m)
            - scalars of turbulence boundary conditions
        mHeight (m)
            - height of atmospheric observations. All observations assumed to come
             from same height at the moment.
        ixDerivMethod (analytical, numerical, or None)
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
    if np.isnan(airTemp) or np.isnan(airPres) or np.isnan(VPair) or np.isnan(groundTemp):
        raise ValueError('Input data includes nan value.')
    # check that measurement height above the ground surface is above the roughness length
    if mHeight < snowDepth+z0Ground:
        raise ValueErrorr('measurement height < snow depth + roughness length')
    # define height above the snow surface
    heightAboveGround  = mHeight - snowDepth

    ########
    # Local variables
    evapSmooth=1.                                  # smoothing parameter for latent heat (W m-2)
    volHeatCapacityAir = iden_air*Cp_air           # volumetric heat capacity of air (J m-3)
    latentHeatConstant = iden_air*w_ratio/airPres  # latent heat constant for (kg m-3 Pa-1)

    ########
    # saturation vapor pressure at the temperature of the ground (Pa)
    (satVP_GroundTemp,empty) = satVapPress(groundTemp - Tfreeze)

    #######
    # Latent Heat - Vaporization or sublimation.
    # NOTE: The physics implied here may be wrong. Point to bring up with Martyn/Bart
    if scalarGroundSnowFraction > 0. and groundTemp < Tfreeze:
        latHeatSubVapGround = LH_sub  # sublimation from snow
    # case when the ground is snow-free
    # evaporation of water in the soil pores, this occurs even if frozen because of super-cooled water
    elif scalarGroundSnowFraction == 0.:
        latHeatSubVapGround = LH_vap
    else:
        latHeatSubVapGround = LH_sub

    ########
    # compute resistances
    derivDesired = 'analytical' in ixDerivMethod or 'numerical' in ixDerivMethod
    stabOut = aStability(
                    derivDesired,               # flag to indicate if analytical derivatives are desired
                    ixStability,                # choice of stability function
                    # input: above-canopy forcing data
                    heightAboveGround,           # measurement height (m)
                    airTemp,                    # air temperature at some height above the surface (K)
                    groundTemp,                 # ground temperature (K)
                    windspd,                    # wind speed at some height above the surface (m s-1)
                    # input: parameters
                    z0Ground,                   # surface roughness length (below canopy/non-vegetated [snow]) (m)
                    )

    # Unpack conductances
    (RiBulkGround,                              # bulk Richardson number (-)
    groundStabilityCorrection,                  # stability correction for turbulent heat fluxes (-)
    conductanceSensible,                        # Conductance parameter for sensible heat exchange
    conductanceLatent,                          # Conductance parameter for latent heat exchange
    _,                                          # derivative in stab. corr. w.r.t. Ri for the ground surface (-)
    _,                                          # derivative in stab. corr. w.r.t. air temperature (K-1)
    _                                           # derivative in stab. corr. w.r.t. sfc temperature (K-1)
    ) = stabOut

    ########
    # compute sensible and latent heat fluxes, and derivatives...
    # (positive downwards)
    senHeatGround      = -volHeatCapacityAir*conductanceSensible*(groundTemp - airTemp)
    latHeatGround      = -latHeatSubVapGround*latentHeatConstant*conductanceLatent * \
                          (satVP_GroundTemp*soilRelHumidity - VPair)


    # compute derivatives
    if ixDerivMethod == 'analytical' and not ixStability == 'moninObukhov':
        # compute derivatives for the ground fluxes w.r.t. ground temperature
        # d(ground sensible heat flux)/d(ground temp)
        dSenHeatGround_dTGround =   (-volHeatCapacityAir*dGroundCondSH_dGroundTemp)*(groundTemp - airTemp) + \
                                    (-volHeatCapacityAir*groundConductanceSH)
        # d(ground latent heat flux)/d(ground temp)
        dLatHeatGround_dTGround =   (-latHeatSubVapGround*latentHeatConstant*dGroundCondLH_dGroundTemp) * \
                                    (satVP_GroundTemp*soilRelHumidity - VPair) + \
                                    (-latHeatSubVapGround*latentHeatConstant*groundConductanceLH) * \
                                    dSVPGround_dGroundTemp*soilRelHumidity

    ########
    # net turbulent flux at the ground surface (W m-2)
    turbFluxGround = senHeatGround + latHeatGround

    # compute derivatives
    if ixDerivMethod == 'analytical':
        # (energy derivatives)
        # derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
        dTurbFluxGround_dTGround = dSenHeatGround_dTGround + dLatHeatGround_dTGround
    else: # (just make sure we return something)
        # (energy derivatives)
        dTurbFluxGround_dTGround = -9999

    return (
        groundConductanceSH,                # ground conductance for sensible heat (m s-1)
        groundConductanceLH,                # ground conductance for latent heat (m s-1)
        # output: fluxes from non-vegetated surfaces
        senHeatGround,                      # sensible heat flux from ground surface (W m-2)
        latHeatGround,                      # latent heat flux from ground surface (W m-2)
        turbFluxGround,                     # net turbulent heat fluxes at the ground surface (W m-2)
        # output: energy flux derivatives
        dTurbFluxGround_dTGround,           # derivative in net turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
        )
###########################################################################################################
