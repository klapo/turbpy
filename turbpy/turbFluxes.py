import numpy as np

def turbFluxes(
        # input: model control
        ixDerivMethod,                  # choice of method used to compute derivative (analytical or numerical)
        ixStability,                    # method for calculating stability
        ixStabParam,                    # Variable holding stability scheme parameters (unitless)
        z0Ground,                       # Surface roughness length for log-layer (m)
        # input: above-canopy forcing data
        airTemp,                        # air temperature at some height above the surface (K)
        airPres,                        # air pressure of the air above the vegetation canopy (Pa)
        VPair,                          # vapor pressure of the air above the vegetation canopy (Pa)
        windspd,                        # wind speed above the canopy (m s-1)
        # input: canopy and ground temperature
        groundTemp,                     # ground temperature (K)
        # input: diagnostic variables
        soilRelHumidity,                # relative humidity in the soil pores [0-1]
        mHeight,                        # height of observations (m)
        snowDepth                       # depth of snow (m)
        ):
# --------------------------------------------------------------------------------------------------------------------
    if np.isnan(airTemp) or np.isnan(airPres) or np.isnan(VPair) or np.isnan(groundTemp):
        raise ValueError('Input data includes nan value.')

    ########
    # Local variables
    evapSmooth=1.                                  # smoothing parameter for latent heat (W m-2)
    volHeatCapacityAir = iden_air*Cp_air           # volumetric heat capacity of air (J m-3)
    latentHeatConstant = iden_air*w_ratio/airPres  # latent heat constant for (kg m-3 Pa-1)
    # soilResistance is an unclear term for snow -- there is no resistance from snow: Sellers (1992)
    # scalarSoilResistance = scalarGroundSnowFraction*1. + (1. - groundSnowFraction)*EXP(8.25 - 4.225*soilEvapFactor)
    soilResistance = 1.

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
    # Unpack stability parameter
    if ~np.isnan(ixParam):
        stabParams[ixStability] = ixParam

    ########
    # compute resistances
    derivDesired = 'analytical' in ixDerivMethod or 'numerical' in ixDerivMethod
    resistOut = aeroResist(
                    derivDesired,               # flag to indicate if analytical derivatives are desired
                    ixStability,                # choice of stability function
                    # input: above-canopy forcing data
                    mHeight,                    # measurement height (m)
                    airTemp,                    # air temperature at some height above the surface (K)
                    windspd,                    # wind speed at some height above the surface (m s-1)
                    # input: temperature (canopy, ground, canopy air space)
                    groundTemp,                 # ground temperature (K)
                    # input: diagnostic variables
                    snowDepth,                  # snow depth (m)
                    # input: parameters
                    z0Ground,                   # surface roughness length (below canopy/non-vegetated [snow]) (m)
                    critRichNumber,             # critical value for the bulk Richardson number (-)
                    Louis79_bparam,             # parameter in Louis (1979) stability function
                    Mahrt87_eScale              # exponential scaling factor in the Mahrt (1987) stability function
                    )

    # Unpack resistances
    (RiBulkGround,                  # bulk Richardson number for the ground surface (-)
    groundStabilityCorrection,      # stability correction for the ground surface (-)
    groundResistance,               # below canopy aerodynamic resistance (s m-1)
    dGroundResistance_dTGround      # derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
    ) = resistOut

    ########
    # compute conductances, and derivatives...
    # NOTE: soilResistance accounts for fractional snow, and =0 when snow cover is 100%
    groundConductanceLH = 1./(groundResistance + soilResistance)
    groundConductanceSH = 1./groundResistance

    ########
    # compute sensible and latent heat fluxes, and derivatives...
    # (positive downwards)
    senHeatGround      = -volHeatCapacityAir*groundConductanceSH*(groundTemp - airTemp)
    latHeatGround      = -latHeatSubVapGround*latentHeatConstant*groundConductanceLH * \
                          (satVP_GroundTemp*soilRelHumidity - VPair)


    # compute derivatives
    if ixDerivMethod == 'analytical':
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
