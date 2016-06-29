import numpy as np

def aeroResist( derivDesired,               # flag to indicate if analytical derivatives are desired
                ixStability,                # choice of stability function
                # input: forcing data
                mHeight,                    # measurement height (m)
                airTemp,                    # air temperature at some height above the surface (K)
                windspd,                    # wind speed at some height above the surface (m s-1)
                # input: diagnostic variables
                groundTemp,                 # ground temperature (K)
                snowDepth,                  # snow depth (m)
                # input: parameters
                z0Ground,                   # surface roughness length (below canopy/non-vegetated [snow]) (m)
                critRichNumber,             # critical value for the bulk Richardson number (-)
                Louis79_bparam,             # parameter in Louis (1979) stability function
                Mahrt87_eScale              # exponential scaling factor in the Mahrt (1987) stability function
                ):
# -------------------------------------------------------------------------------------------------------
# compute aerodynamic resistances
# -------------------------------------------------------------------------------------------------------

    # Constantish
    C_r = 0.3                 # roughness element drag coefficient (-) from Raupach (BLM, 1994)
    C_s = 0.00                # substrate surface drag coefficient (-) from Raupach (BLM, 1994)
    approxDragCoef_max = 0.   # maximum value of the approximate drag coefficient (-) from Raupach (BLM, 1994)
    vkc = 0.4                 # von Karman constant (-)

    ########
    # compute resistances (no canopy)
    # turbulent transfer coefficient under conditions of neutral stability (-)
    groundExNeut = (vkc**2.) / ( np.log((mHeight - snowDepth)/z0Ground)**2.)
    # compute the resistance between the surface and canopy air UNDER NEUTRAL CONDITIONS (s m-1)
    groundResistanceNeutral = 1. / (groundExNeut*windspd)

    # check that measurement height above the ground surface is above the roughness length
    if mHeight < snowDepth+z0Ground:
        raise ValueErrorr('measurement height < snow depth + roughness length')
    # define height above the snow surface
    heightAboveGround  = mHeight - snowDepth

    # check that measurement height above the ground surface is above the roughness length
    if heightAboveGround < z0Ground:
        print(\
                'z0Ground = %d \n' \
                'mHeight  = %d \n' \
                'snowDepth = %d \n' \
                'heightAboveGround = %d \n' \
                , (z0Ground,mHeight,snowDepth,heightAboveGround))
        raise ValueErrorr('Height above ground < roughness length [likely due to snow accumulation]')

    # compute ground stability correction
    aStabilityOut = aStability(
                  # input
                  derivDesired,                                 # logical flag to compute analytical derivatives
                  ixStability,                                  # choice of stability function
                  # input: forcing data, diagnostic and state variables
                  heightAboveGround,                            # measurement height above the ground surface (m)
                  airTemp,                                      # temperature above the ground surface (K)
                  groundTemp,                                   # trial value of surface temperature (K)
                  windspd,                                      # wind speed above the ground surface (m s-1)
                  # input: stability parameters
                  critRichNumber,                               # critical value for the bulk Richardson number (-)
                  Louis79_bparam,                               # parameter in Louis (1979) stability function
                  Mahrt87_eScale,                               # exponential scaling factor in Mahrt stability
                  )

    # Unpack
    (RiBulkGround,                               # bulk Richardson number (-)
    groundStabilityCorrection,                   # stability correction for turbulent heat fluxes (-)
    dGroundStabilityCorrection_dRich,            # derivative in stab. corr. w.r.t. Ri for the ground surface (-)
    dGroundStabilityCorrection_dAirTemp,         # derivative in stab. corr. w.r.t. air temperature (K-1)
    dGroundStabilityCorrection_dSfcTemp)  = aStabilityOut # derivative in stab. corr. w.r.t. surface temperature (K-1)

    # compute the ground resistance (after stability corrections)
    groundResistance = groundResistanceNeutral/groundStabilityCorrection
    if groundResistance < 0.:
        raise ValueErrorr('Ground resistance < 0 [no vegetation]')

    # if analytical derivatives are desired
    if derivDesired:
        # compute derivatives for ground resistance
        dGroundResistance_dTGround = -dGroundStabilityCorrection_dSfcTemp/ \
                                        (windspd*groundExNeut*groundStabilityCorrection**2.)
    else:
        dGroundResistance_dTGround = -9999

    return  (RiBulkGround,              # bulk Richardson number for the ground surface (-)
        groundStabilityCorrection,      # stability correction for the ground surface (-)
        groundResistance,               # below canopy aerodynamic resistance (s m-1)
        dGroundResistance_dTGround)     # derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)

###########################################################################################################
