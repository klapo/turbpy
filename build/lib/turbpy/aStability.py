from multiConst import *
from bulkRichardson import bulkRichardson

def aStability(# input: control
            computeDerivative,              # logical flag to compute analytical derivatives
            ixStability,                    # choice of stability function
            ixStabParam,                    # Variable holding stability scheme parameters (unitless)
            # input: forcing data, diagnostic and state variables
            mHeight,                        # measurement height (m)
            airTemp,                        # air temperature (K)
            sfcTemp,                        # surface temperature (K)
            windspd,                        # wind speed (m s-1)
            # input: stability parameters
            z0Ground,                   # surface roughness length (below canopy/non-vegetated [snow]) (m)
            ):

########
# Sub-functions (stability schemes)
# --------------------------------------------------------------------------------------------------------------------
    #### "standard" stability correction, a la Anderson 1976
    def standard(critRichNumber=stabParams['standard']):
        # compute surface-atmosphere exchange coefficient (-)
        if RiBulk <  critRichNumber:
            stabilityCorrection = (1. - 5.*RiBulk)**2.
        elif RiBulk >= critRichNumber:
            # stabilityCorrection = epsilon(stabilityCorrection)
            stabilityCorrection = machineEpsilon
        else:
            stabilityCorrection = machineEpsilon

        # compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
        if computeDerivative:
            if RiBulk <  critRichNumber:
                dStabilityCorrection_dRich = (-5.) * 2.*(1. - 5.*RiBulk)
            if RiBulk >= critRichNumber:
                dStabilityCorrection_dRich = 0.
        else:
            dStabilityCorrection_dRich = -9999
        return stabilityCorrection,dStabilityCorrection_dRich

    #### Louis 1979
    def louisInversePower(Louis79_bparam=stabParams['louisInversePower']):
        # scale the "b" parameter for stable conditions
        bprime = Louis79_bparam/2.
        # compute surface-atmosphere exchange coefficient (-)
        stabilityCorrection = 1. / ( (1. + bprime*RiBulk)**2. )
        if stabilityCorrection < machineEpsilon:
            stabilityCorrection = machineEpsilon
        # compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
        if computeDerivative:
            dStabilityCorrection_dRich = bprime * (-2.)*(1. + bprime*RiBulk)**(-3.)
        else:
            dStabilityCorrection_dRich = -9999
        return stabilityCorrection,dStabilityCorrection_dRich

    #### Mahrt 1987
    def mahrtExponential(Mahrt87_eScale=stabParams['mahrtExponential']):
        # compute surface-atmosphere exchange coefficient (-)
        stabilityCorrection = np.exp(-Mahrt87_eScale * RiBulk)
        if stabilityCorrection < machineEpsilon:
            stabilityCorrection = machineEpsilon
        # compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
        if computeDerivative:
            dStabilityCorrection_dRich = (-Mahrt87_eScale) * np.exp(-Mahrt87_eScale * RiBulk)
        else:
            dStabilityCorrection_dRich = -9999
        return stabilityCorrection,dStabilityCorrection_dRich

    #### Stability Error Message
    def stabErrMess():
        raise ValueError('Unrecognized stability choice: '+ixStability+'\n',
             'Valid stability options: standard, louisInversePower, mahrtExponential')
# --------------------------------------------------------------------------------------------------------------------

    # compute the bulk Richardson number (-)
    bulkRichardsonOut = bulkRichardson(# input
                                        airTemp,                        # air temperature (K)
                                        sfcTemp,                        # surface temperature (K)
                                        windspd,                        # wind speed (m s-1)
                                        mHeight,                        # measurement height (m)
                                        computeDerivative,              # flag to compute the derivative
                                        )

    # output
    (RiBulk,                             # bulk Richardson number (-)
    dRiBulk_dAirTemp,                    # derivative in the bulk Richardson number w.r.t. air temperature (K-1)
    dRiBulk_dSfcTemp) = bulkRichardsonOut # derivative in the bulk Richardson number w.r.t. surface temperature (K-1)

    # set derivative to one if not computing it
    if not computeDerivative:
        dStabilityCorrection_dRich    = 1.
        dStabilityCorrection_dAirTemp = 1.
        dStabilityCorrection_dSfcTemp = 1.

    #########
    # process unstable cases
    if RiBulk<0.:
        # compute surface-atmosphere exchange coefficient (-)
        stabilityCorrection = (1. - 16.*RiBulk)**0.5
        # compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
        if computeDerivative:
            dStabilityCorrection_dRich    = (-16.) * 0.5 *(1. - 16.*RiBulk)**(-0.5)
            dStabilityCorrection_dAirTemp = dRiBulk_dAirTemp * dStabilityCorrection_dRich
            dStabilityCorrection_dSfcTemp = dRiBulk_dSfcTemp * dStabilityCorrection_dRich

    ########
    # process stable cases
    # Use a dictionary of functions to select stability function
    else:
        stabilityCase = {
            'standard': standard,
            'louisInversePower': louisInversePower,
            'mahrtExponential': mahrtExponential,
        }

        # Get the function from switcher dictionary
        func = stabilityCase.get(ixStability, stabErrMess)
        # Gather stability function parameters
        # if ixStabParam:
        #     stabParams[ixStability] = ixParam
        # Execute the function
        stabilityCorrection,dStabilityCorrection_dRich = func()

        # Derivative of stable cases
        if computeDerivative:
            dStabilityCorrection_dAirTemp = dRiBulk_dAirTemp * dStabilityCorrection_dRich
            dStabilityCorrection_dSfcTemp = dRiBulk_dSfcTemp * dStabilityCorrection_dRich

    return  (RiBulk,                                # bulk Richardson number (-)
            stabilityCorrection,                    # stability correction for turbulent heat fluxes (-)
            dStabilityCorrection_dRich,             # derivative in stab. corr. w.r.t. Ri for the ground surface (-)
            dStabilityCorrection_dAirTemp,          # derivative in stab. corr. w.r.t. air temperature (K-1)
            dStabilityCorrection_dSfcTemp)          # derivative in stab. corr. w.r.t. sfc temperature (K-1)
