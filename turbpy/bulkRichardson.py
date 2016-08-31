import turbpy.multiConst as mc


def bulkRichardson(airTemp,                  # air temperature (K)
                   sfcTemp,                  # surface temperature (K)
                   windspd,                  # wind speed (m s-1)
                   mHeight,                  # measurement height (m)
                   computeDerivative=False,  # flag to compute the derivative
                   ):
    # -------------------------------------------------------------------------------------------------------
    # Local variables
    T_grad = airTemp - sfcTemp
    T_mean = 0.5 * (airTemp + sfcTemp)
    RiMult = (mc.gravity * mHeight) / (windspd**2.)
    # compute the Richardson number
    RiBulk = (T_grad / T_mean) * RiMult

    ########
    # compute the derivative in the Richardson number
    if computeDerivative:
        dRiBulk_dAirTemp = (RiMult / T_mean
                            - (RiMult * T_grad)
                            / (0.5 * ((airTemp + sfcTemp)**2.)))
        dRiBulk_dSfcTemp = (-RiMult / T_mean
                            - RiMult * T_grad
                            / (0.5 * ((airTemp + sfcTemp)**2.)))
    else:
        dRiBulk_dAirTemp = 1.
        dRiBulk_dSfcTemp = 1.

    return (RiBulk,             # bulk Richardson number (-)
            dRiBulk_dAirTemp,   # derivative in the bulk Richardson number w.r.t. air temperature (K-1)
            dRiBulk_dSfcTemp)   # derivative in the bulk Richardson number w.r.t. surface temperature (K-1)
