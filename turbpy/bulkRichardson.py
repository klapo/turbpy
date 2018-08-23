import turbpy.multiConst as mc


def bulkRichardson(airTemp,                  # air temperature (K)
                   sfcTemp,                  # surface temperature (K)
                   windspd,                  # wind speed (m s-1)
                   mHeight,                  # measurement height (m)
                   ):
    # -------------------------------------------------------------------------------------------------------
    # Local variables
    T_grad = airTemp - sfcTemp
    T_mean = 0.5 * (airTemp + sfcTemp)
    RiMult = (mc.gravity * mHeight) / (windspd**2.)
    # compute the Richardson number
    RiBulk = (T_grad / T_mean) * RiMult

    return (RiBulk)             # bulk Richardson number (-)
