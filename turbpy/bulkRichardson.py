import turbpy.multiConst as mc
import numpy as np

def bulkRichardson(airTemp,                  # air temperature (K)
                   sfcTemp,                  # surface temperature (K)
                   windspd,                  # wind speed (m s-1)
                   mHeight,                  # measurement height (m)
                   ):
    # -------------------------------------------------------------------------------------------------------
    # Local variables
    if np.max(airTemp) < 200:
        airTemp = airTemp + 273.15
    if np.max(sfcTemp) < 200:
        sfcTemp = sfcTemp + 273.15

    T_grad = airTemp - sfcTemp
    T_mean = 0.5 * (airTemp + sfcTemp)
    RiMult = (mc.gravity * mHeight) / (windspd**2.)
    # compute the Richardson number
    RiBulk = (T_grad / T_mean) * RiMult

    return (RiBulk)             # bulk Richardson number (-)
