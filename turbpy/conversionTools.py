import numpy as np


def vapPress(q, p):
    # Input
    # q :        specific humidity (g g-1)
    # p :        pressure (Pa)
    # Output
    # vapPress : vapor pressure (Pa)

    w_ratio = 0.622                     # molecular weight ratio of water to dry air (-)
    w = q / (1. - q)                    # mixing ratio (-)
    vapPress = (w / (w + w_ratio)) * p  # vapor pressure (Pa)
    return vapPress


def satVapPress(TC):
    # Uses Teten's formula to compute saturated vapor pressure (Pa)
    # temperature units are degC !!!!
    # input
    # TC : temperature (C)
    # output
    # SVP      : saturation vapor pressure (Pa)
    # dSVP_dT  : d(SVP)/dT

    # Crude check for temperature in celsius
    if np.size(TC) > 1 and min(TC) > 100.:
        TC = TC - 273.
    elif np.size(TC) == 1 and TC > 100.:
        TC = TC - 273.
    X1 = 17.27
    X2 = 237.30
    # Saturation water vapour pressure at 273.16K (Pa)
    SATVPFRZ = 610.8
    # Saturated Vapour Press (Pa)
    SVP = SATVPFRZ * np.exp((X1 * TC) / (X2 + TC))
    dSVP_dT = SVP * (X1 / (X2 + TC) - X1 * TC / (X2 + TC)**2.)

    return (SVP, dSVP_dT)
