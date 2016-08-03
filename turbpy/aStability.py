import numpy as np

from bulkRichardson import bulkRichardson
import multiConst as mc


def aStability(computeDerivative, ixStability, ixStabParam, mHeight, airTemp,
               airVaporPress, sfcTemp, sfcVaporPress, windspd, z0Ground):
    '''
    computeDerivative:      logical flag to compute analytical derivatives
    ixStability             choice of stability function
    ixStabParam             Variable holding stability scheme parameters (unitless)
    mHeight                 measurement height above the surface (m)
    airTemp                 air temperature (K)
    airVaporPress                   vapor pressure of air (Pa)
    sfcTemp                 surface temperature (K)
    VPground                Vapor pressure at the surface (Pa)
    windspd                 wind speed (m s-1)
    z0Ground                surface roughness length (below canopy/non-vegetated [snow]) (m)
    '''
# ------------------------------------------------------------------------------
# Sub-functions (stability schemes)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# "standard" stability correction, a la Anderson 1976
# ------------------------------------------------------------------------------
    def standard(critRichNumber=mc.stabParams['standard']):
        # compute surface-atmosphere exchange coefficient (-)
        if RiBulk < critRichNumber:
            stabilityCorrection = (1. - 5. * RiBulk)**2.
        else:
            stabilityCorrection = mc.machineEpsilon

        # compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
        if computeDerivative:
            if RiBulk < critRichNumber:
                dStabilityCorrection_dRich = -5. * 2. * (1. - 5. * RiBulk)
            if RiBulk >= critRichNumber:
                dStabilityCorrection_dRich = 0.
        else:
            dStabilityCorrection_dRich = np.nan

        # Calculate conductance
        conductance = (conductanceNeutral * windspd * stabilityCorrection)

        # Assume that latent and sensible heat have same conductance
        conductanceSensible = conductance
        conductanceLatent = conductance

        # Collect output
        bulkAerodynamicParameters = {'stabilityCorrection': stabilityCorrection,
                                     }
        bulkAerodynamicDerivatives = {
            'dStabilityCorrection_dRich': dStabilityCorrection_dRich,
            }
        return (bulkAerodynamicParameters,
                bulkAerodynamicDerivatives,
                conductanceSensible,
                conductanceLatent,
                )

# ------------------------------------------------------------------------------
# Louis 1979
# ------------------------------------------------------------------------------
    def louisInversePower(Louis79_bparam=mc.stabParams['louisInversePower']):
        # compute surface-atmosphere exchange coefficient (-)
        bprime = Louis79_bparam / 2.
        stabilityCorrection = 1. / ((1. + bprime * RiBulk)**2.)
        if stabilityCorrection < mc.machineEpsilon:
            stabilityCorrection = mc.machineEpsilon

        # compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
        if computeDerivative:
            dStabilityCorrection_dRich = (bprime * -2.
                                          * (1. + bprime * RiBulk)**(-3.))
        else:
            dStabilityCorrection_dRich = np.nan

        # Calculate conductance
        conductance = (conductanceNeutral * windspd * stabilityCorrection)
        # Assume that latent and sensible heat have same conductance
        conductanceSensible = conductance
        conductanceLatent = conductance

        # Collect output
        bulkAerodynamicParameters = {'stabilityCorrection': stabilityCorrection,
                                     }
        bulkAerodynamicDerivatives = {
            'dStabilityCorrection_dRich': dStabilityCorrection_dRich,
            }
        return (bulkAerodynamicParameters,
                bulkAerodynamicDerivatives,
                conductanceSensible,
                conductanceLatent,
                )

# ------------------------------------------------------------------------------
# Mahrt 1987
# ------------------------------------------------------------------------------
    def mahrtExponential(Mahrt87_eScale=mc.stabParams['mahrtExponential']):
        # compute surface-atmosphere exchange coefficient (-)
        stabilityCorrection = np.exp(-Mahrt87_eScale * RiBulk)
        if stabilityCorrection < mc.machineEpsilon:
            stabilityCorrection = mc.machineEpsilon

        # compute derivative in stability correction w.r.t. temperature (K-1)
        if computeDerivative:
            dStabilityCorrection_dRich = (-Mahrt87_eScale
                                          * np.exp(-Mahrt87_eScale * RiBulk)
                                          )
        else:
            dStabilityCorrection_dRich = np.nan

        # Calculate conductance
        conductance = (conductanceNeutral * windspd * stabilityCorrection)
        # Assume that latent and sensible heat have same conductance
        conductanceSensible = conductance
        conductanceLatent = conductance

        # Collect output
        bulkAerodynamicParameters = {'stabilityCorrection': stabilityCorrection,
                                     }
        bulkAerodynamicDerivatives = {
            'dStabilityCorrection_dRich': dStabilityCorrection_dRich,
            }
        return (bulkAerodynamicParameters,
                bulkAerodynamicDerivatives,
                conductanceSensible,
                conductanceLatent,
                )

    ' Stability Error Message '
    def stabErrMess():
        raise ValueError(
            'Unrecognized stability choice: '
            + ixStability
            + '\n'
            + 'Valid stability options: '
            + stabilityCase.keys
            )

# ------------------------------------------------------------------------------
# Monin-Obukhov stability as implemented in SNTHERM
# ------------------------------------------------------------------------------
    def moninObukhov(z0Ground=z0Ground):
        # Assume that surface drag parameter is equal
        # between momentum, heat, and moisture
        z0Groundw = z0Ground
        z0Groundh = z0Ground
        z0Groundq = z0Ground
        # log gradient for temperature, moisture, and wind
        dlogT = np.log(mHeight / z0Groundh)
        dlogQ = np.log(mHeight / z0Groundq)
        dlogW = np.log(mHeight / z0Groundw)

        # Iteration control
        numMaxIterations = 50       # Number of iterations for convergence
        itTolerance = 1. * 10**(-2)  # Tolerance iterative solution
        # Bounds on zetaT for N-R scheme
        zetaLow = 0.
        zetaHigh = 500.

        # Initial Values: assume neutral stability for first guess
        psiM = 0.					# Stability function for momentum
        psiH = 0.					# Stability function for heat
        psiQ = 0.					# Stability function for water vapor
        deltaT = airTemp - sfcTemp
        deltaE = (airVaporPress - sfcVaporPress) / 100.  # Convert [Pa] -> [hPa]
        Tbar = (airTemp + sfcTemp) / 2.
        zetaT = RiBulk * dlogW      # Stability parameter
        L = mHeight / zetaT 		# Obukhov height

        ########
        # Iterative solution to CH,CD, and L
        for it in np.arange(0, numMaxIterations):
            # Newton-Raphson iterative loop. Loop at least twice.
            if it == 0:
                # Initial pass
                if RiBulk > 1.427:
                    # Exceeds critical Richardson number
                    # Decouple the atmosphere and land surface
                    conductanceSensible = mc.machineEpsilon
                    conductanceLatent = mc.machineEpsilon

                    moninObukhovParameters = {'L': L,
                                              'psiM': psiM,
                                              'psiH': psiH,
                                              'psiQ': psiQ,
                                              'zeta': mHeight / L,
                                              'zetaT': zetaT,
                                              'freelim': 0.,
                                              'freelimv': 0.,
                                              }
                    moninObukhovDerivatives = {}
                    return (moninObukhovParameters, moninObukhovDerivatives,
                            conductanceSensible, conductanceLatent)

            #########
            # a. Compute momentum stability functions, transfer coefficient Cd,
            # and friction veolocity ustar
            zeta = mHeight / L
            if L >= 0.:
                # Stable
                psiM = Dutch(zeta)
            else:
                # Unstable
                psiM = Unstable(L, mHeight, z0Groundh, 1)
            sqrtCd = mc.vkc / max((dlogW - psiM), 0.1)
            ustar = sqrtCd * windspd

            #########
            # b. Estimate scalar roughness lengths (zh, zq,and zm) and
            # stability functions (heat and moisture)
            #
            # Assumes that measurement height for temperature and
            # water vapor are made at the same height as wind. If this
            # assumption is broken, psiQ and psiH need be calculated
            # using the "Dutch" function.
            [z0Groundh, z0Groundq, dlogQ, dlogT] = \
                andreas(ustar, z0Ground, airTemp, mHeight)
            if L >= 0.:
                # Stable
                psiH = psiM
                psiQ = psiH
            else:
                # Unstable
                psiH = Unstable(L, mHeight, z0Groundh, 2)
                psiQ = psiH
            # Scalar lengths for temperature and moisture
            tstar = mc.vkc * deltaT / (dlogT - psiH)
            qstar = (mc.vkc * deltaE * 100
                     / (mc.R_wv * airTemp * mc.iden_air)
                     / (dlogQ - psiQ))

            #########
            # c. Compute Monin-Obukhov length (L)
            if (not tstar == 0.) and (not ustar == 0.):
                # These expressions from SNTHERM are combined below in `L`
                # They seem similar to B6.22, but do not match exactly...
                # L = 1 + 0.61 .* Tbar.*qstar./tstar;
                # dum=L;
                # L = L .* 9.8 .* 0.4 .* tstar./ (Tbar.* wstar.^2);
                # L=1/L;
                B0 = (1. + 0.61 * Tbar * qstar / tstar)
                L = (1. / (B0
                           * (mc.gravity * mc.vkc * tstar)
                           / (Tbar * ustar**2.)
                           )
                     )
            elif (tstar == 0.) and (not ustar == 0.):
                # Extremely stable, decouple the land surface
                L = 1. * 10.**14.

            #########
            # d. Limits for extremely stable and unstable stratification
            freelim = 0
            freelimv = 0
            if zetaT > 500 or mHeight / L > 500:
                # Very stable
                # Decouple the atmosphere and land surface
                conductanceSensible = mc.machineEpsilon
                conductanceLatent = mc.machineEpsilon
                moninObukhovParameters = {'L': L,
                                          'psiM': psiM,
                                          'psiH': psiH,
                                          'psiQ': psiQ,
                                          'zeta': zeta,
                                          'zetaT': zetaT,
                                          'freelim': freelim,
                                          'freelimv': freelimv,
                                          }
                moninObukhovDerivatives = {}
                return (moninObukhovParameters, moninObukhovDerivatives,
                        conductanceSensible, conductanceLatent)

            elif L < 0. and L > -1.0:
                # Very unstable:
                # Compute free convection limits.  Applies generally for L < -0.1
                # Andreas and Cash, Convective heat transfer over wintertime
                # leads and polynyas, JGR 104, C11, 25,721-25,734, 1999.
                # Method is for a smooth surface with unlimited fetch.

                # Difference in bouyancy flux
                deltaB = -(mc.gravity * deltaT / Tbar) * B0
                # Kinematic viscosity of air
                nu = 1.315 * 10**(-5.)
                D = 0.024 / (mc.Cp_air * mc.iden_air)
                Dv = 1.1494 * D

                zscale = (nu * D / deltaB)**(1. / 3.)
                zscalev = (nu * Dv / deltaB)**(1. / 3.)
                freelim = -0.15 * mc.iden_air * mc.Cp_air * D * deltaT / zscale

                # Commented code handles non-snow covered surfaces...
                # if scalarGroundSnowFraction == 1:
                freelimv = (-0.15 * Dv * deltaE * 2.838 * 10**(6.)
                            / (zscalev * mc.R_wv * airTemp * mc.iden_air)
                            )
                # else:
                #   freelimv=(dlogT/dlogQ)*freelimv

            ########
            # e. Compute change in zeta
            y = (mHeight / L) - zetaT
            if (L > 0.) and (y > 0.):
                # Update lower bound
                zetaLow = zetaT
            elif (L > 0.) and (y <= 0.):
                # Update upper bound
                zetaHigh = zetaT

            if (np.abs(y) < itTolerance) and (it >= 2):
                # Within tolerance, leave loop
                conductanceSensible = (mc.vkc**2.
                                       / ((dlogW - psiM) * (dlogT - psiH)))
                conductanceLatent = (mc.vkc**2.
                                     / ((dlogW - psiM) * (dlogQ - psiQ)))
                # Collect MOST related parameters for use outside function
                moninObukhovParameters = {'L': L,
                                          'psiM': psiM,
                                          'psiH': psiH,
                                          'psiQ': psiQ,
                                          'zeta': zeta,
                                          'zetaT': zetaT,
                                          'freelim': freelim,
                                          'freelimv': freelimv,
                                          }
                moninObukhovDerivatives = {}
                return (moninObukhovParameters, moninObukhovDerivatives,
                        conductanceSensible, conductanceLatent)
            if it >= numMaxIterations - 1:
                # NOT within tolerance, leave loop
                conductanceSensible = (mc.vkc**2.
                                       / ((dlogW - psiM) * (dlogT - psiH)))
                conductanceLatent = (mc.vkc**2.
                                     / ((dlogW - psiM) * (dlogQ - psiQ)))
                # Alert the user
                # print('Convergence problem in turbpy.aStability.moninObukhov')
                # Collect MOST related parameters for use outside function
                moninObukhovParameters = {'L': L,
                                          'psiM': psiM,
                                          'psiH': psiH,
                                          'psiQ': psiQ,
                                          'zeta': zeta,
                                          'zetaT': zetaT,
                                          'freelim': freelim,
                                          'freelimv': freelimv,
                                          }
                moninObukhovDerivatives = {}
                return (moninObukhovParameters, moninObukhovDerivatives,
                        conductanceSensible, conductanceLatent)
            #########
            # f. Use Newton-Raphson scheme for stable
            # sequential estimates for unstable
            if L >= 0.:
                # Next is derivative of Psi
                DpsitDzetaT = (-0.70 - 0.75
                               * np.exp(-0.35 * zetaT)
                               * (6. - 0.35 * zetaT)
                               )
                DpsiMDzetaT = (-0.70 - 0.75
                               * np.exp(-0.35 * zeta)
                               * (6. - 0.35 * zeta))
                dumM = 1. / (dlogW - psiM)
                dumT = 1. / (dlogT - psiH)
                # Derivative of y (ignoring humidity effects)
                yprime = ((mHeight / L)
                          * ((DpsitDzetaT * dumT)
                             - (2. * DpsiMDzetaT * dumM)
                             )
                          - 1.
                          )
                zetaT = zetaT - (y / yprime)
                if (zetaT <= zetaLow) or (zetaT >= zetaHigh):
                    # Root out of range. Bisection instead of Newton-Raphson
                    zetaT = (zetaHigh + zetaLow) / 2.
                else:
                    zetaT = mHeight / L
                # Update Obukhov length
                L = mHeight / zetaT

# ------------------------------------------------------------------------------
# Sub-functions to Monin-Obukhov
# ------------------------------------------------------------------------------
    def Unstable(L, Ht, zt, icall):
        # Set an upper limit on function at L = -.1 or  at 200 times
        # the large -Ht/L limit where log(-4*Ht/L) is 0.  This asymptotic
        # limit is from Godfrey and Beljaars (1991).  This is a semi-
        # empirical fix to keep the heat flux from increasing with
        # decreasing -L.
        Llim = min(-.1, -100. * zt)
        zeta = mHeight / min(L, Llim)
        x = (1. - 16. * zeta)**(0.25)
        if icall == 1:
            stab = (np.log((1 + x**2.) / 2.)
                    + 2. * np.log((1. + x) / 2.)
                    - 2. * np.arctan(x) + 1.5707963)
        else:
            stab = 2. * np.log((1. + x**2.) / 2.)
        return stab

# ------------------------------------------------------------------------------
    def Dutch(zeta):
        stab = (-(0.70 * zeta
                + 0.75 * (zeta - 14.28)
                * np.exp(-0.35 * zeta)
                + 10.71))
        return stab

# ------------------------------------------------------------------------------
    def andreas(ustar, z0Ground, airTemp, mHeight):
        #    Compute scaler roughness lengths using procedure
        #    in Andreas, 1987, Boundary-Layer Meteorology 38, 159-184.  Only
        #    use this procedure for snow.
        #
        #     Expression for kinematic viscosity of air is taken from program
        #     of Launiainen and Vihma, 1990, Environmental Software, vol. 5,
        #     No. 3, pp. 113 - 124.

        # Heat
        Reynolds = z0Ground * ustar * 1. * 10.**7. / (.9065 * airTemp - 112.7)
        if Reynolds <= 0.135:
            # Smooth coefficients
            b0 = 1.250
            b1 = 0.
            b2 = 0.
        elif Reynolds < 2.5:
            # Transition coefficients
            b0 = 0.149
            b1 = -0.550
            b2 = 0.
        else:
            # Rough coefficients
            b0 = 0.317
            b1 = -0.565
            b2 = -0.183
        z0Groundh = (z0Ground
                     * np.exp(b0
                              + b1 * np.log(Reynolds)
                              + b2 * np.log(Reynolds)**2.))

        # Moisture
        if Reynolds <= 0.135:
            b0 = 1.61
            b1 = 0.
            b2 = 0.
        elif Reynolds < 2.5:
            b0 = 0.351
            b1 = -0.628
            b2 = 0.
        else:
            b0 = 0.396
            b1 = -0.512
            b2 = -0.180
        z0Groundq = (z0Ground
                     * np.exp(b0
                              + b1 * np.log(Reynolds)
                              + b2 * np.log(Reynolds)**2.
                              )
                     )

        # Log profiles
        dlogT = np.log(mHeight / z0Groundh)
        dlogQ = np.log(mHeight / z0Groundq)

        return (z0Groundh, z0Groundq, dlogQ, dlogT)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
    # Dictionary to stability function calls
    stabilityCase = {
        'standard': standard,
        'louisInversePower': louisInversePower,
        'mahrtExponential': mahrtExponential,
        'moninObukhov': moninObukhov,
        }

    ########
    # Determine atmospheric stability
    # compute the bulk Richardson number (-)
    bulkRichardsonOut = bulkRichardson(
        airTemp,                        # air temperature (K)
        sfcTemp,                        # surface temperature (K)
        windspd,                        # wind speed (m s-1)
        mHeight,                        # measurement height (m)
        computeDerivative,              # flag to compute the derivative
        )

    (RiBulk,                                    # bulk Richardson number (-)
        dRiBulk_dAirTemp,                       # derivative in the bulk Richardson number w.r.t. air temperature (K-1)
        dRiBulk_dSfcTemp) = bulkRichardsonOut   # derivative in the bulk Richardson number w.r.t. surface temperature (K-1)

    # Conductance under conditions of neutral stability (-)
    conductanceNeutral = (mc.vkc**2.) / (np.log((mHeight) / z0Ground)**2.)

    #########
    # Stability Corrections
    if RiBulk < 0. and 'moninObukhov' not in ixStability:
        #########
        # Unstable
        stabilityCorrection = (1. - 16. * RiBulk)**(0.5)
        conductance = (conductanceNeutral * windspd * stabilityCorrection)
        # Assume latent and sensible heat have same conductance.
        conductanceSensible = conductance
        conductanceLatent = conductance
        # Assign to output dictionary
        stabilityCorrectionParameters = {
            'stabilityCorrection': stabilityCorrection}

        # compute derivative in stability  w.r.t. temperature (K-1)
        if computeDerivative:
            dStabilityCorrection_dRich = ((-16.) * 0.5
                                          * (1. - 16. * RiBulk)**(-0.5))
            dStabilityCorrection_dAirTemp = (dRiBulk_dAirTemp
                                             * dStabilityCorrection_dRich)
            dStabilityCorrection_dSfcTemp = (dRiBulk_dSfcTemp
                                             * dStabilityCorrection_dRich)
        else:
            # set derivative to nan if not computing
            dStabilityCorrection_dAirTemp = np.nan
            dStabilityCorrection_dSfcTemp = np.nan
            dStabilityCorrection_dRich = np.nan

        # Collect derivatives for output.
        stabilityCorrectionDerivatives = {
            'dStabilityCorrection_dAirTemp': dStabilityCorrection_dAirTemp,
            'dStabilityCorrection_dSfcTemp': dStabilityCorrection_dSfcTemp,
            'dStabilityCorrection_dRich': dStabilityCorrection_dRich,
            }
    else:
        ########
        # Stable cases
        # Use a dictionary of functions to select stability function.
        func = stabilityCase.get(ixStability, stabErrMess)
        (stabilityCorrectionParameters,
         stabilityCorrectionDerivatives,
         conductanceSensible,
         conductanceLatent) = func()

        ########
        # Derivatives
        if computeDerivative:
            # Derivative of stability functions w.r.t. airTemp and sfcTemp
            dStabilityCorrection_dAirTemp = (dRiBulk_dAirTemp
                                             * dStabilityCorrection_dRich)
            dStabilityCorrection_dSfcTemp = (dRiBulk_dSfcTemp
                                             * dStabilityCorrection_dRich)
        else:
            # set derivative to nan if not computing
            dStabilityCorrection_dAirTemp = np.nan
            dStabilityCorrection_dSfcTemp = np.nan

        # Collect derivatives into dictionary
        stabilityCorrectionDerivatives['dStabilityCorrection_dAirTemp'] = \
            dStabilityCorrection_dAirTemp
        stabilityCorrectionDerivatives['dStabilityCorrection_dSfcTemp'] = \
            dStabilityCorrection_dSfcTemp

    # Add neutral conductance to output dictionary
    stabilityCorrectionParameters['conductanceNeutral'] = conductanceNeutral
    return (RiBulk,
            stabilityCorrectionParameters,
            stabilityCorrectionDerivatives,
            conductanceSensible,
            conductanceLatent,
            )
