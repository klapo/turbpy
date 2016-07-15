def surfFluxCalc_moninObukhov(  # Inputs
                                dlogW,
                                dlogT,
                                psiM,
                                psiH,
                                psiQ,
                                airTemp,
                                windspd,
                                ):
    '''
    The Monin-Obukhov surface turbulent fluxes. Many of these equations and definitions are derived from Chris Bretherton's notes on boundary layer meteorology (atmos.washington.edu/~breth/classes/AS547) or from Andreas 2001. Equation numbers will be noted at either B.XX or A.XX for reference.

    Andreas, E. L. (2001), Parameterizing Scalar Transfer over Snow and Ice: A Review, J. Hydrometeorol., 3, 417–432.
    '''


    # Define constants
    Ck = 0								# Windless exchange coefficient for water vapor [Wm^-2]
    Csk = 2								# Windless exchange coefficient for heat [Wm^-2]

    # bulk aerodynamic transfer coefficients with stability corrections (B.6.18 and 6.19)
    C_D = vkc**2 / ( (dlogW - psiM) * (dlogT - psiH) )
    C_H = vkc**2 / (R_wv*tkair*(dlogW - psiM)*(dlogQ - psiQ))

	sen = Cp_air * rhoair * windspd *
    # Since vapor pressure is in mb, fqlat is multiplied by 100 to
    # convert to Pascals
	lat = windspd * * 10**2 /

    ## Final computation of qsen and qsen.
    # Note qlat includes the windless coefficients
    qlat = Ck + qlat0 * 2.838*10**6
    if L >= 0:
    	qsen=qsen0
    	if (qsen0 < np.max(Csk,.1)) & (qsen0 > 0):
    		# Minimum value
    		qsen = np.max(Csk,.1)

    return qsen, qlat

## Final computation of qsen and qsen.
# Note qlat includes the windless coefficients
# qlat = Ck + qlat0 * 2.838*10**6
# if L >= 0:
# 	qsen=qsen0
# 	if (qsen0 < np.max(Csk,.1)) & (qsen0 > 0):
# 		# Minimum value
# 		qsen = np.max(Csk,.1)
# else:
# 	print('moninObukhov_SNTHERM arrived at a stable Obukhov length')
# 	return None
