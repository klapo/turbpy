#
#
# def surfFluxCalc_moninObukhov(  # Inputs
#     dlogW,
#     dlogT,
#     psiM,
#     psiH,
#     psiQ,
#     airTemp,
#     windspd,
#     ):
#     '''
#     The Monin-Obukhov surface turbulent fluxes. Many of these equations and definitions are derived from Chris Bretherton's notes on boundary layer meteorology (atmos.washington.edu/~breth/classes/AS547) or from Andreas 2001. Equation numbers will be noted at either B.XX or A.XX for reference.
#
#     Andreas, E. L. (2001), Parameterizing Scalar Transfer over Snow and Ice: A Review, J. Hydrometeorol., 3, 417–432.
#     '''
#
#     # Windless exchange coefficients for water vapor, k, and heat, sk [Wm^-2]
#     Ck = 0
#     Csk = 2
#
#     # bulk aerodynamic transfer coefficients with stability corrections (B.6.18 and 6.19)
#     C_D = vkc**2 / ( (dlogW - psiM) * (dlogT - psiH) )
#     C_H = vkc**2 / (R_wv*tkair*(dlogW - psiM)*(dlogQ - psiQ))
#
# 	sen = Cp_air * rhoair * windspd *
#     # Since vapor pressure is in mb, fqlat is multiplied by 100 to
#     # convert to Pascals
# 	lat = windspd * * 10**2 /
#
#     ## Final computation of qsen and qsen.
#     # Note qlat includes the windless coefficients
#     qlat = Ck + qlat0 * 2.838*10**6
#
# #########
# # d. Limits for extremely stable and unstable stratification
# freelim = 0
# freelimv = 0
# if (zetaT > 500) | (mHeight / L > 500):
#     # Very stable, decouple
#     conductanceSensible = machineEpsilon
#     conductanceLatent = machineEpsilon
#     return (psiM, psiH, conductanceSensible, conductanceLatent)
# elif L < 0 & L > -1.0:
#     # Very unstable:
#     # Compute free convection limits.  Applies generally for L < -0.1
#     # Andreas and Cash, Convective heat transfer over wintertime
#     # leads and polynyas, JGR 104, C11, 25,721-25,734, 1999.
#     # Method is for a smooth surface with unlimited fetch.
#     deltaB = (-9.8 * deltaT / Tbar) * \
#         (1 + 0.61 * Tbar * qstar / tstar)
#     nu = 1.315 * 10**(-5) 					# Kinematic viscosity of air
#     D = 0.024 / (cdryair * rhoair) 			# Replace with definitions from multiConst
#     Dv = 1.1494 * D
#     zscale = (nu * D / deltaB)**(1. / 3.)
#     zscalev = (nu * Dv / deltaB)**(1. / 3.)
#     freelim = -0.15 * rhoair * cdryair * D * deltaT / zscale
#     freelimv = -0.15 * Dv * deltaE * con * 2.838 * 10**(6.) / zscalev
#
# # Final computation of qsen and qsen.
# # Note qlat and qsen add the windless coefficients
# qlat = Ck + qlat0 .*2.838.*10.^6;
# if L >= 0
#
# 	qsen=qsen0;
#
# 	if qsen0 < max(Csk,.1) && qsen0 > 0
# 		% Minimum value
# 		qsen = max(Csk,.1);
# 	end
# else
# 	dum = 0;
# 	if deltaT ~= 0
# 		dum=-freelim/deltaT;
# 	end
# 	qsen = max(qsen0,dum);
#
# 	dum = 0;
# 	if deltaE ~= 0
# 		dum = -freelimv/deltaE;
# 	end
# 	qlat = max(qlat,dum);
# end
