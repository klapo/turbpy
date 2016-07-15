def moninObukhov_SNTHERM():

'''
Some assumptions:
	- all observations come from the same height? SNTHERM accomodates different
	observation heights using a separate z and zeta for each observation type

'''

## Constants and derived values
# Assume that surface drag parameter is equal between momentum, heat, and moisture
z0Groundw = z0Ground
z0Groundh = z0Ground
z0Groundq = z0Ground
dlogT = np.log(mHeight/z0Groundh)
dlogQ = np.log(mHeight/z0Groundq)
dlogW = np.log(mHeight/z0Groundw)
Ck = 0										# Windless exchange coefficient for water vapor [Wm^-2]
Csk = 2										# Windless exchange coefficient for heat [Wm^-2]
numMaxIterations = 50						# Maximum number of iterations to search for convergence
itTolerance = 1*10**(-2)                    # Tolerance in stability parameter for iterative solution

# Initial Values: assume neutral stability for first guess
psiM = 0 									# Stability function for momentum
psiH = 0 									# Stability function for heat
psiQ = 0 									# Stability function for water vapor
deltaT = airTemp - sfcTemp
deltaE = Eair - Esurface
Tbar = (airTemp + sfcTemp)/2
zetaT = RiBulk * dlogW 						# Stability parameter
L = mHeight/zetaT 							# Obukhov height
if L < 0:
	raise ValueErrorr('Obukhov length less than zero (unstable) before iteration.')

# Bounds on zetaT for N-R scheme
zetaLow = 0
zetaHigh = 500

########
# Iterative solution to CH,CD, and L
for it in np.arange(0:numMaxIterations)
	# Newton-Raphson iterative loop. Loop at least twice.
	if it == 0:
		# Initial pass
		if RiBulk > 1.427:
			# Exceeds critical Richardson number -- decouple the atmosphere and land surface
			conductanceHeat	  = machineEpsilon
			conductanceLatent = machineEpsilon
			return conductanceSensible, conductanceLatent

	#########
	# a. Compute momentum stability functions, transfer coefficient Cd, and friction veolocity ustar
	zeta = mHeight/L
	if L >= 0:
		# Stable
		psiM = Dutch(zeta)
	else
		# Unstable
		psiM = Unstable(L,H,z0h,1)
	sqrtCd = vkc / np.max( (dlogW - psiM), 0.1)
	ustar = sqrtCd * windspd

	#########
	# b. Estimate scalar roughness lengths (zh, zq,and zm) and stability functions (heat and moisture)
	[z0Groundh,z0Groundq,dlogQ,dlogT] = andreas(ustar,z0Ground,airTemp,mHeight)
	if L >= 0:
		# Stable: there is an assumption here that L can never become negative (unstable)
		psiH = psiM
		psiQ = psiH
	else:
		# Unstable
		Psih = Unstable(L,H,z0h,2)
		Psiq = Psih
	# Scalar lengths for temperature and moisture
	tstar = vkc * deltaT / (dlogT - psiH) # 0.4 = von Karman constant? If so replace it
	qstar = vkc * deltaE * con / (dlogQ - psiQ)

	#########
	# c. Compute Monin-Obukhov length (L)
	if tstar ~= 0 & ustar ~= 0:
		# Need to figure out these defintions...
		# From the matlab port of SNTHERM:
		# L = 1 + 0.61 .* Tbar.*qstar./tstar;
		# dum=L;
		# L = L .* 9.8 .* 0.4 .* tstar./ (Tbar.* wstar.^2);
		# L=1/L;
		# These expressions are combined below, seems similar to B6.22, but not exactly...
		L = 1. / ( (1 + 0.61 .* Tbar*qstar/tstar) * gravity * vkc * tstar / (Tbar * ustar**2) )
	elif (tstar == 0) & (ustar ~= 0):
		# Extremely stable, decouple the land surface
		L=1*10**14

	#########
	# d. Limits for extremely stable and unstable stratification
	freelim = 0
	freelimv = 0
	if (zetaT > 500) | (mHeight/L > 500):
		# Very stable, decouple
		conductanceHeat	  = machineEpsilon
		conductanceLatent = machineEpsilon
		return conductanceHeat,conductanceLatent
	elif L < 0 && L > -1.0:
	    # Very unstable:
	    # Compute free convection limits.  Applies generally for L < -0.1
	    # Andreas and Cash, Convective heat transfer over wintertime
	    # leads and polynyas, JGR 104, C11, 25,721-25,734, 1999.
	    # Method is for a smooth surface with unlimited fetch.
	    deltaB = (-9.8 .* deltaT ./ Tbar) * (1 + 0.61 .* Tbar.*qstar./tstar)
	    nu = 1.315*10**(-5) 							# Kinematic viscosity of air
	    D = 0.024 / (cdryair * rhoair) 					# Replace with definitions from multiConst
	    Dv = 1.1494*D
	    zscale = (nu.*D./deltaB)**(1./3.)
	    zscalev = (nu*Dv/deltaB)**(1./3.)
	    freelim = -0.15  * rhoair * cdryair * D * deltaT /zscale
	    freelimv = -0.15 * Dv * deltaE * con * 2.838*10**(6.) / zscalev

	########
	# e. Compute change in zeta
	y = (mHeight/L) - zetaT
	if L > 0 &  y >  0:
		# Update lower bound
		zetaLow = zetaT
	elif L > 0 & y <= 0:
		# Update upper bound
		zetaHigh = zetaT

	if (np.abs(y) < itTolerance) & (it >= 2):
		# Within tolerance and more than one iteration executed, leave loop
	    conductanceSensible = vkc**2 / ((dlogW - psiM) * (dlogT - psiH))
	    conductanceLatent   = vkc**2 / ((dlogW - psiM) * (dlogQ - psiQ))
		return conductanceSensible, conductanceLatent
	if it > 50:
		# NOT within tolerance and more than 50 iterations executed, leave loop
	    conductanceSensible = vkc**2 / ((dlogW - psiM) * (dlogT - psiH))
	    conductanceLatent   = vkc**2 / ((dlogW - psiM) * (dlogQ - psiQ))
		# Alert the user
		print('Convergence problem in moninObukhov_SNTHERM')
		return conductanceSensible, conductanceLatent

	## f. Use Newton-Raphson scheme for stable, sequential estimates for unstable
	if L >= 0:
		# Next is derivative of Psi
		DpsitDzetaT = -0.70 - 0.75 * exp(-0.35 * zetaT)*(6 - 0.35 * zetaT)
		DpsiMDzetaT = -0.70 - 0.75 * exp(-0.35 * zeta)*(6 - 0.35 * zeta)
		dumM = 1 / (dlogW-psiM)
		dumT = 1 / (dlogT-psiH)	# SNTHERM has dlogt, not dlogT, even though it appears no where else...
		# Next is derivative of y (ignoring humidity effects)
		yprime = (mHeight/L) * (DpsitDzetaT * dumT - 2 * DpsiMDzetaT * dumM) - 1
		zetaT = zetaT - (y/yprime)
		if (zetaT <= zetal) | (zetaT >= zetah):
			# Root out of range. Bisection instead of Newton-Raphson scheme
			zetaT = (zetah+zetal)/2
	else:
		zetaT = mHeight/L
	# Update Obukhov length
	L = mHeight/zetaT

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

########
# Sub-functions (stability schemes)
# --------------------------------------------------------------------------------------------------------------------
def Unstable(L,Ht,zt,icall):
	# Set an upper limit on function at L = -.1 or  at 200 times
	# the large -Ht/L limit where log(-4*Ht/L) is 0.  This asymptotic
	# limit is from Godfrey and Beljaars (1991).  This is a semi-
	# empirical fix to keep the heat flux from increasing with
	# decreasing -L.
	Llim = np.min(-.1,-100.*zt)
	zeta = mHeight / np.min(L,Llim)
	x = (1-16.*zeta)**(0.25)
	if icall == 1:
		stab = np.log((1+x*x)/2)+2 * np.log((1. + x)/2) - 2.*atan(x) + 1.5707963
	else
		stab = 2.*log((1+x*x)/2)

def Dutch(zeta):
	stab = -(0.70.* zeta + 0.75 .*(zeta-14.28) .* exp(-0.35*zeta)  + 10.71);
	return stab

def andreas(ustar,z0Ground,Tkair,mHeight):
	#    Compute scaler roughness lengths using procedure
	#    in Andreas, 1987, Boundary-Layer Meteorology 38, 159-184.  Only
	#    use this procedure for snow.
	#
	#     Expression for kinematic viscosity of air is taken from program
	#     of Launiainen and Vihma, 1990, Environmental Software, vol. 5,
	#     No. 3, pp. 113 - 124.

	# Heat
	Reynolds = z0Ground * ustar * 1*10**7 / (.9065 * airTemp - 112.7)
	if Reynolds <= 0.135:
		# Smooth coefficients
		b0 = 1.250
		b1 = 0
		b2 = 0
	elif Reynolds < 2.5:
		# Transition coefficients
		b0 = 0.149
		b1 = -0.550
		b2 = 0
	else:
		# Rough coefficients
		b0 = 0.317
		b1 = -0.565
		b2 = -0.183
	z0Groundh = z0Ground * np.exp( b0 + b1 * np.log(Reynolds) + b2 * np.log(Reynolds)**2 )

	# Moisture
	if Reynolds <= 0.135:
		b0=1.61
		b1=0
		b2=0
	elif Reynolds < 2.5:
		b0=0.351
		b1=-0.628
		b2=0
	else:
		b0=0.396
		b1=-0.512
		b2=-0.180
	z0Groundq = z0Ground * np.exp(b0 + b1 * np.log(Reynolds) + b2 * np.log(Reynolds)**2 )

	# Log profiles
	dlogT = np.log(mHeight/z0Groundh)
	dlogQ= np.log(mHeight/z0Groundq)

	return (z0Groundh, z0Groundq, dlogQ, dlogT)
