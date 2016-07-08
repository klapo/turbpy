def = moninObukhov_SNTHERM(
							mHeight,          # measurement height (m)
							airTemp,          # air temperature (K)
							sfcTemp,          # surface temperature (K)
							windspd,          # wind speed (m s-1)
							# input: stability parameters
							z0Ground          # surface roughness length (below canopy/non-vegetated [snow]) (m)
							):
# Conversion of SNTHERM's (revision 4) qturb.f returns fluxes (not wind functions)
# qlat and qsen. All input values are assumed to be observed 2m above the snow's surface.
# Fluxes are estimated using the Obukhov parameterization, which requires
# an iterative solution to the stability functions and scalar roughness
# parameters. Surface aerodynamic roughness is not an input parameter as
# the Obukhov length method is not sensitive to that value.
#
# SYNTAX:
#
# INPUT:
#	airTemp		= 1x1 scalar of air temperature [K]
#	sfcTemp	= 1x1 scalar of snow surface temperature [K]
#	Eair		= 1x1 scalar of water vapor pressure [hPa]
#	windspd			= 1x1 scalar of wind speed [ms^-1]
#
# OUTPUT:
#	qsen		= 1x1 scalar of the sensible heat flux [Wm^-2]
#	qlat		= 1x1 scalar of the latent heat flux [Wm^-2 <-? check...]
#	Rb			= 1x1 scalar of the bulk Richardson number [-]

## Definitions from fortran version (cause there's lots). Not all variables needed in conversion:
#
# height(ld)	: height above ground of measured met values
# Tkair			: Air temperature [K]
# sfcTempo		: Past surface temperature [K] when icall = 1
#          		: Actual surface temperature [K] when icall = 2
# Ea			: Vapor pressure of water in air [mb]
# Es			: Past vapor pressure at surface [mb] when icall = 1
#          		: Actual vapor pressure at surface [mb] when icall = 2
# windspd			:  Wind speed [m/s]
# windspdo			: Past wind speed [m/s]
# qlat			: Wind function (coefficient) for latent heat exchange [W-m/kg]
# qlat0			: Same as above, except without windless term and not multipled
#       		 	by the latent heat of sublimation
# qsen			: Wind function (coefficient) for sensible heat exchange [W/(m^2 K)]
# qsen0 		: Same as above, except without windless term
# dliqvol		: Volume fraction of liquid water
# Cp_air		: Specific heat of dry air at 0C (1005.0) [J/kg-k]
# R_wv 			: gas constant for water vapor (461.296) [J/kg-K]
# Ck 			: Windless exchange coefficient for water vapor
# Csk 			: Windless exchange coefficient for heat
# z0Ground 		: Roughness length for turbulent transfer of momentum [m]
# Cdn 			: Drag coefficient at neutral stability
# rhoair 		: Density of air  [kg/m^3]
# bpress		: Barometric pressure (mb)
# dlogT 		: ln(temperature height above surface/roughness length for heat)
# dlogQ 		: ln(humidity height above surface/roughness length for water vapor)
# dlogW 		: ln(wind helght above surface/roughness length for momentum)
# wstar 		: friction velocity  [m/s]
# Rb 			: Bulk Richardson number
# Hw			: Height above snowpack for wind speed, RH, and temperature measurements [m]
# Next 3 variables are the integrated stability corrections, all computed using
#   actual measurement heights - These should end up being equal
# Psim			: Correction for momentum
# Psiq			: Correction for water vapor
# Psih			: Correction for heat
# z0Groundh 			: Roughness length for turbulent transfer of heat [m]
# z0Groundq 			: Roughness length for turbulent transfer of water vapor [m]

## Constants and derived values
z0Groundh = z0Ground						# Assume a water vapor flux characteristic length of .005m
z0Groundq = z0Ground						# Assume a temperature flux characteristic length of .005m
dlogT = log(mHeight/z0Groundh)
dlogQ = log(mHeight/z0Groundq)

# Constants should be imported from multiConst
#R_wv = 461.296								# Gas constant for water vapor [J/kg-K]
Ck = 0										# Windless exchange coefficient for water vapor [Wm^-2]
Csk = 2										# Windless exchange coefficient for heat [Wm^-2]
#Cp_air = 1005								# Specific heat of dy air [J kg^-1 K^-1]

# Computes drag coefficient over snow at neutral stability
dlogW = log(mHeight/z0Ground);

# Initial Values: assume neutral stability for first guess
Psim = 0
Psih = 0
Psiq = 0
deltaT = airTemp - sfcTemp
deltaE = Eair - Esurface
Tbar = (airTemp + sfcTemp)/2
Rbw = 9.8 .* mHeight * deltaT/(Tbar*windspd**2)
zetat = Rbw * dlogW
L = mHeight/zetat
# Gas constant for dry air (287J/kg-K) is from "Handbook of
# Meteorology,"Berry, et al, 1945, McGraw-Hill, p. 353.
rhoair = 100*800/(287*airTemp)				# I'm assuming a higher elevation site...
con = 100./(R_wv*airTemp*rhoair)

## Code
# Compute stability functions using iterative solution.
# Next begins Newton-Raphson iterative loop.  Loop at least twice.
for it in np.arange(0:numMaxIterations)
	if it == 0:
	# Initial pass
		if Rbw > 1.427:
			# Exceeds critical Richardson number
			qsen0= machineEpsilon
			qlat0= machineEpsilon
			return()
 		if L >= 0:
			# Bounds on zetat for N-R scheme
			zetal = 0
			zetah = 500

	## a. Compute momentum stability functions, Cd, and wstar
	zeta = mHeight/L
	if L >= 0:
		# Stable
		Psim = Dutch(zeta)
	else:
		# Unstable
		Psim = Unstable(L,mHeight,z0Groundh,1)

	sqrtCd = 0.4/np.max((dlogW - Psim),0.1)
	wstar = sqrtCd * windspd;

	## b. Estimate scalar roughness lengths and stability functions
	[z0Groundh,z0Groundq,dlogQ,dlogT] = andreas(wstar,z0Ground,airTemp,mHeight);

	if L >= 0:
		# Stable
		Psih=Psim;
		PsimT=Psih;
		Psiq=Psih;
	else:
		# Unstable
		PsimT=Unstable(L,mHeight,z0Groundh,1)
		Psih=Unstable(L,mHeight,z0Groundh,2)
		Psiq=Psih
	tstar = 0.4 * deltaT/(dlogT - Psih) # 0.4 = von Karman constant? If so replace it
	qstar = 0.4 * deltaE*con/ (dlogQ - Psiq)

	## c. Compute Monin-Obukhov length (L)
	if tstar ~= 0 & wstar ~= 0:
		# This code needs desperate help, its hard as fuck to follow
		# Why is the representative length scale for moisture used to calculate the Obukhov length?
		L = 1 + 0.61 .* Tbar*qstar/tstar
		dum = L
		L = L * 9.8 * 0.4 * tstar / (Tbar * wstar**2)
		L = 1/L
	elif (tstar == 0) & (wstar ~= 0):
		L=1*10**14

	## d. Limits for extremely stable and unstable stratification
	freelim = 0
	freelimv = 0
	if (zetat > 500) | (mHeight/L > 500):
		# Very stable
		qsen0 = machineEpsilon
		qlat0 = machineEpsilon
		break

	# What is happening here?
	# elseif L < 0 && L > -1.0

	## e. Compute change in zeta
	y = (mHeight/L) - zetat
	if L > 0 &  y >  0:
		zetal = zetat									# Update lower bound
	elif L > 0 & y <= 0:
		zetah=zetat										# Update upper bound

	if np.abs(y) < 1*10**(-2) & it > 2:
		# Within tolerance and more than one iteration executed, leave loop
		qlat0 = fqlat(dlogW,dlogQ,R_wv,airTemp,Psim,Psiq,windspd);
		qsen0 = fqsen(dlogW,dlogT,Cp_air,rhoair,Psim,Psih,windspd);
		return()
	if it > 50:
		# NOT within tolerance and more than 50 iterations executed, leave loop
		qlat0 = fqlat(dlogW,dlogQ,R_wv,airTemp,Psim,Psiq,windspd);
		qsen0 = fqsen(dlogW,dlogT,Cp_air,rhoair,Psim,Psih,windspd);
		# Raise some sort of error here
		disp('Convergence problem in QTURB')
		return

	## f. Use Newton-Raphson scheme for stable, sequential estimates for unstable
	# Change exist logic here...
	if L >= 0 & ~exist('zetah')
		disp('whoops')

	elif L >= 0:
		# Next is derivative of Psi
		DpsitDzetat = -0.70 - 0.75 * exp(-0.35 * zetat)*(6 - 0.35 * zetat)
		DpsimDzetat = -0.70 - 0.75 * exp(-0.35 * zeta)*(6 - 0.35 * zeta)
		dumM = 1 / (dlogW-Psim)
		dumT = 1 / (dlogT-Psih)	# SNTHERM has dlogt, not dlogT, even though it appears no where else...
		# Next is derivative of y (ignoring humidity effects)
		yprime = (mHeight/L) * (DpsitDzetat * dumT - 2 * DpsimDzetat * dumM) - 1
		zetat = zetat - (y/yprime)
		if (zetat <= zetal) | (zetat >= zetah):
			# Root out of range. Bisection instead of Newton-Raphson scheme
			zetat = (zetah+zetal)/2
	else:
		zetat = mHeight/L
	# Obukhov length
	L = mHeight/zetat

## Final computation of qsen and qsen.
# Note qlat and qsen add the windless coefficients
qlat = Ck + qlat0 * 2.838*10**6
if L >= 0:
	qsen=qsen0
	if (qsen0 < np.max(Csk,.1)) & (qsen0 > 0):
		# Minimum value
		qsen = np.max(Csk,.1)
else:
	dum = 0
	if not (deltaT == 0):
		dum = -freelim/deltaT
	qsen = np.max(qsen0,dum)

	dum = 0
	if not (deltaE == 0):
		dum = -freelimv/deltaE
	qlat = np.max(qlat,dum)

## Compute bulk Richardson number
# First adjust windspeed to height of temperature measurement
if  L > 0:
	PsimT = Psih
else:
	PsimT = Unstable(L,mHeight,z0Groundh,1)
windspdT = (windspd./ (dlogW - Psim)) .* (log(mHeight/z0Ground)-PsimT)
Rb = 9.8 .* mHeight .* deltaT ./ (Tbar*windspdT**2)

########
# Sub-functions (stability schemes)
# --------------------------------------------------------------------------------------------------------------------
def fqsen(dlogW,dlogT,Cp_air,rhoair,Psim,Psih,windspd):
	sen = Cp_air*rhoair*windspd*0.16/((dlogW - Psim)*(dlogT - Psih));
# 	sen = 1/((dlogW - Psim)*(dlogT - Psih));

def fqlat(dlogW,dlogQ,R_wv,tkair,psim,Psiq,windspd):
# Since vapor pressure is in mb, fqlat is multiplied by 100 to
# convert to Pascals
	# Wind coefficient
	lat = windspd * 0.16*10**2 / (R_wv*tkair*(dlogW - psim)*(dlogQ - Psiq))
# 	lat = 1/(R_wv*tkair*(dlogW - psim)*(dlogQ - Psiq))

def Dutch(zeta):
	stab = -(0.70.* zeta + 0.75 .*(zeta-14.28) .* exp(-0.35*zeta)  + 10.71);

def andreas(wstar,z0Ground,Tkair,mHeight):
	#    Compute scaler roughness lengths using procedure
	#    in Andreas, 1987, Boundary-Layer Meteorology 38, 159-184.  Only
	#    use this procedure for snow.

	#     Expression for kinematic viscosity of air is taken from program
	#     of Launiainen and Vihma, 1990, Environmental Software, vol. 5,
	#     No. 3, pp. 113 - 124.

	# Heat
	Reynolds = z0Ground * wstar * 1*10**7 / (.9065 * Tkair - 112.7)
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
	dlogReynolds = np.log(Reynolds)
	dlogsqrd = dlogReynolds * dlogReynolds
	dum= b0 + b1 * dlogReynolds + b2 * dlogsqrd
	z0Groundh = z0Ground * np.exp(dum)

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
	dum = b0 + b1 * dlogReynolds + b2*dlogsqrd;
	z0Groundq = z0Ground*exp(dum);

	# Log profiles
	dlogT = log(mHeight/z0Groundh);
	dlogQ= log(mHeight/z0Groundq);

	return (z0Groundh, z0Groundq, dlogQ, dlogT)
