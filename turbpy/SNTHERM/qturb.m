function [qsen, qlat, Rb] = qturb(Tair,Tsurface,Eair,wsp)
% Conversion of SNTHERM's (revision 4) qturb.f returns fluxes (not wind functions)
% qlat and qsen. All input values are assumed to be observed 2m above the snow's surface.
% Fluxes are estimated using the Obukhov parameterization, which requires
% an iterative solution to the stability functions and scalar roughness
% parameters. Surface aerodynamic roughness is not an input parameter as 
% the Obukhov length method is not sensitive to that value.
%
% SYNTAX:
%	[qsen, qlat, Rb] = TurbulentFluxes_SNTHERM(Tair,Tsurface,Eair,wsp)
%
% INPUT:
%	Tair		= 1x1 scalar of air temperature [K]
%	Tsurface	= 1x1 scalar of snow surface temperature [K]
%	Eair		= 1x1 scalar of water vapor pressure [hPa] 
%	wsp			= 1x1 scalar of wind speed [ms^-1]
%
% OUTPUT:
%	qsen		= 1x1 scalar of the sensible heat flux [Wm^-2]
%	qlat		= 1x1 scalar of the latent heat flux [Wm^-2 <-? check...]
%	Rb			= 1x1 scalar of the bulk Richardson number [-]
%
% DEPENDENCIES:
%	fvapri_SNTHERM.m


%% Definitions from fortran version (cause there's lots). Not all variables needed in conversion:
%
% height(ld)	: height above ground of measured met values
% Tkair			: Air temperature [K]
% Tsurfaceo		: Past surface temperature [K] when icall = 1
%          		: Actual surface temperature [K] when icall = 2
% Ea			: Vapor pressure of water in air [mb]
% Es			: Past vapor pressure at surface [mb] when icall = 1
%          		: Actual vapor pressure at surface [mb] when icall = 2
% wsp			:  Wind speed [m/s]
% wspo			: Past wind speed [m/s]
% qlat			: Wind function (coefficient) for latent heat exchange [W-m/kg]
% qlat0			: Same as above, except without windless term and not multipled 
%       		 	by the latent heat of sublimation
% qsen			: Wind function (coefficient) for sensible heat exchange [W/(m^2 K)]
% qsen0 		: Same as above, except without windless term 
% dliqvol		: Volume fraction of liquid water
% cdryair		: Specific heat of dry air at 0C (1005.0) [J/kg-k]
% Rw 			: gas constant for water vapor (461.296) [J/kg-K]
% Ck 			: Windless exchange coefficient for water vapor
% Csk 			: Windless exchange coefficient for heat
% z0 			: Roughness length for turbulent transfer of momentum [m]
% Cdn 			: Drag coefficient at neutral stability
% rhoair 		: Density of air  [kg/m^3]
% bpress		: Barometric pressure (mb)
% dlogT 		: ln(temperature height above surface/roughness length for heat)
% dlogQ 		: ln(humidity height above surface/roughness length for water vapor)
% dlogW 		: ln(wind helght above surface/roughness length for momentum)
% wstar 		: friction velocity  [m/s]
% icall 		: =1 for call to estimate qsen and qlat before heat balance
%         		  =2 for call after heat balance, to finalize transfer coefficients
% Rb 			: Bulk Richardson number
% Hw			: Height above snowpack for wind speed, RH, and temperature measurements [m]
% Next 3 variables are the integrated stability corrections, all computed using
%   actual measurement heights - These should end up being equal
% Psim			: Correction for momentum
% Psiq			: Correction for water vapor
% Psih			: Correction for heat
% z0h 			: Roughness length for turbulent transfer of heat [m]
% z0q 			: Roughness length for turbulent transfer of water vapor [m]
%

%% Checks
if numel(Tair) > 1 || numel(Tsurface) > 1 || numel(Eair) > 1 || numel(wsp) > 1
	error('This function requires scalar arguments')
end
% Below statement breaks UEB since it searches for the surface temperature
% over a laughably large range.
% if max(Tair) < 200 || max(Tsurface) < 200
% 	error('Air and surface temperature should be in Kelvin')
% end

%% Constants and derived values
z0 = .005;									% Assume a momentum roughness length of .005m
z0h = z0;									% Assume a water vapor flux characteristic length of .005m
z0q = z0;									% Assume a temperature flux characteristic length of .005m
H = 2;										% Assume measurement 2m above snow
dlogT = log(H/z0h);
dlogQ = log(H/z0q);
Rw = 461.296;								% Gas constant for water vapor [J/kg-K]
Ck = 0;										% Windless exchange coefficient for water vapor [Wm^-2]
Csk = 2;									% Windless exchange coefficient for heat [Wm^-2]
cdryair = 1005;								% Specific heat of dy air [J kg^-1 K^-1]

% Water vapor partial pressure (hPa)
Esurface = fvapri_SNTHERM(Tsurface,1);

% Computes drag coefficient over snow at neutral stability
dlogW = log(H./z0);

% Use values at neutral stability for first guess 
Psim = 0;
Psih = 0;
Psiq = 0;
deltaT = Tair - Tsurface;
deltaE = Eair - Esurface;
Tbar = (Tair + Tsurface)./2;
Rbw = 9.8 .* H .* deltaT./(Tbar.*wsp.^2);
zetat = Rbw * dlogW;
L = H/zetat;
% Gas constant for dry air (287J/kg-K) is from "Handbook of 
% Meteorology,"Berry, et al, 1945, McGraw-Hill, p. 353.
rhoair = 100.*800/(287.*Tair);				% I'm assuming a higher elevation site...
con = 100./(Rw.*Tair.*rhoair);

%% Code 
% Compute stability functions using iterative solution.
% Next begins Newton-Raphson iterative loop.  Loop at least twice.
for it = 0:51
	if it == 0
	%% Initial pass
		if Rbw > 1.427
			% Exceeds critical Richardson number
			qsen0= 1.*10.^-12;
			qlat0= 1.*10.^-12;
			break
		end
% 		if L >= 0
			% Bounds on zetat for N-R scheme
			zetal = 0;
			zetah = 500;
% 		end
	end

	%% a. Compute momentum stability functions, Cd, and wstar
	zeta = H/L;
	if L >= 0
		% Stable   
		Psim = Dutch(zeta);
	else
		% Unstable 
		Psim = Unstable(L,H,z0h,1);
	end

	sqrtCd = 0.4./max((dlogW - Psim),0.1);
	wstar = sqrtCd .* wsp;

	%% b. Estimate scalar roughness lengths and stability functions
	[z0h,z0q,dlogQ,dlogT] = andreas(wstar,z0,Tair,H);

	if L >= 0
		% Stable   
		Psih=Psim;
		PsimT=Psih;
		Psiq=Psih;
	else
		% Unstable 
		PsimT=Unstable(L,H,z0h,1);
		Psih=Unstable(L,H,z0h,2);
		Psiq=Psih;
	end 
	tstar = 0.4 .* deltaT./(dlogT - Psih);
	qstar = 0.4 .* deltaE.*con./ (dlogQ - Psiq);

	%% c. Compute Monin-Obukhov length (L)
	if tstar ~= 0 && wstar ~= 0 
		L = 1 + 0.61 .* Tbar.*qstar./tstar;
		dum=L;
		L = L .* 9.8 .* 0.4 .* tstar./ (Tbar.* wstar.^2);
		L=1/L;
	elseif tstar == 0 && wstar ~= 0
		L=1.*10.^14;
	end

	%% d. Limits for extremely stable and unstable stratification 
	freelim = 0;
	freelimv = 0;
	if zetat > 500 || H/L > 500
		% Very stable
		qsen0 = 1.*10.^-12;
		qlat0 = 1.*10.^-12;
		break	

	elseif L < 0 && L > -1.0
		% Very unstable:
		% Compute free convection limits.  Applies generally for L < -0.1
		% Andreas and Cash, Convective heat transfer over wintertime
		% leads and polynyas, JGR 104, C11, 25,721-25,734, 1999.
		% Method is for a smooth surface with unlimited fetch.
		deltaB = (-9.8 .* deltaT ./ Tbar) .* dum;
		nu = 1.315.*10.^-5; 							% Kinematic viscosity of air
		D = 0.024 ./ (cdryair .* rhoair);
		Dv = 1.1494*D;
		zscale = (nu.*D./deltaB).^(1/3);
		zscalev = (nu*Dv/deltaB).^(1/3);
		freelim = -0.15 .*rhoair .* cdryair .* D .* deltaT ./zscale;
		freelimv = -0.15 .* Dv .* deltaE .* con .* 2.838.*10.^6./zscalev;
	end
	%% e. Compute change in zeta
	y = (H./L) - zetat;
	if L > 0 &&  y >  0
		zetal = zetat;									% Update lower bound
	elseif L > 0 && y <= 0
		zetah=zetat;									% Update upper bound
	end

	if abs(y) < 1.*10.^-2 && it > 2
		% Within tolerance and more than one iteration executed, leave loop
		qlat0 = fqlat(dlogW,dlogQ,Rw,Tair,Psim,Psiq,wsp);
		qsen0 = fqsen(dlogW,dlogT,cdryair,rhoair,Psim,Psih,wsp);
		break
	end
	if it > 50
		% NOT within tolerance and more than 50 iterations executed, leave loop
		qlat0 = fqlat(dlogW,dlogQ,Rw,Tair,Psim,Psiq,wsp);
		qsen0 = fqsen(dlogW,dlogT,cdryair,rhoair,Psim,Psih,wsp);
		disp('Convergence problem in QTURB') 
		break
	end

	%% f. Use Newton-Raphson scheme for stable, sequential estimates for unstable
	if L >= 0 && ~exist('zetah')
		disp('whoops')
	
	elseif L >= 0
		% Next is derivative of Psi
		DpsitDzetat = -0.70 - 0.75.* exp(-0.35 .* zetat).*(6 - 0.35 .* zetat);
		DpsimDzetat = -0.70 - 0.75 .* exp(-0.35 .* zeta).*(6 - 0.35 .* zeta);
		dumM = 1./ (dlogW-Psim);
		dumT = 1 ./ (dlogT-Psih);								% SNTHERM has dlogt, not dlogT, even though it appears no where else...
		% Next is derivative of y (ignoring humidity effects)
		yprime = (H./L) .* (DpsitDzetat .* dumT - 2 .* DpsimDzetat .* dumM) - 1;
		zetat = zetat - (y./yprime);
		if zetat <= zetal || zetat >= zetah
			% Root out of range. Bisection instead of Newton-Raphson scheme
			zetat = (zetah+zetal)/2;
		end
	else
		zetat = H/L;
	end
	L = H/zetat;
end

%% Final computation of qsen and qsen. 
% Note qlat and qsen add the windless coefficients
qlat = Ck + qlat0 .*2.838.*10.^6;
if L >= 0
	
	qsen=qsen0;

	if qsen0 < max(Csk,.1) && qsen0 > 0
		% Minimum value
		qsen = max(Csk,.1);
	end
else
	dum = 0;
	if deltaT ~= 0
		dum=-freelim/deltaT;
	end
	qsen = max(qsen0,dum);
	
	dum = 0;
	if deltaE ~= 0
		dum = -freelimv/deltaE;
	end
	qlat = max(qlat,dum);
end

%% Compute bulk Richardson number
% First adjust windspeed to height of temperature measurement
if  L > 0 
	PsimT = Psih;
else
	PsimT = Unstable(L,H,z0h,1);
end
wspT = (wsp./ (dlogW - Psim)) .* (log(H/z0)-PsimT);
Rb = 9.8 .* H .* deltaT ./ (Tbar.*wspT.^2);

%% Convert to fluxes
qsen = qsen .* (Tair - Tsurface);
qlat = qlat .* (Eair - Esurface);
end

%***********************************************************************
%% Sub-funcitons
%***********************************************************************

function sen = fqsen(dlogW,dlogT,cdryair,rhoair,Psim,Psih,wsp)
	sen = cdryair*rhoair*wsp*0.16/((dlogW - Psim)*(dlogT - Psih));
% 	sen = 1/((dlogW - Psim)*(dlogT - Psih));
end  

function lat = fqlat(dlogW,dlogQ,Rw,tkair,psim,Psiq,wsp)
% Since vapor pressure is in mb, fqlat is multiplied by 100 to
% convert to Pascals
	% Wind coefficient
	lat = wsp*0.16.*10.^2/(Rw*tkair*(dlogW - psim)*(dlogQ - Psiq));
% 	lat = 1/(Rw*tkair*(dlogW - psim)*(dlogQ - Psiq));
end 

function stab = Dutch(zeta)
	stab = -(0.70.* zeta + 0.75 .*(zeta-14.28) .* exp(-0.35.*zeta)  + 10.71);
end 

function stab = Unstable(L,Ht,zt,icall)
% Set an upper limit on function at L = -.1 or  at 200 times
% the large -Ht/L limit where log(-4*Ht/L) is 0.  This asymptotic
% limit is from Godfrey and Beljaars (1991).  This is a semi-
% empirical fix to keep the heat flux from increasing with 
% decreasing -L.
Llim = min(-.1,-100.*zt);
zeta = Ht / min(L,Llim);
x = (1-16.*zeta).^0.25;
if icall <= 1      
	stab = log((1+x*x)/2)+2 .* log((1d0+x)/2) - 2.*atan(x) + 1.5707963;
else      
	stab = 2.*log((1+x*x)/2);
end
end      

function [z0h, z0q, dlogQ, dlogT] = andreas(wstar,z0,Tkair,H)
%    Compute scaler roughness lengths using procedure
%    in Andreas, 1987, Boundary-Layer Meteorology 38, 159-184.  Only
%    use this procedure for snow.

%     Expression for kinematic viscosity of air is taken from program
%     of Launiainen and Vihma, 1990, Environmental Software, vol. 5,
%     No. 3, pp. 113 - 124.

% Heat
Reynolds = z0.*wstar.*1.*10.^7./(.9065.*Tkair-112.7);
if Reynolds <= 0.135
	% Smooth coefficients
	b0=1.250;
	b1=0;
	b2=0;
elseif Reynolds < 2.5
	% Transition coefficients
	b0=0.149;
	b1=-0.550;
	b2=0;
else
	% Rough coefficients
	b0=0.317;
	b1=-0.565;
	b2=-0.183;
end
dlogReynolds = log(Reynolds);
dlogsqrd = dlogReynolds .* dlogReynolds;
dum= b0 + b1.*dlogReynolds + b2.*dlogsqrd;
z0h = z0.*exp(dum);

% Moisture 
if Reynolds <= 0.135
	b0=1.61;
	b1=0;
	b2=0;
elseif Reynolds < 2.5
	b0=0.351;
	b1=-0.628;
	b2=0;
else
	b0=0.396;
	b1=-0.512;
	b2=-0.180;
end
dum = b0 + b1 * dlogReynolds + b2*dlogsqrd;
z0q = z0*exp(dum);

% Log profiles
dlogT = log(H/z0h);
dlogQ= log(H/z0q);
end
