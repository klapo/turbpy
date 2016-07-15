import numpy as np

########
# Parameter Values
ave_slp      =  101325.0     		 # mean sea level pressure              (Pa)
vkc          =       0.4     		 # von Karman constant                  (-)
satvpfrz     =     610.8     		 # sat vapour pressure at 273.16K       (Pa)
w_ratio      =       0.622   		 # molecular ratio water to dry air     (-)
R_da         =     287.053   		 # gas constant for dry air             (Pa K-1 m3 kg-1; J kg-1 K-1)
R_wv         =     461.285   		 # gas constant for water vapor         (Pa K-1 m3 kg-1; J kg-1 K-1)
gravity      =       9.80616 		 # acceleration of gravity              (m s-2)
Cp_air       =    1005.      		 # specific heat of air                 (J kg-1 K-1)
Cp_ice       =    2114.      		 # specific heat of ice                 (J kg-1 K-1)
Cp_soil      =     850.      		 # specific heat of soil                (J kg-1 K-1)
Cp_water     =    4181.      		 # specific heat of liquid water        (J kg-1 K-1)
Tfreeze      =     273.16    		 # temperature at freezing              (K)
TriplPt      =     273.16    		 # triple point of water                (K)
LH_fus       =  333700.0     		 # latent heat of fusion                (J kg-1)
LH_vap       = 2501000.0     		 # latent heat of vaporization          (J kg-1)
LH_sub       = 2834700.0     		 # latent heat of sublimation           (J kg-1)
sb           =       5.6705*10**-8   # Stefan Boltzman constant             (W m-2 K-4)
em_sno       =       0.99    		 # emissivity of snow                   (-)
lambda_air   =       0.026   		 # thermal conductivity of air          (W m-1 K-1)
lambda_ice   =       2.50    		 # thermal conductivity of ice          (W m-1 K-1)
lambda_water =       0.60    		 # thermal conductivity of liquid water (W m-1 K-1)
iden_air     =       1.293   		 # intrinsic density of air             (kg m-3)
iden_ice     =     917.0     		 # intrinsic density of ice             (kg m-3)
iden_water   =    1000.0     		 # intrinsic density of liquid water    (kg m-3)
secprday     =   86400.      		 # number of seconds in a day
secprhour    =    3600.      		 # number of seconds in an hour
secprmin     =      60.      		 # number of seconds in a minute
integerMissing = -9999          	 # value for mising integer
machineEpsilon = (np.finfo(float).eps) # zero value for numerical stability
RiBulk_crit  =       1.47			 # Critical bulk Richardson              (unitless)

########
# Dictionaries for switch-case like behavior in python

# Stability parameters linked to parameter names
stabParams = \
{
    'standard'          : .2,   # critical value for the bulk Richardson number (-)
    'louisInversePower' : 9.4,  # parameter in Louis (1979) stability function
    'mahrtExponential'  : 1.,   # exponential scaling factor in the Mahrt (1987) stability function
    'moninObukhov'      : None  # No tunable parameters
}

# stability Functions linked to needed stability parameters
stabMethods = \
{
    'standard'          : 'critRichNumber',
    'louisInversePower' : 'Louis79_bparam',
    'mahrtExponential'  : 'Mahrt87_eScale'
}
