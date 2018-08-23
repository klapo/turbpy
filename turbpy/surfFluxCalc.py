def moninObukhov(airTemp, airVaporPress, sfcTemp, sfcVaporPress,
                 stabilityCorrectionParameters, senHeatGround0, latHeatGround0,
                 conductanceSensible, conductanceLatent,
                 volHeatCapacityAir, latHeatSubVapGround,
                 latentHeatConstant, cap=True):

    # Windless exchange coefficients for water vapor, k, and heat, sk [Wm^-2]
    Ck = 0
    Csk = 2

    # Unpack parameters from aStability.moninObukhov
    L = stabilityCorrectionParameters['L']
    freelim = stabilityCorrectionParameters['freelim']
    freelimv = stabilityCorrectionParameters['freelimv']

    # Gradients in temperature and water vapor
    deltaT = airTemp - sfcTemp
    deltaE = airVaporPress - sfcVaporPress

    ########
    # Windless Exchange -- Calculate sensible and latent heat fluxes.
    if cap == 'windless_exchange':
        # Note that the fluxes include the windless exchange coefficients.
        latHeatGround = Ck + latHeatGround0
        if L >= 0:
            senHeatGround = senHeatGround0
            if senHeatGround0 < max(Csk, .1) and senHeatGround0 > 0:
                # Enforce a minimum value.
                senHeatGround = max(Csk, .1)
        else:
            # Sensible heat maximum for unstable
            dum = 0
            if not deltaT == 0:
                dum = -freelim / deltaT
            senHeatGround = max(senHeatGround0, dum)

            # Latent heat maximum for unstable
            dum = 0
            if not deltaE == 0:
                dum = -freelimv / deltaE
            latHeatGround = max(latHeatGround, dum)

    ########
    # NO Windless Exchange -- Calculate sensible and latent heat fluxes.
    else:
        # Note that the fluxes include the windless exchange coefficients.
        latHeatGround = Ck + latHeatGround0
        if L >= 0:
            senHeatGround = senHeatGround0
        else:
            # Sensible heat maximum for unstable
            dum = 0
            if not deltaT == 0:
                dum = -freelim / deltaT
            senHeatGround = max(senHeatGround0, dum)

            # Latent heat maximum for unstable
            dum = 0
            if not deltaE == 0:
                dum = -freelimv / deltaE
            latHeatGround = max(latHeatGround, dum)

    ########
    # Back out new conductance
    conductanceSensible = senHeatGround / volHeatCapacityAir
    conductanceLatent = latHeatGround / (latentHeatConstant
                                         * conductanceLatent)

    ########
    # Convert to fluxes
    senHeatGround = senHeatGround * deltaT
    latHeatGround = latHeatGround * deltaE

    ########
    # Return turbulent fluxes and newly calculated conductances
    return (senHeatGround, latHeatGround,
            conductanceSensible, conductanceLatent)
