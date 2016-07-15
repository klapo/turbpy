def surfFluxCalc_bulkAero(   # Inputs
                    groundTemp,
                    airTemp,
                    groundConductanceSH,
                    vpGround,
                    vpAir,
                    groundConductanceLH,
                    ):

    senHeatGround      = -volHeatCapacityAir * groundConductanceSH * (groundTemp - airTemp)
    latHeatGround      = -latHeatSubVapGround * latentHeatConstant * groundConductanceLH \
                          * (vpGround - vpAir)
                              
    return (senHeatGround, latHeatGround)
