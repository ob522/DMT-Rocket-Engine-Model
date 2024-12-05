# Physical constants, parameters, or settings

CONSTRAINED_VARIABLES = {
    'exit pressure': 1,    #bar
    'chamber pressure': 30, #bar
    'target thrust': 1,     #kN
    'O/F ratio': 2.3,         
    'contraction ratio': 16,
    'characteristic length': 1,
    # Add other physical parameters as needed
}

NASA_CEA_OUTPUTS = {
    'combustion temperature': 1900,    #K
    'combustion density': 2.9887,    #kg*m-3
    'Cp': 2.1680, #KJ/(KG)(K)
    'gamma': 1.2209,
    # Add other physical parameters as needed
}
"""
NASA_CEA_OUTPUTS = {
    'combustion temperature': 3205.27,    #K
    'combustion density': 2.9018,    #kg*m-3
    'Cp': 4, #KJ/(KG)(K)
    'gamma': 1.14,
    # Add other physical parameters as needed
}
"""
CONICAL_NOZZLE_PARAMETERS = {
    'alpha': 0.26179939,    #rads

}

#FLOW_PARAMETERS = {
#    'viscosity': 0.001,
 #   'density': 1.225,
  #  'inlet_velocity': 15.0,
   # # Add other physical parameters as needed
#}

#TIME_STEP = 0.01  # Time step for the simulation
#TOTAL_TIME = 10.0  # Total simulation time
