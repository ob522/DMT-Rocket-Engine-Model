# Functions to load or generate data

def assign_constants(constants):
    """
    Extracts relevant constants from the provided dictionary structure.
    
    Args:
        constants (module or dictionary): The constants file or dictionary containing the necessary parameters.
    
    Returns:
        A tuple containing the extracted constants in the appropriate units.
    """
    
    # Extract variables from CONSTRAINED_VARIABLES
    Pe = constants.CONSTRAINED_VARIABLES['exit pressure'] * (10**5)    # Exit pressure (Pa)
    P0 = constants.CONSTRAINED_VARIABLES['chamber pressure'] * (10**5)  # Chamber pressure (Pa)
    F = constants.CONSTRAINED_VARIABLES['target thrust'] * (10**3)     # Thrust (N)
    CR = constants.CONSTRAINED_VARIABLES['contraction ratio']          # Contraction ratio
    char_L = constants.CONSTRAINED_VARIABLES['characteristic length']  # Characteristic length (m)

    # Extract variables from NASA_CEA_OUTPUTS
    T0 = constants.NASA_CEA_OUTPUTS['combustion temperature']          # Combustion temperature (K)
    rho0 = constants.NASA_CEA_OUTPUTS['combustion density']            # Combustion density (kg/m^3)
    Cp = constants.NASA_CEA_OUTPUTS['Cp'] * (10**3)                    # Specific heat capacity (J/kg/K)
    gamma = constants.NASA_CEA_OUTPUTS['gamma']                        # Specific heat ratio
    R = Cp - Cp / gamma                                                # Gas constant (J/kg/K)

    # Extract variables from CONICAL_NOZZLE_PARAMETERS
    alpha = constants.CONICAL_NOZZLE_PARAMETERS['alpha']               # Nozzle half-angle (radians)

    return (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha)
