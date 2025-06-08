import math
from scipy.optimize import fsolve

def residual_M1(M1_guess, beta1, beta2, gamma):
    """
    Computes the residual function that should be zero when M1 is correct.
    
    Parameters:
    - M1_guess: Current guess for M1 (scalar)
    - beta1: measured shock angle upstream (in radians)
    - beta2: measured shock angle downstream (in radians)
    - gamma: ratio of specific heats (e.g., 1.4 for air)
    
    Returns:
    - residual: Difference between tan(Theta) computed from Equation 1 and Equation 2
    """
    
    # Equation 1: Compute tan(Theta) based on upstream shock angle beta1 and guess for M1
    # Numerator of Equation 1
    numerator1 = 2 * (1 / math.tan(beta1)) * (M1_guess**2 * math.sin(beta1)**2 - 1)
    
    # Denominator of Equation 1
    denominator1 = M1_guess**2 * (gamma + math.cos(2 * beta1)) + 2
    
    # Calculate tan(Theta) using Equation 1
    tan_theta1 = numerator1 / denominator1
    
    # Calculate Theta (shock deflection angle) in radians
    theta = math.atan(tan_theta1)
    
    # Numerator of the right-hand side of Equation 3:
    numerator_eq3_rhs = 2 + (gamma - 1) * (M1_guess * math.sin(beta1))**2
    
    # Denominator of the right-hand side of Equation 3:
    denominator_eq3_rhs = 2 * gamma * (M1_guess * math.sin(beta1))**2 - (gamma - 1)
    
    # Calculate the right-hand side value:
    # This represents (M2 * sin(beta1 - Theta))^2, which we will use to find M2.
    M2_squared_sin_sq = numerator_eq3_rhs / denominator_eq3_rhs
    
    # Check if result is non-negative before taking the square root
    if M2_squared_sin_sq <= 0:
        return 1e6  # Return a large residual if non-physical solution (avoid math error)
    
    # Calculate M2 using (M2 * sin(beta1 - Theta))^2
    M2 = math.sqrt(M2_squared_sin_sq) / abs(math.sin(beta1 - theta))
    
    # Equation 2: Compute tan(Theta) based on downstream shock angle beta2 and computed M2
    numerator2 = 2 * (1 / math.tan(beta2)) * (M2**2 * math.sin(beta2)**2 - 1)
    denominator2 = M2**2 * (gamma + math.cos(2 * beta2)) + 2
    
    tan_theta2 = numerator2 / denominator2
    
    # Residual: Difference between tan(Theta) values from Equation 1 and Equation 2
    residual = tan_theta1 - tan_theta2
    
    return residual

def solve_M1(beta1_deg, beta2_deg, gamma, M1_guess=3.0):
    """
    Solves for M1 using measured shock angles (beta1, beta2) and gamma.
    
    Parameters:
    - beta1_deg: measured upstream shock angle in degrees
    - beta2_deg: measured downstream shock angle in degrees
    - gamma: ratio of specific heats (e.g., 1.4 for air)
    - M1_guess: initial guess for M1 (default is 2.0)
    
    Returns:
    - M1_solution: value of M1 that satisfies all three equations
    """
    
    # Convert measured shock angles from degrees to radians
    beta1_rad = math.radians(beta1_deg)
    beta2_rad = math.radians(beta2_deg)
    
    # Use scipy.optimize.fsolve to find the root of the residual function
    solution = fsolve(residual_M1, M1_guess, args=(beta1_rad, beta2_rad, gamma))
    
    # Return the first (and only) solution found by fsolve
    return solution[0]

from scipy.optimize import root_scalar

def solve_M1_bracket(beta1_deg, beta2_deg, gamma, M1_lower=1.0, M1_upper=10.0):
    """
    Uses a robust bracketing solver (root_scalar) to find M1.
    Requires a bracket [M1_lower, M1_upper] where the solution is expected.
    """
    beta1_rad = math.radians(beta1_deg)
    beta2_rad = math.radians(beta2_deg)
    
    # Define a lambda that takes a single argument (M1) as required by root_scalar
    func = lambda M1: residual_M1(M1, beta1_rad, beta2_rad, gamma)
    
    # Use Brent’s method (good mix of speed and reliability)
    result = root_scalar(func, bracket=[M1_lower, M1_upper], method='brentq')
    
    if result.converged:
        return result.root
    else:
        raise ValueError("Root finding did not converge.")


# Example usage:
beta1_deg = 80/2#70/2  # Example upstream shock angle in degrees 92
beta2_deg = 70/2 #80/2  # Example downstream shock angle in degrees   67
gamma = 1.26       # Ratio of specific heats for air

# # Solve for M1
# M1 = solve_M1(beta1_deg, beta2_deg, gamma)
# print(f"Solved M1: {M1:.4f}")

try:
    # Solve for M1 using the robust bracketing method
    M1_solution = solve_M1_bracket(beta1_deg, beta2_deg, gamma, M1_lower=1.1, M1_upper=10.0)
    print(f"Solved M1: {M1_solution:.4f}")
except ValueError as e:
    print(f"Error: {e}")
