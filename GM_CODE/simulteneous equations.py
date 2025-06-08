import numpy as np
from coolant_properties import h_coolant

h_g = 7737 #value from RPA simulation
t_w = 0.002 #thickness
k = 20 #thermal conductivity of steel alloy
h_l = h_coolant
T_g = 2057 #temp of gas as found from RPA
T_l = 335 #in kelvin, assumed for the throat

A = np.array([
    [1/h_l,0,-1],  # First equation coefficients
    [1/h_g, 1, 0],  # Second equation coefficients
    [t_w/k,-1, 1],  # Third equation coefficients
])

B = np.array([-T_l, T_g, 0])  # Constant terms

# Solve the equations
try:
    solution = np.linalg.solve(A, B)
    print("Solution:")
    print("q=", solution[0]*10**-3, 'kW/(m^2*K)')
    print("T_wg =", solution[1], 'K')
    print("T_wl =", solution[2], 'K')

except np.linalg.LinAlgError:
    print("The system of equations has no unique solution.")
