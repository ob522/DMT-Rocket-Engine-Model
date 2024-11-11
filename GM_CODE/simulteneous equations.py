import numpy as np
from coolant_properties import h_coolant

h_g = 10000 #value from Bartz Equation solved by OB code
t_w = 0.002 #thickn
k = 125 #thermal conductivity of copper alloy
h_l = h_coolant
T_g = 2200 #temp of gas as found from RPA

A = np.array([
    [1/h_g+t_w/k+1/h_l, 1, 0, 0],  # First equation coefficients
    [1/h_g, 0, 1, 0],  # Second equation coefficients
    [t_w/k, 0, -1, 1],  # Third equation coefficients
    [1/h_l, 1, 0, -1]   # Fourth equation coefficients
])

B = np.array([T_g, T_g, 0, 0])  # Constant terms

# Solve the equations
try:
    solution = np.linalg.solve(A, B)
    print("Solution:")
    print("q=", solution[0])
    print("T_L =", solution[1])
    print("T_wg =", solution[2])
    print("T_wl =", solution[3])

except np.linalg.LinAlgError:
    print("The system of equations has no unique solution.")
