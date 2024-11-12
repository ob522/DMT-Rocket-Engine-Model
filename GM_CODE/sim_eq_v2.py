import numpy as np
from multi_dimension_channels import h_coolant,h,w

# Given constants
h_g = 7737  # heat transfer coefficient from gas side (W/m^2*K)
t_w = 0.003  # wall thickness (m)
k = 15  # thermal conductivity of steel alloy (W/m*K)
T_g = 2057  # temperature of gas (K)
T_l = 335  # temperature of coolant (K)

# Assume h_coolant is now a matrix of coolant heat transfer coefficients (W/m^2*K)


# Initialize an empty list to store T_wg values
T_wg_matrix = np.zeros(np.shape(h_coolant))

# Iterate through each coolant heat transfer coefficient in the h_l matrix
for i in range(0,len(h)):
    for j in range(0,len(w)):
    # Define the coefficient matrix A and constant vector B for each h_l value
        A = np.array([
            [1/h_coolant[i,j], 0, -1],      # First equation coefficients
            [1/h_g, 1, 0],       # Second equation coefficients
            [t_w/k, -1, 1],      # Third equation coefficients
        ])
        B = np.array([-T_l, T_g, 0])  # Constant terms

        # Solve the system of equations
        try:
            solution = np.linalg.solve(A, B)
            # Extract T_wg (hot wall temperature on the gas side) and append it to the results matrix
            T_wg_matrix[i,j] =solution[1]

        except np.linalg.LinAlgError:
            print("The system of equations has no unique solution for h_l =", h_l)
            T_wg_matrix.append(None)  # Append None for cases where the solution doesn't exist

# Convert the list of T_wg values to a numpy array for easy manipulation
T_wg_matrix = np.array(T_wg_matrix)

print("Matrix of T_wg values for each h_l:")
print(T_wg_matrix)
