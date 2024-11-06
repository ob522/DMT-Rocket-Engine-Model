# testing and debugging code

import numpy as np

#x0 = 0
#xi = 10
#dx_0i = 1

#x_0i = np.arange(x0, xi + dx_0i, dx_0i)  # From x0 to xi
#print(x_0i)


'''
# define key points along nozzle length
x0 = 0
xi = 10
xt = 10 + (3-2)/(np.tan(0.26))
xe = xt + (5-2)/(np.tan(0.26))

# set x steps
dx_0i = 0.5   # Step between x0 and xi
dx_it = 0.1   # Step between xi and xt
dx_te = 0.2  # Step between xt and xe

# define x array (with given steps)
x_0i = np.arange(x0, xi, dx_0i)  # From x0 to (xi - dx_0i)
x_it = np.arange(xi, xt, dx_it)  # From xi to (xt - dx_it)
x_te = np.arange(xt, xe, dx_te)  # From xt to (xe - dx_te)

x_values = np.concatenate((x_0i, x_it, x_te))  # Concatenate all the segments into a single array

# Add the final value xe (since arange may not include it)  (no longer needed)
if x_values[-1] < xe:
    x_values = np.append(x_values, xe)


r_values = np.zeros_like(x_values)

# define conical radius function (piecewise function)
r_values[:len(x_0i)] = 1    # For x0 to xi
r_values[len(x_0i):len(x_0i) + len(x_it)] = 2
r_values[len(x_0i) + len(x_it):] = 3

print(x_values)
print(r_values)

'''

def conical_radius_function( Lc, alpha):
    # define key radii along nozzle length
    r0 = 2
    rt = 1
    re = 3
    
    # define key points along nozzle length
    x0 = 0
    xi = Lc
    xt = Lc + (r0-rt)/(np.tan(alpha))
    xe = xt + (re-rt)/(np.tan(alpha))

    # set x steps
    dx_0i = 0.5   # Step between x0 and xi
    dx_it = 0.1   # Step between xi and xt
    dx_te = 0.1   # Step between xt and xe

    # define x array (with given steps)
    x_0i = np.arange(x0, xi, dx_0i)  # From x0 to (xi - dx_0i)
    x_it = np.arange(xi, xt, dx_it)  # From xi to (xt - dx_it)
    x_te = np.arange(xt, xe, dx_te)  # From xt to (xe - dx_te)
    
    # Add the final value xe (since arange may not include it)
    if x_te[-1] < xe:
        x_te = np.append(x_te, xe)

    x_values = np.concatenate((x_0i, x_it, x_te))  # Concatenate all the segments into a single array

    # Initialize an array for radius values
    r_values = np.zeros_like(x_values)

    # define conical radius function (piecewise function)
    r_values[:len(x_0i)] = r0    # For x0 to xi
    r_values[len(x_0i):len(x_0i) + len(x_it)] = -np.tan(alpha)*(x_it-xt)+rt  # For xi to xt
    r_values[len(x_0i) + len(x_it):] = np.tan(alpha)*(x_te-xt)+rt  # For xt to xe

    return (x_values, r_values)

Lc = 5
alpha = 0.26

print(conical_radius_function( Lc, alpha))