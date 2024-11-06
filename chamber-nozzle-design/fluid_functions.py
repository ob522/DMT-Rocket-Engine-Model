# Functions: fluid dynamics calculations

import numpy as np

#####################################################################################
# fluids equations
#####################################################################################

def exit_mach_no(gamma,P0,Pe):
    Me = np.sqrt((2/(gamma-1))*((P0/Pe)**((gamma-1)/gamma)-1))
    return Me

def exit_temp(T0, gamma, Me):
    Te = T0/(1+(gamma-1)*(Me**2)/2)
    return Te

def exit_velocity(Me,gamma, R, Te):
    Ve = Me*np.sqrt(gamma*R*Te)
    return Ve

def mass_flow_rate(F, Ve):
    m_dot = F/Ve
    return m_dot

def throat_area(m_dot, P0, T0, gamma,R ):
    At = m_dot/(P0*(np.sqrt(gamma/(R*T0)))*((2/(gamma+1))**((gamma+1)/(2*(gamma-1)))))
    return At

def exit_area(Me, At, gamma):
    Ae = At*(2*(1+((gamma-1)/2)*Me**2)/(gamma+1))**((gamma+1)/(2*(gamma-1)))/Me
    return Ae

def chamber_sizing(char_L, At, CR):
    Vc = At*char_L
    A0 = CR*At
    Lc = Vc/A0
    return (A0, Lc)

def area_to_radius(A):
    r = np.sqrt(A/np.pi)
    return r

#####################################################################################
# Nozzle geometry
#####################################################################################

def chamber_nozzle_sizing(Pe,P0,F,CR,char_L, T0,gamma,R):
    Me = exit_mach_no(gamma,P0,Pe)
    Te = exit_temp(T0, gamma, Me)
    Ve = exit_velocity(Me,gamma, R, Te)
    m_dot = mass_flow_rate(F, Ve)
    At = throat_area(m_dot, P0, T0, gamma,R )
    Ae = exit_area(Me, At, gamma)
    (A0, Lc) = chamber_sizing(char_L, At, CR)
    print(m_dot)  ############################### Remove
    print("looooooooooooooooooooooook")########## Remove
    return (A0,At,Ae,Lc)

def conical_radius_function(A0, At, Ae, Lc, alpha):
    # define key radii along nozzle length
    r0 = area_to_radius(A0)
    rt = area_to_radius(At)
    re = area_to_radius(Ae)
    
    # define key points along nozzle length
    x0 = 0
    xi = Lc
    xt = Lc + (r0-rt)/(np.tan(alpha))
    xe = xt + (re-rt)/(np.tan(alpha))

    # set x steps
    dx_0i = 0.0002   # Step between x0 and xi
    dx_it = 0.0001   # Step between xi and xt
    dx_te = 0.0001   # Step between xt and xe

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

    throat_index = len(x_0i) + len(x_it)
    return (x_values, r_values, throat_index)

def conical_area_function(A0, At, Ae, Lc, alpha):
    (x_values, r_values, throat_index) = conical_radius_function(A0, At, Ae, Lc, alpha)
    Area_values = (np.pi)*r_values**2
    return (x_values, Area_values, throat_index)

#####################################################################################
# fluids conditions 
#####################################################################################

def mach_no_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma):
    (x_values, Area_values, throat_index) = conical_area_function(A0, At, Ae, Lc, alpha)

    ### Newton-Raphson method
    def func(M,A,At,gamma):
        f = (1 / M) * ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M**2))**((gamma + 1) / (2 * (gamma - 1))) - A/At
        return f

    def myNewton(x0,eps, A,At,gamma):
        # set a dx step
        dx = 0.1
        # set initial guess as solution
        xn = x0
        # set the error large enough, to enter the loop once
        err = 10*eps
        # repeat while the error is too large
        while err>eps:
            # set the current solution as old solution
            xp = xn
            # evaluate the function at xn and xn+dx
            fxn = func(xn, A,At,gamma)
            fxndx = func(xn+dx, A,At,gamma)
            # compute derivative
            dfxn = (fxndx - fxn)/ dx
            # apply NR method
            xn = xp - fxn/dfxn
            # assess the error
            err = abs(xn-xp)
        return xn
    
    #set max error
    eps = 0.0001

    # Subsonic region
    subsonic_values = np.zeros_like(x_values[:throat_index])
    for i in range(len(subsonic_values)):
        subsonic_values[i] = myNewton(eps, eps, Area_values[i], At, gamma)  # M0 = eps for subsonic region

    # Supersonic region
    supersonic_values = np.zeros_like(x_values[throat_index:])
    for i in range(len(supersonic_values)):
        supersonic_values[i] = myNewton(10, eps, Area_values[throat_index + i], At, gamma)  # M0 = 10 for supersonic region

    M_values = np.concatenate((subsonic_values, supersonic_values))

    return (x_values, M_values)

def pressure_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,P0):
    (x_values, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma)
    P_values = P0 * (1 + ((gamma - 1) / 2) * M_values**2) ** (-gamma / (gamma - 1))
    return (x_values,P_values)

def density_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,rho0):
    (x_values, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma)
    rho_values = rho0 * (1 + ((gamma - 1) / 2) * M_values**2) ** (-1 / (gamma - 1))
    return (x_values, rho_values)

def temperature_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,T0):
    (x_values, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma)
    T_values = T0 * (1 + ((gamma - 1) / 2) * M_values**2) ** (-1)
    return (x_values,T_values)

def velocity_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,R,T0):
    (x_values, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma)
    (T_values) = temperature_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,T0)[1]
    c_values = np.sqrt(gamma*R*T_values)
    v_values = M_values*c_values
    return (x_values,v_values)






    
