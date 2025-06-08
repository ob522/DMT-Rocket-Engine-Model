# Functions: nozzle fluid dynamics calculations

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
    return (A0,At,Ae,Lc,m_dot)

def conical_radius_function(A0, At, Ae, Lc, alpha):
    # define key radii along nozzle length
    r0 = area_to_radius(A0)
    rt = area_to_radius(At)
    re = area_to_radius(Ae)
    
    # define key points along nozzle length
    #x0 = 0
    #xi = Lc - (r0-rt)/(np.tan(30*np.pi/180)) #pre nozzle angle of 30 degrees  ## also replaced below
    #xt = Lc 
    #xe = xt + (re-rt)/(np.tan(alpha))

    x0 = 0
    xi = Lc  #pre nozzle angle of 30 degrees  ## also replaced below
    xt = Lc + (r0-rt)/(np.tan(30*np.pi/180))
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
    r_values[len(x_0i):len(x_0i) + len(x_it)] = -np.tan(30*np.pi/180)*(x_it-xt)+rt  # For xi to xt
    r_values[len(x_0i) + len(x_it):] = np.tan(alpha)*(x_te-xt)+rt  # For xt to xe
    throat_index = len(x_0i) + len(x_it)
    chamber_end_index = len(x_0i)
    return (x_values, r_values, throat_index, chamber_end_index)

def parabolic_radius_function(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te):
    # define key radii along nozzle length
    r0 = area_to_radius(A0)
    r1 = r0
    rt = area_to_radius(At)
    re = area_to_radius(Ae)
    
    x0 = 0
    x1 = Lc
    xt = Lt
    xe = L
    
    x2 = x1 + R2*np.sin(b)
    x3 = xt - R1*np.sin(b)
    x4 = xt + Rn*np.sin(Tn)

    # set x steps
    #dx_01 = 0.0001   # Step between x0 and x1
    #dx_12 = 0.0001  
    #dx_23 = 0.0001  
    #dx_3t = 0.0001  
    #dx_t4 = 0.0001 
    #dx_4e = 0.0001 

    dx_01 = 0.001   # Step between x0 and x1
    dx_12 = 0.001  
    dx_23 = 0.001  
    dx_3t = 0.001  
    dx_t4 = 0.001 
    dx_4e = 0.001   

    # define x array (with given steps)

    x_01 = np.arange(x0, x1, dx_01)  # From x0 to (x1 - dx_01)
    x_12 = np.arange(x1, x2, dx_12)
    x_23 = np.arange(x2, x3, dx_23)
    x_3t = np.arange(x3, xt, dx_3t)
    x_t4 = np.arange(xt, x4, dx_t4)
    x_4e = np.arange(x4, xe + dx_4e, dx_4e) # From x4 to xe

    # Add the final value xe (since arange may not include it)
    if x_4e[-1] < xe:
        x_4e = np.append(x_4e, xe)

    x_values = np.concatenate((x_01, x_12, x_23, x_3t, x_t4, x_4e))  # Concatenate all the segments into a single array

    # Initialize an array for radius values
    r_values = np.zeros_like(x_values)

    # define conical radius function (piecewise function)
    r_values[:len(x_01)] = r0    # For x0 to x1
    r_values[len(x_01):len(x_01) + len(x_12)] = r1 - R2 + np.sqrt((R2)**2 - (x_values[len(x_01):len(x_01) + len(x_12)]-x1)**2)  # For x1 to x2
    #r_values[len(x_01) + len(x_12):len(x_01) + len(x_12) + len(x_23)] =  r_values[len(x_01) + len(x_12)-1] - (x_values[len(x_01) + len(x_12):len(x_01) + len(x_12) + len(x_23)] - x2)*np.tan(b)  # For x2 to x3
    # x2 to x3 is done last (prevent discontinuity)
    r_values[len(x_01) + len(x_12) + len(x_23):len(x_01) + len(x_12) + len(x_23) + len(x_3t)] = rt + R1 - np.sqrt((R1)**2 - (xt - x_values[len(x_01) + len(x_12) + len(x_23):len(x_01) + len(x_12) + len(x_23) + len(x_3t)])**2)  # For x3 to xt
    r_values[len(x_01) + len(x_12) + len(x_23) + len(x_3t):len(x_01) + len(x_12) + len(x_23) + len(x_3t) + len(x_t4)] = rt + Rn - np.sqrt(Rn**2 - (x_values[len(x_01) + len(x_12) + len(x_23) + len(x_3t):len(x_01) + len(x_12) + len(x_23) + len(x_3t) + len(x_t4)] - xt)**2) # For xt to x4
    a = (np.tan(Tn)-np.tan(Te))/(2*(x4-xe))
    b = np.tan(Te) - 2*a*xe
    c = re - a*xe**2 - b*xe
    r_values[len(x_01) + len(x_12) + len(x_23) + len(x_3t) + len(x_t4):] = a*(x_values[len(x_01) + len(x_12) + len(x_23) + len(x_3t) + len(x_t4):])**2 + b*(x_values[len(x_01) + len(x_12) + len(x_23) + len(x_3t) + len(x_t4):]) + c # For x4 to xe
    
    grad = (r_values[len(x_01) + len(x_12) + len(x_23)] - r_values[len(x_01) + len(x_12)-1])/(x3-x2)
    intercept = r_values[len(x_01) + len(x_12)-1] - grad*x2
    r_values[len(x_01) + len(x_12):len(x_01) + len(x_12) + len(x_23)] = grad*x_values[len(x_01) + len(x_12):len(x_01) + len(x_12) + len(x_23)] + intercept # For x2 to x3

    throat_index = len(x_01) + len(x_12) + len(x_23) + len(x_3t)
    chamber_end_index = len(x_01)
    return (x_values, r_values, throat_index, chamber_end_index)


def conical_area_function(A0, At, Ae, Lc, alpha):
    (x_values, r_values, throat_index, chamber_end_index) = conical_radius_function(A0, At, Ae, Lc, alpha)
    Area_values = (np.pi)*r_values**2
    return (x_values, Area_values, throat_index, chamber_end_index)

def parabolic_area_function(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te):
    (x_values, r_values, throat_index, chamber_end_index) = parabolic_radius_function(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te)
    Area_values = (np.pi)*r_values**2
    return (x_values, Area_values, throat_index, chamber_end_index)

#####################################################################################
# fluids conditions 
#####################################################################################

def mach_no_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma):
    (x_values, Area_values, throat_index, chamber_end_index) = parabolic_area_function(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te)

    ### Newton-Raphson method
    def func(M,A,At,gamma):
        f = (1 / M) * ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M**2))**((gamma + 1) / (2 * (gamma - 1))) - A/At
        return f

    def myNewton(x0,eps, A,At,gamma):
        # set a dx step
        #dx = 0.1
        dx = 1e-5
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
    eps = 0.00001

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

def pressure_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,P0):
    (x_values, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)
    P_values = P0 * (1 + ((gamma - 1) / 2) * M_values**2) ** (-gamma / (gamma - 1))
    return (x_values,P_values)

def density_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,rho0):
    (x_values, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)
    rho_values = rho0 * (1 + ((gamma - 1) / 2) * M_values**2) ** (-1 / (gamma - 1))
    return (x_values, rho_values)

def temperature_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,T0):
    (x_values, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)
    T_values = T0 * (1 + ((gamma - 1) / 2) * M_values**2) ** (-1)
    return (x_values,T_values)

def velocity_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,T0):
    (x_values, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)
    (T_values) = temperature_as_function_of_distance(A0, At, Ae, Lc, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,T0)[1]
    c_values = np.sqrt(gamma*R*T_values)
    v_values = M_values*c_values
    return (x_values,v_values)






    
