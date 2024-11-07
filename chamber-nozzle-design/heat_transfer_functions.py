# Functions: Heat Transfer Functions

import numpy as np

from fluid_functions import conical_area_function, mach_no_as_function_of_distance

#####################################################################################
# Bertz' equation (hg, convective ht tran coef as function of distance)
#####################################################################################

def Prandlt_estimate(gamma):
    return (4*gamma/(9*gamma-5))

def hg_as_function_of_distance(P0,A0, At, Ae, Lc, alpha, gamma,R,Cp,T0,Tw):
    Pr = Prandlt_estimate(gamma)
    mu = 7.02257 *10**(-5)  ## Temporary## ## need to replace with function for calculation ############
    c_star = np.sqrt(R*T0)/(np.sqrt(2*gamma/(gamma-1))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))) ## may require update ##### based on research ############
    rc = 0.003 #throat radius of curvature

    Dt = np.sqrt((4*At/np.pi))

    (x_values, Area_values, throat_index) = conical_area_function(A0, At, Ae, Lc, alpha)
    (x_values, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma)

    sigma = 1/((0.5*(Tw/T0)*(1+(gamma-1)*(M_values**2)/2)+0.5)**(0.68)*(1+(gamma-1)*(M_values**2)/2)**(0.12))  #correction factor
    hg_values = ((0.026/Dt**0.2)*((mu**0.2)*Cp/(Pr**0.6))*(P0/c_star)**(0.8)*(Dt/rc)**(0.1))*(At/Area_values)**(0.9)*sigma
    return (x_values,hg_values)

def calculate_hg(P0, At, gamma,R,Cp,T0,Tw, Area, M):
    Pr = Prandlt_estimate(gamma)
    #Pr = 0.61
    mu = 7.02257 *10**(-5)  ## Temporary## ## need to replace with function for calculation ############
    #mu = 6.5*10**(-5)
    #c_star = np.sqrt(R*T0)/(np.sqrt(2*gamma/(gamma-1))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))) ## may require update ##### based on research ############
    c_star = np.sqrt(gamma*R*T0)/(gamma*np.sqrt((2/(gamma+1))**((gamma+1)/(gamma-1))))
    rc = 0.003 #throat radius of curvature

    Dt = np.sqrt((4*At/np.pi))

    sigma = 1/((0.5*(Tw/T0)*(1+(gamma-1)*(M**2)/2)+0.5)**(0.68)*(1+(gamma-1)*(M**2)/2)**(0.12))  #correction factor
    hg = ((0.026/Dt**0.2)*((mu**0.2)*Cp/(Pr**0.6))*(P0/c_star)**(0.8)*(Dt/rc)**(0.1))*(At/Area)**(0.9)*sigma
    return (hg)

#####################################################################################
# Cooling convection (hl, convective ht tran coef calc)
#####################################################################################

def constant_hl_value(Cp, mu, k, m_dot, area, D):
    hl = 0.023*Cp*(m_dot/area)*(D*m_dot/(area*mu))**(-0.2)*(mu*Cp/k)**(-2/3)
    return hl



#####################################################################################
# Steady State equation solver
#####################################################################################

def SS_HT_analysis_equation_solver(hg,T,tw,Kw,hl):
    Const = 1 / ((1 / hg) + (tw / Kw) + (1 / hl))

    A = np.array([
        [1, Const, 0, 0],
        [1, 0, hg, 0],
        [1, 0, -Kw / tw, Kw / tw],
        [1, hl, 0, -hl]
    ])

    B = np.array([
        Const * T,
        hg * T,
        0,
        0
    ])

    q, Tl, Twg, Twl = np.linalg.solve(A, B)
    return (q, Tl, Twg, Twl)


# solve along chamber length WITHOUT ITTERATION !!! (chaotic results)
def SS_HT_analysis_along_chamber(x_values, hg_values, T_values, tw, Kw, hl): 
    # Initialize lists to collect results for each element
    q_values = []
    Tl_values = []
    Twg_values = []
    Twl_values = []
    
    # Loop over each element
    for x, hg, T in zip(x_values, hg_values, T_values):
        Const = 1 / ((1 / hg) + (tw / Kw) + (1 / hl))
        
        A = np.array([
            [1, Const, 0, 0],
            [1, 0, hg, 0],
            [1, 0, -Kw / tw, Kw / tw],
            [1, hl, 0, -hl]
        ])
        
        B = np.array([
            Const * T,
            hg * T,
            0,
            0
        ])
        
        # Solve the system of equations for this element
        q, Tl, Twg, Twl = np.linalg.solve(A, B)
        
        # Append results to lists
        q_values.append(q)
        Tl_values.append(Tl)
        Twg_values.append(Twg)
        Twl_values.append(Twl)
    
    # Convert lists to arrays for output
    return np.array(q_values), np.array(Tl_values), np.array(Twg_values), np.array(Twl_values)

#####################################################################################
# Steady State itterative heat transfer solver
#####################################################################################

def SS_HT_itterative_analysis_along_chamber( A_values,M_values,Tg_values,tw,Kw,hl,   P0,At,gamma,R,Cp,T0,Twg_start_val):
    #start with const Twg values - will change with itteration
    ##Twg_values = [Twg_start_val] * len(Tg_values)
    Twg_values = np.full(Tg_values.shape, Twg_start_val)

    # Initialize lists to collect results for each element
    q_values = []
    Tl_values = []
    Twl_values = []

    # Loop over each element
    for idx, (A, M, Tg, Twg) in enumerate(zip(A_values, M_values, Tg_values, Twg_values)):
        itt_complete = False

        while (not itt_complete):
            # Calculate hg
            hg = calculate_hg(P0, At, gamma,R,Cp,T0,Twg, A, M)

            # Solve the equations to find q, Tl, Twg_new, and Twl
            (q, Tl, Twg_new, Twl) = SS_HT_analysis_equation_solver(hg,Tg,tw,Kw,hl)
            print(hg)

            # Check for convergence of Twg
            if ((abs(Twg_new - Twg) / abs(Twg_new)) <= 0.1): #10% tolerence
                Twg_values[idx] = Twg_new
                q_values.append(q)
                Tl_values.append(Tl)
                Twl_values.append(Twl)
                itt_complete = True
            elif((abs(Twg_new - Twg) / abs(Twg_new)) >= 10): #stability check
                raise ValueError("Error: itteration is not converging")
            else:
                Twg = Twg_new # Update Twg for the next iteration

    # Convert lists to arrays for output
    return np.array(q_values), np.array(Tl_values), Twg_values, np.array(Twl_values)
