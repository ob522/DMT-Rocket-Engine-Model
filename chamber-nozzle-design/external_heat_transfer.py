import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve

def get_channel_CHT_coef(q_dot_values, Twl_values, Tl_values):
    # ensure we have NumPy arrays
    q_dot = np.asarray(q_dot_values, dtype=float)
    Twl   = np.asarray(Twl_values, dtype=float)
    Tl    = np.asarray(Tl_values, dtype=float)

    # now subtraction and division work element‐wise
    hl_values = q_dot / (Twl - Tl)
    return hl_values

def calc_free_CHT_coef(Ta, Ts, D):
    g, Pa = 9.81, 101325
    fluid = 'Air'

    rho = CP.PropsSI('Dmass', 'T', Ta, 'P', Pa, fluid)
    k   = CP.PropsSI('conductivity', 'T', Ta, 'P', Pa, fluid)
    cp  = CP.PropsSI('Cpmass', 'T', Ta, 'P', Pa, fluid)
    mu  = CP.PropsSI('viscosity', 'T', Ta, 'P', Pa, fluid)

    alpha = k / (rho * cp)
    nu    = mu / rho
    Pr    = nu/alpha

    Tf   = (Ts + Ta)/2
    beta = 1/Tf
    Ra   = g*beta*(Ts - Ta)*D**3/(nu*alpha)
    Nu   = 0.6 + (0.387*Ra**(1/6))/((1 + (0.559/Pr)**(9/16))**(8/27))
    h = k*Nu/D
    print('free CHTC:',h)
    return h

def calc_external_wall_temp(q_dot_values,
                            Twl_values,
                            Tl_values,
                            r_values,
                            Ta, Pa,
                            k_w,    # wall thermal conductivity, W/m·K
                            t_w):   # wall thickness, m

    # cast inputs
    q_dot = np.asarray(q_dot_values, dtype=float)
    Twl   = np.asarray(Twl_values, dtype=float)
    Tl    = np.asarray(Tl_values, dtype=float)
    r     = np.asarray(r_values, dtype=float)

    # 1) channel‐side HTC
    h_l = get_channel_CHT_coef(q_dot, Twl, Tl)

    # 2) hydraulic diameter for each element
    D = (r + 0.0003 + 0.002 + 0.003) * 2

    Ts_values = np.zeros_like(q_dot)

    for i in range(len(q_dot)):
        hl = h_l[i]
        Di = D[i]
        qd = q_dot[i]
        Tli = Tl[i]

        # unknowns: Ts, Ti, q
        def residuals(vars):
            Ts, Ti, q = vars
            he = calc_free_CHT_coef(Ta, Ts, Di)
            # 1) conduction through wall
            eq1 = q - (k_w/t_w)*(Ti - Ts)
            # 2) free convection outside
            eq2 = q - he*(Ts - Ta)
            # 3) channel‐side convection
            eq3 = q - hl*(Tli - Ti)
            return [eq1, eq2, eq3]

        # initial guesses
        Ts0 = (Twl[i] + Ta)/2
        Ti0 = (Twl[i] + Tli)/2
        q0  = qd

        sol = fsolve(residuals, x0=[Ts0, Ti0, q0])
        Ts_values[i] = sol[0]

    return Ts_values


######################################
# Transient heat analysis
######################################

def lumped_capacitance1(h, A, rho, c_p, V, T_i, T_inf, t_step=1.0):
    """
    Calculate transient temperature using the lumped capacitance model (in Kelvin).
    
    Parameters:
    - h: convective heat transfer coefficient (W/m²·K)
    - A: surface area of the object (m²)
    - rho: density of the object (kg/m³)
    - c_p: specific heat capacity of the object (J/kg·K)
    - V: volume of the object (m³)
    - T_i: initial temperature of the object (K)
    - T_inf: ambient temperature (K)
    - t_step: time step for the simulation (s)
    
    Returns:
    - times: array of time values (s)
    - temps: array of corresponding temperature values (K)
    """
    
    tau = (rho * c_p * V) / (h * A)  # time constant (s)
    
    # approximate steady state: within 1% of ambient
    target_temp = T_inf + 0.01 * (T_i - T_inf)
    
    times = []
    temps = []
    
    # t = 0.0
    # T = T_i
    
    # while abs(T - T_inf) > abs(0.01 * (T_i - T_inf)):
    #     times.append(t)
    #     temps.append(T)
    #     T = T_inf + (T_i - T_inf) * np.exp(-t / tau)
    #     t += t_step

    t = 0.0
    T = T_i

    while t <= 70.0:
    #while abs(T - T_inf) > abs(0.01 * (T_i - T_inf)):
        times.append(t)
        temps.append(T)
        t += t_step
        T = T_inf + (T_i - T_inf) * np.exp(-t / tau)
        
    # add final time step at steady state
    times.append(t)
    temps.append(T_inf)
    
    return np.array(times), np.array(temps)

def lumped_capacitance(h, A, rho, c_p, V, T_i, T_inf, t_step=1.0):
    """
    Calculate transient temperature using the lumped capacitance model (in Kelvin).
    
    Parameters:
    - h: convective heat transfer coefficient (W/m²·K)
    - A: surface area of the object (m²)
    - rho: density of the object (kg/m³)
    - c_p: specific heat capacity of the object (J/kg·K)
    - V: volume of the object (m³)
    - T_i: initial temperature of the object (K)
    - T_inf: ambient temperature (K)
    - t_step: time step for the simulation (s)
    
    Returns:
    - times: array of time values (s)
    - temps: array of corresponding temperature values (K)
    """
    
    tau = (rho * c_p * V) / (h * A)  # time constant (s)
    
    times = []
    temps = []
    
    t = 0.0
    T = T_i

    while t <= 70.0:
        times.append(t)
        temps.append(T)
        t += t_step
        T = T_inf + (T_i - T_inf) * np.exp(-t / tau)
        
    times.append(t)
    temps.append(T_inf)
    
    return np.array(times), np.array(temps)

# def transient_conduction_1D(k, rho, c, L, Ti, Tinf):
#     C1 = 1.1191
#     zeta = 0.8603
#     alpha = k/(rho*c)
#     Fo = alpha*t/(L**2)
#     theta_star  = C1*np.exp(-Fo*zeta**2)

#     T = theta_star*(Ti-Tinf) + Tinf


def transient_conduction_1D(k, rho, c, L, Ti, Tinf, t_step=1.0):
    """
    Calculate transient surface temperature using the 1D transient conduction model (in Kelvin).
    
    Parameters:
    - k: thermal conductivity (W/m·K)
    - rho: density (kg/m³)
    - c: specific heat capacity (J/kg·K)
    - L: wall thickness (m)
    - Ti: initial temperature (K)
    - Tinf: ambient temperature (K)
    - t_step: time step for the simulation (s)
    
    Returns:
    - times: array of time values (s)
    - temps: array of corresponding surface temperature values (K)
    """
    
    C1 = 1.1191
    zeta = 0.8603
    alpha = k / (rho * c)
    
    times = []
    temps = []
    
    t = 0.0
    
    while t <= 70.0:
        Fo = alpha * t / (L ** 2)
        theta_star = C1 * np.exp(-Fo * zeta ** 2)
        T = theta_star * (Ti - Tinf) + Tinf
        
        times.append(t)
        temps.append(T)
        
        t += t_step
    
    times.append(t)
    temps.append(Tinf)
    
    return np.array(times), np.array(temps)

