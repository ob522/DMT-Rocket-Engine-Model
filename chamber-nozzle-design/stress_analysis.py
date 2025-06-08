import numpy as np

def calc_hoop_stress(q_dot, Pl, Pg, w, tw, E,v, alpha,k):
    sigma_hoop = ((Pl-Pg)*(w/tw)**2)/2 + E*alpha*q_dot*tw/(2*(1-v)*k)
    return sigma_hoop

def calc_radial_stress(Twg, Twc, E, alpha):
    sigma_r = E*alpha*(Twg-Twc)
    return sigma_r

def calc_max_axial_stress(Pl_max, Al_max, Pg_max, Ac, At,rt,tw,hc, t_outer ):
    F_max = Pl_max*Al_max + Pg_max*(Ac-At)
    A_min = tw*2*np.pi*rt + t_outer*2*np.pi*(rt+tw+hc)
    sigma_axial_max = F_max/A_min
    return sigma_axial_max

def calc_Von_Mises_equivalent_stress(sigma_hoop, sigma_r, sigma_axial):
    sigma_VM = np.sqrt((sigma_hoop-sigma_r)**2+(sigma_r-sigma_axial)**2+(sigma_axial-sigma_hoop)**2)
    return sigma_VM

def stress_analysis_calculations(Twg_values, Twc_values, q_dot_values, Pl_values, Pg_values, w_values, tw, E,v, alpha,k, Pl_max, Al_max, Pg_max, Ac, At,rt,hc, t_outer):
    sigma_hoop_values = np.zeros_like(Twg_values)
    sigma_r_values =np.zeros_like(Twg_values)
    sigma_axial_max_values = np.zeros_like(Twg_values)
    sigma_VM_values = np.zeros_like(Twg_values)

    for i in range(len(sigma_hoop_values)):
        sigma_hoop_values[i] = calc_hoop_stress(q_dot_values[i], Pl_values[i], Pg_values[i], w_values[i], tw, E,v, alpha,k)
    
    for i in range(len(sigma_r_values)):
        sigma_r_values[i] = calc_radial_stress(Twg_values[i], Twc_values[i], E, alpha)

    for i in range(len(sigma_axial_max_values)):
        sigma_axial_max_values[i] = calc_max_axial_stress(Pl_max, Al_max, Pg_max, Ac, At,rt,tw,hc, t_outer)

    for i in range(len(sigma_VM_values)):
        sigma_VM_values[i] = calc_Von_Mises_equivalent_stress(sigma_hoop_values[i], sigma_r_values[i], sigma_axial_max_values[i])

    return sigma_hoop_values, sigma_r_values, sigma_axial_max_values, sigma_VM_values

