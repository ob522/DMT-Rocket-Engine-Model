# For tesing code snippets

import numpy as np
import matplotlib.pyplot as plt

import CoolProp.CoolProp as CP

def calc_hl_value2(m_dot, A, Dh, Cp, mu, Kl, Twl, Tl):
    Re = m_dot*Dh/(A*mu)
    Pr = mu*Cp/Kl
    hl = 0.023*(Kl/Dh)*Re**(0.8)*Pr**(0.4)*(Twl/Tl)**(-0.3)
    #hl = 0.023*(Cp*m_dot/A)*(Dh*m_dot/A*mu)**(-0.2)*(mu*Cp/Kl)**(-2/3)
    return hl

T1 = 300
T2 = 400
P = 10**5
Cp1 = CP.PropsSI('Cpmass', 'T', T1, 'P', P, 'Ethanol') #Specific Heat at Constant Pressure (Cp) in J/(kg·K)
mu1 = CP.PropsSI('VISCOSITY', 'T', T1, 'P', P, 'Ethanol') #Dynamic Viscosity (μ) in Pa·s
Kl1 = CP.PropsSI('CONDUCTIVITY', 'T', T1, 'P', P, 'Ethanol') #Thermal Conductivity (k) in W/(m·K)
#Dh = calc_hydraulic_diameter(Ac_values[i], Wc_values[i])
Cp2 = CP.PropsSI('Cpmass', 'T', T2, 'P', P, 'Ethanol') #Specific Heat at Constant Pressure (Cp) in J/(kg·K)
mu2 = CP.PropsSI('VISCOSITY', 'T', T2, 'P', P, 'Ethanol') #Dynamic Viscosity (μ) in Pa·s
Kl2 = CP.PropsSI('CONDUCTIVITY', 'T', T2, 'P', P, 'Ethanol') #Thermal Conductivity (k) in W/(m·K)

m_dot = 1
A = 1
Dh = 1
Twl =1
Tl =1 

h1 = calc_hl_value2(m_dot, A, Dh, Cp1, mu1, Kl1, Twl, Tl)
h2 = calc_hl_value2(m_dot, A, Dh, Cp1, mu1, Kl1, Twl, Tl)

print(h1,h2)