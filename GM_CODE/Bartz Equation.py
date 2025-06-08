#Bartz Equation
import math
import numpy as np

D_t = 19*10**(-3) #throat diameter found from RPA

p_c = 2*10**6 #stagnation pressure (make sure is correct)

R_U = 0.06 #wall radius of curvature at throat

A = 0.007 #cross sectional area at point x (varies throught chamber/nozzle)
A_t = np.pi*D_t**2/4 #throat cross sectional area
T_s = 300 #Stagnation Temperature ??idk what is is
T_hw = 1000 #hot wall temperature ?? must make assumption

gamma = 1.2 #combustion gas specific heat  ??dk what to put
M = 2 #Mach Number
#convective heat transfer coefficient


T_free = #free stream temperature
T = T_free*(1+0.032*M**2+0.58*(T_hw/T_free-1))

a = [0,1,2,3,4] #constants for solving C_p: https://www.cryo-rocket.com/equation-of-state/2.4-combustion-products/#2.4.1

C_p = R_U*(a[0]+a[1]*T+a[2]*T**2+a[3]*T**3+a[4]*T**4)#combustion gas specific heat

eta = #combustion gas viscosity: https://www.cryo-rocket.com/equation-of-state/2.6-transport-properties/#2.6.12
lamda =  #combustion gas thermal conductivity: https://www.cryo-rocket.com/equation-of-state/2.6-transport-properties/#2.6.12
Pr = C_p*eta/lamda #combustion gas Prandtl number


R_s = #specific gas constant
c = (np.sqrt(gamma*R_s*T))/(gamma*np.sqrt((2/(gamma+1))**((gamma+1)/(gamma-1)))) #characteristic velocity: https://www.cryo-rocket.com/combustion-model/3.4-mixture-ratio-optimization/#3.4.1

sigma = 1/(0.5*T_hw/T_s*(1+(gamma-1)/2*M**2+0.5)**0.68*(1+(gamma-1)/2*M**2))**0.12
h_gas = 0.026/(D_t)**2*(eta**0.2*C_p/(Pr**0.6))*(p_c/c)**0.8*(D_t/R_U)**0.1*(A_t/A)**0.9*sigma
