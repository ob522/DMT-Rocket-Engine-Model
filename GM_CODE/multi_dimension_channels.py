import numpy as np
import math


from coolant_properties import rho,conductivity,mu,C_p, mu_w

m_tot = 0.158 # from rpa
d = 17.33 #[mm] diameter of throat
circ = np.pi*d
min_rib = 1 #[mm] minimum width of ribs as should be calculated by stress analysis

h = np.arange(0.5,5.5,0.5) #height of channel at throat
w =np.arange(0.5,5.5,0.5) #width of channel at throat

N = np.zeros(len(w))
m_channel = np.zeros(len(w))

dc = 69.26 #[mm] diameter of chamber
circ_c = np.pi*dc

de = 36.08 #[mm] diameter of exit
circ_e = np.pi*de

w_c = np.arange(len(w))
w_e = np.arange(len(w))




for i in range(0,len(w)-1):
    N[i] = math.floor(circ/(w[i]+min_rib)) #rounds down to the nearest number of channels
    m_channel[i]= m_tot/N[i]
    w_c[i] = circ_c/N[i]-min_rib
    w_e[i] = circ_e/N[i]-min_rib
    

Area= h[:, np.newaxis] * w #matrix of areas of channels (in mm^2 for now)
v = np.zeros(np.shape(Area))

for i in range(0,len(h)):
    for j in range(0,len(w)):
        v[i,j]=m_channel[i]/(rho*Area[i,j]*10**-6) #output of v is in m/s

Pr = mu*C_p/conductivity #pradtl number

d = np.zeros(np.shape(Area)) #effective diameter of channels 
Re = np.zeros(np.shape(Area)) #Reynolds number effecitvely
h_coolant = np.zeros(np.shape(Area)) #convective heat transfer coefficient
h_sieder = np.zeros(np.shape(Area)) #convective heat transfer coefficient


for i in range(0,len(h)):
    for j in range(0,len(w)):
        d[i,j] = 1.30*(h[i]*w[j])**0.625/(h[i]+w[j])**0.25
        Re[i,j] = rho*v[i,j]*d[i,j]/mu*10**-3 #make d into m rather than mm
        h_coolant[i,j] = 0.023*conductivity/d[i,j]*Re[i,j]**0.8*Pr**0.4*10**3 #dittus boelter
        h_sieder[i,j] = 0.027*conductivity/d[i,j]*Re[i,j]**0.8*Pr**0.33*10**3*(mu/mu_w)**0.14 #sider tate


#print(Re) # a lot of the values are under 4000 which is weird
print(h_coolant)
print(N,w_c,w_e)

h_coolant_value_0_5 = h_coolant[1, 3]  # h=1 mm and w=0.6 mm





