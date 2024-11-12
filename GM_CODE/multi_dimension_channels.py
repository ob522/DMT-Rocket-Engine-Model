import numpy as np
import math


from coolant_properties import rho,conductivity,mu,C_p

m_tot = 0.135
d = 17.33 #[mm] diameter of throat
circ = np.pi*d
min_rib = 0.4 #[mm] minimum width of ribs as should be calculated by stress analysis

h = np.arange(0.5,5.5,0.5) #height of channel at throat
w =np.arange(0.2,5.2,0.1) #width of channel at throat

N = np.zeros(len(w))
m_channel = np.zeros(len(w))


for i in range(0,len(w)-1):
    N[i] = math.floor(circ/(w[i]+min_rib)) #rounds down to the nearest number of channels
    m_channel[i]= m_tot/N[i]
    

Area= h[:, np.newaxis] * w #matrix of areas of channels (in mm^2 for now)
v = np.zeros(np.shape(Area))

for i in range(0,len(h)):
    for j in range(0,len(w)):
        v[i,j]=m_channel[i]/(rho*Area[i,j]*10**-6) #output of v is in m/s

Pr = mu*C_p/conductivity #pradtl number

d = np.zeros(np.shape(Area)) #effective diameter of channels 
Re = np.zeros(np.shape(Area)) #Reynolds number effecitvely
h_coolant = np.zeros(np.shape(Area)) #convective heat transfer coefficient


for i in range(0,len(h)):
    for j in range(0,len(w)):
        d[i,j] = 1.30*(h[i]*w[j])**0.625/(h[i]+w[j])**0.25
        Re[i,j] = rho*v[i,j]*d[i,j]/mu*10**-3 #make d into m rather than mm
        h_coolant[i,j] = 0.023*conductivity/d[i,j]*Re[i,j]**0.8*Pr**0.4*10**3


print(Re) # a lot of the values are under 4000 which is weird
print(h_coolant)


Re_value_0_5 = Re[1, 3]        # h=1 mm and w=0.6 mm

h_coolant_value_0_5 = h_coolant[1, 3]  # h=1 mm and w=0.6 mm

print(Re_value_0_5,h_coolant_value_0_5)





