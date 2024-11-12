import numpy as np
import math


from coolant_properties import rho,conductivity,mu,C_p

m_tot = 0.135
d = 17.33 #[mm] diameter of throat
circ = np.pi*d
min_rib = 0.4 #[mm] minimum width of ribs as should be calculated by stress analysis

h = np.arange(0.5,5.5,0.5) #height of channel at throat
w =np.arange(0.2,5.2,0.2) #width of channel at throat

N = np.zeros(len(w))
m_channel = np.zeros(len(w))


for i in range(0,len(w)-1):
    N[i] = math.floor(circ/(w[i]+0.4)) #rounds down to the nearest number of channels
    m_channel[i]= m_tot/N[i]
    

Area= h[:, np.newaxis] * w #matrix of areas of channels (in mm^2 for now)
v = np.zeros(np.shape(Area))

for i in range(0,len(h)):
    for j in range(0,len(w)):
        v[i,j]=m_channel[i]/(rho*Area[i,j]*10**-6) #output of v is in m/s

d = np.zeros(np.shape(Area))

for i in range(0,len(h)):
    for j in range(0,len(w)):
        d = 1.30*(h[i]*w[j])**0.625/(h[i]+w[j])**0.25

print(d)


