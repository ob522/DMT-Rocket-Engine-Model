import numpy as np
import math
import matplotlib.pyplot as plt

T = np.linspace(500,2600,100) #temperature
index = [0,37,54,100] #four indeces
temp_index = (T,index) #expecting the tuble outout to be like this. if it isn't code below may not work

#all constants below found from the following website: https://ntrs.nasa.gov/api/citations/19940013151/downloads/19940013151.pdf
con1000 = [ #constants from 1000 to 5000 K
    [0.65060585, 28.517449, -16690.236, 1.5223271], #N2 constants [A,B,C,D]
    [0.50714993, -689.66913, 87454.75, 3.0285155], #H2O constants
    [0.65060585, 28.517449, -16690.236, 1.5223271], #CO constants
    [0.65318879, 51.738759, -62834.882, 1.5227045] #CO2
]
con300 = [ #constants from 300 to 1000K
    [0.60443938, -43.632704, -884.41949, 1.897215], #N2 [A,B,C,D]
    [0.7838778, -382.60408, 49040.158, 0.85222785], #H2O
    [0.60443938, -43.632704, -884.41949, 1.897215], #CO
    [-0.54330318, -188.23898, 8872.6567, 2.4499362] #CO2
]

mass_fraction = [ #mass fraction of each gas found from RPA
    [0.4628923, 0.4628923, 0.4629015, 0.4629033], #N2
    [0.1590833, 0.1590833, 0.1557418, 0.1231363], #H2O
    [0.2821263, 0.2821263, 0.2767468, 0.2259908], #CO
    [0.0777402, 0.0777402, 0.0861982, 0.1659354] #CO2
]

#mass fraction linear interpolation done below (dw about why it works, because it does)
mass_frac_interp = np.zeros((4,temp_index[1][-1]))
print(mass_frac_interp)
for i in range(0,np.size(temp_index[1])):
    mass_frac_interp[i] = np.interp(np.arange(0,temp_index[1][-1],1),temp_index[1],mass_fraction[i])
    plt.plot(np.arange(0,temp_index[1][-1],1),mass_frac_interp[i])
plt.show()

viscocity = np.zeros(temp_index[1][-1])


for i in range(0,np.size(temp_index[0])):
    visc = 0
    mass_sum = 0
    for j in range(0,np.size(con300[1])):
        if T[i]<1000:
            visc +=mass_frac_interp[j][i]*math.exp(con300[j][0]*math.log(temp_index[0][i])+con300[j][1]/temp_index[0][i]+con300[j][2]/temp_index[0][i]**2+con300[j][3])
        elif T[i]>=1000:
            visc +=mass_frac_interp[j][i]*math.exp(con1000[j][0]*math.log(temp_index[0][i])+con1000[j][1]/temp_index[0][i]+con1000[j][2]/temp_index[0][i]**2+con1000[j][3])
        mass_sum += mass_frac_interp[j][i]
    viscocity[i]= visc/mass_sum

plt.plot(np.arange(0,temp_index[1][-1],1),viscocity)
plt.show()


print(viscocity)
#print(temp_index)
