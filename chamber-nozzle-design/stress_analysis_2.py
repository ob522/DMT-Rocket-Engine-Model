# Plz Delet all the titles for the text file before using this program
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker

#file_path = '0.3mm-film from face2.txt' # this file is wrong
file_path = 'final2.txt'

#t_w = 0.001
t_w = 0.0003
Nc_main = 27
Nc_chamber = Nc_main*2
channel_transition = 0.08 # m

max_pressure_diff = 40e5 # Pa

#with open('data-for-ICLR2.txt', 'r') as input_file:
#    lines = input_file.readlines()

with open(file_path, 'r') as input_file:
    lines = input_file.readlines()

#lines = lines[8:]
lines = lines[8:] if any('#' in line for line in lines) else lines

filtered_lines = [line for line in lines if not line.startswith('#')]

# Define a regular expression to match numbers
number_pattern = re.compile(r'[-+]?\d*\.\d+|\d+')

filtered_lines = []
for line in lines:
    # Extract numerical values from the line using the regular expression
    numbers = number_pattern.findall(line)
    # Join the numbers into a comma-separated string and add it to the filtered lines
    filtered_line = (' '.join(numbers) + '\n')
    if filtered_line.strip():
        filtered_lines.append(filtered_line)

#with open('data-for-ICLR2.txt', 'w') as output_file:
#    output_file.writelines(filtered_lines)

with open(file_path, 'w') as output_file:
    output_file.writelines(filtered_lines)

x = np.array([24, 90, 150, 200, 260, 320, 370, 430, 480, 540, 590, 650, 700, 760, 820]) #this is temperature in C
y = np.array([19.9, 19.8, 19.3, 18.9, 18.5, 18.0, 17.5, 16.7, 16.5, 16.0, 15.6, 15.1, 14.6, 14.0, 13.4]) * 9.807* 1e9  # Convert to GPa * 1e9 #this is elastic modulus (Pa) # data from: https://bssa.org.uk/bssa_articles/elevated-temperature-physical-properties-of-stainless-steels/
x2 = np.array([293.15, 473.15, 673.15, 873.15, 923.15])-273.15  #this is temperature in C
y2 = np.array([450, 380, 350, 290, 240]) * 1e6 #this is yield stress (Pa)

# Fit a quadratic polynomials (degree 2) to data
coefficients = np.polyfit(x, y, 2)
coefficients2 = np.polyfit(x2, y2, 2)

# Extract and convert quadratic coefficients
A = float(coefficients[0])
B = float(coefficients[1])
C = float(coefficients[2])
A2 = float(coefficients2[0])
B2 = float(coefficients2[1])
C2 = float(coefficients2[2])

#Thermal Expansion coefficient (Kelvin^(-1))
a = 16 * 10 ** -6
#Thermal Conductivity (W/m*k)
k = 20
#Poissons Ratio
v = 0.29
#Wall Thickness (m)
#t_w = 0.001


#file_path = 'data-for-ICLR2.txt' # name of the text file filtered
#file_path = '0.3mm.txt' # name of the text file filtered



pos = []
rad = []
twc = []
twg = []
tc = []
yieldstress = []
tempstress_t = []
tempstress_l = []
tempstress_p = []
von_mises = []
a_channel = []
safety_factor = []

qtotal = []

data_arrays = []

with open(file_path, 'r') as file:
    for line in file:

        line_array = line.split()
        line_array = [float(item) for item in line_array]
        data_arrays.append(line_array)

        # check position to determine Nc
        if ((float(line_array[0]) / 1000)<channel_transition):
            Nc = Nc_chamber
        else:
            Nc = Nc_main

        T = float(line_array[6] - 273.15)
        if 0 < T < 1600:
            E = A * (T ** 2) + B * T + C
            Ys = A2 * (T ** 2) + B2 * T + C2
        else:
            E = 0
            Ys = 0
        q = float(line_array[5]) * 1000
        stress_t = (E * a * q * t_w) / (2 * (1 - v) * k)
        stress_t2 = E * a * (float(line_array[6]) - float(line_array[8]))
        #stress_p = 35e5 * 0.5 * (0.00475 * float(line_array[1]) / 47.13 / t_w) ** 2
        #stress_p = 50e5 * 0.5 * (0.00475 * float(line_array[1]) / 47.13 / t_w) ** 2
        stress_p = max_pressure_diff * 0.5 * (((2*np.pi*(float(line_array[1]))/(Nc))-1)*10**(-3) / t_w) ** 2


        # a_channel.append(0.00475 * float(line_array[1]) / 47.13)             
       
        qtotal.append(float(line_array[5]))
        twg.append(float(line_array[6]))
        twc.append(float(line_array[8]))
        tc.append(float(line_array[9]))
        #print(float(line_array[6]) - float(line_array[8]))
        line_array.append(stress_t)
       
        pos.append(float(line_array[0]) / 1000)
        rad.append(float(line_array[1]) / 1000)
       
        yieldstress.append(Ys * 1e-6)
        tempstress_t.append(stress_t * 1e-6)
        tempstress_l.append(stress_t2 * 1e-6)
        tempstress_p.append(stress_p* 1e-6)
        s1 = stress_t + stress_p
        s2 = stress_t2
        vm = np.sqrt(0.5 * ((s1 - s2)**2 + (s2)**2 + (s1)**2)) * 1e-6
        von_mises.append(vm)
        safety_factor.append((Ys * 1e-6)/vm)

scale = 32/27
pos = np.array(pos)*scale*1000
"""

fig, ax = plt.subplots(1, 3, figsize=(15, 5))

ax2 = ax[0]
#ax2.plot(pos, yieldstress, color="tab:green"    , label="Yield stress")
ax2.plot(pos, tempstress_t, color="tab:pink"    , label="Tangential Thermal")
ax2.plot(pos, tempstress_l, color="tab:purple"  , label="Longitudinal Thermal")
ax2.plot(pos, tempstress_p, color="tab:orange"  , label="Tangential Pressure")
ax2.plot(pos, von_mises, color="tab:red"        , label="Von-Mises")
ax2.set_ylabel("Stress (MPA)")
ax2.set_xlabel("Axial Distance")
ax[0].grid()
ax2.legend()


ax3 = ax[1]
ax3.plot(pos, twg, label="twg")
ax3.plot(pos, twc, label='twc')
ax3.plot(pos, tc, label='tc')
ax3.set_ylabel("Temperature (K)")
ax3.set_xlabel("Axial Distance")
ax3.legend()
ax3.grid()


ax4 = ax[2]
ax4.plot(pos, safety_factor, label="Safety Factor")
ax4.set_ylabel("Safety Factor")
ax4.set_xlabel("Axial Distance")
ax4.set_ylim(bottom=0)
ax4.set_ylim(top=5)
ax4.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax4.legend()
ax4.grid()
plt.show()

"""
fig, ax = plt.subplots(figsize=(6, 5))

# Plotting the safety factor
ax.plot(pos, safety_factor, label="Safety Factor")
ax.set_xlabel("Axial Distance", fontsize=14)
ax.set_ylabel("Safety Factor", fontsize=14)
ax.set_ylabel("Safety Factor")
ax.set_xlabel("Axial Distance")
ax.set_ylim(bottom=0)
ax.set_ylim(top=5)
#ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.4))
#ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5)) 
#ax.legend()
ax.grid()
ax.set_title("Steady State Safety Factor vs. Axial Distance")

plt.show()


#total_power = 0
#for i, q_i in enumerate(qtotal):
#    total_power += q_i * 2 * rad[i] * np.pi * (pos[1] - pos[0]) * 1000
#print(total_power)
#print((tc[0] - tc[-1]) * 2530)

# ax4 = ax[3]
# ax4.plot(pos, stress_p, label="p stress")
# ax4.set_ylabel("p stress")
# ax4.set_xlabel("Axial Distance")
# #ax4.set_ylim(bottom=0)
# #ax4.set_ylim(top=5)
# ax4.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
# ax4.legend()
# ax4.grid()
# plt.show()

# plt.plot(stress_p)
# plt.show()

#print(tempstress_p)
