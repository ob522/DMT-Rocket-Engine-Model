import CoolProp as CP
import numpy as np

from pyfluids import Fluid, Mixture, FluidsList, Input

ethanol = Fluid(FluidsList.Ethanol).with_state(
    Input.pressure(4.7e6), Input.temperature(25) #at 4.5 MPa and 335 K (will enter nozzle at 300 K, assume there will be 35 K increase on temp up to throat)
)

rho = ethanol.density 
conductivity = ethanol.conductivity
mu = ethanol.dynamic_viscosity
C_p = ethanol.specific_heat

m_tot = 0.158 #from RPA software [kg/s]
max_pressure = 5*10**6

v = 10 #must find velocity or can't find dimensions of channel (maybe relate velocity to pressure and temp?)

A = m_tot/(rho*v)

d = np.sqrt(4*A/np.pi) 

print('Main Channel diameter[mm]',d*10**(3))

safety_factor = 4

yield_strength = 150*10**6 #change later

t = max_pressure*d*safety_factor/(4*yield_strength) #assuming thin wall cylinder

print('Main Channel thickness[mm]',t*10**(3)) #result rn is way less than a mm so can make it a mm to for manufacturing purposes
