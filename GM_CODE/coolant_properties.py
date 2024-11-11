import CoolProp as CP

from pyfluids import Fluid, Mixture, FluidsList, Input

ethanol = Fluid(FluidsList.Ethanol).with_state(
    Input.pressure(4e6), Input.temperature(75) #at 4 MPa and 350 K
)

rho = ethanol.density 
conductivity = ethanol.conductivity
mu = ethanol.dynamic_viscosity
C_p = ethanol.specific_heat

m_tot = 0.13 #kg/s as estimated by rpa
N= 50 # number of channels (assumed)
m_channel = m_tot/N
d = 0.002 #assumed
A = d**2
v = m_channel/(rho*A)

Re = rho*v*d/mu
print(Re)
Pr = mu*C_p/conductivity
h_coolant = 0.023*conductivity/d*Re**0.8*Pr**0.4
print(h_coolant)

