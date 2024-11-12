import CoolProp as CP

from pyfluids import Fluid, Mixture, FluidsList, Input

ethanol = Fluid(FluidsList.Ethanol).with_state(
    Input.pressure(4.5e6), Input.temperature(75) #at 4.5 MPa and 335 K (will enter nozzle at 300 K, assume there will be 35 K increase on temp up to throat)
)

rho = ethanol.density 
conductivity = ethanol.conductivity
mu = ethanol.dynamic_viscosity
C_p = ethanol.specific_heat

m_tot = 0.135 #kg/s as estimated by rpa
N= 50 # number of channels (assumed)
m_channel = m_tot/N

d_t = 17 #diameter of throat in mm
h = 0.001 #assumed
w = 0.0005 #width at throat must be very small as it has small region
A = h*w
v = m_channel/(rho*A)

d = 1.30*(h*w)**0.625/(h+w)**0.25 #https://www.engineeringtoolbox.com/equivalent-diameter-d_205.html
print(d)
Re = rho*v*d/mu #as long as Re>4000 we can use the equation below
print(Re)
Pr = mu*C_p/conductivity
h_coolant = 0.023*conductivity/d*Re**0.8*Pr**0.4
print(h_coolant)

