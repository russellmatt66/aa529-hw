"""
Matt Russell
AA529: Space Propulsion
HW1 Problem 2
Fuel-rich operation of a bi-propellant rocket
"""
deltaH_H2O = -241.93 # (g) [kJ/mol]
deltaH_H2 = -7.03 # (l) [kJ/mol]
deltaH_O2 = -9.42 # (l) [kJ/mol]

""" Function Declarations """
def FlameTemp(Tc_initial,Q_avail,):
"""
Tc_initial - Assumed value for Tc = Tcprime
Q_avail = Available Heat
"""
Q_req =

""" Bulk """
import numpy as np
import matplotlib.pyplot as plt

# eps = fuel-rich operation mode: 1,2,3,4
# Going to linearly interpolate enthalpy for Q_req
T0 = 298.0 # deg K
Tf = 4000.0 # deg K
H0 = 0.0 # kJ/mol
Hf_H2 = 126.220 # kJ/mol
Hf_H2O = 184.226 # kJ/mol
Hf_O2 = 138.778 #  kJ/mol
dHdT_H2 = (Hf_H2 - H0)/(Tf - T0)
dHdT_H2O = (Hf_H2O - H0)/(Tf - T0)
dHdT_O2 = (Hf_O2 - H0)/(Tf - T0) 

Q_avail = 2.0*deltaH_H2O + 2.0*(float(eps) - 1.0)*deltaH_H2 - 2.0*float(eps)*deltaH_H2 \
- deltaH_O2
