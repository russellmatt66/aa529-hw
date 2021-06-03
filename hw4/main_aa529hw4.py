"""
Matt Russell
AA529: Space Propulsion, HW4
University of Washington Department of Aeronautics & Astronautics
6/4/21
"""
import numpy as np
import matplotlib.pyplot as plt

""" Constants """
amu_to_kg = 1.67e-27 # 1 [amu] = this many [kg]
Joules_to_eV = (1.0/1.6)*10.0**19 # 1 [J] = this man [eV]
m_Xe = 131.0*amu_to_kg # [kg]<-[amu]
g0 = 9.82 # [m/s^2], gravitational acceleration at sea level

""" Problem 1 """
eta_T1 = 0.6 # [dimless], thrust efficiency
Isp_1 = 1200.0 # [s]
KE_ion = (m_Xe/2.0)*(g0*Isp_1)**2 # [J]
E_cost = Joules_to_eV*KE_ion*(1.0 - eta_T1)/eta_T1 # [eV/ion]<-[J/ion]

print("The cost for a Hall/ion thruster to operate at this level is %f [eV] per ion" %E_cost)

""" Problem 2 """
