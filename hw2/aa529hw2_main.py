"""
Matt Russell
University of Washington Department of Aeronautics & Astronautics
AA529: Space Propulsion
5/5/21
Code for HW2
"""
import numpy as np
import matplotlib.pyplot as plt

""" Problem 1 """
mass_xenon = 131.0*1.67*10**(-27) # [kg]
U_I = 12.13 # [eV]

v_1 = np.sqrt((2.0*U_I*1.6*10**(-19))/mass_xenon) # [m/s]

print("The velocity at which the ion kinetic energy is equal to the first ionization energy is %f" %v_1)

LiqH2O2rxn_deltaH = 460.38 * 10**(3) * (1.0/1.6) * 10**(19) # [kJ]->[J]->[eV]

print("The energy per reaction for liquid H2/O2 is %.2E [eV]" %LiqH2O2rxn_deltaH)

""" Problem 2 """
n0 = 1.0*10**(18) # [1/m^{3}]
Te = 10.0 # [eV]
B0 = 0.02 # [T]
E = 300.0 # [V/cm]
me = 9.11*10**(-31) # [kg]
e = 1.6*10**(-19) # [C]
