"""
Matt Russell
University of Washington
Department of Aeronautics & Astronautics
6/9/2021
Module for final
"""
import numpy as np

k_B = 1.38e-23 # [J/degK]
e = 1.6e-19 # [C]
m_e = 9.11e-31 # [kg]

U_ex1 = 11.6 # [eV] Xenon
U_I = 12.13 # [eV] Xenon

def Rex1(Te): # Pass Te as eV, [m^3 s^-1]
    result_ex1 = (1.931e-19)*np.sqrt((8.0*e)/(np.pi*m_e))*np.exp(-U_ex1/Te)
    return result_ex1

def Rion(Te): # Pass Te as eV [m^3 s^-1]
    result_ion = (1.0e-20)*((-1.031e-4)*Te**2 + 6.386*np.exp(-U_I/Te))*np.sqrt((8.0*e*Te)/(np.pi*m_e))
    return result_ion

def UIstar(Te): # Pass Te as eV
    result_star = U_I + (Rex1(Te)/Rion(Te))*U_ex1
    return result_star
