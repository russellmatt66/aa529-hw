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

print("The velocity at which the ion kinetic energy is equal to the first ionization energy is %f [m/s]" %v_1)

LiqH2O2rxn_deltaH = 460.38 * 10**(3) * (1.0/1.6) * 10**(19) # [kJ]->[J]->[eV]

print("The energy per reaction for liquid H2/O2 is %.2E [eV]" %LiqH2O2rxn_deltaH)

""" Problem 2 """
n0 = 1.0*10**(18) # [1/m^{3}]
Te = 10.0 # [eV]
B0 = 0.02 # [T]
E = 300.0 # [V/cm]
me = 9.11*10**(-31) # [kg]
e = 1.6*10**(-19) # [C]
epsilon_0 = 8.854 * 10**(-12) # [F/m]

""" Problem 4 """
eta_thrust = 0.65 # thrust efficiency
t_firing = 40000.0*60.0*60.0 # [sec]
g0 = 9.81 # [m/s^{2}]

alpha = 25.0 * 10**(-3) # [kg/W] - from midterm ...
mass_xenon = 131.0 * 1.67 * 10**(-27) # [kg]
x_a = 0.01 # [m], represents grid spacing - selected as free parameter
d_grid = 50.0e-2 # [m], diameter of grid - selected as free parameter
A_grid = np.pi*(d_grid/2.0)**(2) # [m^2] - grid area

Isp_optimal = (1.0/g0)*np.sqrt((2.0*eta_thrust*t_firing)/alpha) # calculate from dm/dIsp = 0 with m = mp + mps

print("The optimal Isp for the Problem 4 mission is %f [sec]" %Isp_optimal)

V0 = ((g0*Isp_optimal)**(2)*mass_xenon)/(2.0*e)

print("The discharge voltage that fits this number is %f [volts]" %V0)

J_ChildLangmuir = (4.0*epsilon_0/9.0)*np.sqrt((2.0*e)/mass_xenon)*(V0**(3/2)/x_a**(2))

print("The current is space-charge limited to %f [A/m^2]" %J_ChildLangmuir)

Thrust = (8.0 * epsilon_0 * A_grid)/(9.0) * (V0/x_a)**(2)

print("The thrust of the propulsion system is %f [N]" %Thrust)

P_total = 0.5 * (Thrust * g0 * Isp_optimal)/eta_thrust

print("The total power required for this is %f [W]" %P_total)

m_powersupply = alpha*P_total

print("The power supply necessary to operate this system would weigh %f [kg]" %m_powersupply)

u_exhaust = np.sqrt(2.0*e*V0/mass_xenon)

print("The exhaust velocity is %f [m/s]" %u_exhaust)

""" Problem 5 """
Isp_5 = 1000.0 # [sec]
P_5 = 1.0 # [atm]

print("The exhaust velocity required for the given Isp in problem 5 is %f [m/s]" %(g0*Isp_5))

IsobaricSpecificHeat_Carbon = 0.71 # [J/g K], room temperature - couldn't find anything else
IsobaricSpecificHeat_Lithium = 3.60 # [J/g K], " " " " " "
IsobaricSpecificHeat_Helium = 5.19 # [J/g K], " " " " " "

# IsobaricSpecificHeat_List = [IsobaricSpecificHeat_Carbon, IsobaricSpecificHeat_Lithium, IsobaricSpecificHeat_Helium]
#
# for SpecificHeat in IsobaricSpecificHeat_List:
#     Tc_idealideal = 10**(-3)*((g0*Isp_5)**(2))/(2.0*SpecificHeat) # mass conversion factor
#     print("The chamber temperature in the ideal-ideal case is %f [K]" %Tc_idealideal)
