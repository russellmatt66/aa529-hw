"""
Matt Russell
University of Washington Department of Aeronautics & Astronautics
AA529: Space Propulsion
5/5/21
Software for the midterm
"""
import numpy as np
import sys

def findIndex(T,T_instance):
    # Binary search
    # T - array of linearly spaced and monotonically increasing values
    # T_instance - number to search array and find index for
    i = 0 # the brakes
    i_limit = int(np.sqrt(np.size(T))) # when to pump them
    jlow = 0
    jhigh = T.shape[0]-1
    jguess = int(np.floor((jlow + jhigh)/2))
    deltaT = float((T[jhigh] - T[jlow])/np.size(T))
    # eps = deltaT*10**(-1)
    eps = 10.0
    while(np.abs(T[jguess] - T_instance) > eps):
        if(T[jguess] > T_instance):
            jhigh = jguess
            jguess = int(np.floor((jlow + jhigh)/2))
        if(T[jguess] < T_instance):
            jlow = jguess
            jguess = int(np.floor((jlow + jhigh)/2))
        i += 1
        if(i > i_limit):
            print("Search halted manually, index not found, returning -1")
            return -1
    return jguess
"""
Problem 1
Combustion of Methane
"""
a = 1
b = 2
HeatForm_CH4l = -75.00 # product, [kJ/mol]
HeatForm_O2l = -9.42 # product, [kJ/mol]
HeatForm_CO2g = -393.50 # reactant, [kJ/mol]
HeatForm_H2Og = -241.93 # reactant, [kJ/mol]
Q_avail = HeatForm_CH4l + float(2*HeatForm_O2l) + float(a*HeatForm_CO2g) + float(b*HeatForm_H2Og)

print("The reaction supplies %6.3f [kJ] to the surroundings" %np.abs(Q_avail))

Tc_low = 300.0
Tc_high = 6000.0
num_T = 10000
T = np.linspace(Tc_low,Tc_high,num_T)

A = np.array([[1500.0, 1.0],[300.0, 1.0]]) # isobaric specific heats found using linear interpolation
kg_to_mol_CO2 = (1.0/44.009) * 10**(3) # 1 kg CO2 = this many mols CO2
kg_to_mol_H2O = (1.0/18.015) * 10**(3) # 1 kg H2O = this many mols H2O
CisobaricEmpirical_CO2 = (1.0/kg_to_mol_CO2)*np.array([[1.476],[0.846]]) # [kJ/(kg deg K)] - > [kJ/(mol degK)]
CisobaricEmpirical_H2O = (1.0/kg_to_mol_H2O)*np.array([[3.350],[1.864]]) # [kJ/(kg deg K)] - > [kJ/(mol degK)]
Cisobariccoeff_CO2 = np.dot(np.linalg.inv(A),CisobaricEmpirical_CO2) # m and b for linear interlpolation using empirical data for T = (300,1500) deg K
Cisobariccoeff_H2O = np.dot(np.linalg.inv(A),CisobaricEmpirical_H2O)

CisobaricLinInterp_CO2 = Cisobariccoeff_CO2[0]*T + Cisobariccoeff_CO2[1]
CisobaricLinInterp_H2O = Cisobariccoeff_H2O[0]*T + Cisobariccoeff_H2O[1]

Q_req = 0.0
Tc_prime = float((Tc_low + Tc_high)/2) # guess
eps = 1.0 # search tolerance
i = 0
i_limit = 25
while(np.abs(Q_req - np.abs(Q_avail)) > eps): # Q_req for these products is also a monotonically increasing function of temperature, quadratic
    j_prime = findIndex(T,Tc_prime) # index of guess, i.e, Tc_prime
    print(j_prime)
    print(Tc_prime)
    Q_req = np.trapz(CisobaricLinInterp_CO2[0:j_prime],T[0:j_prime]) + float(2*np.trapz(CisobaricLinInterp_H2O[0:j_prime],T[0:j_prime]))
    print(Q_req)
    if(Q_req > np.abs(Q_avail)): # Tc_prime too high
        Tc_high = Tc_prime
        Tc_prime = float((Tc_low + Tc_high)/2)
    if(Q_req < np.abs(Q_avail)):
        Tc_low = Tc_prime
        Tc_prime = float((Tc_low + Tc_high)/2)
    i += 1
    if(i > i_limit):
        sys.exit("Method of Available Heat failed, halting ...")

print("The chamber temperature is calculated to be %7.3f [deg K]" %Tc_prime)

nmols_CO2 = 1
nmols_H2O = 2
molmass_CO2 = 44.009e-3 # [kg/mol]
molmass_H2O = 18.015e-3 # [kg/mol]
mhat_prod = float((nmols_CO2*molmass_CO2 + nmols_H2O*molmass_H2O)/(nmols_CO2 + nmols_H2O))

print("The average molecular mass of the products is %f [kg]" %mhat_prod)

Isp_estimate = np.sqrt(Tc_prime/mhat_prod)

print("The specific impulse is estimated as %7.3f [s]" %Isp_estimate)
