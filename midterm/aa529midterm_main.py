"""
Matt Russell
University of Washington Department of Aeronautics & Astronautics
AA529: Space Propulsion
5/5/21
Code for the midterm
"""
import numpy as np
import matplotlib.pyplot as plt
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

Isp_Problem1 = np.sqrt(Tc_prime/mhat_prod)

print("The specific impulse for I.(v) is estimated as %7.3f [s]" %Isp_Problem1)

"""
Problem 2
A specific impulse calculation
"""
g0 = 9.81 # [m/s^{2}]
eta_Thrust = 0.60
t_firing = 50.0*86400.0 # [days] -> [s]
alpha_ii = 25.0/(10.0**(3)) # [kg/kW] -> [kg/W]
alpha_iii = alpha_ii/10.0 # Improved specific impulse from II.(iii)

Isp_Problem2ii = (1.0/g0)*np.sqrt((2.0*eta_Thrust*t_firing)/alpha_ii)

print("The optimal specific impulse for the maneuever described in II.(ii) is %7.3f" %Isp_Problem2ii)

Isp_Problem2iii = (1.0/g0)*np.sqrt((2.0*eta_Thrust*t_firing)/alpha_iii)

print("The optimal specific impulse with the improvement to power supply specific mass described in II.(iii) is %7.3f" %Isp_Problem2iii)

"""
Problem 3
Examining gridded-ion thruster performance scaling using a basic model for the effect
of electron-impact ionization on operation
"""
""" Basic Parameters """
V0 = 1300.0 # Thruster voltage [V]
d_grid = 20.0 * 10**(-2) # grid diameter [m]
d_chamber = 20.0 * 10**(-2) # chamber diameter [m]
l_chamber = 15.0 * 10**(-2) # chamber length [m]
theta_grid = 0.80 # Grid transparency
theta_div = 10.0 # beam divergence [degrees]
xi_c = 0.15 # ratio of doubly- to singly-ionized xenon ions
eta_mass = 0.95 # proportion of the mass that is utilized for thrust
mass_xenon = 131.0 * 1.67 * 10**(-27) # mass of a xenon atom [kg]
U_I = 12.13 # first ionization energy [eV]
U_ex1 = 11.60 # energy to excite atom of xenon [eV]
Aa_over_Agrid = np.array([0.1, 1.0])

""" Computed """
A_grid = np.pi*(d_grid/2.0)**(2)
V_chamber = np.pi*(d_chamber/2.0)**(2)*l_chamber
alpha_mass = (1.0 + xi_c/2.0)/(1.0 + xi_c) # correction due to mass utilization inefficiency
alpha_thrust = (1.0 + xi_c/np.sqrt(2.0))/(1.0 + xi_c) # thrust correction to account for ionization

""" Electrons """
Te = np.linspace(1.0,20.0,1000)
me = 9.11*10**(-31) # [kg]
e = 1.6*10**(-19) # [C]

""" (i) """
# R_ion = np.empty(Te.shape[0])
# R_ex1 = np.empty(Te.shape[0])
# U_Ieff = np.empty(Te.shape[0])
#
# for tidx in np.arange(Te.shape[0]):
#     R_ion[tidx] = 10.0**(-20)*(-(1.031*10**(-4))*Te[tidx]**(2) + 6.386*np.exp(-(U_I)/Te[tidx])*np.sqrt((8.0*e*Te[tidx])/(np.pi*me)))
#     R_ex1[tidx] = 1.931*10.0**(-19)*(np.exp(-U_ex1/Te[tidx])/np.sqrt(Te[tidx]))*np.sqrt((8.0*e*Te[tidx])/(np.pi*me))
#     U_Ieff[tidx] = U_I + (R_ex1[tidx]/R_ion[tidx])*U_ex1

R_ion = 10.0**(-20)*(-(1.031*10**(-4))*Te**(2) + 6.386*np.exp(-(U_I)/Te)*np.sqrt((8.0*e*Te)/(np.pi*me)))
R_ex1 = 1.931*10.0**(-19)*(np.exp(-U_ex1/Te)/np.sqrt(Te))*np.sqrt((8.0*e*Te)/(np.pi*me))
U_Ieff = U_I + (R_ex1/R_ion)*U_ex1

U_IeffFig = plt.figure()
plt.plot(Te,U_Ieff)
plt.xlim((0.0,Te[Te.shape[0]-1]))
plt.ylim(0.0,70.0)
plt.axhline(y = U_I,linestyle = 'dashed')
plt.title('Effective ionization energy for a Xenon gridded-ion thruster')
plt.xlabel('Electron temperature [eV]')
plt.ylabel('$U_{I}^{\\ast}$ [eV]')

""" (ii) - (v) """
VlossFig = plt.figure()
eta_electricalFig = plt.figure()

for value in Aa_over_Agrid:
    VlossV0 = (1.0/(e*V0*1.6*10**(19)))*(U_Ieff + (5.0/2.0 - 2.0*theta_grid + (2.0 - theta_grid)*np.log(value*np.sqrt((2.0*mass_xenon)/(np.pi*me))))*Te)
    fig = plt.figure(VlossFig.number)
    plt.plot(Te,VlossV0,label="$\\frac{A_{a}}{A} = %3.2f$" %value)
    eta_electrical = (1.0/(1.0 + (1.0/theta_grid)*VlossV0))
    fig = plt.figure(eta_electricalFig.number)
    plt.plot(Te,eta_electrical,label="$\\frac{A_{a}}{A} = %3.2f$" %value)
    Temax_index = np.where(eta_electrical == np.max(eta_electrical))
    Temax = Te[Temax_index[0]]
    neutraldensity = (A_grid/V_chamber)*(1.0/(2.0*R_ion[Temax_index[0]]))*np.sqrt(Temax*1.6*10**(-19)/mass_xenon)
    print("The neutral density at which the efficiency is maximized for area ratio %3.2f is %.1E" %(value,neutraldensity))
    if(value == 0.1):
        ni = neutraldensity
        ionsoundspeed = np.sqrt(1.6*10**(-19)*Temax/mass_xenon) # [eV] -> [J] gives [m/s]
        I_ion = e*A_grid*0.5*ni*ionsoundspeed
        I_beam = theta_grid*I_ion
        thrust_ideal = I_beam*np.sqrt(2.0*mass_xenon*V0/e)
        mdot = (alpha_mass * I_beam * mass_xenon)/(e * eta_mass)
        Isp_v = eta_mass * alpha_thrust * np.cos(theta_div*(np.pi/180.0)) * (thrust_ideal/(mdot*g0))
        eta_thrust = eta_electrical[Temax_index[0]]*eta_mass*alpha_thrust**(2)*np.cos(theta_div*(np.pi/180.0))**(2)
        print("The specific impulse for (v) is %f [sec]" %Isp_v)
        print("The thrust efficiency for (v) is %f" %eta_thrust)

Isp_ideal = (1.0/g0)*np.sqrt(2.0*e*V0/mass_xenon)
print("The ideal specific impulse is %f [sec]" %Isp_ideal)

fig = plt.figure(VlossFig.number)
# plt.ylim((0.0,0.02))
plt.xlim((Te[0],Te[Te.shape[0]-1]))
plt.legend()
plt.title('Power loss for a Xenon gridded-ion thruster')
plt.xlabel('Electron temperature [eV]')
plt.ylabel('$\\frac{V_{loss}}{V_{0}}$')

fig = plt.figure(eta_electricalFig.number)
plt.xlim((Te[0],Te[Te.shape[0]-1]))
plt.ylim(top=1.0)
plt.legend()
plt.title('Electrical efficiency of a Xenon gridded-ion thruster')
plt.xlabel('Electron temperature [eV]')
plt.ylabel('$\\eta_{e}$')

plt.show()
