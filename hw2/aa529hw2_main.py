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
species = ["electron", "xenon"]
n0 = 1.0*10**(18) # [1/m^{3}]
Te = 10.0 # [eV]
B0 = 0.02 # [T]
E = 300.0 # [V/cm]
me = 9.11*10**(-31) # [kg]
mass_xenon = 131 * 1.67 * 10**(-27) # [kg]
e = 1.6*10**(-19) # [C]
Z = 1 # ion charge state
epsilon_0 = 8.854 * 10**(-12) # [F/m]

mass_species = [me, mass_xenon]
mu = me/mass_xenon # electron-ion mass ratio

lambda_ie = 20.0 # Coulomb logarithm for ion-electron collisions

# Part (a): Ion and Electron Larmor Radii
# Part (b): Ion and Electron Cycloctron frequency
i = 0
for mass in mass_species:
    r_larmor = 10**(2)*mass*(E/(e*B0**(2))) # [cm] -> [m]
    omega_c = e*B0/mass # [rad/s]
    print("The %s Larmor radius is %f [m]" %(species[i],r_larmor))
    print("The %s cyclotron frequency is %f [rad/s]" %(species[i],omega_c))
    i += 1

# Part (c): ExB drift velocity
v_EB = (E*B0)/B0**(2) # |ExB| = |E||B|, i.e, E and B are orthogonal
print("The ExB drift velocity is %f [m/s]" %(10**(2)*v_EB)) # E: [V/cm] -> [V/m]

# Part (d): Debye Length
lambda_De = 6.9*np.sqrt((11600.0)*Te/(10**(-6)*n0)) # Te: [eV] -> [K] and n0: [1/m^3] -> [1/cm^3]
print("The Debye Length is %f [cm]" %lambda_De)

# Part (e): Plasma Frequency
omega_p = np.sqrt((n0*e**(2))/(me*epsilon_0)) # [rad/s]
print("The plasma frequency is %f [rad/s]" %omega_p)

# Part (f): Ion-electron collision frequency
nu_ie = (10**(-6)*n0)*Z**(2)*lambda_ie*(1.6*10**(-9)*np.sqrt(mu)*Te**(-3/2)) # n0: [1/m^3] -> [1/cm^3], assume quasineutrality
print("The ion-electron collision frequency is %f [cm^3/s]" %nu_ie)

""" Problem 3 """
# Part (a): Symbolic math
# Part (b): Compute components of electron conductivity tensor
nu_e = 2.9*10**(-6)*(10**(-6))*n0*lambda_ie*Te**(-3/2) # [1/sec]
Omega_H = (e*B0/me)/nu_e # recomputing omega_ce due to lack of foresight
sigma_parallel_electrons = (n0*e**(2))/(me*nu_e) # [S/m]
sigma_perp_electrons = sigma_parallel_electrons/(1.0 + Omega_H**(2)) # [S/m]
sigma_Pederson_electrons = sigma_perp_electrons*Omega_H # [S/m]
# Part(c): Compare magnetized conductivities to conductivity of copper
print("The magnetized electron parallel conductivity is %f [Siemens/m]" %sigma_parallel_electrons)
print("The magnetized electron perpendicular conductivity is %f [Siemens/m]" %sigma_perp_electrons)
print("The magnetized electron Pederson conductivity is %f [Siemens/m]" %sigma_Pederson_electrons)

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
k_b = 1.38 * 10**(-23) # [m^2 kg s^-2 K^-1]
h = 6.62 * 10**(-34) # planck's constant, [m^2 kg s^-1]
c = 3.0 * 10**(8) # speed of light, [m s^-1]
Joules_to_eV = 1.6 * 10**(19) # conversion factor: 1 [J] = 1.6*10**(19) [eV]

# BASIC INFORMATION
# Assume quasineutrality and only atoms, i.e, no molecules form (technically untrue but complicated)
alpha_1 = 1.0
alpha_1p = 1.0
alpha_e = 1.0
alpha_p = alpha_1 + alpha_1p + alpha_e
Isp_5 = 1000.0 # [sec]

# Atomic masses, i.e 1/N0
He_AtomicMass = 4.00 * 1.67 * 10**(-27) # amu -> kg, He-4
Li_AtomicMass = 7.02 * 1.67 * 10**(-27) # amu -> kg, Li-7
C_AtomicMass = 12.00 * 1.67 * 10**(-27) # amu -> kg, Carbon-12

# Ground State ionization energy
He_GroundIonization = 24.59 # [eV]
# He_GroundIonization = 0.0 # [eV]
Li_GroundIonization = 5.39 # [eV]
C_GroundIonization = 11.26 # [eV]

# Energy levels are given in units of [1/cm] from NIST - multiply by h*c to turn into energy
# These are the \beta_{k} from Jahn Ch.6
# He_Neutral_Excitation = 159856.0 # 1s^2 -> 1s2s and J = 0 -> 1
# Li_Neutral_Excitation = 14903.6 # 2s -> 2p and J = 1/2 -> 1/2
# C_Neutral_Excitation = 33735.2 # 2s^2 2p^2 -> 2s2p^3 and J = 0 -> 2
He_Neutral_Excitation = 0.0 # 1s^2 -> 1s2s and J = 0 -> 1
Li_Neutral_Excitation = 0.0 # 2s -> 2p and J = 1/2 -> 1/2
C_Neutral_Excitation = 0.0 # 2s^2 2p^2 -> 2s2p^3 and J = 0 -> 2

# There are the \beta_{m} from Jahn Ch. 6
# He_SingIon_Excitation = 329179.3 # 1s -> 2p and J = 1/2 -> 1/2
# Li_SingIon_Excitation = 476034.2 # 1s^2 -> 1s2s and J = 0 -> 0
# C_SingIon_Excitation =  43003.3 # 2s^2 2p -> 2s 2p^2 and J = 1/2 -> 1/2
He_SingIon_Excitation = 0.0 # 1s -> 2p and J = 1/2 -> 1/2
Li_SingIon_Excitation = 0.0 # 1s^2 -> 1s2s and J = 0 -> 0
C_SingIon_Excitation =  0.0 # 2s^2 2p -> 2s 2p^2 and J = 1/2 -> 1/2

# Collect numbers together to do computation in a for loop
Elements = ["Helium", "Lithium", "Carbon"]
N0 = [1.0/He_AtomicMass, 1.0/Li_AtomicMass, 1.0/C_AtomicMass] # [1/kg]
NeutralExcitation = [He_Neutral_Excitation, Li_Neutral_Excitation, C_Neutral_Excitation]
SingIonExcitation = [He_SingIon_Excitation, Li_SingIon_Excitation, C_SingIon_Excitation]
GroundIonization = [He_GroundIonization, Li_GroundIonization, C_GroundIonization]

xi = np.linspace(0.0,1.0,25) # What proportion of the flow's internal energy is lost
Tc = np.empty((xi.shape[0],len(N0))) # Chamber temperature for a given xi and particular atom

TcProblem5Fig = plt.figure()

for eidx in np.arange(len(Elements)):
    NeutralExcitation[eidx] = 10**(2)*h*c*NeutralExcitation[eidx] # [1/cm]-> [J]
    SingIonExcitation[eidx] = 10**(2)*h*c*SingIonExcitation[eidx] # [1/cm] -> [J]
    GroundIonization[eidx] = (1.0/Joules_to_eV)*GroundIonization[eidx] # [eV] -> [J]
    for xidx in np.arange(xi.shape[0]):
        e_Frozen = alpha_1*NeutralExcitation[eidx] + alpha_1p*SingIonExcitation[eidx] + alpha_1p*GroundIonization[eidx]
        Tc[xidx,eidx] = (2.0/(3.0*k_b*alpha_p))*((g0*Isp_5)**(2)/(2.0*N0[eidx]) - (1.0 - xi[xidx])*e_Frozen)
    plt.plot(xi,Tc[:,eidx],label=Elements[eidx])
    print("For %s at xi = %f, Tc = %f [K]" %(Elements[eidx],xi[0],Tc[0,eidx]))
    print("For %s at xi = %f, Tc = %f [K]" %(Elements[eidx],xi[xi.shape[0]-1],Tc[xi.shape[0]-1,eidx]))

plt.xlabel('Frozen Flow loss parameter $\\xi$')
plt.ylabel('Chamber Temperature, $T_{c}$ [K]')
plt.title('$T_{c}$ vs. $\\xi$')
plt.legend()

plt.show()
# Isp_5 = 1000.0 # [sec]
# P_5 = 1.0 # [atm]
#
# print("The exhaust velocity required for the given Isp in problem 5 is %f [m/s] " %(g0*Isp_5))
#
# IsobaricSpecificHeat_Carbon = 0.71 # [J/g K], room temperature - couldn't find anything else
# IsobaricSpecificHeat_Lithium = 3.60 # [J/g K], " " " " " "
# IsobaricSpecificHeat_Helium = 5.19 # [J/g K], " " " " " "

# IsobaricSpecificHeat_List = [IsobaricSpecificHeat_Carbon, IsobaricSpecificHeat_Lithium, IsobaricSpecificHeat_Helium]
#
# for SpecificHeat in IsobaricSpecificHeat_List:
#     Tc_idealideal = 10**(-3)*((g0*Isp_5)**(2))/(2.0*SpecificHeat) # mass conversion factor
#     print("The chamber temperature in the ideal-ideal case is %f [K]" %Tc_idealideal)
