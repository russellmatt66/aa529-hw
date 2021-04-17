"""
Matt Russell
AA529: Space Propulsion
HW1 Problem 2
Fuel-rich operation of a bi-propellant rocket
"""
k_b = 1.38e-23 # [J/K]
R = 8.314 # [J/deg K* mol]
Cp_H2O = 54.0 # [J/mol deg K]
Cp_H2 = 37.0 # [J/mol deg K]

def mdot_piece(gamma,mp,Tc):
    return np.sqrt((gamma*mp)/(k_b*Tc))*((gamma + 1.0)/2.0)\
    **(-(gamma + 1.0)/(2.0*(gamma - 1.0)))

def ExpRat_expression(Me,gamma):
    return ((gamma + 1.0) / 2.0)**(-(gamma + 1.0) / (2.0 * (gamma - 1.0))) \
* (1.0 / Me)*(1.0 + ((gamma - 1.0) / 2.0)*Me**2) \
**((gamma + 1.0) / (2.0 * (gamma - 1.0)))

def computeExhaustMachNumber(ExpRat,gamma):
    Me_low = 1.0 # Assuming that the number is between these Me_low and Me_high
    Me_high = 100.0
    Me_guess = 0.5*(Me_low + Me_high)
    delta = 1e-13
    ExpRat_check = ExpRat_expression(Me_guess,gamma)
    diff = ExpRat_check - float(ExpRat)
    i = 0
    limit = 1000
    while np.abs(diff) > delta:
        if i > limit:
            break
        if diff > 0.0: # Guess too large
            Me_high = Me_guess
            Me_guess = 0.5*(Me_low + Me_high)
        if diff < 0.0: # Guess too small
            Me_low = Me_guess
            Me_guess = 0.5*(Me_low + Me_high)
        i+=1
        ExpRat_check = ExpRat_expression(Me_guess,gamma)
        diff = ExpRat_check - ExpRat
    print("Computed Expansion Ratio is %5.3f" %ExpRat_check)
    print("Real Expansion Ratio is %i" %ExpRat)
    return Me_guess

""" Bulk """
import numpy as np
import matplotlib.pyplot as plt
p_c = 100.0*6894.757 # [Pa]
ExpRat = 50
d_t = 2e-3 # throat diameter [m]
A_t = np.pi*(d_t/2.0)**2
A_e = ExpRat*A_t
print(A_e)
gamma_H2O = 1.33 # (g)
gamma_H2 = 1.41 # (g)
gamma_avg = (gamma_H2O + gamma_H2)/2.0


# Going to linearly interpolate enthalpy for Q_req
T0 = 298.0 # inlet temperature, deg K
Tf = 4000.0 # deg K
H0 = 0.0 # kJ/mol
Hf_H2 = 126.220 # kJ/mol from enthalpy table given with assignment
Hf_H2O = 184.226 # kJ/mol
dHdT_H2 = (Hf_H2 - H0)/(Tf - T0)
dHdT_H2O = (Hf_H2O - H0)/(Tf - T0)

n_H2Oprod = 2.0 # mols of H2O produced

# Enthalpies of Formation
deltaH_H2O = -241.93 # (g) [kJ/mol]
deltaH_H2 = -7.03 # (l) [kJ/mol]
deltaH_O2 = -9.42 # (l) [kJ/mol]

# Atomic Masses
amutokg = 1.66e-27
m_H2 = 2.016*amutokg # atomic mass -> [kg/particle]
m_H2O = m_H2 + 15.999*amutokg # atomic mass
m_O2 = 2*(m_H2O - m_H2) # atomic mass

mp_H2O = n_H2Oprod*m_H2O # [mol]*[kg/particle]

eps = np.array([1,2,3,4]) #fuel-rich operation mode: 1,2,3,4
Tcprime = np.array([T0,T0,T0,T0])
mp_average = np.array([0.0,0.0,0.0,0.0])
Isp = np.array([0.0,0.0,0.0,0.0])

for i in np.arange(np.size(eps)):
    print("Epsilon = %i" %eps[i])
    n_H2prod = 2.0*(float(eps[i]) - 1.0) # mols of H2 produced
    mp_H2 = n_H2prod*m_H2 # changes with eps
    mp = n_H2Oprod*m_H2O + n_H2prod*m_H2 # n_H2prod changes with eps
    mp_average[i] = mp/(n_H2Oprod + n_H2prod)
    print("Average Molecular Mass %e" %mp_average[i])
    Cp_mix = (n_H2prod*Cp_H2 + n_H2Oprod*Cp_H2O)/(n_H2prod + n_H2Oprod)
    gamma_mix = Cp_mix/(Cp_mix - R)
    print("Mixture gamma = %f" %gamma_mix)
    Q_avail = np.abs(2.0*deltaH_H2O + 2.0*(float(eps[i]) - 1.0)*deltaH_H2 \
    - 2.0*float(eps[i])*deltaH_H2 - deltaH_O2)
    Tcprime[i] = T0 + Q_avail/(n_H2Oprod*dHdT_H2O + n_H2prod*dHdT_H2)
    m_dot = A_t*p_c*mdot_piece(gamma_mix,mp_average[i],Tcprime[i])
    print("Mass flow rate %f" %m_dot)
    exhMachNumber = computeExhaustMachNumber(ExpRat,gamma_mix)
    print("Exit Mach Number %f" %exhMachNumber)
    p_e = p_c*(1.0+(gamma_mix - 1.0)/(2.0)*exhMachNumber**2)**(-gamma_mix/(gamma_mix-1.0))
    print("Exhaust pressure %f" %p_e)
    c_e = np.sqrt(((k_b*Tcprime[i])/mp_average[i]) * (2.0*gamma_mix)/(gamma_mix - 1.0)\
    *(1.0 - (p_e/p_c)**((gamma_mix - 1.0)/gamma_mix)))
    print("Exhaust speed %f" %c_e)
    F = m_dot*c_e + (p_e - 0.0)*A_e # Thrust
    print("Thrust %f" %F)
    Isp[i] = F/(m_dot*9.8)

print(Tcprime)
print(mp_average)
print(Isp)

""" Plotting """
trendfig, trendax = plt.subplots(ncols=3)
trendax[0].plot(eps,Tcprime)
trendax[1].plot(eps,mp_average)
trendax[2].plot(eps,Isp)

plt.show()
