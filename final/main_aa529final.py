"""
Matt Russell
University of Washington
Department of Aeronautics & Astronautics
6/9/2021
Main program file for final
"""
import numpy as np
import matplotlib.pyplot as plt
""" Problem 2 """
phiS_minus_phiP = np.asarray([-50.0,-100.0,-150.0,-200.0]) # [V] - electrostatic potential difference between sheath (S) and plasma (P)
Te = np.asarray([12.5, 16.5, 12.5, 10.0]) # [eV] - determined empirically

k_B = 1.38e-23 # [J/degK]
e = 1.6e-19

eV_to_degK = 11600.0 # 1 eV = this many degK
Joules_to_eV = (1.0/1.6)*(1.0e19)
amu_to_kg = 1.67e-27
cm_to_m = 1.0e-2

m_e = 9.11e-31 # [kg]
m_i = 131.0*amu_to_kg # [kg]<-[amu] Xenon atomic mass
n_s = 0.25e12/(cm_to_m**3) # [m^-3]<-[cm^-3] Quasineutral sheath density

r_out = 7.0*cm_to_m # [m]<-[cm]
r_in = 3.0*cm_to_m # [m]<-[cm]
L = 4.0*cm_to_m # [m]<-[cm] Length of the plasma channel
A_w = 2.0*np.pi*L*(r_in + r_out) # [m] Total area of the wall

eta_current = 0.75
eta_voltage = 0.90
I_d = 10.0 # [A]
V_d = 300.0 # [V]
I_b = eta_current*I_d
V_b = eta_voltage*V_d
P_d = I_d * V_b

P_we = np.zeros(phiS_minus_phiP.size)
P_wi = np.zeros(phiS_minus_phiP.size)
P_w = np.zeros(phiS_minus_phiP.size)
E_ions = np.zeros(phiS_minus_phiP.size)

idx = 0
Te_vs_Phi = plt.figure()
Pwe_vs_Phi = plt.figure()
Pwi_vs_Phi = plt.figure()
Pw_vs_Phi = plt.figure()
Eion_vs_Phi = plt.figure()
for phi in phiS_minus_phiP:
    v_is = np.sqrt((k_B*eV_to_degK*Te[idx])/m_i)
    P_we[idx] = 0.25*np.sqrt((8.0*k_B*eV_to_degK*Te[idx])/(np.pi*m_e))*n_s*A_w \
    *np.exp((e*phi)/(k_B*eV_to_degK*Te[idx]))*2.0*k_B*eV_to_degK*Te[idx]
    P_wi[idx] = n_s*v_is*A_w*(0.5*m_i*v_is**2 - e*phi)
    P_we[idx] = P_we[idx]/P_d
    P_wi[idx] = P_wi[idx]/P_d
    P_w[idx] = P_we[idx] + P_wi[idx]
    E_ions[idx] = Joules_to_eV*k_B*eV_to_degK*Te[idx]*(0.5 + np.log(np.sqrt((2.0*m_i)/(np.pi*m_e))))

plt.figure(Te_vs_Phi.number)
plt.scatter(phiS_minus_phiP,Te)
plt.xlabel('$\phi_{s} - \phi_{p}$ [V]')
plt.ylabel('Electron Temperature [eV]')
plt.title('Electron Temperature')

plt.figure(Pwe_vs_Phi.number)
plt.scatter(phiS_minus_phiP,P_we)
plt.xlabel('$\phi_{s} - \phi_{p}$ [V]')
plt.ylabel('$\\frac{P_{w,e}}{P_{d}}$')
plt.title('Wall loss - electrons')

plt.figure(Pwi_vs_Phi.number)
plt.scatter(phiS_minus_phiP,P_wi)
plt.xlabel('$\phi_{s} - \phi_{p}$ [V]')
plt.ylabel('$\\frac{P_{w,i}}{P_{d}}$')
plt.title('Wall loss - ions')

plt.figure(Pw_vs_Phi.number)
plt.scatter(phiS_minus_phiP,P_w)
plt.xlabel('$\phi_{s} - \phi_{p}$ [V]')
plt.ylabel('$\\frac{P_{w}}{P_{d}}$')
plt.title('Total Wall losses')

print(E_ions)
plt.figure(Eion_vs_Phi.number)
plt.scatter(phiS_minus_phiP,E_ions)
plt.xlabel('$\phi_{s} - \phi_{p}$ [V]')
plt.ylabel('Ion energy [eV]')
plt.title('Energy of wall bombardment')

""" Problem 3 """

plt.show()
