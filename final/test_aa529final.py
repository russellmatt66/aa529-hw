"""
Matt Russell
University of Washington
Department of Aeronautics & Astronautics
6/9/2021
Testbed for final
"""
import numpy as np
import matplotlib.pyplot as plt
import module_aa529final as fmod

k_B = 1.38e-23 # [J/degK]
e = 1.6e-19 # [C]

eV_to_degK = 11600.0 # 1 eV = this many degK
eV_to_Joules = 1.6e-19 # 1 eV = this many Joules
cm_to_m = 1.0e-2 # 1 cm = this many meters
amu_to_kg = 1.67e-27 # 1 amu = this many kg

# phiS_minus_phiP = np.asarray([-50.0,-100.0,-150.0,-200.0]) # [V]
phiS_minus_phiP = np.asarray([-200.0]) # <- Debugging why the total power explodes
Te = np.linspace(2.0,50.0,1000) # [eV] thermal energy of the electrons

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

P_d = V_d * I_d # [W] Total Discharge Power
P_b = V_b * I_b # [W] Beam Power
P_beam = (P_b/P_d)*np.ones(Te.size) # Vectorized

P_a = 2.0*I_d*(k_B*eV_to_degK*Te)/e # [W] = [J s^-1] = [A*J*degK^-1*degK*C^-1] Power lost to anode surface
v_is = np.sqrt(k_B*eV_to_degK*Te/m_i) # [m/s] = sqrt([J*kg^-1])
I_iw = n_s*v_is*e*A_w # [A] = [C*s^-1] = [m^-3*m*s^-1*C*m^2]
UI_Star = fmod.UIstar(Te) # [eV] Effective ionization energy
# UIStarFig = plt.figure()
# PR_plus_Pion = (1.0/e)*UI_Star*(I_b + I_iw) # Power losses due to radiation and ionization
PR_plus_Pion = np.empty(Te.size)

for Tidx in np.arange(Te.size):
    P_a[Tidx] = P_a[Tidx]/P_d
    PR_plus_Pion[Tidx] = (1.0/e)*fmod.UIstar(Te[Tidx])*(I_b + I_iw[Tidx])*eV_to_Joules
    PR_plus_Pion[Tidx] = PR_plus_Pion[Tidx]/P_d

idx = 0
ElectronTempFig = plt.figure()
for phi in phiS_minus_phiP:
    # idx = idx + 1
    # P_w = 0.25 * np.sqrt((8.0*eV_to_degK*Te)/(np.pi*m_e)) * n_s*A_w*2.0*eV_to_degK*Te \
    # * np.exp(e*phi/(eV_to_degK*Te)) + n_s*v_is*A_w*(0.5*m_i*v_is**2 - e*phi)
    # P_total = (1.0/P_d)*(P_b + P_w + P_a + PR_plus_Pion)
    P_w = np.empty(Te.size)
    P_we = np.empty(Te.size)
    P_wi = np.empty(Te.size)
    P_total = np.empty(Te.size)
    for Tidx in np.arange(Te.size):
        P_we[Tidx] = 0.25 * np.sqrt((8.0*k_B*eV_to_degK*Te[Tidx])/(np.pi*m_e)) * n_s*A_w*2.0*k_B*eV_to_degK*Te[Tidx] \
        * np.exp(e*phi/(k_B*eV_to_degK*Te[Tidx]))
        P_wi[Tidx] = n_s*v_is[Tidx]*A_w*(0.5*m_i*v_is[Tidx]**2 - e*phi)
        P_we[Tidx] = P_we[Tidx]/P_d
        P_wi[Tidx] = P_wi[Tidx]/P_d
        P_w[Tidx] = P_we[Tidx] + P_wi[Tidx]
        # P_w[Tidx] = 0.25 * np.sqrt((8.0*k_B*eV_to_degK*Te[Tidx])/(np.pi*m_e)) * n_s*A_w*2.0*k_B*eV_to_degK*Te[Tidx] \
        # * np.exp(e*phi/(k_B*eV_to_degK*Te[Tidx])) + n_s*v_is[Tidx]*A_w*(0.5*m_i*v_is[Tidx]**2 - e*phi)
        P_w[Tidx] = P_w[Tidx]
        P_total[Tidx] = P_beam[Tidx] + P_w[Tidx] + P_a[Tidx] + PR_plus_Pion[Tidx]
    plt.figure(ElectronTempFig.number)
    plt.plot(Te,P_a,label='anode')
    plt.plot(Te,P_beam,label='beam')
    # plt.plot(Te,P_w,label='wall')
    plt.plot(Te,P_we,label='wall-electrons')
    plt.plot(Te,P_wi,label='wall-ions')
    plt.plot(Te,PR_plus_Pion,label='collisional-radiative')
    plt.plot(Te,P_total,label='total power loss')
    plt.title('phi = %.1f V' %phi)
    # plt.legend()
    idx += 1

plt.figure(ElectronTempFig.number)
plt.plot(Te,np.ones(Te.size),label='Normalized Discharge Power')
plt.xlabel('Electron temperature [eV]')
plt.ylabel('$\\frac{P}{P_{d}}$')
plt.legend()

# plt.figure(UIStarFig.number)
# plt.plot(Te,UI_Star,label='effective ionization energy Xenon')
# plt.legend()

plt.show()
