"""
Matt Russell
University of Washington Aeronautics & Astronautics
AA529: Space Propulsion HW3
5/21/21
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import special as sp

epsilon0 = 8.854e-12 # [F/m] - vacuum permittivity
k_b = 1.38e-23 # [J/deg K]
amutokg = 1.67e-27 # 1 amu = this many kg

""" Problem 1 """
V0 = 1300.0 # [V] - Grid voltage
nn = 1.0e18 # [m^-3] - Neutral density intragrid
xa = 0.5e-3 # [m]<-[mm] - Grid spacing
ra = 0.5e-3 # [m]<-[mm] - Aperture radius
wa = 0.5e-3 # [m]<-[mm] - Grid thickness
Iap = 2.5e-4 # [A] - Current per aperture
Y = 1.0 # [atom/ion] - Sputter yield
IB = 12.2 # [eV] - First ionization energy of atomic Xe
U_ex1 = 11.60 # [eV] excitation energy of xenon
e = 1.6e-19 # [C] - fundamental unit of electric charge

print("Problem 1 Computations")

r_ABmax = (10**(19)/1.6)*3.0*(e**2)/(4.0*np.pi*epsilon0*IB) # [m] - maximum distance between interacting Xe
sigma_cex = np.pi * r_ABmax**2 # [m^2] - classical cross section for resonant charge exchange
I_cex = Iap*nn*sigma_cex*xa # [A] - Current generated via charge-exchange collisions per aperture
print("The total current generated via charge-exchange collisions related to the aperture is %.2E [A]" %I_cex)

nsputtering_dot = (I_cex*Y)/e # [atom/(ion*s)] - Sputtering rate
print("The sputtering rate is %.2E [particles/s]" %nsputtering_dot)

mg = 95.9*amutokg # [kg] - grid atomic mass
rhog = 10.28*10**(-3)*10**(6) # [kg/m^3]<-[g/cm^3] - grid density
Va_dot = nsputtering_dot/(rhog/mg) # [m^3/s] - Volumetric loss rate of grid material
print("The volumetric loss rate is %.2E [m^3/s]" %Va_dot)

ra_dot = Va_dot/(2.0*np.pi*ra*wa) # [m/s] - Rate of increase of grid aperture radius
print("The rate of increase of grid aperture radius is %.2E [m/s]" %ra_dot)

deltat = (1.0/3600.0)*0.5*ra/ra_dot # [hours]<-[s] - Lifetime of thruster based on assumption of resistance to electron backstreaming
print("The lifetime of the thruster is %f [hours]" %deltat)

""" Problem 2 """
eta_b = 0.75 # Current efficiency
eta_V = 0.9 # Voltage efficiency
Id = 10.0 # [A] - Discharge current
Vd = 300.0 # [V] - Discharge voltage
rout = 7.0 # [cm] - Outer wall radius
rin = 3.0 # [cm] - Inner wall radius
L = 4.0 # [cm] - Plasma length
m_Xe = 131.0*amutokg # [kg] - Atomic mass of Xe
me = 9.11e-31 # [kg] - electron mass
a = 0.123 # 2ndary e^{-} emission coefficients for BNSiO_{2}
b = 0.528 # "
A_in = 2.0*np.pi*(rin*10**(-2))*(L*10**(-2)) # [m^2]
A_out = 2.0*np.pi*(rout*10**(-2))*(L*10**(-2)) # [m^2]

print("Problem 2 Computations")

Te = np.linspace(1.0,50.0,100) # [eV]
GammaProd = sp.gamma(2.0 + b)*a # part of the secondary electron emission yield parameter
v_is = np.sqrt(1.6*10**(-19)*Te/m_Xe) # [m/s], Bohm velocity of ions

Te_trans = ((1.0/GammaProd)*(1.0-np.exp(1.02)*np.sqrt(me/(2.0*m_Xe))))**(1.0/b) # [eV] ?
print("The temperature at which the sheath transitions to a space-charge limited region is %f [eV]" %Te_trans)

R_ion = 10.0**(-20)*(-(1.031*10**(-4))*Te**(2) + 6.386*np.exp(-(IB)/Te))*np.sqrt((8.0*e*Te)/(np.pi*me)) # []
R_ex1 = 1.931*10.0**(-19)*(np.exp(-U_ex1/Te)/np.sqrt(Te))*np.sqrt((8.0*e*Te)/(np.pi*me)) # []
U_Ieff = IB + (R_ex1/R_ion)*U_ex1 # [eV]

Ac = np.pi*(rout**2 - rin**2)*10**(-4) # [m^2] - Channel area
n0 = (eta_b*Id)/(e*Ac*np.sqrt(2.0*e*eta_V*Vd/m_Xe)) # [m^-3] - Center channel density
ns = n0/2.0 # [m^-3]
print("The center channel density is %.2E [m^-3]" %n0)

As = A_in + A_out # [m^2] Total area where sheath is interacting with wall
I_iw = ns*e*v_is*As # [A] ion current to wall

Pd = Vd*Id*np.ones(Te.shape[0]) # [W]
Pd_norm = Pd/Pd # [dimless]
Pb_norm = eta_b*eta_V*Id*Vd/Pd # [dimless]
Pa_norm = (1.0/Pd)*2.0*Id*(1.0/e)*(1.0/6.0)*1.6*10**(-19)*Te # [dimless]
PrPion_norm = (1.0/Pd)*(1.0/e)*1.6*10**(-19)*U_Ieff*(eta_b*Id+I_iw) # [dimless]

phi_rel = np.zeros(Te.shape[0])
# Pwe_norm = np.zeros(Te.shape[0])
# Pwi_norm = np.zeros(Te.shape[0])

for pidx in np.arange(phi_rel.shape[0]):
    if Te[pidx] <= Te_trans: # Not yet at space-charge limited region
        phi_rel[pidx] = -((1.6*10**(-19)*Te[pidx])/e)*np.log((1.0-GammaProd*(Te[pidx])**b)*np.sqrt(2.0*m_Xe/me)) # []
    elif Te[pidx] > Te_trans:
        phi_rel[pidx] = -1.02*1.6*10**(-19)*(Te[pidx])/e # []

Pwe_norm = (1.0/Pd)*(1.0/4.0)*np.sqrt((8.0*1.6*10**(-19)*Te)/(np.pi*me))*ns*As*np.exp(e*phi_rel/(1.6*10**(-19)*Te))*2.0*1.6*10**(-19)*Te # []
Pwi_norm = (1.0/Pd)*ns*v_is*As*(0.5*m_Xe*v_is**2 - e*phi_rel) # []

powerFig = plt.figure()
plt.semilogy(Te,Pd_norm,label='Discharge')
plt.semilogy(Te,Pb_norm,label='Beam')
plt.semilogy(Te,Pa_norm,label='Anode')
plt.semilogy(Te,PrPion_norm,label='Radiation+Ionization')
plt.plot(Te,Pwe_norm,label='Wall loss - Electron')
plt.semilogy(Te,Pwi_norm,label='Wall loss - Xenon')
plt.semilogy(Te,Pb_norm+Pa_norm+PrPion_norm+Pwi_norm+Pwe_norm,label='Sum')
plt.legend()
plt.xlabel('Electron Temperature [eV]')
plt.ylabel('Normalized Power')
plt.title('Power losses in a Xenon Hall Thruster')
plt.xlim((Te[0],Te[Te.shape[0]-1]))

""" Problem 3 """

print("Problem 3 Computations")

""" Problem 4 """

print("Problem 4 Computations")

# """ Problem 5 """
#
# print("Problem 5 Computations")
plt.show()
