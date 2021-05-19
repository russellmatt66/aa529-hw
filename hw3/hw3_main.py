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
amutokg = 1.67e-27 # 1 amu = this many kg

""" Problem 1 """
V0 = 1300.0 # [V] - Grid voltage
nn = 10.0e18 # [m^3] - Neutral density intragrid
xa = 0.5 # [mm] - Grid spacing
ra = 0.5 # [mm] - Aperture radius
wa = 0.5 # [mm] - Grid thickness
Iap = 2.5e-4 # [A] - Current per aperture
Y = 1.0 # [atom/ion] - Sputter yield
IB = 12.2 # [eV] - First ionization energy of atomic Xe
e = 1.6e-19 # [C] - fundamental unit of electric charge

print("Problem 1 Computations")

r_ABmax = 3.0*(e**2)/(4.0*np.pi*epsilon0*IB) # [] - maximum distance between interacting Xe
sigma_cex = np.pi * r_ABmax**2 # [] - classical cross section for resonant charge exchange
I_cex = Iap*nn*sigma_cex*xa # [] - Current generated via charge-exchange collisions per aperture
print("The total current generated via charge-exchange collisions related to the aperture is %f [A]" %I_cex)

nsputtering_dot = (I_cex*Y)/e # [] - Sputtering rate
print("The sputtering rate is %f [particles/s]" %nsputtering_dot)

mg = 95.9*amutokg # [kg] - grid atomic mass
rhog = # [] - grid density
Va_dot = nsputtering_dot/(rhog/mg) # [] - Volumetric loss rate of grid material
print("The volumetric loss rate is %f [m^3/s]" %Va_dot)

ra_dot = Va_dot/(2.0*np.pi*ra*wa) # [] - Rate of increase of grid aperture radius
print("The rate of increase of grid aperture radius is %f [m/s]" %ra_dot)

deltat = 0.5*ra/ra_dot # [] - Lifetime of thruster based on assumption of resistance to electron backstreaming
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

print("Problem 2 Computations")

Te = np.linspace(0.1,100.0,1000)
GammaProd = sp.gamma(2.0 + b)*a # part of the secondary electron emission yield parameter
Te_trans = ((1.0/GammaProd)*(1.0-np.exp(1.02)*np.sqrt(me/(2.0*m_Xe))))**(1.0/b)
print("The temperature at which the sheath transitions to a space-charge limited region is %f []" %Te_trans)

Ac = np.pi*(rout - rin)**2 # [] - Channel area
n0 = (eta_b*Id)/(e*Ac*np.sqrt(2.0*e*eta_V*Vd/mi)) # [] - Center channel density
print("The center channel density is %f []" %n0)

phi_rel = np.empty(Te.shape[0])

""" Problem 3 """

print("Problem 3 Computations")

""" Problem 4 """

print("Problem 4 Computations")

""" Problem 5 """

print("Problem 5 Computations")
