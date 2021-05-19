"""
Matt Russell
University of Washington Aeronautics & Astronautics
AA529: Space Propulsion HW3
5/21/21
"""
import numpy as np
import matplotlib.pyplot as plt

epsilon0 = 8.854e-12 # [F/m] - vacuum permittivity

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

n_sputtering = (I_cex*Y)/e # [] - Sputtering rate
print("The sputtering rate is %f [particles/s]" %n_sputtering)

mg = # [] - grid atomic mass
rhog = # [] - grid density
Va = n_sputtering/(rhog/mg) # [] - Volumetric loss rate of grid material
print("The volumetric ")

""" Problem 2 """
eta_b = 0.75 # Current efficiency
eta_V = 0.9 # Voltage efficiency
Id = 10.0 # [A] - Discharge current
Vd = 300.0 # [V] - Discharge voltage
rout = 7.0 # [cm] - Outer wall radius
rin = 3.0 # [cm] - Inner wall radius
L = 4.0 # [cm] - Plasma length
a = 0.123 # 2ndary e^{-} emission coefficients for BNSiO_{2}
b = 0.528 # "

print("Problem 2 Computations")

""" Problem 3 """

print("Problem 3 Computations")

""" Problem 4 """

print("Problem 4 Computations")

""" Problem 5 """

print("Problem 5 Computations")
