import numpy as np
import matplotlib.pyplot as plt

me = 9.11*10**(-31) # [kg]
IB = 12.2 # [eV] - First ionization energy of atomic Xe
U_ex1 = 11.60 # [eV] excitation energy of xenon
e = 1.6e-19 # [C] - fundamental unit of electric charge
Te = np.linspace(1.0,50.0,100) # [eV]
R_ion = 10.0**(-20)*(-(1.031*10**(-4))*Te**(2) + 6.386*np.exp(-(IB)/Te))*np.sqrt((8.0*e*Te)/(np.pi*me)) # []
R_ex1 = 1.931*10.0**(-19)*(np.exp(-U_ex1/Te)/np.sqrt(Te))*np.sqrt((8.0*e*Te)/(np.pi*me)) # []
U_Ieff = IB + (R_ex1/R_ion)*U_ex1 # [eV]

plt.plot(Te,U_Ieff)
plt.show()
