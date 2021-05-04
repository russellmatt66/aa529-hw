"""
Matt Russell
University of Washington Department of Aeronautics & Astronautics
AA529: Space Propulsion
5/5/21
Software for the midterm
"""
import numpy as np
import matplotlib.pyplot as plt
# from aa529midterm_main import findIndex

Tc_low = 300.0
Tc_high = 1500.0
num_T = 100
T = np.linspace(Tc_low,Tc_high,num_T)

""" jguess' type - [FIXED]"""
# jlow = 0
# jhigh = T.shape[0]-1
# jguess = int(np.floor((jlow + jhigh)/2))
# print(type(jguess))
# print(T[jguess])

""" findIndex debugging """
# T_instance =

# plt.plot(T)
# plt.show()
