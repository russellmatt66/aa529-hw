"""
Matt Russell
AA529: Space Propulsion
HW1 Problem 4
Brief analysis of Xenon Hall Thruster
"""
import math
""" Part (a): firing time calculation """
deltav = 3.0 # km/s
m_object = 1000.0 # kg
g_0 = 9.81 # [m/s^{2}]
I_sp = 2000.0 # [s]
eta_thrust = 0.6 # Thrust Efficiency of EP device
P_device = 5.0 # [kW]

t_firing = 1.0e-3 * m_object * (g_0 * I_sp)**2/(2.0 * eta_thrust * P_device) # [s]
mdot_p = m_object/t_firing # [kg/s]
t_firing = (1.0/3600.0)*(t_firing * (1.0 - math.exp(-(1.0e3*deltav/(g_0 * I_sp))))) # [hours]
print("The thruster would need to fire for %f hours" %t_firing)

""" Part (b): Mass of propellant """
m_p = mdot_p * 3600.0 * t_firing # [kg]
print("The propellant mass needed is %f kg" %m_p)
