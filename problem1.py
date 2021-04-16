"""
Matt Russell
AA529: Space Propulsion
HW1 Problem1
Performance Estimate for Cold Gas Thruster
"""
""" Function Declarations """
def ExpRat_expression(Me):
    return ((gamma_N2 + 1.0) / 2.0)**(-(gamma_N2 + 1.0) / (2.0 * (gamma_N2 - 1.0))) \
* (1.0 / Me)*(1.0 + ((gamma_N2 - 1.0) / 2.0)*Me**2) \
**((gamma_N2 + 1.0) / (2.0 * (gamma_N2 - 1.0)))

""" Bulk """
import numpy as np
import matplotlib.pyplot as plt

""" Known Parameters """
k_b = 1.38e-23 # Boltzmann's constant [J/K]
ExpRat = 40.0 # Expansion Ratio: Ae/At []
p_c = 260.0 * 6894.757 # Inlet (Chamber) Pressure: [psia]->[Pa]
T_c = 20.0 + 293.150 # Operating Temperature: [deg C]->[deg K]
d_t = 1.6e-3 # Throat Diameter: [mm]->[m]
A_t = np.pi*(d_t/2.0)**(2) # Throat Area: [m]
A_e = ExpRat*A_t
gamma_N2 = 1.40 # Ratio of Specific Heat of Nitrogen gas at 300 deg K
m_p = 2.0*14.007*1.66e-27 # Atomic mass of N2
g_0 = 9.8 # [m/s^2]

""" Calculate Mach Number, M_{e} = Me """
# num_guesses = 500
# Me_guesses = np.linspace(1.0,10.0,num_guesses)
# ExpRat_expression = ((gamma_N2 + 1.0) / 2.0)**(-(gamma_N2 + 1.0) / (2.0 * (gamma_N2 - 1.0))) \
# * (1.0 / Me_guesses)*(1.0 + ((gamma_N2 - 1.0) / 2.0)*np.multiply(Me_guesses,Me_guesses)) \
# **((gamma_N2 + 1.0) / (2.0 * (gamma_N2 - 1.0)))
# Me = 5.625 # Determined empirically
# ExpRat_guessed = ((gamma_N2 + 1.0) / 2.0)**(-(gamma_N2 + 1.0) / (2.0 * (gamma_N2 - 1.0))) \
# * (1.0 / Me)*(1.0 + ((gamma_N2 - 1.0) / 2.0)*Me**2) \
# **((gamma_N2 + 1.0) / (2.0 * (gamma_N2 - 1.0)))
Me_low = 1.0
Me_high = 10.0
Me_guess = 0.5*(Me_low + Me_high)
eps = 1e-13
ExpRat_check = ExpRat_expression(Me_guess)
diff = ExpRat_check - ExpRat
i = 0
limit = 1000
while np.abs(diff) > eps:
    print("i is %i" %i)
    print("diff is %16.14f" %diff)
    print("Me_low is %f" %Me_low)
    print("Me_high is %f" %Me_high)
    print("Me_guess is %16.14f" %Me_guess)
    if i > limit:
        break
    if diff > 0.0: # Guess too large
        Me_high = Me_guess
        Me_guess = 0.5*(Me_low + Me_high)
    if diff < 0.0: # Guess too small
        Me_low = Me_guess
        Me_guess = 0.5*(Me_low + Me_high)
    i+=1
    ExpRat_check = ExpRat_expression(Me_guess)
    diff = ExpRat_check - ExpRat

Me = Me_guess

print("Part (a) The estimated exhaust Mach Number is %16.14f" %Me)
print("The Expansion Ratio using the estimated Mach number is %16.14f" %ExpRat_check)

""" Calculate exhaust Pressure, p_{e} in [Pa]"""
p_e = p_c*(1.0 + ((gamma_N2 - 1.0) / 2.0)*Me**2)**(-gamma_N2/(gamma_N2-1.0))

print("Part (b) The exhaust pressure is %7.3f [Pa]" %p_e)

""" Calculate exhaust speed, c_{e} in [m/s]"""
c_e = np.sqrt(((k_b*T_c)/m_p)*(2.0*gamma_N2)/(gamma_N2 - 1.0)*(1.0 - (p_e/p_c)**((gamma_N2-1.0)/gamma_N2)))

print("Part (c) The exhaust speed is %6.3f [m/s]" %c_e)

""" Calculate Mass flow rate, \dot{m} [kg/s] """
m_dot = A_t*p_c*np.sqrt(gamma_N2*m_p/(k_b*T_c))*((gamma_N2 + 1.0)/2.0)**(-(gamma_N2+1.0)/(2.0*(gamma_N2 - 1.0)))

print("Part (d) The mass flow rate is %4.3f [kg/s]" %m_dot)

""" Calculate Thrust """
F = m_dot*c_e + (p_e - 0.0)*A_e # "0.0" represents vacuum pressure

print("Part (e) The thrust is %4.3f [N]" %F)

""" Calculate Specific Impulse """
I_sp = F/(m_dot*g_0)

print("Part (f) The specific impulse is %5.3f [s]" %I_sp)

""" Plotting """
# Part (a): exhaust Mach number
# Me_fig = plt.figure()
# plt.plot(Me_guesses,ExpRat_expression,'r',label='Expansion Ratio($M_{e}$)')
# plt.axhline(y=ExpRat,label='$\\frac{A_{e}}{A_{t}} = %3.1f$' %ExpRat)
# plt.title('$\\frac{A_{e}}{A_{t}}$ vs. $M_{e}$')
# plt.xlabel('$M_{e}$')
# plt.ylabel('Expansion Ratio')
# plt.legend()

# plt.show()
