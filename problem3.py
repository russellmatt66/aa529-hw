"""
Matt Russell
AA529: Space Propulsion
HW1 Problem 3
Feasibility analysis of Electric Propulsed (EP) Manned Maneuvering Unit (MMU)
"""
import math
""" Required Power """
a_MMU = 0.04 # m/s^{2}
eta_T = 0.6 # Thrust efficiency
g_0 = 9.8 # m/s^{2}
I_sp = 3000.0 # [s]
m_T = 213.63 # mass of astronaut+MMU [kg]
P_req = (a_MMU*g_0*I_sp*m_T)/(2.0*eta_T) # [Watts]
print("The necessary power is %9.2f Watts" %P_req)

""" Solar Array Calculations """
# Array numbers based on ISS solar panels
P_avg = 100.0 # average power from entire array [kW]
N_cells = 262400 # Number of solar cells that constitute entire array
A_total = 2500.0 # Total area of solar array [m^{2}]
P_cellavg = P_avg/float(N_cells) # Average power from a single solar cell
A_cell = A_total/float(N_cells)

""" MMU Power system analysis """
N_req = math.ceil(1e-3*P_req/P_cellavg) # P_req in [Watts]
print("The number of solar arrays required is %i" %N_req)
A_req = float(A_cell*N_req)
print("This system would occupy an area of %f meters squared" %A_req)
