"""
Matt Russell
AA529: Space Propulsion, HW4
University of Washington Department of Aeronautics & Astronautics
6/4/21
"""
import numpy as np
import matplotlib.pyplot as plt

""" Functions """
def CalculateBeta(Rc,Rp,C,L0):
    beta = (Rc + Rp) * np.sqrt(C/L0)
    return beta

def CalculateAlpha(C,V0,m0,xdot,Ldot,Rc,Rp,beta):
    alpha = 0.5 * (0.5 * C * V0**2)/(0.5 * m0 * xdot**2) \
        *(Ldot/(Rc + Rp))**2 * beta**2
    return alpha

""" Constants """
amu_to_kg = 1.667e-27 # 1 [amu] = this many [kg]
Joules_to_eV = (1.0/1.6)*10.0**19 # 1 [J] = this man [eV]
Pascal_to_Atm = 9.869e-6 # 1 [Pa] = this many [atm]
lbf_to_Newton = 4.448 # 1 [lbf] = this many [N]
degF_to_degK = 255.927 # 1 [degF] = this many [degK]
m_Xe = 131.0*amu_to_kg # [kg]<-[amu]
g0 = 9.819 # [m/s^2], gravitational acceleration at sea level
mu_0 = 4*np.pi*10**(-7) # [H/m]

""" Problem 1 """
eta_T1 = 0.6 # [dimless], thrust efficiency
Isp_1 = 1200.0 # [s]
KE_ion = (m_Xe/2.0)*(g0*Isp_1)**2 # [J]
E_cost = Joules_to_eV*KE_ion*(1.0 - eta_T1)/eta_T1 # [eV/ion]<-[J/ion]

print("P1: The cost for a Hall/ion thruster to operate at this level is %f [eV] per ion" %E_cost)

""" Problem 2 """
rc = 0.02 # [m]
ra = 0.07 # [m]
J = 20.0e3 # [A]
lenr = 100

P2Fig_Btheta = plt.figure()
P2Fig_pdiff = plt.figure()

r2 = np.linspace(0.01,0.1,lenr) # [m]
r2c = np.linspace(1.0e-3*rc,rc,lenr) # [m]
Btheta_ofr = np.empty((lenr)) # Piecewise function
pdiff = np.empty((lenr)) # Pressure differential, p - p_{0}

for ridx in np.arange(lenr): # Calculate azimuthal magnetic field and pressure gradient
    pdiff[ridx] = 9.869e-6*(mu_0*J**2)/(4.0*np.pi**2*rc**2)*(1.0 - r2c[ridx]**2/rc**2) # [Pa]<-[atm]
    if r2[ridx] < rc:
        Btheta_ofr[ridx] = (mu_0*J)/(2.0*np.pi*rc**2)*r2[ridx] # [T]
    elif r2[ridx] > rc:
        Btheta_ofr[ridx] = (mu_0*J)/(2.0*np.pi*r2[ridx]) # [T]

CT = np.log(ra/rc) + 0.75 # Maeker Equation
print("P2a: The thrust coefficient is %f" %CT)

Bmax_find = np.asarray(Btheta_ofr == np.amax(Btheta_ofr)).nonzero()
Bmax_idx = Bmax_find[0] # Btheta has a global maxima
Btheta_max = Btheta_ofr[Bmax_idx]
print("P2b: The maximum field is %f [T]" %Btheta_max)
print("P2b: The maximum field strength is located at r = %f [m]" %r2[Bmax_idx])

pdiffmax_find = np.asarray(pdiff == np.amax(pdiff)).nonzero()
pdiff_idx = pdiffmax_find[0] # pdiff has a global maxima
pdiff_max = pdiff[pdiff_idx]
print("P2c: The maximum pressure gradient is %f [atm]" %pdiff_max)
print("P2c: The location of maximum pressure gradient is r = %3.2E [m]" %r2c[pdiff_idx])

F_blowing = (mu_0*J**2)/(4.0*np.pi)*(np.log(ra/rc) + 1.0) # [N]
F_pumping = (mu_0*J**2)/(8.0*np.pi) # [N]
Thrust = F_blowing + F_pumping # [N]

print("P2d: The blowing force is %f [N]" %F_blowing)
print("P2e: The pumping force is %f [N]" %F_pumping)
print("P2f: The thrust is %f [N]" %Thrust)

""" Problem 3 """
# [Thruster A, Thruster B, Thruster C]
xdot = 5e4 # [m/s], Exhaust Speed
Rc = np.array([5.0, 5.0, 50.0])*10**(-3) # [Ohm]<-[mOhm], Circuit Resistance
Rp = np.array([5.0, 5.0, 5.0])*10**(-3) # [Ohm]<-[mOhm], Plasma Resistance
Lc = np.array([10.0, 10.0, 10.0])*10**(-9) # [H]<-[nH], Stray Inductance
C = np.array([10.0,10.0,10.0])*10**(-6) # [F]<-[uF], Capacitance
V0 = np.array([3.0, 3.0, 7.0])*10**(3) # [V]<-[kV], Discharge Voltage
ChannelAspectRatio = np.array([0.2, 0.05, 0.05]) #
ChannelLength = np.array([20.0, 20.0, 20.0])*10**(-2) # [m]<-[cm]
SlugMass = np.array([10.0, 10.0, 10.0])*10**(-6) # [g]<-[ug]
alphas = np.empty((3))
betas = np.empty((3))

for idx in np.arange(alphas.size): # Calculate for A,B, and C
    Lpdot = mu_0*ChannelAspectRatio[idx]*xdot
    DeltaLp = mu_0*ChannelAspectRatio[idx]*ChannelLength[idx]
    Lovberg = DeltaLp/Lc[idx]
    betas[idx] = CalculateBeta(Rc[idx],Rp[idx],C[idx],Lc[idx]) # Lp(0) = 0
    alphas[idx] = CalculateAlpha(C[idx],V0[idx],SlugMass[idx],xdot,Lpdot,Rc[idx],Rp[idx],betas[idx])
    print("P3a: Thruster %i alpha = %f" %(idx+1,alphas[idx]))
    print("P3a: Thruster %i beta = %f" %(idx+1,betas[idx]))
    print("P3c: Thruster %i Lovberg parameter is %f" %(idx+1,Lovberg))

""" Problem 4 """
P = 90.0 # [W], Power input to 1965 Vela Resistojet
T = 0.042*lbf_to_Newton # [N]<-[lbf], Thrust
eta_T4 = 0.65 # Ad hoc thrust efficiency
Temp_Heater = 1000.0*degF_to_degK # [degK]<-[degF], Heater is Resistance Element
Cp_N2 = 1.32e3 # [J kg^-1 degK^-1] Specific Heat Constant Pressure  ~3000 deg K
ue = np.sqrt(2.0*Cp_N2*Temp_Heater) # [m/s], exhaust speed
# Comparing two different ways of calculating
Isp1 = ue/g0
print("P4: Estimated Specific Impulse Nozzle Flow Way is %f [s]" %Isp1)
Isp2 = 2.0*(eta_T4/g0)*(P/T)
print("P4: Estimated Specific Impulse Thrust Efficiency is %f [s]" %Isp2)

""" Plotting """
plt.figure(P2Fig_Btheta.number)
plt.plot(r2,Btheta_ofr)
plt.axvline(x = r2[Bmax_idx], ls = '--')
plt.xlabel('Chamber radius [m]')
plt.ylabel('Magnetic Induction [T]')
plt.title('Field Strength in an MPDT')
plt.xlim((r2[0],r2[r2.size-1]))
plt.ylim((0.0,Btheta_max))

plt.figure(P2Fig_pdiff.number)
plt.plot(r2c,pdiff)
plt.axvline(x = r2c[pdiff_idx], ls = '--')
plt.xlabel('Cathode radius [m]')
plt.ylabel('Pressure Differential [atm]')
plt.title('MPDT Chamber Pressure Profile')
plt.xlim((r2c[0],r2c[r2c.size-1]))
plt.ylim((0.0,pdiff_max))

plt.show()
