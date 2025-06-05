import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


Total_mass = 1200
g_axial = 8.5
g_lateral = 3
g = 9.80665
M_fuel_tank = 200       #Propellant tank mass (fueld) in kg
M = 100                 #Mass supported by the side panles at launch
E = 114 *10**9          #Panel elastic module in Pa
t_p = 0.005               #Panel thickness in m
w_1 = 1                 #Panel 1 width in m
L_1 = 1                 #Panel 1 length in m (height)
w_2 = 1                 #Panel 2 width in m
L_2 = 1                 #Panel 2 length in m (height)

rho_panels = 4429
r_outer_tanks = 0.3
t_tanks = 0.01

r_inner = r_outer_tanks - t_tanks #tank inner radius m
A_axial = 2*w_1*t_p + 2*w_2*t_p + 4* (2*np.pi*r_outer_tanks - 2*np.pi*r_inner) #area that holds axial loads
#print(A_axial)
A_support = 0 #area of the support beam
A_lateral = 2*w_2*t_p + 2*L_1*t_p + 2*A_support  #area that holds lateral loads

M_t_full = 150 #ull fuel tank mass kg
M_axial = w_1*w_2*t_p*rho_panels + rho_panels*(2*w_1*L_1*t_p + 2*w_2*L_1*t_p) + 4*M_t_full   #Mass carried in the axial direction kg
M_lateral = 2*M_t_full + L_1*w_2*t_p*rho_panels +  rho_panels*(2*w_1*w_2*t_p + 2*w_1*L_1*t_p)     #Mass carried in the lateral direction kg


K_panel_1 = E*w_1*t_p/L_1     #Single panel stiffness in N/m
K_panel_2 = E*w_2*t_p/L_2     #Single panel stiffness in N/m
K_panels = 2*K_panel_1 + 2*K_panel_2    #4 Panel stiffness in N/m
K_Tanks = E*(np.pi*(r_outer_tanks**2-(r_outer_tanks-t_tanks)**2))/L_1*4
K_total = K_panels + K_Tanks  #Total stiffness in N/m
omega_total = (1/(2*np.pi))*np.sqrt(K_total/Total_mass)  #Total natural frequency in Hz
print (f'Total natural frequency: {omega_total:.2f} Hz')



# panels_submasses = np.linspace(0, M_panels, 100)
# tanks_submasses = np.linspace(0, M_fuel_tank, 100)
# for panel_submass in panels_submasses:
#     for tank_submass in tanks_submasses:
A = g
omega = 2*np.pi*100     

# ODE system: convert 2nd-order ODE to two 1st-order ODEs
def mass_spring_forced(t, y):
    x, v = y
    F_t = A * np.sin(omega * t)    # external force
    dxdt = v
    dvdt = (F_t - K_total * x) / M_axial
    return [dxdt, dvdt]

# Initial conditions
y0 = [0.0, 0.0]  # x(0) = 0, v(0) = 0

# Time settings
t_span = (0, 10)
t_eval = np.linspace(*t_span, 10000)

# Solve ODE
sol = solve_ivp(mass_spring_forced, t_span, y0, t_eval=t_eval)

# Plot results
plt.figure(figsize=(10, 5))
plt.plot(sol.t, sol.y[0], label='Displacement x(t)')
plt.plot(sol.t, A * np.sin(omega * sol.t), '--', label='Forcing Function F(t)', alpha=0.6)
plt.xlabel('Time [s]')
plt.ylabel('Response')
plt.title('Forced Mass-Spring System')
plt.grid(True)
plt.legend()
plt.show()