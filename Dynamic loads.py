import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

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
print(M_axial)

K_panel_1 = E*w_1*t_p/L_1     #Single panel stiffness in N/m
K_panel_2 = E*w_2*t_p/L_2     #Single panel stiffness in N/m
K_panel_3 = E*w_1*t_p/w_2
K_panel_4 = E*L_1*t_p/w_2
K_panels_axial = 2*K_panel_1 + 2*K_panel_2    #4 Panel stiffness in N/m
K_panels_lateral = 2*K_panel_3 + 2*K_panel_4
K_Tanks = E*(np.pi*(r_outer_tanks**2-(r_outer_tanks-t_tanks)**2))/L_1*4
K_total_axial = K_panels_axial + K_Tanks  #Total stiffness in N/m
K_total_lateral = K_panels_lateral

# === Damping ===
zeta = 0.01  # 1% damping ratio
omega_n_axial = np.sqrt(K_total_axial / M_axial)  # rad/s
C_axial= 2 * zeta * np.sqrt(K_total_axial * M_axial)  # Ns/m
omega_n_lateral = np.sqrt(K_total_lateral/ M_lateral)  # rad/s
C_lateral= 2 * zeta * np.sqrt(K_total_lateral * M_lateral)  # Ns/m

# === Forcing: 1g sinusoidal acceleration at 100 Hz ===
f_drive = 100  # Hz
omega_drive = 2 * np.pi * f_drive  # rad/s
A_force = M_axial * 9.81  # N

# === Time Domain Simulation ===
t_span = (0, 0.05)  # simulate for 50 ms
t_eval = np.linspace(t_span[0], t_span[1], 10000)  # high-res output

def systemaxial(t, y):
    x, v = y  # displacement and velocity
    a_t = A_force * np.sin(omega_drive * t)
    dxdt = v
    dvdt = (a_t - C_axial * v - K_total_axial * x) / M_axial
    return [dxdt, dvdt]

def systemlateral(t, y):
    x, v = y  # displacement and velocity
    a_t = A_force * np.sin(omega_drive * t)
    dxdt = v
    dvdt = (a_t - C_lateral* v - K_total_lateral * x) / M_lateral
    return [dxdt, dvdt]

# Initial conditions: [displacement, velocity]
y0 = [0, 0]

# Integrate
sol_axial = solve_ivp(systemaxial, t_span, y0, t_eval=t_eval, method='RK45')
sol_lateral = solve_ivp(systemlateral, t_span, y0, t_eval=t_eval, method='RK45')

# === Plot Results Axial ===
plt.figure(figsize=(10, 5))
plt.plot(sol_axial.t*1000, sol_axial.y[0], label='Displacement (m)', color='blue')
plt.xlabel('Time (ms)')
plt.ylabel('Axial Displacement (m)')
plt.title('Axial Displacement Response to 1g Sine Acceleration at 100 Hz')
plt.grid(True)
plt.tight_layout()
plt.show()


# === Plot Results Lateral ===
plt.figure(figsize=(10, 5))
plt.plot(sol_lateral.t*1000, sol_lateral.y[0], label='Displacement (m)', color='blue')
plt.xlabel('Time (ms)')
plt.ylabel('Lateral Displacement (m)')
plt.title('Lateral Displacement Response to 1g Sine Acceleration at 100 Hz')
plt.grid(True)
plt.tight_layout()
plt.show()