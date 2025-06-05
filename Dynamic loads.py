import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# === System Parameters ===
Total_mass = 1200  # kg
E = 114e9  # Pa
t_p = 0.005  # m
w_1 = w_2 = 1.0  # m
L_1 = L_2 = 1.0  # m
rho_panels = 4429  # kg/m^3
r_outer_tanks = 0.3
t_tanks = 0.01
r_inner = r_outer_tanks - t_tanks
M_t_full = 150  # full tank mass (each) in kg

# === Stiffness Calculations ===
K_panel_1 = E * w_1 * t_p / L_1
K_panel_2 = E * w_2 * t_p / L_2
K_panels = 2 * K_panel_1 + 2 * K_panel_2

A_tank = np.pi * (r_outer_tanks**2 - r_inner**2)
K_tanks = E * A_tank / L_1 * 4
K_total = K_panels + K_tanks  # N/m

# === Damping ===
zeta = 0.01  # 1% damping ratio
omega_n = np.sqrt(K_total / Total_mass)  # rad/s
C = 2 * zeta * np.sqrt(K_total * Total_mass)  # Ns/m

# === Forcing: 1g sinusoidal acceleration at 100 Hz ===
f_drive = 100  # Hz
omega_drive = 2 * np.pi * f_drive  # rad/s
A_force = Total_mass * 9.81  # N

# === Time Domain Simulation ===
t_span = (0, 0.05)  # simulate for 50 ms
t_eval = np.linspace(t_span[0], t_span[1], 10000)  # high-res output

def system(t, y):
    x, v = y  # displacement and velocity
    a_t = A_force * np.sin(omega_drive * t)
    dxdt = v
    dvdt = (a_t - C * v - K_total * x) / Total_mass
    return [dxdt, dvdt]

# Initial conditions: [displacement, velocity]
y0 = [0, 0]

# Integrate
sol = solve_ivp(system, t_span, y0, t_eval=t_eval, method='RK45')

# === Plot Results ===
plt.figure(figsize=(10, 5))
plt.plot(sol.t * 1000, sol.y[0], label='Displacement (m)', color='blue')
plt.xlabel('Time (ms)')
plt.ylabel('Displacement (m)')
plt.title('Displacement Response to 1g Sine Acceleration at 100 Hz')
plt.grid(True)
plt.tight_layout()
plt.show()