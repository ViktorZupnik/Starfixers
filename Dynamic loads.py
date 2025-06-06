import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


g = 9.80665

def bending_stress_at_x(
    x,            # Position along tank [m]
    total_mass,   # Tank total mass [kg]
    beam_length,  # Tank height [m]
    r,            # Tank radius [m]
    t,            # Tank thickness [m]
    g_lateral     # Lateral acceleration  [g] (1g = 9.80665 m/s^2)
):
    w = (g_lateral * total_mass * g) / beam_length  # distributed load [N/m]

    # Bending moment at position x
    Mx = (w * beam_length / 12) * (6 * x - beam_length) - (w * x**2) / 2

    # Section properties
    I = np.pi*(r**4-(r-t)**4)/4

    # Bending stress at x
    stress = (Mx * r) / I  # in Pascals
    return stress


print (bending_stress_at_x(0, 150, 1.065, 0.533/2, 0.002, 3))

g_axial = 8.5
g_lateral = 3
M_fuel_tank = 200       #Propellant tank mass (fueld) in kg
M = 100                 #Mass supported by the side panles at launch
E = 114 *10**9          #Panel elastic module in Pa
t_p = 0.002              #Panel thickness in m
w_1 = 1.13                 #Panel 1 width in m
L_1 = 1.065                 #Panel 1 length in m (height)
w_2 = 1.13                 #Panel 2 width in m
L_2 = L_1                 #Panel 2 length in m (height)
r_outer_rod = 0.02
t_rod = 0.001

rho_panels = 2810
r_outer_tanks = 0.533/2
t_tanks = 0.002

r_inner = r_outer_tanks - t_tanks #tank inner radius m
A_axial = 2*w_1*t_p + 2*w_2*t_p + 4* (2*np.pi*r_outer_tanks - 2*np.pi*r_inner) #area that holds axial loads
#print(A_axial)
A_support =  np.pi*(r_outer_rod**2-(r_outer_rod-t_rod)**2) #area of the support beam
A_lateral = 2*w_2*t_p + 2*L_1*t_p + 2*A_support  #area that holds lateral loads

M_t_full = 150 #ull fuel tank mass kg
M_rod = A_support*(w_1-4*r_outer_tanks)*rho_panels
print(M_rod)
M_axial = w_1*w_2*t_p*rho_panels + rho_panels*(2*w_1*L_1*t_p + 2*w_2*L_1*t_p) + 4*M_t_full   #Mass carried in the axial direction kg
M_lateral = 2*M_t_full + L_1*w_2*t_p*rho_panels +  rho_panels*(2*w_1*w_2*t_p + 2*w_1*L_1*t_p) + 2*M_rod    #Mass carried in the lateral direction kg
print(M_axial)


K_panel_1 = E*w_1*t_p/L_1     #Single panel stiffness in N/m
K_panel_2 = E*w_2*t_p/L_2     #Single panel stiffness in N/m
K_panel_3 = E*w_1*t_p/w_2
K_panel_4 = E*L_1*t_p/w_2
K_panels_axial = 2*K_panel_1 + 2*K_panel_2    #4 Panel stiffness in N/m
K_panels_lateral = 2*K_panel_3 + 2*K_panel_4
K_Tanks = E*(np.pi*(r_outer_tanks**2-(r_outer_tanks-t_tanks)**2))/L_1*4
K_rod = E*(np.pi*(r_outer_rod**2-(r_outer_rod-t_rod)**2))/L_1*4
K_total_axial = K_panels_axial + K_Tanks  #Total stiffness in N/m
K_total_lateral = K_panels_lateral +2*K_rod

# === Damping ===
zeta = 0.01  # 1% damping ratio
omega_n_axial = np.sqrt(K_total_axial / M_axial)  # rad/s
C_axial= 2 * zeta * np.sqrt(K_total_axial * M_axial)  # Ns/m
omega_n_lateral = np.sqrt(K_total_lateral/ M_lateral)  # rad/s
C_lateral= 2 * zeta * np.sqrt(K_total_lateral * M_lateral)  # Ns/m

# === Forcing: 1g sinusoidal acceleration at 100 Hz ===
f_drive_axial = 925#100  # Hz
f_drive_lateral = 925#100  # Hz
omega_drive_axial = 2 * np.pi * f_drive_axial  # rad/s
omega_drive_lateral = 2 * np.pi * f_drive_lateral  # rad/s
A_force_axial = M_axial * 9.81*5.13  # N
A_force_lateral = M_lateral * 9.81*5.13  # N

# === Time Domain Simulation ===
t_span = (0, 0.05)  # simulate for 50 ms
t_eval = np.linspace(t_span[0], t_span[1], 10000)  # high-res output

def systemaxial(t, y):
    x, v = y  # displacement and velocity
    a_t = A_force_axial * np.sin(omega_drive_axial * t)
    dxdt = v
    dvdt = (a_t - C_axial * v - K_total_axial * x) / M_axial
    return [dxdt, dvdt]

def systemlateral(t, y):
    x, v = y  # displacement and velocity
    a_t = A_force_lateral * np.sin(omega_drive_lateral * t)
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


max_disp_axial = np.max(np.abs(sol_axial.y[0]))
max_disp_lateral = np.max(np.abs(sol_lateral.y[0]))

print(f"Maximum Axial Displacement: {max_disp_axial:.6e} m")
print(f"Maximum Lateral Displacement: {max_disp_lateral:.6e} m")

# === Compute Strain and Stress ===
# For axial: L = L_1 (height of panel)
# For lateral: L = w_2 (width of panel)
strain_axial = max_disp_axial / L_1
strain_lateral = max_disp_lateral / w_1

stress_axial = E * strain_axial
stress_lateral = E * strain_lateral

print(f"Axial Panel Stress: {stress_axial:.2f} Pa")
print(f"Lateral Panel Stress: {stress_lateral:.2f} Pa")

#------ ACOUSTIC------
import pandas as pd

# Reference pressure
p_ref = 2e-5  # Pa
p_rms = p_ref * 10**(137.2/20)  # Convert dB to Pa
p_peak = p_rms * np.sqrt(2)  # Convert RMS to peak pressure
print(p_rms,p_peak)  # Convert dB to Pa

#------------------adding stiffeners to the panels - hat stiffeners--------------
alpha = 0.8
n = 0.6
E = 71.7*10**9  # Elastic module in Pa
sigma_yield = 503 * 10**6  # Yield strength in Pa
sigma_buck_sheet = 0.83 *10**6
C15 = 0.425      # buckling coefficient for stiffeners 
C234 = 4         # buckling coefficient for stiffeners
v = 0.334         # Poisson's ratio
#omega stiffener with whwhw (1,2,3,4,5)
h = 0.035
w = 0.02
t = 0.002
A24 = h*t
A135 = w*t #Area of the stiffener in m^2
A_tot = A24*2+A135*3
print(A_tot)
# Critical stress for different parts of the stiffeners
sigma_crippling15 = alpha*(C15/sigma_yield*E*np.pi**2/(12*(1-v**2))*(t/w)**2)**(1-n)*sigma_yield
sigma_crippling24 = alpha*(C234/sigma_yield*E*np.pi**2/(12*(1-v**2))*(t/h)**2)**(1-n)*sigma_yield
sigma_crippling3 = alpha*(C234/sigma_yield*E*np.pi**2/(12*(1-v**2))*(t/w)**2)**(1-n)*sigma_yield
# Total crippling stress for the stiffeners
sigma_crippling_tot = (2*sigma_crippling15*A135 + 2*sigma_crippling24*A24 + sigma_crippling3 *A135 )/ A_tot
print(f"Total crippling stress: {sigma_crippling_tot:.2f} Pa")
b = w_2 /2 # stiffener spacing in m 
print(b)
scr  = 4*E*np.pi**2/(12*(1-v**2))*(t_p/b)**2    #panel 

s_tot = (scr*b*t_p+sigma_crippling_tot*A_tot)/(A_tot + b*t_p)
print(f"Total stress with stiffeners: {s_tot:.2f} Pa")


# # Your 1/3 octave band center frequencies (Hz)
# frequencies = np.array([
#     31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500,
#     630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000
# ])

# # Example: SPL values (from the table, WEST RANGE, WITH blankets)
# spl_db = np.array([
#     119.8, 120.0, 120.0, 120.0, 119.8, 120.5, 121.5, 122.0, 121.5, 120.5, 119.0,
#     117.0, 115.0, 113.0, 111.0, 109.5, 108.0, 107.0, 106.0, 105.0, 104.0,
#     103.0, 102.0, 101.0, 100.0, 99.0
# ])

# # Convert SPL (dB) to RMS pressure (Pa)
# p_rms = p_ref * 10**(spl_db / 20)
# plt.plot(frequencies,p_rms)
# plt.show()
# # Estimate bandwidth for 1/3 octave
# factor = 2**(1/6)-2**(-1/6)
# factor_to_db = -10*np.log10(factor)
# delta_f = frequencies / factor_to_db

# # Compute pressure PSD in Pa^2/Hz
# psd = p_rms**2 / delta_f

# # Use trapezoidal integration over the frequency array
# integrated_power = np.trapz(psd, x=frequencies)  # Pa^2

# # Calculate total RMS stress from integrated power
# total_rms_pressure = np.sqrt(integrated_power)

# print(f"Total RMS Pressure for acoustic vibrations: {total_rms_pressure:.2f} Pa")
# #--------------Random vibration analysis-------------------
# frequency_rnd = np.array([20, 100, 300, 700, 800, 925, 2000])  # Hz
# psd_rnd = np.array([0.0044, 0.0044, 0.01, 0.01, 0.03, 0.03, 0.00644])  # g^2/Hz
# # Use trapezoidal integration over the frequency array
# integrated_power_rnd = np.trapz(psd_rnd, x=frequency_rnd)  # Pa^2

# # Calculate total RMS stress from integrated power
# total_rms_pressure_rnd = np.sqrt(integrated_power_rnd)

# print(f"Total RMS Pressure for random vibration: {total_rms_pressure_rnd:.2f} Pa")



