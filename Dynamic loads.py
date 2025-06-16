import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt





g = 9.80665
E = 71.7*10**9  # Elastic module in Pa
def bending_stress_at_x( #https://www.engineeringtoolbox.com/beams-fixed-both-ends-support-loads-deflection-d_809.html
    x,            # Position along tank [m]
    total_mass,   # Tank total mass [kg]
    beam_length,  # Tank height [m]
    r,            # Tank radius [m]
    t,            # Tank thickness [m]
    g_lateral):    # Lateral acceleration  [g] (1g = 9.80665 m/s^2)
#     A_rod,          # Support rod area m**2
#     E_rod            # Support rod elastic modulus

# ):
    w = (g_lateral * total_mass * g) / beam_length  # distributed load [N/m]

    # Bending moment at position x
    Mx = (w * beam_length / 12) * (6 * x - beam_length) - (w * x**2) / 2

    # Section properties
    I = np.pi*(r**4-(r-t)**4)/4

    # Bending stress at x
    stress = (Mx * r) / I  # in Pascals
    displacement = 1/(E*I)*((-w *beam_length**2 * x**2) / 24 + (w * beam_length * x**3) / 12- (w * x**4) / 24)


    return stress, displacement


print ("bending stress in the middle of a tank: ", bending_stress_at_x(1.065/2, 150, 1.065, 0.533/2, 0.002, 3))

g_axial = 8.5
g_lateral = 3
M_fuel_tank = 95.474+15.711       #Propellant tank mass (fueld) in kg
alpha = 0.8
v = 0.334
n = 0.6
E = 71.7*10**9  # Elastic module in Pa
sigma_yield = 503 * 10**6  # Yield strength in Pa
M = 100                 #Mass supported by the side panles at launch
t_p = 0.003              #Panel thickness in m
w_1 = 1.00               #Panel 1 width in m
L_1 = 0.823                #Panel 1 length in m (height)
w_2 = 1.00                 #Panel 2 width in m
L_2 = L_1                 #Panel 2 length in m (height)
r_outer_rod = 0.02          #Outer radius of the support rod in m
t_rod = 0.002
print(L_1*w_1*4+w_1*w_2*2)
rho_panels = 2810
r_outer_tanks = 0.412/2
t_tanks = 0.003

#stiffener dimensions
h_stiff = 0.01
w_stiff = 0.01
t_stiff = 0.001

def omega_stringer(h, w, t,C15=0.425, C234=4):
    # Areas of stiffener segments (based on "whwhw" configuration)
    A24 = h * t        # vertical parts
    A135 = w * t       # horizontal parts
    A_stiff = (2 * A24 + 3 * A135)

    # Crippling stresses
    sigma_crippling15 = alpha * (C15 / sigma_yield * E * np.pi**2 / (12 * (1 - v**2)) * (t / w)**2)**(1 - n) * sigma_yield
    sigma_crippling24 = alpha * (C234 / sigma_yield * E * np.pi**2 / (12 * (1 - v**2)) * (t / h)**2)**(1 - n) * sigma_yield
    sigma_crippling3  = alpha * (C234 / sigma_yield * E * np.pi**2 / (12 * (1 - v**2)) * (t / w)**2)**(1 - n) * sigma_yield

    # Combined crippling stress
    sigma_stiffener= (2 * sigma_crippling15 * A135 + 2 * sigma_crippling24 * A24 + sigma_crippling3 * A135) / A_stiff
    # Stiffener spacing (assuming 1 stiffener in the center, edges supported)
    b_side= L_1/3
    b_topbottom = w_2/2
    # Buckling of panel width b
    sigma_newsheet_side = 4 * E * np.pi**2 / (12 * (1 - v**2)) * (t_p / b_side)**2
    sigma_newsheet_topbottom = 4 * E * np.pi**2 / (12 * (1 - v**2)) * (t_p / b_topbottom)**2
    # Total effective stress of panel + stiffener
    sigma_with_stiff_side = (sigma_newsheet_side * b_side * t_p + sigma_stiffener * A_stiff) / (A_stiff + b_side * t_p)
    sigma_with_stiff_topbottom = (sigma_newsheet_topbottom * b_topbottom * t_p + sigma_stiffener * A_stiff) / (A_stiff + b_topbottom * t_p)

    return A_stiff, sigma_stiffener, sigma_with_stiff_side, sigma_with_stiff_topbottom
print("omega stiffener new sheet sides (2 stiffeners)", omega_stringer(h_stiff, w_stiff, t_stiff)[2]/1e6, " MPa"
      "omega stiffener new sheet top/bottom (1stiffener both directions)", omega_stringer(h_stiff, w_stiff, t_stiff)[3]/1e6, " MPa")
print(omega_stringer(h_stiff, w_stiff, t_stiff)[0], " m^2", omega_stringer(h_stiff, w_stiff, t_stiff)[1])  # Area of the stiffeners
# def halfpipe_stringer(r_outer_tanks, r_outer_rod, t_tanks, sigma_yield, E,C=0.366 ):
#     A_stiff = np.pi * (r_outer_tanks**2 - (r_outer_tanks - t_tanks)**2)/2  
#     Fcy = sigma_yield*A_stiff
#     a = np.sqrt(r_outer_tanks**2-r_outer_rod**2)
#     tan_alpha = r_outer_rod/(r_outer_tanks-a)
#     h = tan_alpha*2*r_outer_tanks
#     b = np.sqrt((h-r_outer_rod)**2 + (r_outer_tanks+a)**2)
#     b_prime = (h-b)/2
#     sigma_cr = (C*np.sqrt(Fcy*E)/(b_prime/t_tanks)**0.75)*A_stiff
#     # Buckling of panel width b
#     sigma_newsheet = 4 * E * np.pi**2 / (12 * (1 - v**2)) * (t_p / b)**2
#     return A_stiff, sigma_cr
def halfpipe_stringer(r_outer_tanks, t_tanks, sigma_yield, E): 
    A_stiff = np.pi * (r_outer_tanks**2 - (r_outer_tanks - t_tanks)**2)   
    sigma_stiffener = alpha * (2 / sigma_yield * E *0.605*t_tanks/r_outer_tanks)**(1 - n) * sigma_yield
    b = r_outer_tanks  # Effective width of the panel
    # Buckling of panel width b
    sigma_newsheet = 4 * E * np.pi**2 / (12 * (1 - v**2)) * (t_p / b)**2 #reduced width of the panel
    # Total effective stress of panel + stiffener
    sigma_with_stiff = (sigma_newsheet * b * t_p + sigma_stiffener * A_stiff) / (A_stiff + b * t_p)

    return sigma_newsheet, sigma_with_stiff

print("small wall buckling stress: ", halfpipe_stringer(r_outer_tanks, t_tanks, sigma_yield, E)[0])

#Area calculations of panels, tanks, rod and stiffeners
r_inner = r_outer_tanks - t_tanks #tank inner radius m
A_axial = 2*w_1*t_p + 2*w_2*t_p + 4* (2*np.pi*r_outer_tanks - 2*np.pi*r_inner) #area that holds axial loads
A_support =  np.pi*(r_outer_rod**2-(r_outer_rod-t_rod)**2) #area of the support beam
A_lateral = 2*w_2*t_p + 2*L_1*t_p + 2*A_support  #area that holds lateral loads
A_stiff = omega_stringer(h_stiff, w_stiff, t_stiff)[0]  #area of the stiffeners

#defining spring stiffnesses
K_panel_1 = E*w_1*t_p/L_1     #Single panel stiffness in N/m
K_panel_2 = E*w_2*t_p/L_2     #Single panel stiffness in N/m
K_panel_3 = E*w_1*t_p/w_2
K_panel_4 = E*L_1*t_p/w_2
K_panels_axial = 2*K_panel_1 + 2*K_panel_2    #4 Panel stiffness in N/m
K_panels_lateral = 2*K_panel_3 + 2*K_panel_4
K_panels_critical = ((3*E*w_1*t_p**3/12)/L_1**3)+((3*E*w_2**3*t_p/12)/L_1**3) #Critical buckling stiffness in N/m
K_Tanks = E*(np.pi*(r_outer_tanks**2-(r_outer_tanks-t_tanks)**2))/L_1*4
K_Tanks_critical = (3*E*(np.pi/4*(r_outer_tanks**4-(r_outer_tanks-t_tanks)**4))/L_1**3)
K_rod = E*(np.pi*(r_outer_rod**2-(r_outer_rod-t_rod)**2))/ (L_1)*4
K_stiff = E*A_stiff/w_1
K_total_axial = K_panels_axial + K_Tanks  #Total stiffness in N/m
K_total_lateral = K_panels_lateral +2*K_rod +6*K_stiff  #Total stiffness in N/m
#K_total_lateral = 2* K_panels_critical + 4*K_Tanks_critical

M_t_full = 150 #full fuel tank mass kg
M_rod = A_support*(w_1-4*r_outer_tanks)*rho_panels
M_stiff = A_stiff * rho_panels *w_1  #Mass of the stiffeners kg
#print(M_rod)
M_axial = w_1*w_2*t_p*rho_panels + rho_panels*(2*w_1*L_1*t_p + 2*w_2*L_1*t_p) + 4*M_t_full+10*M_stiff +4*M_rod  #Mass carried in the axial direction kg
M_lateral = 2*M_t_full + L_1*w_2*t_p*rho_panels +  rho_panels*(2*w_1*w_2*t_p + 2*w_1*L_1*t_p) + 3*M_rod +10*M_stiff   #Mass carried in the lateral direction kg
M_total = 2*w_1*w_2*t_p*rho_panels + rho_panels*(2*w_1*L_1*t_p + 2*w_2*L_1*t_p) + 12*M_stiff + 4*M_rod  #Total mass of the structure kg
#===== Static Loads=======
print(M_total)
Axial_static_stress = g_axial*g*M_axial/A_axial
Lateral_static_stress = g_lateral*g*M_lateral/A_lateral
critical_stress_lat = (M_axial*3*g*L_1*w_2/2)/((w_2**3*t_p/12)/L_1**3)
print("The axial static stress is: ", Axial_static_stress , " Pa")
print("The lateral static stress is: ", Lateral_static_stress , " Pa")
print("The critical stress in lateral direction is: ", critical_stress_lat, " Pa")


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
A_force_axial = M_axial * g*5.13  # N
A_force_lateral = M_lateral * g*5.13  # N

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
    dvdt = (a_t - C_lateral* v - K_total_lateral * x) / M_axial
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
plt.title('Axial Displacement Response to 5.13g Acceleration at 925 Hz')
plt.grid(True)
plt.tight_layout()
plt.show()

# === Plot Results Lateral ===
plt.figure(figsize=(10, 5))
plt.plot(sol_lateral.t*1000, sol_lateral.y[0], label='Displacement (m)', color='blue')
plt.xlabel('Time (ms)')
plt.ylabel('Lateral Displacement (m)')
plt.title('Lateral Displacement Response to 5.13g Acceleration at 925 Hz - Bottom Clamped')
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
max_root_stress = (3*E*max_disp_lateral*w_2/2)/L_1**2
strain_axial = max_disp_axial / L_1
strain_lateral = max_disp_lateral / w_1

stress_axial = E * strain_axial
stress_lateral = E * strain_lateral

print(f"Axial Panel Stress: {stress_axial:.2f} Pa")
print(f"Lateral Panel Stress: {stress_lateral:.2f} Pa")
print(f"Max Root Stress: {max_root_stress:.2f} Pa")



walls_axial = {
    "Front": (w_1, L_1),
    "Back": (w_1, L_1),
    "Left": (w_2, L_1),
    "Right": (w_2, L_1)
}

walls_lateral = {
    "Top": (w_1, w_2),
    "Bottom": (w_1, w_2),
    "Front": (w_1, L_1),
    "Back": (w_1, L_1),
}

def sigma_cr_plate(E, nu, t_p, b, k, SF):

    return (1/SF)*(k * (np.pi**2) * E / (12 * (1 - nu**2))) * (t_p / b)**2

def cylinder_buckling(E, t_t, r_outer, Kd):

    sigma_cr_ideal = 0.605 * E * t_t / (r_outer - t_t/2)
    sigma_cr_real = Kd * sigma_cr_ideal
    return sigma_cr_ideal, sigma_cr_real

Kd = 0.5        #Knockdown factor
ideal, real = cylinder_buckling(E, t_tanks, r_outer_tanks, Kd)

print(f"Ideal buckling stress for cylinder: {ideal/1e6:.4f} MPa")
print(f"Realistic (with knockdown) for cylinder: {real/1e6:.4f} MPa")
k = 4 #bucklin fcator for totally clamped plate
nu= 0.334 # poisson's ratio
SF = 1.1

print("Critical Buckling Stress for Cube Walls in axial direction:")
for name, (a, b) in walls_axial.items():
    b_eff = min(a, b)  # Buckling depends on shorter dimension
    sigma_cr = sigma_cr_plate(E, nu, t_p, b_eff, k, SF)
    print(f"{name} wall: {sigma_cr/1e6:.2f} MPa")

print("\n\nCritical Buckling Stress for Cube Walls in lateral direction:")
for name, (a, b) in walls_lateral.items():
    b_eff = min(a, b)  # Buckling depends on shorter dimension
    sigma_cr = sigma_cr_plate(E, nu, t_p, b_eff, k, SF)
    print(f"{name} wall: {sigma_cr / 1e6:.2f} MPa")
print(f"Buckling stress of a wall with 1 stiffener:{omega_stringer(h_stiff, w_stiff, t_stiff)[2]/1e6:.2f} MPa")

#------ ACOUSTIC------
# Reference pressure
p_ref = 2e-5  # Pa
p_rms = p_ref * 10**(137.9/20)  # Convert dB to Pa
p_peak = p_rms * np.sqrt(2)  # Convert RMS to peak pressure
print(f"rms pressure and peak pressure for acoustic loads in Pa:{p_rms,p_peak}")  # Convert dB to Pa
















'''leave this code here pls'''
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



