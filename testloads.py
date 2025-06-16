import numpy as np
import pytest
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
    b = w_2-2*r_outer_tanks  # Effective width of the panel
    # Buckling of panel width b
    sigma_newsheet = 4 * E * np.pi**2 / (12 * (1 - v**2)) * (t_p / b)**2 #reduced width of the panel
    # Total effective stress of panel + stiffener
    sigma_with_stiff = (sigma_newsheet * b * t_p + sigma_stiffener * A_stiff) / (A_stiff + b * t_p)

    return sigma_newsheet, sigma_with_stiff

print("pipe stringer crippling stress: ", halfpipe_stringer(r_outer_tanks, t_tanks, sigma_yield, E)[0])

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
K_Tanks = E*(np.pi*(r_outer_tanks**2-(r_outer_tanks-t_tanks)**2))/L_1*4
K_rod = E*(np.pi*(r_outer_rod**2-(r_outer_rod-t_rod)**2))/ (L_1)*4
K_stiff = E*A_stiff/w_1
K_total_axial = K_panels_axial + K_Tanks  #Total stiffness in N/m
K_total_lateral = K_panels_lateral +2*K_rod +6*K_stiff  #Total stiffness in N/m

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
print("The axial static stress is: ", Axial_static_stress , " Pa")
print("The lateral static stress is: ", Lateral_static_stress , " Pa")
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
plt.title('Axial Displacement Response to 5.13g Acceleration at 925 Hz')
plt.grid(True)
plt.tight_layout()
plt.show()

# === Plot Results Lateral ===
plt.figure(figsize=(10, 5))
plt.plot(sol_lateral.t*1000, sol_lateral.y[0], label='Displacement (m)', color='blue')
plt.xlabel('Time (ms)')
plt.ylabel('Lateral Displacement (m)')
plt.title('Lateral Displacement Response to 5.13g Acceleration at 925 Hz')
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




# === Dummy implementations — replace with real ones ===
def bending_stress_at_x(x, g_lateral=9.81, total_mass=100, L=10):
    """Stress from distributed load due to lateral gravity"""
    # Simple parabolic approximation (e.g. beam with pinned ends)
    return (total_mass * g_lateral * x * (L - x) / L**2,)

def compute_moment_of_inertia(r, t):
    """Moment of inertia for thin-walled cylinder approximation"""
    return (np.pi / 4) * (r**4 - (r - t)**4)

def axial_natural_frequency(K, M):
    """Natural frequency: ω = sqrt(K/M)"""
    return np.sqrt(K / M)

def compute_stiffness(area, E, L):
    """Axial stiffness: K = EA/L"""
    return E * area / L

def sigma_with_stiff(sigma_panel, sigma_stiff, b, a):
    """Composite stress from weighted average over areas"""
    return (sigma_panel * b + sigma_stiff * a) / (b + a)

# === Tests ===

def test_bending_stress_zero_load():
    assert bending_stress_at_x(3, g_lateral=0, total_mass=0)[0] == 0

def test_bending_stress_symmetry():
    L = 10
    delta = 1e-4
    left = bending_stress_at_x(L / 2 - delta, L=L)[0]
    right = bending_stress_at_x(L / 2 + delta, L=L)[0]
    assert np.isclose(left, right, rtol=1e-5)

def test_bending_stress_maximum_at_center():
    L = 10
    center_stress = bending_stress_at_x(L / 2, L=L)[0]
    edge_stress = bending_stress_at_x(0, L=L)[0]
    assert center_stress > edge_stress

def test_moment_of_inertia_limit_thick_equals_solid():
    r = 0.5
    t = 0.5
    I_hollow = compute_moment_of_inertia(r, t)
    I_solid = (np.pi / 4) * r**4
    assert np.isclose(I_hollow, I_solid, rtol=1e-5)

def test_moment_of_inertia_very_thin_wall():
    r = 0.5
    I_thin = compute_moment_of_inertia(r, 1e-6)
    I_thicker = compute_moment_of_inertia(r, 0.01)
    assert I_thin < I_thicker * 1e-3

def test_axial_natural_frequency_scaling():
    K = 1000
    M = 100
    omega_1 = axial_natural_frequency(K, M)
    omega_2 = axial_natural_frequency(K, M * 4)
    assert np.isclose(omega_2, omega_1 / 2, rtol=1e-3)

def test_stiffness_vs_length():
    A = 0.01
    E = 70e9
    L1 = 1
    L2 = 2
    K1 = compute_stiffness(A, E, L1)
    K2 = compute_stiffness(A, E, L2)
    assert np.isclose(K2, K1 / 2, rtol=1e-5)

def test_sigma_with_stiff_bounds():
    s1 = 100
    s2 = 300
    a, b = 0.03, 0.01
    combined = sigma_with_stiff(s1, s2, b, a)
    assert s1 < combined < s2

def test_sigma_with_stiff_area_limits():
    s1 = 150
    s2 = 350
    assert np.isclose(sigma_with_stiff(s1, s2, 0.01, 0), s1)
    assert np.isclose(sigma_with_stiff(s1, s2, 0, 0.01), s2)

def test_sigma_with_stiff_monotonicity():
    s1, s2, b = 100, 300, 0.01
    values = [sigma_with_stiff(s1, s2, b, a) for a in np.linspace(0.001, 0.1, 10)]
    assert all(v2 > v1 for v1, v2 in zip(values, values[1:]))

def test_unit_consistency_stress():
    stress = bending_stress_at_x(5, g_lateral=9.81, total_mass=10)[0]
    assert 0 < stress < 1e6  # in Pascals, rough upper bound

# To run in interactive mode (e.g., Jupyter)
def run_all_tests():
    pytest.main(["-v", "--disable-warnings"])

# Uncomment this line to run in script
# if __name__ == "__main__":
#     run_all_tests()


# test_forced_vibration.py

import numpy as np
from scipy.integrate import solve_ivp

# === Constants and Parameters ===
g = 9.81
zeta = 0.01

# Example physical parameters — replace with your actual ones
K_total_axial = 1e5
K_total_lateral = 1e4
M_axial = 2.0
M_lateral = 2.0

f_drive_axial = 925  # Hz
f_drive_lateral = 925  # Hz

omega_drive_axial = 2 * np.pi * f_drive_axial
omega_drive_lateral = 2 * np.pi * f_drive_lateral

A_force_axial = M_axial * g * 5.13
A_force_lateral = M_lateral * g * 5.13

omega_n_axial = np.sqrt(K_total_axial / M_axial)
omega_n_lateral = np.sqrt(K_total_lateral / M_lateral)

C_axial = 2 * zeta * np.sqrt(K_total_axial * M_axial)
C_lateral = 2 * zeta * np.sqrt(K_total_lateral * M_lateral)

t_span = (0, 0.05)
t_eval = np.linspace(t_span[0], t_span[1], 10000)
y0 = [0, 0]

# === Systems to Test ===
def systemaxial(t, y):
    x, v = y
    a_t = A_force_axial * np.sin(omega_drive_axial * t)
    dxdt = v
    dvdt = (a_t - C_axial * v - K_total_axial * x) / M_axial
    return [dxdt, dvdt]

def systemlateral(t, y):
    x, v = y
    a_t = A_force_lateral * np.sin(omega_drive_lateral * t)
    dxdt = v
    dvdt = (a_t - C_lateral * v - K_total_lateral * x) / M_lateral
    return [dxdt, dvdt]

# === Helper Function ===
def base_excitation_amplitude(acc_peak, k, m, zeta, omega_drive):
    omega_n = np.sqrt(k / m)
    denom = np.sqrt((omega_n**2 - omega_drive**2)**2 + (2 * zeta * omega_n * omega_drive)**2)
    return acc_peak / denom

# === Tests ===

def test_axial_steady_state_matches_analytical():
    sol = solve_ivp(systemaxial, t_span, y0, t_eval=t_eval)
    x = sol.y[0]
    x_peak = np.max(np.abs(x[-1000:]))  # final section = steady state
    X_theoretical = base_excitation_amplitude(
    acc_peak=5.13 * g,
    k=K_total_axial,
    m=M_axial,
    zeta=zeta,
    omega_drive=omega_drive_axial)
    assert np.isclose(x_peak, X_theoretical, rtol=0.05), f"{x_peak} vs {X_theoretical}"

def test_lateral_steady_state_matches_analytical():
    sol = solve_ivp(systemlateral, t_span, y0, t_eval=t_eval)
    x = sol.y[0]
    x_peak = np.max(np.abs(x[-1000:]))
    X_theoretical = analytical_amplitude(
        A_force_lateral, M_lateral, K_total_lateral, C_lateral, omega_drive_lateral, zeta
    )
    assert np.isclose(x_peak, X_theoretical, rtol=0.05), f"{x_peak} vs {X_theoretical}"

def test_axial_response_peak_near_resonance():
    # Test response amplitude near ω_n vs far from ω_n
    omega_drive_near = omega_n_axial
    omega_drive_far = 2 * omega_n_axial

    def system_near(t, y):
        a_t = A_force_axial * np.sin(omega_drive_near * t)
        x, v = y
        dxdt = v
        dvdt = (a_t - C_axial * v - K_total_axial * x) / M_axial
        return [dxdt, dvdt]

    def system_far(t, y):
        a_t = A_force_axial * np.sin(omega_drive_far * t)
        x, v = y
        dxdt = v
        dvdt = (a_t - C_axial * v - K_total_axial * x) / M_axial
        return [dxdt, dvdt]

    sol_near = solve_ivp(system_near, t_span, y0, t_eval=t_eval)
    sol_far = solve_ivp(system_far, t_span, y0, t_eval=t_eval)

    amp_near = np.max(np.abs(sol_near.y[0][-1000:]))
    amp_far = np.max(np.abs(sol_far.y[0][-1000:]))

    assert amp_near > 3 * amp_far, f"Resonant amp: {amp_near} vs off-resonant: {amp_far}"
