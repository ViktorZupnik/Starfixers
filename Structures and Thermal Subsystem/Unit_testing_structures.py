"Run the code to see the tests that were performed on the Loads.py code"
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

def omega_stringer(h, w, t,sidestiff, topbottomstiff,C15=0.425, C234=4):
    if t > 0.0:
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
    
        b_side= L_1/(sidestiff+1)  #2 stiffeners on the side
        b_topbottom = w_2/(topbottomstiff+1) #1 stiffener on the top and bottom
        # Buckling of panel width b
        sigma_newsheet_side = 4 * E * np.pi**2 / (12 * (1 - v**2)) * (t_p / b_side)**2
        sigma_newsheet_topbottom = 4 * E * np.pi**2 / (12 * (1 - v**2)) * (t_p / b_topbottom)**2
        # Total effective stress of panel + stiffener
        sigma_with_stiff_side = (sigma_newsheet_side * b_side * t_p + sigma_stiffener * A_stiff) / (A_stiff + b_side * t_p)
        sigma_with_stiff_topbottom = (sigma_newsheet_topbottom * b_topbottom * t_p + sigma_stiffener * A_stiff) / (A_stiff + b_topbottom * t_p)
    else: 
        b_side = L_1
        b_topbottom = w_2
        # Total effective stress of panel + stiffener
        sigma_newsheet_side = 1/1.1*4 * E * np.pi**2 / (12 * (1 - v**2)) * (t_p / b_side)**2
        sigma_with_stiff_side = sigma_newsheet_side
        sigma_newsheet_topbottom = 1/1.1*4 * E * np.pi**2 / (12 * (1 - v**2)) * (t_p / b_topbottom)**2
        sigma_with_stiff_topbottom = sigma_newsheet_topbottom
        A_stiff = 0.0
        sigma_stiffener = 0.0

    return A_stiff, sigma_stiffener, sigma_with_stiff_side, sigma_with_stiff_topbottom
def pipe_stringer(r_outer_tanks, t_tanks, sigma_yield, E): 
    A_stiff = np.pi * (r_outer_tanks**2 - (r_outer_tanks - t_tanks)**2)   
    sigma_stiffener = alpha * (2 / sigma_yield * E *0.605*t_tanks/r_outer_tanks)**(1 - n) * sigma_yield
    b = w_2-2*r_outer_tanks  # Effective width of the panel
    # Buckling of panel width b
    sigma_newsheet = 1/1.1 * 4 * E * np.pi**2 / (12 * (1 - v**2)) * (t_p / b)**2 #reduced width of the panel
    # Total effective stress of panel + stiffener
    sigma_with_stiff = (sigma_newsheet * b * t_p + sigma_stiffener * A_stiff) / (A_stiff + b * t_p)

    return sigma_newsheet, sigma_with_stiff

print("pipe stringer crippling stress: ", pipe_stringer(r_outer_tanks, t_tanks, sigma_yield, E)[0])

# === Tests ===

def test_bending_stress_zerodisplacement_atedge():
    total_mass =100
    beam_length =10
    r = 0.2
    t = 0.002
    x =0
    
    g_lateral = 9.81
    assert bending_stress_at_x(x,total_mass,beam_length,r,t,g_lateral)[1] == 0

#check no stiffness increase if stiffness thickness is zero
def test_omegastiff():
    t = 0.0
    w = 0.1
    h = 0.1 
    expected = 3.21   #buckling stress without stiffeners
    assert abs(omega_stringer(h, w, t,3,2)[2]*10**(-6) == expected) < 1e-2
# check that small tanks have no stiffening effect on sheet width
def test_smalltankssstiff():
    r_outer_tanks = 0.01
    t_tanks = 0.003
    E = 71.7*10**9
    sigma_yield = 503 * 10**6  # Yield strength in Pa
    sheet = pipe_stringer(r_outer_tanks, t_tanks, sigma_yield, E)
    expected = 2.17
    assert abs(sheet[0]*10**(-6)- expected) < 10  #buckling stress without stiffeners
def test_amplitude_g():
    K_total_axial = 1e6  # N/m
    K_total_lateral = 1e6  # N/m    
    M_axial = 100  # kg
    M_lateral = 100  # kg
    zeta = 0.01  # 1% damping ratiodef systemaxial(t, y):
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
    sol_axial1 = solve_ivp(systemaxial, t_span, y0, t_eval=t_eval, method='RK45')
    sol_lateral1 = solve_ivp(systemlateral, t_span, y0, t_eval=t_eval, method='RK45')
    max_disp_axial1 = np.max(np.abs(sol_axial1.y[0]))
    max_disp_lateral1 = np.max(np.abs(sol_lateral1.y[0]))
    A_force_axial= M_axial * g*10  # N
    A_force_lateral = M_lateral * g*10 # N
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
    sol_axial2 = solve_ivp(systemaxial, t_span, y0, t_eval=t_eval, method='RK45')
    sol_lateral2 = solve_ivp(systemlateral, t_span, y0, t_eval=t_eval, method='RK45')
    max_disp_axial2 = np.max(np.abs(sol_axial2.y[0]))
    max_disp_lateral2 = np.max(np.abs(sol_lateral2.y[0]))
    assert max_disp_axial2 > max_disp_axial1
    assert max_disp_lateral2 > max_disp_lateral1


   
    

if __name__ == "__main__":
    import sys
    import pytest
    sys.exit(pytest.main([__file__]))