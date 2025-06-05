import numpy as np
import matplotlib.pyplot as plt
Total_mass = 1200
g_axial = 8.5
g_lateral = 3
M_fuel_tank = 200       #Propellant tank mass (fueld) in kg
M = 100                 #Mass supported by the side panles at launch
E = 114 *10**9          #Panel elastic module in Pa
t = 0.005               #Panel thickness in m
w_1 = 1                 #Panel 1 width in m
L_1 = 1                 #Panel 1 length in m (height)
w_2 = 1                 #Panel 2 width in m
L_2 = 1                 #Panel 2 length in m (height)





rho_panels = 4429
r_outer_tanks = 0.3
t_tanks = 0.01

r_inner = r_outer_tanks - t_tanks #tank inner radius m
A_axial = 2*w_1*t + 2*w_2*t + 4* (2*np.pi*r_outer_tanks - 2*np.pi*r_inner) #area that holds axial loads
#print(A_axial)
A_support = 0 #area of the support beam
A_lateral = 2*w_2*t + 2*L_1*t + 2*A_support  #area that holds lateral loads

M_t_full = 150 #ull fuel tank mass kg
M_axial = w_1*w_2*t*rho_panels + rho_panels*(2*w_1*L_1*t + 2*w_2*L_1*t) + 4*M_t_full   #Mass carried in the axial direction kg
M_lateral = 2*M_t_full + L_1*w_2*t*rho_panels +  rho_panels*(2*w_1*w_2*t + 2*w_1*L_1*t)     #Mass carried in the lateral direction kg


K_panel_1 = E*w_1*t/L_1     #Single panel stiffness in N/m
K_panel_2 = E*w_2*t/L_2     #Single panel stiffness in N/m
K_panels = 2*K_panel_1 + 2*K_panel_2    #4 Panel stiffness in N/m
K_Tanks = E*(np.pi*(r_outer_tanks**2-(r_outer_tanks-t_tanks)**2))/L_1*4
K_total = K_panels + K_Tanks  #Total stiffness in N/m
omega_total = (1/(2*np.pi))*np.sqrt(K_total/Total_mass)  #Total natural frequency in Hz
print (f'Total natural frequency: {omega_total:.2f} Hz')

# panels_submasses = np.linspace(0, M_panels, 100)
# tanks_submasses = np.linspace(0, M_fuel_tank, 100)
# for panel_submass in panels_submasses:
#     for tank_submass in tanks_submasses:
        