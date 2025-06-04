import numpy as np
import matplotlib.pyplot as plt

g_axial = 8.5
g_lateral = 3
E = 114 * 10**9 #Elastic module Pa
rho = 4429 #Density kg/m**3

M_fuel_tank = 200       #Propellant tank mass (fueld) in kg
M = 100                 #Mass supported by the side panles at launch
E = 200 *10**9          #Panel elastic module in Pa
t = 0.005               #Panel thickness in m
w_1 = 1                 #Panel 1 width in m
L_1 = 1                 #Panel 1 length in m (height)
w_2 = 1                 #Panel 2 width in m
L_2 = L_1                #Panel 2 length in m (height)
K_panel_1 = E*w_1*t/L_1     #Single panel stiffness in N/m
K_panel_2 = E*w_2*t/L_2     #Single panel stiffness in N/m
K_panels = 2*K_panel_1 + 2*K_panel_2    #4 Panel stiffness in N/m
#print(K_panels)
omega_bus = (1/(2*np.pi))*np.sqrt(K_panels/M)   #Bus natual frequency in Hz using a slightly bs method
#print(omega_bus)