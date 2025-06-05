import numpy as np
import matplotlib.pyplot as plt
g = 9.80665
g_axial = 8.5 * g
g_lateral = 3 * g
SF = 1.1    #Safety factor

# ---------Ti6Al4V-----------
# E = 114 * 10**9 #Elastic module Pa
# rho = 4429 #Density kg/m**3
# s_yeild_t = 880 * 10**6 # yeals strength tensile Pa
# s_ult_t = 950 * 10**6 # ult strength tensile Pa
# s_yeild_c = 970 * 10**6 # yeals strength compressive Pa

# ----------Al7075-T6-----------
E = 71.7 * 10**9 #Elastic module Pa
rho = 2810 #Density kg/m**3
s_yeild_t = 503 * 10**6 # yeals strength tensile Pa
s_ult_t = 572 * 10**6 # ult strength tensile Pa
#s_yeild_c = 970 * 10**6 # yeals strength compressive Pa

r_outer_rod = 0.04
t_rod = 0.004

w_1 = 1     # Width of 2 sides m
w_2 = 1     # Width of 2 sides m
w_3 = 1     # height of a side plate m
t = np.arange(0.0001, 0.005, 0.0001)  #bus wall thickness m
r = 0.533/2  #tank outer radius m
t_t = 0.002 #tank wall thickness m
y1=[]
y2=[]
for i in t:
    r_inner = r - t_t #tank inner radius m
    A_axial = 2*w_1*i + 2*w_2*i + 4* (2*np.pi*r - 2*np.pi*r_inner) #area that holds axial loads
    #print(A_axial)
    A_support = np.pi*(r_outer_rod**2-(r_outer_rod-t_rod)**2) #area of the support beam
    A_lateral = 2*w_2*i + 2*w_3*i + 2*A_support  #area that holds lateral loads

    M_t_full = 150 #ull fuel tank mass kg
    M_axial = w_1*w_2*i*rho + rho*(2*w_1*w_3*i + 2*w_2*w_3*i) + 4*M_t_full   #Mass carried in the axial direction kg
    M_lateral = 2*M_t_full + w_3*w_2*i*rho      #Mass carried in the lateral direction kg


    stress_axial = g_axial*M_axial/A_axial
    stress_lateral = g_lateral*M_lateral/A_lateral

    y1.append(stress_axial)
    y2.append(stress_lateral)

plt.plot(t*1000, y1, 'r', label="Axial Stress" )
plt.plot(t*1000, y2, 'g', label="Lateral Stress")
plt.ylabel("stress [Pa]")
plt.xlabel("Bus wall thickness [mm]")
plt.legend()
plt.grid(True)
plt.show()




t0 = 0.002
r_inner = r - t_t #tank inner radius m
A_axial = 2*w_1*t0 + 2*w_2*t0 + 4* (2*np.pi*r - 2*np.pi*r_inner) #area that holds axial loads
#print(A_axial)
A_support = np.pi*(r_outer_rod**2-(r_outer_rod-t_rod)**2) #area of the support beam
A_lateral = 2*w_2*t0 + 2*w_3*t0 + 2*A_support  #area that holds lateral loads

M_t_full = 150 #ull fuel tank mass kg
M_axial = w_1*w_2*t0*rho + rho*(2*w_1*w_3*t0 + 2*w_2*w_3*t0)    #Mass carried in the axial direction kg
M_lateral = 2*M_t_full + w_3*w_2*t0*rho      #Mass carried in the lateral direction kg


stress_axial = g_axial*M_axial/A_axial
stress_lateral = g_lateral*M_lateral/A_lateral 
print("axial stress: ", stress_axial)
print("lateral stress: ", stress_lateral)
