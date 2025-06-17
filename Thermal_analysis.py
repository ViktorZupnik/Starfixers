import numpy as np
from jinja2.nodes import Break
from matplotlib import pyplot as plt

Js = 1361       #w/m*m solar flux at earth
sigma = 5.67*10**(-8)       # Boltzmann constant W/m**2K**4
a = 0.14                #absoption coef
epsilon_IR = 0.09      #IR absorption coef
epsilon = 0.09           #emission coef
As = 0.823                  #sun exposed area
Ae = 0.823                 #Earth exposed area
A_emitted = 5.292      #emitting area
A_r = 0                #no radiators 
A_tot = A_emitted + A_r  # Total area for radiation balance
epsilon_r = 0.92
Jir = 237              #IR heat flux
T = 300                 #s/c temperature
albedo =0.3            #albedo

Re = 6371
R01 = 6971
R02 = 6752
f1 = (1-np.sqrt(1-(Re/R01)**2))/2 #Earth fill factor at 600km
f2 = (1-np.sqrt(1-(Re/R02)**2))/2 #Earth fill factor at 381km

Q_sun = a * Js * As              #Solar radiation W
Q_albedo = a * Js * albedo * Ae *f2             #albedo radiation W
Q_IR = epsilon_IR * Jir * Ae *f2 #Earth Ir radiation W
Qemitted = sigma* T**4  *(epsilon* A_emitted + epsilon_r* A_r)  #Emitted radiation W
Q_internal = 127    # Internal heat W (example) +RANGE OF NUMBERS


def Eclipse_temperature(Qearth, P_dis, epsilon, s, A_emitted):
    T = ((Qearth + P_dis)/(epsilon*s*A_emitted))**0.25
    return T

print (f"min eclipse temp is:{Eclipse_temperature(Q_IR, Q_internal, epsilon, sigma, A_emitted)}")
def Sunside_temperature(Qsun, Qalb, Qearth, P_dis, epsilon, s, A_emitted):
    T = ((Qsun + Qalb + Qearth + P_dis)/(epsilon*s*A_emitted))**0.25
    return T
print(f"max sunside temp is:{Sunside_temperature(Q_sun, Q_albedo, Q_IR,Q_internal, epsilon, sigma, A_emitted)}")

sigma = 5.67e-8  # W/m²K⁴
mass = 1000       # kg (aluminum)
cp = 900          # J/kg·K (aluminum spec. heat)
ext_area = A_tot   # m² (cube side ~1.15 m for 1000 kg Al)
U_MLI = 0.5          # W/m²K (effective heat transfer coeff)



T_inner11 = Sunside_temperature(Q_sun, Q_albedo, Q_IR, Q_internal, epsilon, sigma, A_emitted)- Q_internal / (U_MLI * A_emitted)
print(f"Inner surface temperature in sunlight: {T_inner11:.2f} K")
T_inner22 = Eclipse_temperature(Q_IR, Q_internal, epsilon, sigma, A_emitted) + Q_internal / (U_MLI * A_emitted)
print(f"Inner surface temperature in eclipse: {T_inner22:.2f} K")

# Environment (600 km LEO)
T_eclipse = Eclipse_temperature(Q_IR, Q_internal, epsilon, sigma, A_emitted)       # K
T_sunlight = Sunside_temperature(Q_sun, Q_albedo, Q_IR, Q_internal, epsilon, sigma, A_emitted)        # K
solar_flux = 1361        # W/m²
orbit_period = 96 * 60   # seconds (96 min)
eclipse_frac = 0.375    # 35% of orbit in eclipse


# # Time-stepping simulation (1 orbit)
# time = np.linspace(0, orbit_period, 1000)
# T_ext = np.where(time < eclipse_frac * orbit_period, T_eclipse, T_sunlight)
# T_internal = np.ones_like(time) * 290  # Initial guess

# for i in range(1, len(time)):
#     dt = time[i] - time[i-1]
#     Q_MLI = U_MLI * ext_area * (T_ext[i] - T_internal[i-1])
#     dT = (Q_MLI + Q_internal) * dt / (mass * cp)
#     T_internal[i] = T_internal[i-1] + dT

# # Results
# min_T = np.min(T_internal)
# max_T = np.max(T_internal)
# #print(f"Internal temp range: {min_T:.1f} K to {max_T:.1f} K")

a1 = [0.65, 0.16, 0.2,0.35,0.27]
e1 = [0.82, 0.03, 0.15,0.79,0.82]
alpha_epsilon_pairs = list(zip(a1, e1))

for alpha, epsilon in alpha_epsilon_pairs:
    Q_sun = alpha * Js * As
    Q_albedo = alpha * Js * albedo * Ae * f2
    Q_IR = epsilon_IR * Jir * Ae * f2

    T_eclipse = Eclipse_temperature(Q_IR, Q_internal, epsilon, sigma, A_emitted)
    T_sun = Sunside_temperature(Q_sun, Q_albedo, Q_IR, Q_internal, epsilon, sigma, A_emitted)

    print(f"{alpha:5.2f} {epsilon:5.2f} {T_eclipse:15.2f} {T_sun:15.2f}")












