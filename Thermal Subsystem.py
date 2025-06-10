import numpy as np
from matplotlib import pyplot as plt
Js = 1361       #w/m*m solar flux at earth
sigma = 5.67*10**(-8)       # Boltzmann constant W/m**2K**4
a = 0.1                  #absoption coef
epsilon_IR = 0.05         #IR absorption coef
epsilon = 0.05            #emission coef
As = 1                  #sun exposed area
Ae = 1                  #Earth exposed area
A_emitted = 4           #emitting area
Jir = 237              #IR heat flux
T = 300                 #s/c temperature
albedo =0.3             #albedo


Re = 6371
R01 = 6971
R02 = 6752
f1 = (1-np.sqrt(1-(Re/R01)**2))/2 #Earth fill factor at 600km
f2 = (1-np.sqrt(1-(Re/R02)**2))/2 #Earth fill factor at 381km

Q_sun = a * Js * As              #Solar radiation W
Q_albedo = a * Js * albedo * Ae *f2             #albedo radiation W
Q_IR = epsilon_IR * Jir * Ae *f2 #Earth Ir radiation W
Qemitted = epsilon * sigma * T**4 * A_emitted   #Emitted radiation W
Q_internal = 0                       #DIssipated energy W


def Eclipse_temperature(Qearth, P_dis, epsilon, s, A_emitted):
    T = ((Qearth + P_dis)/(epsilon*s*A_emitted))**0.25
    return T

#print ( Eclipse_temperature(Qearth, P_dis, epsilon, s, A_emitted))
def Sunside_temperature(Qsun, Qalb, Qearth, P_dis, epsilon, s, A_emitted):
    T = ((Qsun + Qalb + Qearth + P_dis)/(epsilon*s*A_emitted))**0.25
    return T
#print(Sunside_temperature(Qsun, Qalb, Qearth, P_dis, epsilon, s, A_emitted))


# Spacecraft properties
mass = 500  # kg
c = 900    # J/kg·K (e.g., aluminum)

# Orbit properties
period = 5760  # sec (90 min)
eclipse_fraction = 0.375
dt = 1  # time step in seconds
t_end = 3 * period  # simulate 3 orbits

# Time array
time = np.arange(0, t_end, dt)
T = np.zeros_like(time)
T[0] = 300  # initial temp in K

# Heat input function
def Q_in(t):
    phase = t % period
    in_eclipse = phase < (eclipse_fraction * period)
    if in_eclipse:
        return Q_IR + Q_internal
    else:
        return Q_sun + Q_albedo + Q_IR + Q_internal

# Time integration using Euler method
for i in range(1, len(time)):
    Qin = Q_in(time[i])
    Qout = epsilon * sigma * A_emitted * T[i-1]**4
    dTdt = (Qin - Qout) / (mass * c)
    T[i] = T[i-1] + dTdt * dt

# Plot
plt.figure(figsize=(10, 5))
plt.plot(time / 60, T, label='Temperature (K)')
plt.xlabel('Time (minutes)')
plt.ylabel('Temperature (K)')
plt.title('Transient Thermal Response Over Orbit')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()