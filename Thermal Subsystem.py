import numpy as np
Js = 1361       #w/m*m solar flux at earth
s = 5.67*10**(-8)       # Boltzmann constant W/m**2K**4
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

Qsun = a * Js * As              #Solar radiation W
Qalb = a * Js * albedo * Ae *f2             #albedo radiation W
Qearth = epsilon_IR * Jir * Ae *f2 #Earth Ir radiation W
Qemitted = epsilon * s * T**4 * A_emitted   #Emitted radiation W
P_dis = 0                       #DIssipated energy W


def Eclipse_temperature(Qearth, P_dis, epsilon, s, A_emitted):
    T = ((Qearth + P_dis)/(epsilon*s*A_emitted))**0.25
    return T

print ( Eclipse_temperature(Qearth, P_dis, epsilon, s, A_emitted))
def Sunside_temperature(Qsun, Qalb, Qearth, P_dis, epsilon, s, A_emitted):
    T = ((Qsun + Qalb + Qearth + P_dis)/(epsilon*s*A_emitted))**0.25
    return T
print(Sunside_temperature(Qsun, Qalb, Qearth, P_dis, epsilon, s, A_emitted))