import numpy as np

#Spacecraft constants
I_x = 201.9317163 # [kg*m2]
I_y = 131.4344193 # [kg*m2]
M = 646 # [kg]
l = 1 # [m]

#Thruster constants
T_cg = np.array([0.12, 0.9, 3.6, 20]) # [N]
Isp_cg = np.array([65, 70, 57, 242])
Temp_gas = 293 # [K]
Pres_gas = np.array([6.9, 90, 15.7, 100])*10**5 # [Pa]

#Debris constants
v_rel = 4.55 # [m/s]
critical_freedom = 10 # [m]
approach_distance = 1500 # [m]

#Allround constants
g = 9.80665 # [m/s2]
R = 8.314 # [-]
m_molN2 = 0.028 # [kg/mol]


#ADCS constants
comp_time = 120 # [s]
num_rendezvous = 1
num_slew = 1
saturation_factor = 1   # 0-1
n_desat = 10
H_sat = 12 # [Nms]


#We burst at right moment to slow down to catch up with debris
#Then we have to rotate 180deg to point thruster toward debris again
#If thruster did not fire, we have to abort mission and get into another orbit using our cold gas thrusters



def exit_check(T_cg, v_rel, critical_freedom, M, approach_distance):
    a = 2 * T_cg / M  # [m/s^2]

    # Case 1: Full burn to clear the danger zone
    t_full = (2 * critical_freedom / a)**0.5
    safe_distance = v_rel * t_full

    # Case 2: Burn then coast
    T = approach_distance / v_rel
    # Solve: (a/2) * t_b^2 - a * t_window * t_b + safe_distance = 0
    A = 0.5 * a
    B = -a * T
    C = critical_freedom

    discriminant = B**2 - 4 * A * C
    
    t_short = (-B - np.sqrt(discriminant)) / (2 * A)

    return t_full, safe_distance, t_short




def cgt(burn_time, T, Isp, g):
    v_e = Isp * g
    m_flow = T / v_e

    mass = m_flow * burn_time
    
    density = Pres_gas * m_molN2 / (R * Temp_gas)

    volume = mass / density * 10**3


    return mass, volume



def desaturate(l, H_sat, T):
    torque = 2*T*l
    time = H_sat / torque

    return time








t_full, safe_distance,  t_short = exit_check(T_cg, v_rel, critical_freedom, M, approach_distance)
m_gas_exit, V_gas_exit = cgt(t_short, T_cg, Isp_cg, g)
t_burn_desat =  desaturate(l, H_sat, T_cg)
m_gas_desat, V_gas_desat = cgt(t_burn_desat, T_cg, Isp_cg, g)


#print("Full-burn escape time: ",t_full,"s")
#print("Safe distance reached:",safe_distance, "m")
    
print("Partial-burn time for safe escape:", t_short, "s")
print("Total burn time for", n_desat,"desaturations:", n_desat * t_burn_desat, "s")
print("Total volume of gas:", m_gas_exit, V_gas_exit + 2 * V_gas_desat, "L")





slew_rate = 0.5 * np.pi / 200 # [rad/s]
momentum = slew_rate * I_y # [Nm*s]

total_momentum = momentum * num_slew * num_rendezvous * saturation_factor # [Nm*s]





print(f"Momentum for 1 slew maneuver: {momentum} [Nm*s]")
print(f"Total momentum per wheel required for slew maneuvers: {total_momentum} [Nm*s]")








