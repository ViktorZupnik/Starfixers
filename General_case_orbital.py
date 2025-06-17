"version 2"
import numpy as np
import matplotlib.pyplot as plt
from Efficiency import calculate_efficiency

# Constants
#Ts = np.linspace(100, 1000, 50)  # Test multiple thrusts for thrust optimization 
debris_amount = 10  # Number of debris objects
md = [260]*debris_amount
#520kg debris, 1340kg fuel, 470kg dry
#260kg debris, 778kg fuel, 470kg dry

fs = {                      #Fart settings
    'or': None,  # Radius of the debris object in meters, set to None for rectangular objects
    'ow': 3.1,     # Width of the debris object in meters, set to None for circular objects
    'ol': 10.9,   # Length of the debris object in meters, set to None for circular objects
    'hca': np.radians(15),  # Half cone angle in radians
    'res': 50,      # Resolution for the efficiency calculation
    'dif': 0.5,        # Diffusivity factor, can be adjusted based on exhaust characteristics
    }
M = 646.2 # Initial mass of the spacecraft in kg
Mi=M
Isp = 342
g0 = 9.80665
Vros = np.arange(3.5, 8, 0.1)   # Wider range for Vro
t = np.linspace(0.1, 15, 150)     # Time vector (start from 0.1 to avoid log(0))
T = 465
op_dist = 9  
Sro = -op_dist
minimum_distance_range = (-4, -3) #Need some sources on this
M_dry = 264

# Function to compute debris delta-V
def DebrisDeltaV(T, eta, t, md):
    return T * eta * t / md

# Function to compute optimal Vro for a given thrust
def OptVro(t, T, md, M, Sro, Vros, Isp, half_cone_angle=fs['hca'], resolution=fs['res'], diffuse=fs['dif'], object_radius=fs['or'], object_width=fs['ow'], object_length=fs['ol']):
    found = False
    step = Vros[1] - Vros[0]
    start = Vros[0]
    stop = Vros[-1]
    while not found:    
        for Vro in Vros:
            st = []
            for i in t:
                m_dot = T / (Isp * g0)
                m_i = M - m_dot * i
                if m_i <= 0:
                    break
                term1 = Sro + Vro * i
                if not st:
                    s = Sro
                eta = calculate_efficiency(half_cone_angle, s, resolution, diffuse, object_radius=object_radius, object_width=object_width, object_length=object_length)
                term2 = -T * eta / (2 * md) * i**2
                term3 = -Isp * g0 * (
                    (m_i * np.log(m_i) - M * np.log(M)) / (T / (Isp * g0)) + i + np.log(M) * i
                )
                s = term1 + term2 + term3
                st.append(s)
            if st and minimum_distance_range[0] <= np.max(st) <= minimum_distance_range[1]:
                found = True
                return Vro
        step /= 2
        Vros = np.arange(start, stop, step)
        if step < 0.001:
            print("No suitable Vro found within the specified range.")
            break  # Prevent infinite loop if no suitable Vro is found
    return None

# Function to compute s_r(t)
def sr(t, T, md, M, Sro, Vro, Isp, half_cone_angle=fs['hca'], resolution=fs['res'], diffuse=fs['dif'], object_radius=fs['or'], object_width=fs['ow'], object_length=fs['ol']):
    m_dot = T / (Isp * g0)
    sr = []
    for i in range(len(t)):
        m_i = M - m_dot * t[i]
        term1 = Sro + Vro * t[i]
        if not sr:
            s = Sro
        eta = calculate_efficiency(half_cone_angle, s, resolution, diffuse, object_radius=object_radius, object_width=object_width, object_length=object_length)
        term2 = -T * eta / (2 * md) * t[i]**2
        term3 = -Isp * g0 * (
            (m_i * np.log(m_i) - M * np.log(M)) / (T / (Isp * g0)) + t[i] + np.log(M) * t[i]
        )
        s = term1 + term2 + term3
        sr.append(s)
    return np.array(sr)

# Time when s_r(t) crosses -5 again
def TimeUnderdistm(srt, t, dist):
    for i in range(1, len(t)):
        if srt[i - 1] > -dist and srt[i] < -dist:
            return t[i]
    return None

def calculate_DV_debris(T, t, md, t_under_10, srt, half_cone_angle=fs['hca'], resolution=fs['res'], diffuse=fs['dif'], object_radius=fs['or'], object_width=fs['ow'], object_length=fs['ol']):
    """
    Calculate the delta-V applied to debris.
    """
    DV = 0
    srt = srt[t <= t_under_10][:-1]  # Filter srt to only include times up to t_under_10

    for s in srt:
        eta = calculate_efficiency(half_cone_angle, s, resolution, diffuse, object_radius=object_radius, object_width=object_width, object_length=object_length)
        DV += T * eta * (t[1]-t[0]) / md
    return DV


#------------------------------------------------------MAIN DELTA V SIMULATION--------------------------------------


'''---------explanation of naming convention for speeds ---------------------
Vd = velocity of debris (starting in circular orbit)
Vm = Velocity of spacecraft starting in Vd - relative velocity 
D_Vbd = The burn delta V applied to the debris after each burn 
D_Vbdtot = total delta V applied to debris after X rendezvous
D_Vbm = change of delta V applied to ourselves after a burn 
D_Vm = delta V for next rendez vous
D_Vtot =  total delta V required, taking into account delta V for next rdv and delta V applied during rdv
V_trans = velocity for transfer to next debris'''

Vd = np.sqrt(3.986*10**14/((600+6371)*1000))
Vro = OptVro(t, T, md[0], M, Sro, Vros, Isp) 
Vm = Vd - Vro
D_Vtot = np.abs(Vro)
M = M/(np.exp(np.abs(Vro)/(Isp*g0)))
mu = 3.986*10**14 #in SI units

#find delta V required to meet after two orbits
def twoorbit(Vd, Vm):                    
    ad = 1/(2/(600000+6371000)-Vd**2/mu)                          #calculate new semi major axis of debris
    Td = 2*np.pi*np.sqrt(ad**3/mu)                          #calculate new orbital period of debris
    am = 1/(2/(600000+6371000)-Vm**2/mu)                          #calculate new semi major axis of spacecraft
    Tm = 2*np.pi*np.sqrt(am**3/mu)                          #calculate new orbital period of spacecraft
    T_desired = 2*Td-Tm                                     #Desired new orbital period: debris does two periods in same orbit while we change orbit
    a_desired = ((T_desired/(2*np.pi))**(2/3))*mu**(1/3)    #new semi maor axis of the spacecraft
    V_desired = np.sqrt(mu*(2/(600000+6371000)-1/a_desired))      #desired velocity at apogee to meet on time
    dvm = Vm-V_desired                                      #delta V required
   
    return dvm

def SmaandE(vm):
    sma = 1/(2/(600000+6371000)-vm**2/mu)                          #calculate new semi major axis of spacecraft
    rp = sma*2-600000-6371000
    e = (600000+6371000-rp)/(600000+rp+6371000*2)
    return sma, e

bs = 0 
for i in range(debris_amount):   #10 debris 

    D_Vbdtot= 0
    Vd = np.sqrt(mu/((600+6371)*1000))  
    b = 0                                                                   
    while D_Vbdtot <= 60.58:                                                   #stop the while loop when delta V applied to debris is enough to deorbit
        b +=1                                                              #number of rdv per debris
        Vro = OptVro(t, T, md[i], M, Sro, Vros, Isp)   
        print(Vro)                   #update Vro with new mass
        srt = sr(t, T, md[i], M, Sro, Vro, Isp)
        t_under_10 = TimeUnderdistm(srt, t, op_dist)                                   #update time between 2-x m
        D_Vbd = calculate_DV_debris(T, t, md[i], t_under_10, srt)               #Delta V applied to debris for this rdv
        Vd = Vd - D_Vbd  
        D_Vbdtot += D_Vbd                                                      #update total debris velocity change                                              
        D_Vbm = Isp*g0*np.log(M/(M-t_under_10*T/(Isp*g0)))                      #Delta V applied to ourselves during burn
        Vm = Vm + D_Vbm                                                        #update our spacecraft velocity
        M = M/(np.exp(np.abs(D_Vbm)/(Isp*g0)))                                         #update our mass after momentum transfer   
        if M < M_dry:  
            print("problem")                                              #check if we have enough fuel left
        #check if it was the last burn 
        if D_Vbdtot >= 60.58:
           print('last burn')
           break                                      
        D_Vm = twoorbit(Vd,Vm)                                                 #see function explanation above
        Vm -= D_Vm                                                             #update velocity
        M = M/(np.exp(np.abs(D_Vm)/(Isp*g0))) 
        if M < M_dry:  
            print("problem")                                           #update mass
        Vro = OptVro(t, T, md[i], M, Sro, Vros, Isp)                    #update Vro
        print(Vro)
        D_V_corr2 = (Vd-Vm) - Vro                                              #correction to achieve desired relative velocity
        Vm += D_V_corr2 
        M = M/(np.exp(np.abs(D_V_corr2)/(Isp*g0))) 
        if M < M_dry:  
            print("problem")                              #update Vm again to prepare for new momentum transfer
        D_Vtot = D_Vtot + D_Vm + D_Vbm + np.abs(D_V_corr2)
        
    bs += b

    #print(f'number of rdv for debris{i+1}: {b}')

    if i < debris_amount - 1:  # If not the last debris, calculate transfer velocity 
        Vro = OptVro(t, T, md[i], M, Sro, Vros, Isp)                 #not take extra transfer into account for last debris (EOL)
        print(Vro)
        D_V_trans = np.sqrt(3.986*10**14/((600+6371)*1000)) -Vm -Vro
        Vm += D_V_trans              #add transfer velocity to rdv with new debris 
        D_Vtot = D_Vtot + np.abs(D_V_trans)
        M = M/(np.exp(np.abs(D_V_trans)/(Isp*g0)))  
        if M < M_dry:  
            print("problem")  
     
        
    print(f'delta-V {i+1}: {D_Vtot}', 'md', md[i])                         # total delta V for all debris and all manoeuvres

fuel_mass = Mi - M                         
print(f'fuel mass:{fuel_mass}', f'and dry mass is {M}' )
print (bs)

#----------------------Run optimization------------------------------
#valid_Ts, deltaVs, valid_v = OptimizeThrust(t, Ts, eta, md, M, Sro, Vros, Isp)
'''
# Plot
plt.plot( deltaVs,valid_Ts, marker='o')
plt.xlabel('Debris ΔV (m/s)')
plt.ylabel('Thrust (N)')
plt.title('Thrust vs Optimal Debris ΔV')
plt.grid(True)
plt.show()

plt.plot(valid_v,valid_Ts, marker='o')
plt.xlabel('vros')
plt.ylabel('Thrust (N)')
plt.title('Thrust vs vros')
plt.grid(True)
plt.show()'''