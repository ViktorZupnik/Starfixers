"version 2"
import numpy as np
import matplotlib.pyplot as plt
import itertools
# Constants
Isp = 342
g0 = 9.80665  # Standard gravity in m/s^2
mu = 3.986 * 10**14  # Gravitational parameter for Earth in m^3/s^2
# Function to compute debris delta-V
def DebrisDeltaV(T, eta, t, md):
    return T * eta * t / md

# Function to compute optimal Vro for a given thrust
def OptVro(t, T, eta, md, M, Sro, Vros, Isp):
    for Vro in Vros:
        st = []
        for i in t:
            m_dot = T / (Isp * g0)
            m_i = M - m_dot * i
            if m_i <= 0:
                break
            term1 = Sro + Vro * i
            term2 = -T * eta / (2 * md) * i**2
            term3 = -Isp * g0 * (
                (m_i * np.log(m_i) - M * np.log(M)) / (T / (Isp * g0)) + i + np.log(M) * i
            )
            s = term1 + term2 + term3
            st.append(s)
        if st and np.max(st) >= -2:
            return Vro
    return None

# Function to compute s_r(t)
def sr(t, T, eta, md, M, Sro, Vro, Isp):
    m_dot = T / (Isp * g0)
    m_i = M - m_dot * t
    s = np.full_like(t, np.nan)
    valid = m_i > 0
    term1 = Sro + Vro * t
    term2 = -T * eta / (2 * md) * t**2
    term3 = -Isp * g0 * (
        (m_i * np.log(m_i) - M * np.log(M)) / (T / (Isp * g0)) + t + np.log(M) * t
    )
    s[valid] = term1[valid] + term2[valid] + term3[valid]
    return s

# Time when s_r(t) crosses -5 again
def TimeUnder5m(srt, t):
    for i in range(1, len(t)):
        if srt[i - 1] > -5 and srt[i] < -5:
            return t[i]
    return None


#------------------------------------------------------MAIN DELTA V SIMULATION--------------------------------------


'''---------explanation of naming convention for speeds ---------------------
Vd = velocity of debris (starting in circular orbit)
Vma = Velocity of spacecraft starting in Vd - relative velocity 
D_Vbd = The burn delta V applied to the debris after each burn 
D_Vbdtot = total delta V applied to debris after X rendezvous
D_Vbm = change of delta V applied to ourselves after a burn 
D_Vma = delta V for next rendez vous
D_Vtot =  total delta V required, taking into account delta V for next rdv and delta V applied during rdv
V_trans = velocity for transfer to next debris'''



#find delta V required to meet after two orbits
def twoorbit(Vd, Vma,hd):                    
    ad = 1/(2/((hd*1000)+6371000)-Vd**2/mu)                          #calculate new semi major axis of debris
    Td = 2*np.pi*np.sqrt(ad**3/mu)                          #calculate new orbital period of debris
    am = 1/(2/((hd*1000)+6371000)-Vma**2/mu)                          #calculate new semi major axis of spacecraft
    Tm = 2*np.pi*np.sqrt(am**3/mu)                          #calculate new orbital period of spacecraft
    T_desired = 2*Td-Tm                                     #Desired new orbital period: debris does two periods in same orbit while we change orbit
    a_desired = ((T_desired/(2*np.pi))**(2/3))*mu**(1/3)    #new semi maor axis of the spacecraft
    V_desired = np.sqrt(mu*(2/((hd*1000)+6371000)-1/a_desired))      #desired velocity at apogee to meet on time
    dVma = V_desired - Vma                                      #delta V required

    return dVma

def dvbtot(h):
    # h in km
    r0 = (h+6371) * 1000  # Initial radius in meters
    r1 = (381 + 6371) * 1000  # Final radius in meters
# Velocity at final radius for elliptical orbit
    dv = np.sqrt(mu / r0) - np.sqrt((2/r0-1/((r0 + r1) / 2)) * mu)  # Delta V required to change from circular orbit to elliptical orbit'
    return dv

def period (ha, va):
    # ha in km, va in m/s
    ha *= 1000  # Convert altitude from km to m
    sma = 1/(2/(ha+6371000)-va**2/mu)
    return 2 * np.pi * np.sqrt(sma**3 / mu)  # Period in seconds

def phi_change(h,t):
    if (h + 6371) <= 0:
        raise ValueError(f"Invalid altitude h={h}: must be >= -6371")
    T = 2*np.pi*np.sqrt((((h+6371)*1000)**3)/mu)
    rest = t%T
    d_phi = rest/T*360
    return d_phi

def catchup_phi (a_adr, h_debris):
    # h_debris in km, a_adr in m
    h_debris *= 1000  # Convert altitude from km to m
    T_debris = 2 * np.pi * np.sqrt((h_debris + 6371000)**3 / mu)  # Orbital period of debris
    T_adr = 2 * np.pi * np.sqrt(a_adr**3 / mu)  # Orbital period of our spacecraft
    delta_T = T_adr - T_debris # Time difference in seconds
    d_phi = (delta_T / T_debris) * 360  # Phase angle change in degrees
    return d_phi

eta = 0.2
#for 10 debris 
'''md = [250, 250, 250, 250, 250, 500, 500, 500, 500, 500]  # Mass of debris in kg
h = [600, 600, 600, 600, 600, 600, 600, 600, 600, 600]  # Altitude of debris in km
phi = np.array([0, 15, 25, 35, 45, 55, 65, 75, 85, 95]) '''
#for 3 debris 
md = [500,300, 250]  # Mass of debris in kg
h = [530, 555, 600] # Altitude of debris in km
phi = np.array([20, 15, 25])# Phase angle of debris in degrees
minimum_dv = 10000000
sequence = []
'''k = [2,0,1]
phi_copy = (phi-phi[k][0]) %360
print (phi_copy)'''
debris_array = np.column_stack((h, md, phi)).astype(float)  # shape (10, 3)

indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] 
indices = [0, 1, 2]# Indices for debris
indicess = []
for perm in itertools.permutations(indices):
    indicess.append(list(perm))

for k in range(len(indicess)):
    phi_copy = ((phi - phi[indicess[k][0]]) % 360).astype(float) 
    time = 0 
    M = 1494.3 
    Mi=M
    Sro = -5
    Vros = np.linspace(0, 3, 100)   # Wider range for Vro
    t = np.linspace(0.1, 15, 150)     # Time vector (start from 0.1 to avoid log(0))
    T = 465  
    Vd = np.sqrt(mu/((debris_array[indicess[k][0],0]+6371)*1000))
    Vro = OptVro(t, T, eta, debris_array[indicess[k][0],1], M, Sro, Vros, Isp) 
    Vma = Vd - Vro
    D_Vtot = np.abs(Vro)
    M = M/(np.exp(np.abs(Vro)/(Isp*g0)))
    bs = 0
    for j in range(len(debris_array)):   #10 debris 
        i = indicess[k][j]
        
        D_Vbdtot= 0
        Vd = np.sqrt(mu/((debris_array[i,0]+6371)*1000))                                            #debris velocity in circular orbit
        b = 0                                                                   
        while D_Vbdtot <= dvbtot(debris_array[i,0]):                                                   #stop the while loop when delta V applied to debris is enough to deorbit
            b +=1                                                                  #number of rdv per debris
            # Vro = OptVro(t, T, eta, debris_array[i,1], M, Sro, Vros, Isp)                      #update Vro with new mass
            Vro = 2.5
            srt = sr(t, T, eta, debris_array[i,1], M, Sro, Vro, Isp)
            t_under_5 = TimeUnder5m(srt, t)                                        #update time between 2-5m
            D_Vbd = T * eta *t_under_5/debris_array[i,1]                           #Delta V applied to debris for this rdv
            Vd = Vd - D_Vbd  
            D_Vbdtot += D_Vbd                                                      #update total debris velocity change                                              
            D_Vbm = Isp*g0*np.log(M/(M-t_under_5*T/(Isp*g0)))                      #Delta V applied to ourselves during burn
            Vma = Vma + D_Vbm                                                        #update our spacecraft velocity
            M = M/(np.exp(np.abs(D_Vbm)/(Isp*g0)))                                         #update our mass after momentum transfer   
                                             #update time after rdv
            #check if it was the last burn 
            if D_Vbdtot >= dvbtot(debris_array[i,0]):
                time += period(debris_array[i,0], Vma)/2
                break                        
            time += period(debris_array[i,0], Vma)              
            D_Vma = twoorbit(Vd,Vma,debris_array[i,0])                                                 #see function explanation above
            Vma += D_Vma                                                             #update velocity
            M = M/(np.exp(np.abs(D_Vma)/(Isp*g0)))  
            time += period(debris_array[i,0], Vma)                                        #update mass
            # Vro = OptVro(t, T, eta, debris_array[i,1], M, Sro, Vros, Isp)                      #update Vro
            Vro = 2.5
            D_V_corr2 = (Vd-Vma) - Vro                                              #correction to achieve desired relative velocity
            Vma += D_V_corr2 
            M = M/(np.exp(np.abs(D_V_corr2)/(Isp*g0)))                             #update Vma again to prepare for new momentum transfer
            D_Vtot = D_Vtot + np.abs(D_Vma) + np.abs(D_Vbm) + np.abs(D_V_corr2)
            
        bs += b

        if j < len(debris_array)-1:  
                      
            i_next = indicess[k][j+1]         
            # Vro = OptVro(t, T, eta, debris_array[i,1], M, Sro, Vros, Isp)                                                             #not take extra transfer into account for last debris (EOL)
            Vro = 2.5
            a1 = 1/(2/(debris_array[i, 0]*1000+6371000)-Vma**2/mu) #m
            hp = 2*a1-(debris_array[i, 0]*1000+6371000)-6371000 #periapsis of our orbit m
            a2 = (hp + 6371000*2 + debris_array[i_next , 0]*1000)/2  #semi major axis of next debris m 
            vp1 = np.sqrt(mu*(2/(hp+6371000)-1/(a1))) # m/s 
            vp2 = np.sqrt(mu*(2/(hp+6371000)-1/(a2))) # m/s
            D_V_Peri = vp2 - vp1
            Vma = np.sqrt(mu*(2/((debris_array[i_next,0]+6371)*1000)-1/(a2)))   #find Vma at apogee m/s
            D_Vtot += np.abs(D_V_Peri)
            M = M/(np.exp(np.abs(D_V_Peri)/(Isp*g0)))  #update mass after transfer velocity
            time += period(debris_array[i_next,0], Vma)/2  #update time after transfer velocity

            for w in range(j+1,len(debris_array)):
                z = indicess[k][w]
                phi_copy[z] = (phi_copy[z] + phi_change(debris_array[z,0],time))%360
                
            time = 0
            while phi_copy[i_next ] > np.abs(catchup_phi(a2, debris_array[i_next ,0])):
                delta_phi = catchup_phi(a2, debris_array[i_next ,0])  #calculate phase angle change for next rendezvous
                phi_copy[i_next ]= (phi_copy[i_next] + delta_phi)  #update phase angle for next rendezvous
                time += period(debris_array[i_next ,0], Vma)  #update time after phase angle change
            Vc = np.sqrt(mu/((debris_array[i_next ,0]+6371)*1000))  #velocity of debris in circular orbit
            max_delta_phi = (period(debris_array[i_next ,0], Vc) - period(debris_array[i_next ,0], Vc-5))/(period(debris_array[i_next ,0], Vc))*360  #(0.7)maximum phase angle change for next rendezvous
            N = 1
            while np.abs(phi_copy[i_next ]/N) > np.abs(max_delta_phi):
                N += 1
            delta_phi = -phi_copy[i_next ]/N  #calculate phase angle change for next rendezvous
            T_adr = delta_phi*period(debris_array[i_next ,0], Vc)/360+ period(debris_array[i_next,0], Vc)
            a_adr = ((T_adr/(2*np.pi))**(2/3))*mu**(1/3)  #semi major axis of our spacecraft after phase angle change
            D_V_trans = np.sqrt(mu*(2/((debris_array[i_next ,0]+6371)*1000)-1/a_adr))-Vma  #update Vma for next rendezvous
            D_Vtot += np.abs(D_V_trans)
            Vma += D_V_trans
            M = M/(np.exp(np.abs(D_V_trans)/(Isp*g0)))  #update mass after transfer velocity
            time += N*period(debris_array[i_next ,0], Vma)  #update time after transfer velocity
            

            D_V_corr = np.sqrt(mu/((debris_array[i_next,0] +6371)*1000)) -Vma - Vro
            Vma += D_V_corr              #add transfer velocity to rdv with new debris 
            D_Vtot = D_Vtot + np.abs(D_V_corr)
            M = M/(np.exp(np.abs(D_V_corr)/(Isp*g0)))  
        
              

    fuel_mass = Mi - Mi/(np.exp(D_Vtot/(Isp*g0)))
    # print(f'fuel mass:{fuel_mass}', f'and dry mass is {M}' )
    # print (bs)
    # print(indicess[k])
    if D_Vtot < minimum_dv:
        minimum_dv = D_Vtot
        sequence = indicess[k]
    print(D_Vtot)
    print (k)
print (f'minimum delta V: {minimum_dv} for sequence {sequence}')

# Calculate correct delta-V with correct dry mass

