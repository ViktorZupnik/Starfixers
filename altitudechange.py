"version 2"
import numpy as np
import matplotlib.pyplot as plt

# Constants

eta = 0.2
md = np.linspace(260,260, 10)
#520kg debris, 1340kg fuel, 470kg dry
#260kg debris, 778kg fuel, 470kg dry
h = np.linspace(600,600,10)
phi = np.linspace(0,360,10)
debris_array = np.column_stack((h, md, phi))  # shape (10, 3)

M = 1200 
Mi=M
Isp = 342
g0 = 9.80665
Sro = -5
Vros = np.linspace(0, 3, 100)   # Wider range for Vro
t = np.linspace(0.1, 15, 150)     # Time vector (start from 0.1 to avoid log(0))
T = 465  


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
Vm = Velocity of spacecraft starting in Vd - relative velocity 
D_Vbd = The burn delta V applied to the debris after each burn 
D_Vbdtot = total delta V applied to debris after X rendezvous
D_Vbm = change of delta V applied to ourselves after a burn 
D_Vm = delta V for next rendez vous
D_Vtot =  total delta V required, taking into account delta V for next rdv and delta V applied during rdv
V_trans = velocity for transfer to next debris'''

Vd = np.sqrt(3.986*10**14/((debris_array[0][0]+6371)*1000))
Vro = OptVro(t, T, eta, debris_array[0][1], M, Sro, Vros, Isp) 
Vm = Vd - Vro
D_Vtot = np.abs(Vro)
M = M/(np.exp(np.abs(Vro)/(Isp*g0)))
mu = 3.986*10**14 #in SI units

#find delta V required to meet after two orbits
def twoorbit(Vd, Vm):                    
    ad = 1/(2/(debris_array[i,0]*1000+6371000)-Vd**2/mu)                          #calculate new semi major axis of debris
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
def dvbtot(h):
    # h in km
    mu = 3.986 * 10**14  # Gravitational parameter for Earth in m^3/s^2
    r0 = (h+6371) * 1000  # Initial radius in meters
    r1 = (381 + 6371) * 1000  # Final radius in meters
  # Velocity at final radius for elliptical orbit
    dv = np.sqrt(mu / r0) - np.sqrt((2/r0-1/((r0 + r1) / 2)) * mu)  # Delta V required to change from circular orbit to elliptical orbit'
    return dv


bs = 0 
for i in range(len(debris_array)):   #10 debris 

    D_Vbdtot= 0
    Vd = np.sqrt(mu/((debris_array[i,0]+6371)*1000))                                            #debris velocity in circular orbit

    if i == 0:
        print (SmaandE(Vm))
    b = 0                                                                   
    while D_Vbdtot <= dvbtot(debris_array[i,0]):                                                   #stop the while loop when delta V applied to debris is enough to deorbit
        b +=1                                                                  #number of rdv per debris
        Vro = OptVro(t, T, eta, debris_array[i,1], M, Sro, Vros, Isp)                      #update Vro with new mass
        srt = sr(t, T, eta, debris_array[i,1], M, Sro, Vro, Isp)
        t_under_5 = TimeUnder5m(srt, t)                                        #update time between 2-5m
        D_Vbd = T * eta *t_under_5/debris_array[i,1]                           #Delta V applied to debris for this rdv
        Vd = Vd - D_Vbd  
        D_Vbdtot += D_Vbd                                                      #update total debris velocity change                                              
        D_Vbm = Isp*g0*np.log(M/(M-t_under_5*T/(Isp*g0)))                      #Delta V applied to ourselves during burn
        Vm = Vm + D_Vbm                                                        #update our spacecraft velocity
        M = M/(np.exp(np.abs(D_Vbm)/(Isp*g0)))                                         #update our mass after momentum transfer   

        #check if it was the last burn 
        if D_Vbdtot >= dvbtot(debris_array[i,0]):
           break                                      
        D_Vm = twoorbit(Vd,Vm)                                                 #see function explanation above
        Vm -= D_Vm                                                             #update velocity
        M = M/(np.exp(np.abs(D_Vm)/(Isp*g0)))                                          #update mass
        Vro = OptVro(t, T, eta, debris_array[i,1], M, Sro, Vros, Isp)                      #update Vro
        D_V_corr2 = (Vd-Vm) - Vro                                              #correction to achieve desired relative velocity
        Vm += D_V_corr2 
        M = M/(np.exp(np.abs(D_V_corr2)/(Isp*g0)))                             #update Vm again to prepare for new momentum transfer
        D_Vtot = D_Vtot + D_Vm + D_Vbm + np.abs(D_V_corr2)
        
    bs += b

    #print(f'number of rdv for debris{i+1}: {b}')

    if i < 9:           
        Vro = OptVro(t, T, eta, debris_array[i,1], M, Sro, Vros, Isp)                                                             #not take extra transfer into account for last debris (EOL)
        D_V_trans = np.sqrt(3.986*10**14/((600+6371)*1000)) -Vm -Vro
        Vm += D_V_trans              #add transfer velocity to rdv with new debris 
        D_Vtot = D_Vtot + np.abs(D_V_trans)
        M = M/(np.exp(np.abs(D_V_trans)/(Isp*g0)))  
     
        
    print(f'delta-V {i+1}: {D_Vtot}', 'md', md[i])                         # total delta V for all debris and all manoeuvres

fuel_mass = Mi - Mi/(np.exp(D_Vtot/(Isp*g0)))
print(f'fuel mass:{fuel_mass}', f'and dry mass is {M}' )
print (bs)

