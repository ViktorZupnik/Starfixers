"version 2"
import numpy as np
import matplotlib.pyplot as plt

# Constants
#Ts = np.linspace(100, 1000, 50)  # Test multiple thrusts for thrust optimization 
eta = 0.2
md = 260
M = 1300
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

# Main optimizer: loop over thrusts
#def OptimizeThrust(t, Ts, eta, md, M, Sro, Vros, Isp):
    deltaVsDebris = []
    valid_Ts = []
    valid_v = []

    for T in Ts:
        Vro = OptVro(t, T, eta, md, M, Sro, Vros, Isp)
        if Vro is None:
            continue
        srt = sr(t, T, eta, md, M, Sro, Vro, Isp)
        t_under5 = TimeUnder5m(srt, t)
        if t_under5 is None:
            continue
        dV = DebrisDeltaV(T, eta, t_under5, md)
        deltaVsDebris.append(dV)
        valid_Ts.append(T)
        valid_v.append(Vro)

    return valid_Ts, deltaVsDebris, valid_v


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

Vro = OptVro(t, 465, eta, md, M, Sro, Vros, Isp)
D_Vtot =0
print(Vro)
Vd = np.sqrt(3.986*10**14/((600+6371)*1000))
Vm = Vd-Vro

#calc delta V 

for i in range(10):   #10 debris 

    D_Vbdtot= 0
    Vd = np.sqrt(3.986*10**14/((600+6371)*1000))

    b = 0                                                                    
    while D_Vbdtot <= 60.58:                                                 #stop the while loop when delta V applied to debris is enough to deorbit
        b +=1                                                                #number of rdv per debris
        Vro = OptVro(t, 465, eta, md, M, Sro, Vros, Isp)                     #update Vro with new mass
        srt = sr(t, 465, eta, md, M, Sro, Vro, Isp)
        t_under_5 = TimeUnder5m(srt, t)                                      #update time between 2-5m
        D_Vbd = 465 * eta *t_under_5/md                                      #Delta V applied to debris for this rdv
        Vd = Vd - D_Vbd                                                      #update debris velocity
        D_Vbdtot += D_Vbd                                                    #update total debris velocity change
        D_Vbm = Isp*g0*np.log(M/(M-t_under_5*465/(Isp*g0)))                  #Delta V applied to ourselves during burn
        Vm = Vm + D_Vbm                                                      #update our velocity
        D_Vm = (Vm-Vd)+Vro                                                   # extra delta V required for next rdv
        Vm -= D_Vm                                                           #update our velocity (decelerate for next rdv)
        D_Vtot = D_Vtot + D_Vm + D_Vbm                                       #update Total Delta V required 
        M = M - t_under_5*465/(Isp*g0) - M + M/(np.exp(D_Vm/(Isp*g0)))       # update mass by taking the 2 burns into account (shooting+rdv)
    print(f'number of rdv for debris{i+1}: {b}')
    if i < 9:                                                                #not take extra transfer into account for last debris (EOL)
        V_trans = np.sqrt(3.986*10**14/((600+6371)*1000)) -Vm -Vro               #add transfer velocity to rdv with new debris 
        D_Vtot = D_Vtot + V_trans
        M = M- M + M/(np.exp(V_trans/(Isp*g0)))   
        
    print(f'delta-V {i+1}: {D_Vtot}')                         # total delta V for all debris and all manoeuvres

fuel_mass = 2100 - 2100/(np.exp(D_Vtot/(Isp*g0)))
print(f'fuel mass:{fuel_mass}', f'and dry mass is {M}' )
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