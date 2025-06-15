import numpy as np 
g0 = 9.80665
def sr(t, T, md, M, Sro, Vro, Isp):
    m_dot = T / (Isp * g0)
    sr = []
    for i in range(len(t)):
        m_i = M - m_dot * t[i]
        term1 = Sro + Vro * t[i]
        if not sr:
            s = Sro
        eta = 0.2
        term2 = -T * eta / (2 * md) * t[i]**2
        term3 = -Isp * g0 * (
            (m_i * np.log(m_i) - M * np.log(M)) / (T / (Isp * g0)) + t[i] + np.log(M) * t[i]
        )
        s = term1 + term2 + term3
        sr.append(s)
    return np.array(sr)

def TimeUnderdistm(srt, t, dist):
    for i in range(1, len(t)):
        if srt[i - 1] > -dist and srt[i] < -dist:
            return t[i]
    return None
mu = 3.986e14  # Gravitational parameter for Earth in m^3/s^2
def twoorbit(Vd, Vm):                    
    ad = 1/(2/(600000+6371000)-Vd**2/mu)                         
    Td = 2*np.pi*np.sqrt(ad**3/mu)                         
    am = 1/(2/(600000+6371000)-Vm**2/mu)                          
    print(am) #check correct semi major axis
    Tm = 2*np.pi*np.sqrt(am**3/mu)                          
    print(Tm)
    T_desired = 2*Td-Tm                                     
    a_desired = ((T_desired/(2*np.pi))**(2/3))*mu**(1/3)    
    V_desired = np.sqrt(mu*(2/(600000+6371000)-1/a_desired))      
    dvm = Vm-V_desired                                     
    print("Delta V required to meet debris:", dvm)          #check correct sign of delta V
    return dvm
twoorbit(7200, 7503.16)  # Example call to the twoorbit function

def calculate_DV_debris(T, t, md, t_under_10, srt, eta=0.2):
    """
    Calculate the delta-V applied to debris.
    """
    DV = 0
    srt = srt[t <= t_under_10][:-1]  # Filter srt to only include times up to t_under_10

    for s in srt:
        DV += T * eta * (t[1]-t[0]) / md
    return DV
def test_sr_basic():
    t = np.array([0,2,4,6,7])
    T = 465
    md = 260
    M = 646.2
    Sro = -9
    Vro = 3.3
    Isp = 342

    result = sr(t, T, md, M, Sro, Vro, Isp)
    print(result)
    # Expected value calculated manually or from a reliable source
    expected_value = np.array([-9.0, -4.55478, -4.41992,-8.59667,-12.30228])  #values from the plot on desmos
    assert np.allclose(result, expected_value)  # Check first few values for correctness
#test passes if the function returns the time when the underdistance occurs
def test_TimeUnderdistm_basic():
    t = np.linspace(0.1, 15, 2000)
    # srt = np.array([-9, -8,-6,-10,-11])
    srt = sr(t, 465, 260, 646.2, -9, 3.3, 342)
    dist = 9
    result = TimeUnderdistm(srt, t, dist)
    expected = 6.12473    #plot on desmos
    assert abs(result - expected) < 1e-2  # Use a tolerance for float comparison
 

#test passes if the function returns None when no crossing occurs - debris never enters the range at which we can shoot at it                       
def test_TimeUnderdistm_no_crossing():
    t = np.array([0, 1, 2])
    srt = np.array([3, 2, 1])
    dist = 2
    result = TimeUnderdistm(srt, t, dist)
    assert result is None   
              
def test_calculate_DV_debris_constant_eta():
    T = 10
    md = 2
    t = np.array([0, 1, 2, 3, 4])
    srt = np.array([0.5, 0.4, 0.3, 0.2, 0.1])
    t_under_10 = 3

    result = calculate_DV_debris(T, t, md, t_under_10, srt,0.2)

    print("Calculated DV:", result)
    expected_value = 3    #hand calculation 
    assert abs(result - expected_value) < 1e-6  # Use a tolerance for float comparison

print(calculate_DV_debris(10, np.array([0, 1, 2, 3,4]), 2, 3, np.array([0.5, 0.4, 0.3, 0.2,0.1]),0.2))


#check case if there is not time under t_under_10 and that DV is 0
def test_calculate_DV_debris_no_time_in_range():
    T = 10
    md = 2
    t = np.array([0, 1, 2])
    srt = np.array([-0.5, -0.4, - 0.3])
    t_under_10 = -1  # before first t
    DV = calculate_DV_debris(
        T, t, md, t_under_10, srt,0.2)
    assert np.isclose(DV, 0.0)

def test_calculate_DV_debris_zero_efficiency(monkeypatch):
    eta = 0
    T = 10
    md = 2
    t = np.array([0, 1, 2, 3])
    srt = np.array([0.5, 0.4, 0.3, 0.2])
    t_under_10 = 3

    DV = calculate_DV_debris(T, t, md, t_under_10, srt,eta)
    assert np.isclose(DV, 0.0)
#test with hand calculations and internet tools
def test_twoorbit():
    Vd = 7200
    Vm = 7503.16
    expected_dvm = 690.2811422079149  # Expected delta V to meet debris
    dvm = twoorbit(Vd, Vm)
    assert np.isclose(dvm, expected_dvm, atol=1e-2)  # Allow a small tolerance for floating point comparison


"version 2"
import numpy as np
import matplotlib.pyplot as plt
from efficiency import calculate_efficiency

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



# Function to compute debris delta-V
def DebrisDeltaV(T, eta, t, md):
    return T * eta * t / md

# Function to compute optimal Vro for a given thrust
def OptVro(t, T, md, M, Sro, Vros, Isp):
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
                eta = 0.2
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
M = 650 # Initial mass of the spacecraft in kg
Mi=M
Isp = 342
g0 = 9.80665
Vros = np.arange(3.5, 8, 0.1)   # Wider range for Vro
t = np.linspace(0.1, 15, 150)     # Time vector (start from 0.1 to avoid log(0))
T = 465
op_dist = 9  
Sro = -op_dist
minimum_distance_range = (-4, -3) #Need some sources on this
bs = 0 
Vd = np.sqrt(3.986*10**14/((600+6371)*1000))
Vro = OptVro(t, T, md[0], M, Sro, Vros, Isp) 
Vm = Vd - Vro
D_Vtot = np.abs(Vro)
M = M/(np.exp(np.abs(Vro)/(Isp*g0)))
mu = 3.986*10**14 #in SI units
# Function to compute s_r(t)
def sr(t, T, md, M, Sro, Vro, Isp):
    m_dot = T / (Isp * g0)
    sr = []
    for i in range(len(t)):
        m_i = M - m_dot * t[i]
        term1 = Sro + Vro * t[i]
        if not sr:
            s = Sro
        eta = 0.2
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
#check operation
def test_operation():
    M_dry = 266
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
    bs = 0 
    Vd = np.sqrt(3.986*10**14/((600+6371)*1000))
    Vro = OptVro(t, T, md[0], M, Sro, Vros, Isp) 
    Vm = Vd - Vro
    D_Vtot = np.abs(Vro)
    M = M/(np.exp(np.abs(Vro)/(Isp*g0)))
    mu = 3.986*10**14 #in SI units
    for i in range(debris_amount):   #10 debris 

        D_Vbdtot= 0
        Vd = np.sqrt(mu/((600+6371)*1000))  
        b = 0                                                                   
        while D_Vbdtot <= 60.58:                                                   #stop the while loop when delta V applied to debris is enough to deorbit
            b +=1                                                              #number of rdv per debris
            Vro = OptVro(t, T, md[i], M, Sro, Vros, Isp)                      #update Vro with new mass
            srt = sr(t, T, md[i], M, Sro, Vro, Isp)
            t_under_10 = TimeUnderdistm(srt, t, op_dist)                                   #update time between 2-x m
            D_Vbd = calculate_DV_debris(T, t, md[i], t_under_10, srt,0.2)               #Delta V applied to debris for this rdv
            #assert abs(D_Vbd - 12.555)< 10e-5      #checked with hand calculation
            Vd = Vd - D_Vbd 
            rrg = D_Vbdtot 
            D_Vbdtot += D_Vbd                                                                                               
            D_Vbm = Isp*g0*np.log(M/(M-t_under_10*T/(Isp*g0)))                      
            Vm = Vm + D_Vbm
            ghj = M                                                        
            M = M/(np.exp(np.abs(D_Vbm)/(Isp*g0))) 
            assert M<ghj and M > M_dry
                                                             

            #check if it was the last burn 
            if D_Vbdtot >= 60.58:
                print('last burn')
                break                                      
            D_Vm = twoorbit(Vd,Vm)                                       
            Vm -= D_Vm 
            am = 1/(2/(600000+6371000)-Vm**2/mu)  
            assert am>6861       #check if semi major axis is larger than 6721km (100+600+6371*2)/2 to ensure we are in orbit and don't crash into earth 
            ghj = M                                                           
            M = M/(np.exp(np.abs(D_Vm)/(Isp*g0)))                                            
            assert M<ghj and M > (M_dry)
          
            Vro = OptVro(t, T, md[i], M, Sro, Vros, Isp)                    
            D_V_corr2 = (Vd-Vm) - Vro                                              
            Vm += D_V_corr2 
            ghj = M
            M = M/(np.exp(np.abs(D_V_corr2)/(Isp*g0)))   
            assert M<ghj and M > (M_dry)
                                                          
            assert D_Vm >0
            assert D_Vbd > 0
            assert D_Vbm > 0
            D_Vtot = D_Vtot + D_Vm + D_Vbm + np.abs(D_V_corr2)
            
        bs += b

        #print(f'number of rdv for debris{i+1}: {b}')

        if i < debris_amount - 1:  # If not the last debris, calculate transfer velocity 
            Vro = OptVro(t, T, md[i], M, Sro, Vros, Isp)                 #not take extra transfer into account for last debris (EOL)
            D_V_trans = np.sqrt(3.986*10**14/((600+6371)*1000)) -Vm -Vro
            Vm += D_V_trans              #add transfer velocity to rdv with new debris 
            D_Vtot = D_Vtot + np.abs(D_V_trans)
            ghj = M
            M = M/(np.exp(np.abs(D_V_trans)/(Isp*g0)))  
            assert M<ghj and M > (M_dry)
        
            
        print(f'delta-V {i+1}: {D_Vtot}', 'md', md[i])                         # total delta V for all debris and all manoeuvres

    fuel_mass = Mi - M                         
    print(f'fuel mass:{fuel_mass}', f'and dry mass is {M}' )
    print (bs)
