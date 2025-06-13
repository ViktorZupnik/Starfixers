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
def calculate_DV_debris(T, t, md, t_under_10, srt,eta = 0.2):
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

    result = calculate_DV_debris(T, t, md, t_under_10, srt,eta = 0.2)

    print("Calculated DV:", result)
    expected_value = 3    #hand calculation 
    assert abs(result - expected_value) < 1e-6  # Use a tolerance for float comparison

print(calculate_DV_debris(10, np.array([0, 1, 2, 3,4]), 2, 3, np.array([0.5, 0.4, 0.3, 0.2,0.1])))


#check case if there is not time under t_under_10 and that DV is 0
def test_calculate_DV_debris_no_time_in_range():
    T = 10
    md = 2
    t = np.array([0, 1, 2])
    srt = np.array([-0.5, -0.4, - 0.3])
    t_under_10 = -1  # before first t
    DV = calculate_DV_debris(
        T, t, md, t_under_10, srt,eta = 0.2)
    assert np.isclose(DV, 0.0)

def test_calculate_DV_debris_zero_efficiency(monkeypatch):
    eta = 0
    T = 10
    md = 2
    t = np.array([0, 1, 2, 3])
    srt = np.array([0.5, 0.4, 0.3, 0.2])
    t_under_10 = 3

    DV = calculate_DV_debris(T, t, md, t_under_10, srt,eta=0)
    assert np.isclose(DV, 0.0)
#test with hand calculations and internet tools
def test_twoorbit():
    Vd = 7200
    Vm = 7503.16
    expected_dvm = 690.2811422079149  # Expected delta V to meet debris
    dvm = twoorbit(Vd, Vm)
    assert np.isclose(dvm, expected_dvm, atol=1e-2)  # Allow a small tolerance for floating point comparison
