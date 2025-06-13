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

def calculate_DV_debris(T, t, md, t_under_10, srt,eta = 0.2):
    """
    Calculate the delta-V applied to debris.
    """
    DV = 0
    srt = srt[t <= t_under_10][:-1]  # Filter srt to only include times up to t_under_10

    for s in srt:
        DV += T * eta * (t[1]-t[0]) / md
    return DV

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