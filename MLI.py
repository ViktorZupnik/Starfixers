import numpy as np

# Constants
sigma = 5.67e-8  # W/m²K⁴
mass = 1000       # kg (aluminum)
cp = 900          # J/kg·K (aluminum spec. heat)
ext_area = 6.0    # m² (cube side ~1.15 m for 1000 kg Al)

# MLI properties (conservative)
epsilon_MLI = 0.03      # Effective emissivity
alpha_MLI = 0.1         # Solar absorptivity
U_MLI = 0.5             # W/m²K (effective heat transfer coeff)

# Environment (600 km LEO)
T_eclipse = 136.7          # K
T_sunlight = 386.7        # K
solar_flux = 1361        # W/m²
orbit_period = 96 * 60   # seconds (96 min)
eclipse_frac = 0.375    # 35% of orbit in eclipse

# Internal heat (adjust based on payload)
Q_internal = 80          # W (example)

# Time-stepping simulation (1 orbit)
time = np.linspace(0, orbit_period, 1000)
T_ext = np.where(time < eclipse_frac * orbit_period, T_eclipse, T_sunlight)
T_internal = np.ones_like(time) * 290  # Initial guess

for i in range(1, len(time)):
    dt = time[i] - time[i-1]
    Q_MLI = U_MLI * ext_area * (T_ext[i] - T_internal[i-1])
    dT = (Q_MLI + Q_internal) * dt / (mass * cp)
    T_internal[i] = T_internal[i-1] + dT

# Results
min_T = np.min(T_internal)
max_T = np.max(T_internal)
print(f"Internal temp range: {min_T:.1f} K to {max_T:.1f} K")