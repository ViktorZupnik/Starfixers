import numpy as np
import matplotlib.pyplot as plt

# ==== 1. CONFIGURATION & PHYSICAL CONSTANTS ==== #

# Starlink‐V1 “flat‐puck” dimensions (m)
#a, b, c = 8.86, 3.7, 1.5

# Masses
m_bus = 128
m_panel = 132
m=260

# Bus (stowed) dims [m]
Lx_bus, Ly_bus, Lz_bus = 0.1, 3.7, 1.5

# Panel dims [m]: span along Z, depth along Y, negligible thickness in X
Lx_pan, Ly_pan = 8.86, 3.7
Lz_pan = 0.01  # small thickness for inertia calc

com = np.array([2.32, 0.00, 0.37])

# Bus center in global‐COM frame
r_bus = np.array([0.00, 0.00, 0.00])  
d_bus = r_bus - com  # offset from global COM

r_panel = np.array([Lx_pan/2, Ly_pan/2, 0.005 ])  # ≈ [0, 0, +4.48]
d_panel = r_panel - com

# ==== 2. INERTIA ABOUT EACH PART’S OWN CENTER ==== #

# Bus block about its own center
Ibx0 = (1/12) * m_bus * (Ly_bus**2 + Lz_bus**2)  # about X
Iby0 = (1/12) * m_bus * (Lx_bus**2 + Lz_bus**2)  # about Y
Ibz0 = (1/12) * m_bus * (Lx_bus**2 + Ly_bus**2)  # about Z

# Panel (thin plate) about its own center
Ipx0 = (1/12) * m_panel * (Ly_pan**2 + Lz_pan**2)   # about X‐axis of panel
Ipy0 = (1/12) * m_panel * (Lx_pan**2 + Lz_pan**2)   # about Y‐axis of panel
Ipz0 = (1/12) * m_panel * (Ly_pan**2 + Lx_pan**2)  # about Z‐axis (normal)

# ==== 3. PARALLEL‐AXIS SHIFTS TO GLOBAL COM ==== #

# Bus shifted
Ix_bus = Ibx0 + m_bus * (d_bus[1]**2 + d_bus[2]**2)
Iy_bus = Iby0 + m_bus * (d_bus[0]**2 + d_bus[2]**2)
Iz_bus = Ibz0 + m_bus * (d_bus[0]**2 + d_bus[1]**2)

# Panel shifted
Ix_pan = Ipx0 + m_panel * (d_panel[1]**2 + d_panel[2]**2)
Iy_pan = Ipy0 + m_panel * (d_panel[0]**2 + d_panel[2]**2)
Iz_pan = Ipz0 + m_panel * (d_panel[0]**2 + d_panel[1]**2)

# ==== 4. TOTAL INERTIA ABOUT GLOBAL COM ==== #
Ix = Ix_bus + Ix_pan
Iy = Iy_bus + Iy_pan
Iz = Iz_bus + Iz_pan

I = np.array([Ix, Iy, Iz])

# Initial tumble rate (same for X, Y, Z cases)
tumble_rpm = 3.0
omega_mag  = tumble_rpm * 2 * np.pi / 60  # rad/s ≈ 0.314 rad/s

# Thruster parameters (475 N thruster)
thrust_force    = 475.0       # [N]
efficiency      = 1        # 60% of that 475 N gives usable torque
      
#Lever arm calculations       # [m] nominal distance from COM to “hit point”
lever_armx=8.86-2.32
lever_army=3.7/2


effective_force = thrust_force * efficiency  # [N]
torque_magx      = effective_force * lever_armx  # [N·m]
torque_magy      = effective_force * lever_army  # [N·m]

# Per‐axis angular accelerations: αᵢ = τ / Iᵢ (rad/s²)
alphax_scalar = torque_magx / I[0]
alphay_scalar = torque_magy / I[1]

# “Detumble” threshold (rad/s)
omega_threshold = 0.005  # if |ω| falls below this, we call it “detumbled”

# Simulation time parameters
dt                 = 0.01                     # [s] timestep
mission_duration_s = 10.0                     # [s]
N_steps            = int(mission_duration_s / dt)  # total timesteps


# ==== 2. DETUMBLE SIMULATION FUNCTION ==== #
def simulate_tumble(axis_index):
    """
    Simulate a non‐contact “remote detumble” about a single principal axis:
      axis_index = 0 for X, 1 for Y, 2 for Z.

    Returns:
      - t_det: time (s) when |ω| ≤ omega_threshold
      - dv_used: cumulative ΔV (m/s) used up to that moment
      - time_hist: array of times at each step
      - omega_hist: array of ω_component (rad/s) at each step
      - dv_hist:    array of cumulative ΔV at each step
    """
    # Initialize body‐rate: tumble only about the chosen axis
    omega_body = np.zeros(3)
    omega_body[axis_index] = omega_mag

    # Storage for this axis
    time_hist  = []
    omega_hist = []
    dv_hist    = []

    cumulative_dv = 0.0

    for i in range(N_steps):
        t = i * dt
        time_hist.append(t)

        # Record the current tumble rate on that axis
        omega_comp = omega_body[axis_index]
        omega_hist.append(omega_comp)

        # Check if magnitude has dropped below threshold
        if abs(omega_comp) <= omega_threshold:
            dv_hist.append(cumulative_dv)
            return t, cumulative_dv, np.array(time_hist), np.array(omega_hist), np.array(dv_hist)

        # Apply a “torque chunk” on that same axis, using the correct α for X, Y or Z
        if axis_index == 0:
            domega_control = alphax_scalar * dt
        elif axis_index == 1:
            domega_control = alphay_scalar * dt
       
        if abs(omega_comp) <= domega_control:
            omega_body[axis_index] = 0.0
        else:
            omega_body[axis_index] -= np.sign(omega_comp) * domega_control

        # Accumulate ΔV from the thruster over this dt
        dv_step       = (thrust_force * dt) / m
        cumulative_dv += dv_step
        dv_hist.append(cumulative_dv)

    # If loop finishes without dropping below threshold
    return mission_duration_s, cumulative_dv, np.array(time_hist), np.array(omega_hist), np.array(dv_hist)


# ==== 3. RUN ALL THREE AXES & PRINT RESULTS ==== #
axes = ['X', 'Y']
results = {}

print("\n===== STARLINK‐V1 DETUMBLE SIMULATION =====\n")
print(f"Moments of Inertia (kg·m²):  Iₓ = {Ix:.2f},  Iᵧ = {Iy:.2f},  I_z = {Iz:.2f}")
print(f"Thruster force: {thrust_force:.1f} N  |  Efficiency: {efficiency*100:.0f}%  |  Lever arm: {lever_armx:.1f} m | Lever arm: {lever_army:.1f} m")
print(f"Resulting torque around x (τ) = {torque_magx:.1f} N·m")
print(f"Resulting torque around y (τ) = {torque_magy:.1f} N·m")
print(f"Angular accelerations: αₓ = {alphax_scalar:.3f},  αᵧ = {alphay_scalar:.3f}  (rad/s²)")
print(f"Initial tumble rate: {tumble_rpm:.2f} rpm  = {omega_mag:.3f} rad/s")
print(f"Detumble threshold: |ω| ≤ {omega_threshold:.3f} rad/s\n")
print(" Axis   |  Time to detumble (s)  |  ΔV used (m/s) ")
print("--------+-------------------------+----------------")

for idx, ax in enumerate(axes):
    t_det, dv_used, time_hist, omega_hist, dv_hist = simulate_tumble(idx)
    results[ax] = (t_det, dv_used, time_hist, omega_hist, dv_hist)
    print(f"   {ax}    |       {t_det:>6.2f} s        |     {dv_used:>6.3f}   m/s")


# ==== 4. CHART & DIAGRAMS ==== #

# 4.2 – Plot each axis's |ω| (which equals |ω_i| since only that component is nonzero)
plt.figure(figsize=(10, 4))
for idx, ax in enumerate(axes):
    time_hist, omega_hist = results[ax][2], results[ax][3]
    plt.plot(time_hist, np.abs(omega_hist), label=f"|ω_{ax}|(t)")
plt.axhline(omega_threshold, color='red', linestyle='--', label="ω_threshold")
plt.xlabel("Time (s)")
plt.ylabel("|ω| (rad/s)")
plt.title("Magnitude of ω vs Time During Detumble (Single‐Axis Cases)")
plt.grid(True)
plt.legend()
plt.tight_layout()

# 4.4 – Plot ΔV vs time for each axis
plt.figure(figsize=(10, 4))
for idx, ax in enumerate(axes):
    time_hist, dv_hist = results[ax][2], results[ax][4]
    plt.plot(time_hist, dv_hist, label=f"ΔV until detumble about {ax}")
plt.xlabel("Time (s)")
plt.ylabel("Cumulative ΔV (m/s)")
plt.title("ΔV Consumption vs Time for X, Y  Detumble")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.show()
