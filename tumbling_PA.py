import numpy as np
import matplotlib.pyplot as plt

# ==== 1. CONFIGURATION & PHYSICAL CONSTANTS ==== #

# Starlink‐V1 “flat‐puck” dimensions (m)
a, b, c = 8.86, 0.1, 3.7    # length along Z, thickness along X, width along Y
m = 227.0                   # kg (debris mass for moment of inertia)

# Principal moments of inertia (kg·m²)
Ix = (1/12) * m * (b**2 + c**2)
Iy = (1/12) * m * (a**2 + c**2)
Iz = (1/12) * m * (a**2 + b**2)
I = np.array([Ix, Iy, Iz])  # [Iₓ, Iᵧ, I_z]

# Initial tumble rate (same for X, Y, Z cases)
tumble_rpm = 3.0
omega_mag  = tumble_rpm * 2 * np.pi / 60  # rad/s ≈ 0.314 rad/s

# Thruster parameters (475 N thruster)
thrust_force    = 475.0       # [N]
efficiency      = 0.60        # 60% of that 475 N gives usable torque
lever_arm       = 5.0         # [m] nominal distance from COM to “hit point”
effective_force = thrust_force * efficiency  # [N]
torque_mag      = effective_force * lever_arm  # [N·m]

# Per‐axis angular accelerations: αᵢ = τ / Iᵢ (rad/s²)
alpha = torque_mag / I  # array of length 3

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

        # Apply a “torque chunk” on that same axis
        domega_control = alpha[axis_index] * dt
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
axes = ['X', 'Y', 'Z']
results = {}

print("\n===== STARLINK‐V1 DETUMBLE SIMULATION =====\n")
print(f"Moments of Inertia (kg·m²):  Iₓ = {Ix:.2f},  Iᵧ = {Iy:.2f},  I_z = {Iz:.2f}")
print(f"Thruster force: {thrust_force:.1f} N  |  Efficiency: {efficiency*100:.0f}%  |  Lever arm: {lever_arm:.1f} m")
print(f"Resulting torque (τ) = {torque_mag:.1f} N·m")
print(f"Angular accelerations: αₓ = {alpha[0]:.3f},  αᵧ = {alpha[1]:.3f},  α_z = {alpha[2]:.3f}  (rad/s²)")
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
plt.title("ΔV Consumption vs Time for X, Y, Z Detumble")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.show()
