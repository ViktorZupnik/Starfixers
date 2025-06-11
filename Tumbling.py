# import numpy as np
# import matplotlib.pyplot as plt

# # ---- 1. CONFIGURATION ----
# # Satellite dimensions & mass (Starlink V1 example, but we only need Iz)
# a, b, c = 8.86, 0.1, 3.7   # [m]
# m = 227.0                  # [kg]

# # Compute Moments of Inertia (kg·m^2)
# Ix = (1/12) * m * (b**2 + c**2)
# Iy = (1/12) * m * (a**2 + c**2)
# Iz = (1/12) * m * (a**2 + b**2)
# I = np.array([Ix, Iy, Iz])

# # Choose an initial tumble rate about Z (rpm → rad/s)
# tumble_rpm = 3.0
# omega_mag = tumble_rpm * 2 * np.pi / 60    # rad/s

# # Initialize body‐rate so that it is purely about Z:
# omega_body = np.array([0.0, 0.0, omega_mag])  # [ωx, ωy, ωz]
# #omega_body = np.array([0.0, omega_mag, 0.0])  # [ωx, ωy, ωz]
# #omega_body = np.array([omega_mag, 0.0, 0.0])  # [ωx, ωy, ωz]

# # Initial quaternion (identity: no rotation)
# q = np.array([1.0, 0.0, 0.0, 0.0])

# # Thruster parameters (actual 475 N thruster)
# thrust_force = 475          # [N]
# efficiency  = 0.6              # 60% of force goes into usable torque
# effective_force = thrust_force * efficiency  # [N] for torque calculation
# lever_arm = 5                # [m] distance from COM to thruster

# # Angular‐rate threshold (rad/s) below which we say “detumbled”
# omega_threshold = 0.005       # rad/s

# # Simulation Time
# dt = 0.1                      # [s] timestep
# mission_duration_sec = 10.0  # [s] (for example)
# N_steps = int(mission_duration_sec / dt)

# # Storage arrays
# time_arr      = np.zeros(N_steps)
# omegaz_arr    = np.zeros(N_steps)   # just ωz each step
# omega_mag_arr = np.zeros(N_steps)
# dv_arr        = np.zeros(N_steps)

# def euler_derivatives(omega, I, damping_coeff=1e-3):
#     """
#     Compute dω/dt from Euler's equations (small linear damping added).
#     ω = [ωx, ωy, ωz], I = [Ix, Iy, Iz].
#     """
#     wx, wy, wz = omega
#     dω = np.zeros(3)
#     dω[0] = ((I[1] - I[2]) * wy * wz) / I[0] - damping_coeff * wx
#     dω[1] = ((I[2] - I[0]) * wz * wx) / I[1] - damping_coeff * wy
#     dω[2] = ((I[0] - I[1]) * wx * wy) / I[2] - damping_coeff * wz
#     return dω

# def quat_derivative(q, omega):
#     """
#     dq/dt = 0.5 * Ω(ω) * q, where q=[w,x,y,z], ω=[ox,oy,oz].
#     """
#     w, x, y, z = q
#     ox, oy, oz = omega
#     dqdt = 0.5 * np.array([
#         -x*ox - y*oy - z*oz,
#          w*ox + y*oz - z*oy,
#          w*oy - x*oz + z*ox,
#          w*oz + x*oy - y*ox
#     ])
#     return dqdt

# def normalize_quaternion(q):
#     norm = np.linalg.norm(q)
#     if norm < 1e-8:
#         return np.array([1.0, 0.0, 0.0, 0.0])
#     return q / norm

# # ---- 2. MAIN SIMULATION LOOP ----
# cumulative_dv = 0.0
# detumbled = False

# # Pre‐compute constant torque/acceleration about Z:
# torque_z    = effective_force * lever_arm          # [N·m]
# alpha_z     = torque_z / Iz                       # [rad/s^2]

# for i in range(N_steps):
#     t = i * dt
#     time_arr[i] = t

#      # 3.1 – Natural dynamics: Euler's eqns + small damping
#     dω_natural = euler_derivatives(omega_body, I)
#     omega_body += dω_natural * dt

#     # 3.2 – Integrate quaternion (attitude) if you want to track actual orientation
#     dqdt = quat_derivative(q, omega_body)
#     q    = normalize_quaternion(q + dqdt * dt)

#     # 3.3 – Current |ω| magnitude
#     omega_norm = np.linalg.norm(omega_body)
#     omega_mag_arr[i] = omega_norm

#     #only rotate about Z, we ignore Euler‐derivatives coupling for X,Y:
#     #include a tiny damping on ωz so it slowly drifts 
#     damping_coeff = 1e-3
#     omega_z = omega_body[2]
#     domegadz_dt = -damping_coeff * omega_z       # simple linear damping on w_z
#     omega_z += domegadz_dt * dt

#     # Check detumble condition (|ωz| <= threshold)
#     if abs(omega_z) <= omega_threshold:
#         #effectively detumbled about Z
#         omega_z = 0.0
#         omega_body = np.array([0.0, 0.0, omega_z])
#         omegaz_arr[i] = omega_z
#         dv_arr[i] = cumulative_dv
#         detumbled = True
#         break

#     # Otherwise, apply thruster torque about Z to reduce ωz
#     # Δωz = α_z * dt (opposite sign of current ωz)
#     domega_control = alpha_z * dt
#     if abs(omega_z) <= domega_control:
#         # If the thrust impulse would “flip” ωz, just zero it out
#         omega_z = 0.0
#     else:
#         # Subtract domega_control in the direction opposite to ωz
#         omega_z -= np.sign(omega_z) * domega_control

#     omega_body = np.array([0.0, 0.0, omega_z])
#     omegaz_arr[i] = omega_z

#     # Accumulate ΔV from firing the thruster continuously
#     # dv = (thrust_force / m) * dt
#     dv_step = (thrust_force / m) * dt
#     cumulative_dv += dv_step
#     dv_arr[i] = cumulative_dv

# # If we never broke out (i.e. never “detumbled”), record final ΔV:
# if not detumbled:
#     dv_arr[N_steps - 1] = cumulative_dv

# # ---- 3. PRINT RESULTS ----
# if detumbled:
#     print(f"Detumble complete at t = {time_arr[i]:.2f} s   (|ωz| ≤ {omega_threshold:.3f} rad/s)")
#     print(f"Total ΔV used for detumble = {cumulative_dv:.3f} m/s")
# else:
#     print("Did not achieve detumble within allowed mission time.")
#     print(f"Final ΔV consumed = {cumulative_dv:.3f} m/s over {mission_duration_sec:.1f} s")


# print(f"Satellite Inertia (kg·m^2):   Ix={Ix:.2f}, Iy={Iy:.2f}, Iz={Iz:.2f}")
# print(f"Initial Tumble Rate about Z:  {tumble_rpm:.2f} rpm ({omega_mag:.3f} rad/s)")
# print(f"Chosen ω_threshold (rad/s):   {omega_threshold:.3f}")

# # ---- 4. PLOT RESULTS ----
# # Truncate arrays to the actual final index `i`
# time_plot = time_arr[: i+1]
# omega_plot = omegaz_arr[: i+1]
# dv_plot    = dv_arr[: i+1]

# # (a) |ωz| vs. time
# plt.figure(figsize=(10, 4))
# plt.plot(time_plot, np.abs(omega_plot), lw=1.2, label="|ωz|(t)")
# plt.axhline(omega_threshold, color='r', linestyle='--',
#             label=f"ω_threshold = {omega_threshold:.3f} rad/s")
# plt.xlabel("Time (s)")
# plt.ylabel("|ωz| (rad/s)")
# plt.title("Rotation‐Rate about Z vs. Time During Detumble")
# plt.grid(True)
# plt.legend()

# # (b) Cumulative ΔV vs. time
# plt.figure(figsize=(10, 4))
# plt.plot(time_plot, dv_plot, lw=1.2, color="orange")
# plt.xlabel("Time (s)")
# plt.ylabel("Cumulative ΔV (m/s)")
# plt.title("ΔV Consumption During Detumble")
# plt.grid(True)

# plt.show()


















# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.spatial.transform import Rotation as R

# # ---- 1. CONFIGURATION ----
# # Starlink V1 satellite dimensions (meters)
# a, b, c = 8.86, 3.7, 0.1
# m = 227  # kg

# # Moments of Inertia (kg·m^2)
# Ix = (1/12) * m * (b**2 + c**2)
# Iy = (1/12) * m * (a**2 + c**2)
# Iz = (1/12) * m * (a**2 + b**2)
# I = np.array([Ix, Iy, Iz])

# # Choose a tumbling rate (rpm → rad/s)
# tumble_rpm = 1.0
# omega_mag = tumble_rpm * 2 * np.pi / 60  # rad/s

# # Simulate tumbling around a random initial axis
# np.random.seed(42)  # For reproducibility
# initial_axis = np.random.randn(3)
# initial_axis /= np.linalg.norm(initial_axis)
# omega_body = omega_mag * initial_axis  # angular velocity in body frame (rad/s)

# # Initial quaternion (identity: no rotation)
# q = np.array([1.0, 0.0, 0.0, 0.0])

# # Thruster / detumble parameters
# thrust_force = 475      # [N]  — cold gas thruster force magnitude
# efficiency = 0.6        # 50% effective momentum transfer for torque
# effective_force = thrust_force * efficiency  # [N] used for torque calculation
# lever_arm = 5.0         # [m] distance from COM to where thruster applies force

# # Decide on an angular‐rate threshold (rad/s) below which we consider “detumbled”
# omega_threshold = 0.005  # rad/s (≈ 0.86 deg/s).  Adjust as needed.

# # Simulation Time
# dt = 0.1                                # [s] time step
# mission_duration_sec = 110        # total available contact time (s)
# N_steps = int(mission_duration_sec / dt)  # number of timesteps

# # Preallocate arrays for recording
# time_arr       = np.zeros(N_steps)
# omega_mag_arr  = np.zeros(N_steps)
# dv_arr         = np.zeros(N_steps)

# # ---- 2. DYNAMICS FUNCTIONS ----
# def euler_derivatives(omega, I, damping_coeff=1e-3):
#     """
#     Compute natural time derivative of body‐rates using Euler's equations
#     (with a small linear damping term).
    
#     ω = [ωx, ωy, ωz]
#     I = [Ix, Iy, Iz]
#     """
#     wx, wy, wz = omega
#     d_omega = np.zeros(3)
#     d_omega[0] = ((I[1] - I[2]) * wy * wz) / I[0] - damping_coeff * wx
#     d_omega[1] = ((I[2] - I[0]) * wz * wx) / I[1] - damping_coeff * wy
#     d_omega[2] = ((I[0] - I[1]) * wx * wy) / I[2] - damping_coeff * wz
#     return d_omega

# def quat_derivative(q, omega):
#     """
#     Given quaternion q = [w, x, y, z] and body‐rate ω = [ox, oy, oz],
#     return dq/dt = 0.5 * Ω(ω) * q.
#     """
#     w, x, y, z = q
#     ox, oy, oz = omega
#     dqdt = 0.5 * np.array([
#         -x*ox - y*oy - z*oz,
#          w*ox + y*oz - z*oy,
#          w*oy - x*oz + z*ox,
#          w*oz + x*oy - y*ox
#     ])
#     return dqdt

# def normalize_quaternion(q):
#     norm = np.linalg.norm(q)
#     if norm < 1e-8:
#         return np.array([1.0, 0.0, 0.0, 0.0])
#     return q / norm

# # ---- 3. MAIN SIMULATION LOOP ----
# cumulative_dv = 0.0
# aligned = False

# for i in range(N_steps):
#     t = i * dt
#     time_arr[i] = t

#     # 3.1  Update body‐rates with Euler's equations + damping
#     domega_dt = euler_derivatives(omega_body, I)
#     omega_body += domega_dt * dt

#     # 3.2  Integrate quaternion
#     dqdt = quat_derivative(q, omega_body)
#     q += dqdt * dt
#     q = normalize_quaternion(q)

#     # 3.3  Compute current angular‐rate magnitude
#     omega_norm = np.linalg.norm(omega_body)
#     omega_mag_arr[i] = omega_norm

#     # 3.4  Check “detumble‐complete” condition (ω below threshold)
#     if omega_norm <= omega_threshold:
#         aligned = True
#         dv_arr[i] = cumulative_dv
#         break

#     # 3.5  If still tumbling (ω > threshold), apply detumble torque
#     #        apply a thruster to produce torque opposing ω:
#     #       - torque magnitude = (effective_force * lever_arm)
#     #       - direction = opposite to ω vector
#     #       => domega_control = (torque / I_max) * dt, applied along (-ω_unit)
#     I_max = np.max(I)  # worst‐case principal inertia

#     torque_mag = effective_force * lever_arm  # [N·m]
#     alpha_scalar = torque_mag / I_max         # [rad/s^2]
#     domega_control = alpha_scalar * dt        # how much we reduce ω‐magnitude each step

#     # Avoid dividing by zero:
#     if omega_norm > 0:
#         omega_body -= (omega_body / omega_norm) * domega_control

#     # 3.6  Accumulate ΔV from firing the thruster (for detumble)
#     #       We assume continuous thrust at thrust_force [N], so
#     #       dv = F/m * dt (the thruster is on continuously while detumbling)
#     dv_step = (thrust_force * dt) / m
#     cumulative_dv += dv_step
#     dv_arr[i] = cumulative_dv

# # If never broke out of the loop, record final dv:
# if not aligned:
#     dv_arr[N_steps-1] = cumulative_dv

# # ---- 4. POST‐PROCESS & PRINT RESULTS ----
# if aligned:
#     print(f"Detumble complete at t = {time_arr[i]:.1f} s  (ω dropped below {omega_threshold:.3f} rad/s)")
#     print(f"Total ΔV used for detumble = {cumulative_dv:.3f} m/s")
# else:
#     print("Did not achieve detumble within allowed mission time.")
#     print(f"Final ΔV consumed = {cumulative_dv:.3f} m/s in {mission_duration_sec/60:.1f} min")

# # Print summary of parameters
# ideal_dv_per_sec = thrust_force / m
# ideal_total_dv = ideal_dv_per_sec * mission_duration_sec

# # Compute what fraction of time we spent firing vs. ideal
# time_used = time_arr[i] if aligned else mission_duration_sec
# effective_time_needed_sec = (ideal_total_dv / cumulative_dv) * time_used if cumulative_dv > 0 else np.inf
# #extra_time_percent = (effective_time_needed_sec - time_used) / time_used * 100 if cumulative_dv > 0 else np.nan

# print(f"Satellite Inertia (kg·m^2):  Ix={Ix:.2f}, Iy={Iy:.2f}, Iz={Iz:.2f}")
# print(f"Initial Tumbling Rate: {tumble_rpm:.2f} rpm ({omega_mag:.3f} rad/s)")
# print(f"Chosen ω_threshold (rad/s): {omega_threshold:.3f}")
# print(f"Cumulative ΔV: {cumulative_dv:.2f} m/s")
# print(f"Ideal ΔV (no tumbling): {ideal_total_dv:.2f} m/s")
# #print(f"Extra Mission Time Needed: {extra_time_percent:.2f}%\n")

# # ---- 5. PLOTS ----
# time_plot = time_arr[:i+1]
# omega_plot = omega_mag_arr[:i+1]
# dv_plot    = dv_arr[:i+1]

# # (a) Plot angular‐speed vs. time
# plt.figure(figsize=(10, 4))
# plt.plot(time_plot, omega_plot, lw=1.0)
# plt.axhline(omega_threshold, color='r', linestyle='--', label=f'ω_threshold = {omega_threshold:.3f} rad/s')
# plt.xlabel('Time (s)')
# plt.ylabel('|ω| (rad/s)')
# plt.title('Angular Speed |ω| vs Time During Detumble')
# plt.grid(True)
# plt.legend()

# # (b) Plot cumulative ΔV vs. time
# plt.figure(figsize=(10, 4))
# plt.plot(time_plot, dv_plot, lw=1.0)
# plt.xlabel('Time (s)')
# plt.ylabel('Cumulative ΔV (m/s)')
# plt.title('ΔV Consumption During Detumble')
# plt.grid(True)

# plt.show()

import numpy as np
import matplotlib.pyplot as plt

# ==== 1. CONFIGURATION & PHYSICAL CONSTANTS ==== #

m_bus   = 128.0
m_panel = 132.0
m       = m_bus + m_panel

# Bus dims: 0.1 m in X, 3.7 m in Y, 1.5 m in Z
Lx_bus, Ly_bus, Lz_bus = 0.1, 3.7, 1.5

# Panel dims: 8.86 m in X, 3.7 m in Y, 0.1 m in Z
Lx_pan, Ly_pan, Lz_pan = 8.86, 3.7, 0.1

# Global COM
com = np.array([2.32, 0.00, 0.37])

# Offsets from part centers to global COM
d_bus   = np.array([0.00, 0.00, 0.00]) - com
r_panel = np.array([Lx_pan/2, Ly_pan/2, Lz_pan/2])
d_panel = r_panel - com

# ==== 2. INERTIA ABOUT EACH PART’S OWN CENTER ==== ####
# Bus block about its own center
Ibx0 = (1/12)*m_bus*(Ly_bus**2 + Lz_bus**2)
Iby0 = (1/12)*m_bus*(Lx_bus**2 + Lz_bus**2)
Ibz0 = (1/12)*m_bus*(Lx_bus**2 + Ly_bus**2)

# Panel (thin plate) about its own center
Ipx0 = (1/12)*m_panel*(Ly_pan**2 + Lz_pan**2)
Ipy0 = (1/12)*m_panel*(Lx_pan**2 + Lz_pan**2)
Ipz0 = (1/12)*m_panel*(Ly_pan**2 + Lx_pan**2)

# ==== 3. PARALLEL‐AXIS SHIFTS TO GLOBAL COM ==== ####
Ix_bus = Ibx0 + m_bus*(d_bus[1]**2 + d_bus[2]**2)
Iy_bus = Iby0 + m_bus*(d_bus[0]**2 + d_bus[2]**2)
Iz_bus = Ibz0 + m_bus*(d_bus[0]**2 + d_bus[1]**2)

Ix_pan = Ipx0 + m_panel*(d_panel[1]**2 + d_panel[2]**2)
Iy_pan = Ipy0 + m_panel*(d_panel[0]**2 + d_panel[2]**2)
Iz_pan = Ipz0 + m_panel*(d_panel[0]**2 + d_panel[1]**2)

# ==== 4. TOTAL INERTIA ABOUT GLOBAL COM ==== ####
Ix, Iy, Iz = Ix_bus + Ix_pan, Iy_bus + Iy_pan, Iz_bus + Iz_pan
I = np.array([Ix, Iy, Iz])

# ==== 5. DETUMBLE SETTINGS ==== ####
tumble_rpm  = 3.0
omega_mag   = tumble_rpm * 2*np.pi/60  # rad/s
thrust_force= 475.0                    # N
efficiency  = 0.60

# Lever arms (just reuse your values or recompute if needed)
lever_armx = 8.86 - 2.32
lever_army = 3.7/2

tau_x = thrust_force*efficiency*lever_armx
tau_y = thrust_force*efficiency*lever_army

alphax = tau_x / Ix
alphay = tau_y / Iy

omega_threshold = 0.005  # rad/s

dt = 0.01
T  = 10.0
N  = int(T/dt)

def simulate_tumble(axis):
    omega_body = np.zeros(3)
    omega_body[axis] = omega_mag
    dv = 0.0
    for i in range(N):
        w = omega_body[axis]
        if abs(w) <= omega_threshold:
            return i*dt, dv
        dω = (alphax if axis==0 else alphay) * dt
        if abs(w) <= dω:
            omega_body[axis] = 0.0
        else:
            omega_body[axis] -= np.sign(w)*dω
        dv += thrust_force*dt/m
    return T, dv

# ==== 6. RUN X & Y ONLY & PLOT ==== ####
axes = ['X','Y']
print("Axis | t_det (s) | ΔV (m/s)")
print("-------------------------")
results = {}
for i,ax in enumerate(axes):
    t, dv = simulate_tumble(i)
    results[ax] = (t, dv)
    print(f"  {ax}   |   {t:6.3f}   |  {dv:6.3f}")

# ω vs time
plt.figure()
for i,ax in enumerate(axes):
    hist=[]
    w_vec=np.zeros(3); w_vec[i]=omega_mag
    for _ in range(N):
        hist.append(abs(w_vec[i]))
        dω=(alphax if i==0 else alphay)*dt
        if abs(w_vec[i])<=dω: break
        w_vec[i]-=np.sign(w_vec[i])*dω
    plt.plot(np.arange(len(hist))*dt, hist, label=f"|ω_{ax}|")
plt.axhline(omega_threshold,color='r',ls='--')
plt.legend(); plt.grid(); plt.title("ω vs Time"); plt.show()

# ΔV vs time
plt.figure()
for i,ax in enumerate(axes):
    dvh=[]
    w_vec=np.zeros(3); w_vec[i]=omega_mag
    dv=0.0
    for _ in range(N):
        dvh.append(dv)
        dω = (alphax if i==0 else alphay)*dt
        if abs(w_vec[i])<=dω: break
        w_vec[i]-=np.sign(w_vec[i])*dω
        dv+=thrust_force*dt/m
    plt.plot(np.arange(len(dvh))*dt, dvh, label=f"ΔV_{ax}")
plt.legend(); plt.grid(); plt.title("ΔV vs Time"); plt.show()
