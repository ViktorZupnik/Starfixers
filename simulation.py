import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# === Kepler's Equation Solver ===
def solve_kepler(M, e, tol=1e-8):
    E = M.copy()
    for _ in range(100):
        delta = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= delta
        if np.max(np.abs(delta)) < tol:
            break
    return E

# === Position computation for interpolated perigee ===
def compute_position_from_rp(r_p, elapsed_t):
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)
    T = 2 * np.pi * np.sqrt(a**3 / mu)
    M = 2 * np.pi * elapsed_t / T
    E = solve_kepler(np.array([M]), e)[0]
    theta = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2),
                           np.sqrt(1 - e) * np.cos(E / 2))
    r = a * (1 - e**2) / (1 + e * np.cos(theta))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

# === Inputs ===
perigees = input("Enter perigee radii (comma-separated, km): ")
perigees = np.array([float(r.strip()) for r in perigees.split(",")])  # array of perigees
r_a = float(input("Enter fixed apogee radius (km): "))                # fixed apogee

# Constants
mu = 398600.4418  # Earth's gravitational parameter [km^3/s^2]

# Settings
num_points_per_orbit = 1000  # frames per orbit segment
frame_interval = 10          # ms between frames
speed_factor = 1             # speed multiplier (1 = real time)
transition_frames = 50       # frames for smooth transition

# Precompute orbits for each perigee
orbits = []
periods = []
for r_p in perigees:
    a = (r_p + r_a) / 2
    e = (r_a - r_p) / (r_a + r_p)
    T = 2 * np.pi * np.sqrt(a**3 / mu)
    periods.append(T)
    
    t = np.linspace(0, T, num_points_per_orbit)
    M = 2 * np.pi * t / T
    E = solve_kepler(M, e)
    theta = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2),
                           np.sqrt(1 - e) * np.cos(E / 2))
    r = a * (1 - e**2) / (1 + e * np.cos(theta))
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # Rotate arrays so orbit starts and ends at apogee (theta=pi)
    idx_apogee = np.argmin(np.abs(theta - np.pi))
    x = np.roll(x, -idx_apogee)
    y = np.roll(y, -idx_apogee)
    theta = np.roll(theta, -idx_apogee)

    orbits.append((x, y, theta))

# Circular orbit parameters for the second satellite
r_circ = r_a
T_circ = 2 * np.pi * np.sqrt(r_circ**3 / mu)

# Total animation time and frames
total_time = sum(periods) / speed_factor
fps = 1000 / frame_interval
total_transition_time = transition_frames / fps
total_frames = int(total_time * fps + total_transition_time * (len(orbits) - 1))

# Fix broadcasting issue here:
all_x = np.hstack([orbit[0] for orbit in orbits] + [np.array([r_circ, -r_circ])])
all_y = np.hstack([orbit[1] for orbit in orbits] + [np.array([r_circ, -r_circ])])
margin = 0.1 * max(np.max(all_x) - np.min(all_x), np.max(all_y) - np.min(all_y))

# Cumulative periods for orbit switching
cumulative_periods = np.cumsum(periods) / speed_factor

# Fix cumulative times addition to same length arrays
n = len(orbits)
extended_cumulative_times = np.insert(cumulative_periods, 0, 0) + np.concatenate([np.arange(n) * total_transition_time, [0]])

# --- Precompute circular orbit path ---
theta_circ_path = np.linspace(0, 2*np.pi, 1000)
x_circ_path = r_circ * np.cos(theta_circ_path)
y_circ_path = r_circ * np.sin(theta_circ_path)

# === Plot Setup ===
fig, ax = plt.subplots()
ax.set_aspect('equal')
central_body = plt.Circle((0, 0), 200, color='gold', zorder=5, label='Central Body')
ax.add_patch(central_body)

orbit_line, = ax.plot([], [], 'b--', alpha=0.5, label='Elliptical Orbit')
circular_orbit_line, = ax.plot(x_circ_path, y_circ_path, 'g--', alpha=0.3, label='Circular Orbit')  # circular orbit path
satellite_ellip, = ax.plot([], [], 'ro', label='Elliptical Satellite')
satellite_circ, = ax.plot([], [], 'go', label='Circular Satellite')

ax.set_xlim(np.min(all_x) - margin, np.max(all_x) + margin)
ax.set_ylim(np.min(all_y) - margin, np.max(all_y) + margin)
ax.set_title("Elliptical and Circular Satellites (Different Periods)")
ax.legend()

def update(frame):
    elapsed_t = frame / fps*100  # simulation time in seconds

    # Find which orbit segment we're in
    for i in range(len(extended_cumulative_times) - 1):
        if extended_cumulative_times[i] <= elapsed_t < extended_cumulative_times[i + 1]:
            seg_idx = i
            break
    else:
        seg_idx = 0
        elapsed_t = elapsed_t % extended_cumulative_times[-1]

    t_start = extended_cumulative_times[seg_idx]
    t_in_segment = elapsed_t - t_start
    orbit_duration = periods[seg_idx] / speed_factor
    transition_duration = total_transition_time

    if seg_idx == len(orbits) - 1 or t_in_segment <= orbit_duration:
        # Normal orbit phase
        local_t = min(t_in_segment, orbit_duration)
        x_e, y_e, _ = orbits[seg_idx]
        T_e = periods[seg_idx] / speed_factor
        frame_float = (local_t / T_e) * (len(x_e) - 1)
        i_low = int(np.floor(frame_float))
        i_high = (i_low + 1) % len(x_e)
        alpha = frame_float - i_low
        x_ellip = (1 - alpha) * x_e[i_low] + alpha * x_e[i_high]
        y_ellip = (1 - alpha) * y_e[i_low] + alpha * y_e[i_high]
        satellite_ellip.set_data(x_ellip, y_ellip)
        orbit_line.set_data(x_e, y_e)
    else:
        # Transition phase between orbits
        alpha = (t_in_segment - orbit_duration) / transition_duration
        r_p_old = perigees[seg_idx]
        r_p_new = perigees[(seg_idx + 1) % len(perigees)]
        r_p_interp = (1 - alpha) * r_p_old + alpha * r_p_new
        elapsed_transition_t = (t_in_segment - orbit_duration)
        x_ellip, y_ellip = compute_position_from_rp(r_p_interp, elapsed_transition_t)
        satellite_ellip.set_data(x_ellip, y_ellip)
        orbit_line.set_data([], [])  # Hide orbit line during transition

    # Circular satellite
    theta_circ = 2 * np.pi * (elapsed_t % T_circ) / T_circ
    x_circ = r_circ * np.cos(theta_circ)
    y_circ = r_circ * np.sin(theta_circ)
    satellite_circ.set_data(x_circ, y_circ)

    return satellite_ellip, orbit_line, satellite_circ, circular_orbit_line

ani = FuncAnimation(fig, update, frames=total_frames, interval=frame_interval, blit=True)
plt.show()
