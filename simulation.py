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

# === Position computation for given orbital parameters ===
def compute_position(a, e, elapsed_t):
    # Ensure a is positive for a valid orbit and e is within valid range [0, 1) for ellipse
    if a <= 0:
        return np.array([0]), np.array([0]), 0, 0 # Return arrays for x, y
    
    # Cap eccentricity to prevent issues with hyperbolic/parabolic orbits, keep it elliptical
    e = np.clip(e, 0, 0.99999) 

    T = 2 * np.pi * np.sqrt(a**3 / mu)
    
    # Avoid division by zero if T is effectively zero (e.g., very small a)
    if T == 0:
        return np.array([0]), np.array([0]), 0, 0

    # Ensure elapsed_t is within one period to avoid large M values unnecessarily
    elapsed_t = elapsed_t % T 

    M = 2 * np.pi * elapsed_t / T
    E = solve_kepler(np.array([M]), e)[0]
    
    theta = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2),
                             np.sqrt(1 - e) * np.cos(E / 2))
    r = a * (1 - e**2) / (1 + e * np.cos(theta))
    
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return np.array([x]), np.array([y]), theta, T # Return arrays for x, y

# Constants
mu = 398600.4418  # Earth's gravitational parameter [km^3/s^2]

# Settings
num_points_per_orbit_path = 500  # Points to draw the full orbit path
frame_interval = 10          # ms between frames
speed_factor = 2000            # Speed multiplier (higher for faster animation)
transition_frames = 100      # Frames for smooth visual transition

# Define an apogee tolerance for triggering transitions (in radians)
APOGEE_TOLERANCE_ANGLE = np.deg2rad(5) # within 5 degrees of apogee (180 deg true anomaly)

# === Inputs for Elliptical Satellite 1 ===
print("--- Elliptical Satellite 1 ---")
perigees_1_input = input("Enter perigee radii (comma-separated, km) for Satellite 1: ")
perigees_1 = np.array([float(r.strip()) for r in perigees_1_input.split(",")])
apogees_1_input = input("Enter apogee radii (comma-separated, km) for Satellite 1: ")
apogees_1 = np.array([float(r.strip()) for r in apogees_1_input.split(",")])

if len(perigees_1) != len(apogees_1):
    raise ValueError("Number of perigees and apogees for Satellite 1 must be the same.")

# === Inputs for Elliptical Satellite 2 ===
print("\n--- Elliptical Satellite 2 ---")
perigees_2_input = input("Enter perigee radii (comma-separated, km) for Satellite 2: ")
perigees_2 = np.array([float(r.strip()) for r in perigees_2_input.split(",")])
apogees_2_input = input("Enter apogee radii (comma-separated, km) for Satellite 2: ")
apogees_2 = np.array([float(r.strip()) for r in apogees_2_input.split(",")])

if len(perigees_2) != len(apogees_2):
    raise ValueError("Number of perigees and apogees for Satellite 2 must be the same.")

# --- Initial states for dynamic orbital parameters ---
# Satellite 1
current_perigee_1 = perigees_1[0]
current_apogee_1 = apogees_1[0]
current_orbit_idx_1 = 0
time_in_current_orbit_1 = 0.0 # Time counter since the start of the current *defined* orbit
transitioning_1 = False
transition_progress_1 = 0.0
apogee_passage_detected_1 = False # Flag to prevent multiple triggers during one apogee passage

# Satellite 2
current_perigee_2 = perigees_2[0]
current_apogee_2 = apogees_2[0]
current_orbit_idx_2 = 0
time_in_current_orbit_2 = 0.0 # Time counter since the start of the current *defined* orbit
transitioning_2 = False
transition_progress_2 = 0.0
apogee_passage_detected_2 = False # Flag to prevent multiple triggers during one apogee passage

# --- Circular Satellite ---
r_circ = 600 # Example fixed radius for circular satellite
T_circ = 2 * np.pi * np.sqrt(r_circ**3 / mu)

# --- Determine the overall maximum extent for plotting ---
all_possible_radii = np.concatenate([perigees_1, apogees_1, perigees_2, apogees_2, [r_circ]])
max_r = np.max(all_possible_radii)
margin = 0.1 * max_r

# --- Precompute circular orbit path (fixed) ---
theta_circ_path = np.linspace(0, 2*np.pi, num_points_per_orbit_path)
x_circ_path = r_circ * np.cos(theta_circ_path)
y_circ_path = r_circ * np.sin(theta_circ_path)

# --- Define fps and total_frames BEFORE update function ---
fps = 1000 / frame_interval
total_frames = 30000 # Increased frames for longer animation

# === Plot Setup ===
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect('equal')
central_body = plt.Circle((0, 0), 200, color='gold', zorder=5, label='Central Body')
ax.add_patch(central_body)

# Elliptical Satellite 1 elements
orbit_line_1, = ax.plot([], [], 'b--', alpha=0.5, label='Elliptical Orbit 1')
satellite_ellip_1, = ax.plot([], [], 'ro', markersize=8, label='Elliptical Satellite 1')
perigee_marker_1, = ax.plot([], [], 'rx', markersize=10, mew=2, label='Perigee 1') # Perigee marker
apogee_marker_1, = ax.plot([], [], 'r+', markersize=10, mew=2, label='Apogee 1') # Apogee marker

# Elliptical Satellite 2 elements
orbit_line_2, = ax.plot([], [], 'm--', alpha=0.5, label='Elliptical Orbit 2') # Magenta for Sat 2
satellite_ellip_2, = ax.plot([], [], 'co', markersize=8, label='Elliptical Satellite 2') # Cyan for Sat 2
perigee_marker_2, = ax.plot([], [], 'cx', markersize=10, mew=2, label='Perigee 2') # Perigee marker
apogee_marker_2, = ax.plot([], [], 'c+', markersize=10, mew=2, label='Apogee 2') # Apogee marker

# Circular Satellite element
circular_orbit_line, = ax.plot(x_circ_path, y_circ_path, 'g--', alpha=0.3, label='Circular Orbit')
satellite_circ, = ax.plot([], [], 'go', markersize=8, label='Circular Satellite')

ax.set_xlim(-max_r - margin, max_r + margin)
ax.set_ylim(-max_r - margin, max_r + margin)
ax.set_title("Multiple Satellite Orbits (Perigee changes at Apogee)")
ax.legend(loc='upper right')
ax.grid(True, linestyle=':', alpha=0.7)

# Store state variables for easy access in update
state = {
    'current_perigee_1': current_perigee_1,
    'current_apogee_1': current_apogee_1,
    'current_orbit_idx_1': current_orbit_idx_1,
    'time_in_current_orbit_1': time_in_current_orbit_1,
    'transitioning_1': transitioning_1,
    'transition_progress_1': 0.0, # Reset to 0 when transition starts
    'apogee_passage_detected_1': apogee_passage_detected_1,

    'r_p_start_1': 0.0, 'r_a_start_1': 0.0, 'r_p_end_1': 0.0, 'r_a_end_1': 0.0, # For transition interpolation

    'current_perigee_2': current_perigee_2,
    'current_apogee_2': current_apogee_2,
    'current_orbit_idx_2': current_orbit_idx_2,
    'time_in_current_orbit_2': time_in_current_orbit_2,
    'transitioning_2': transitioning_2,
    'transition_progress_2': 0.0, # Reset to 0 when transition starts
    'apogee_passage_detected_2': apogee_passage_detected_2,

    'r_p_start_2': 0.0, 'r_a_start_2': 0.0, 'r_p_end_2': 0.0, 'r_a_end_2': 0.0, # For transition interpolation
}

def update(frame):
    dt = (frame_interval / 1000.0) * speed_factor # Time step in seconds

    # --- Update Elliptical Satellite 1 ---
    a_1_calc = (state['current_perigee_1'] + state['current_apogee_1']) / 2
    e_1_calc = (state['current_apogee_1'] - state['current_perigee_1']) / (state['current_apogee_1'] + state['current_perigee_1'])

    if state['transitioning_1']:
        state['transition_progress_1'] += dt / (transition_frames * (frame_interval / 1000.0))
        alpha = np.clip(state['transition_progress_1'], 0, 1)

        r_p_interp_1 = (1 - alpha) * state['r_p_start_1'] + alpha * state['r_p_end_1']
        r_a_interp_1 = (1 - alpha) * state['r_a_start_1'] + alpha * state['r_a_end_1']
        
        a_1_calc = (r_p_interp_1 + r_a_interp_1) / 2
        e_1_calc = (r_a_interp_1 - r_p_interp_1) / (r_a_interp_1 + r_p_interp_1)

        if alpha >= 1.0: # Transition complete
            state['transitioning_1'] = False
            state['current_orbit_idx_1'] = (state['current_orbit_idx_1'] + 1) % len(perigees_1)
            state['current_perigee_1'] = perigees_1[state['current_orbit_idx_1']]
            state['current_apogee_1'] = apogees_1[state['current_orbit_idx_1']]
            
            a_new_1 = (state['current_perigee_1'] + state['current_apogee_1']) / 2
            e_new_1 = (state['current_apogee_1'] - state['current_perigee_1']) / (state['current_apogee_1'] + state['current_perigee_1'])
            T_new_1 = 2 * np.pi * np.sqrt(a_new_1**3 / mu)
            
            state['time_in_current_orbit_1'] = T_new_1 / 2.0 
            
            print(f"Sat 1 Transition COMPLETE. New Perigee: {state['current_perigee_1']:.2f} km, New Apogee: {state['current_apogee_1']:.2f} km")
            
        # Update the orbit path line and markers using interpolated values during transition
        theta_path_1 = np.linspace(0, 2*np.pi, num_points_per_orbit_path)
        r_path_1 = a_1_calc * (1 - e_1_calc**2) / (1 + e_1_calc * np.cos(theta_path_1))
        x_path_1 = r_path_1 * np.cos(theta_path_1)
        y_path_1 = r_path_1 * np.sin(theta_path_1)
        idx_apogee_path_1 = np.argmin(np.abs(theta_path_1 - np.pi))
        orbit_line_1.set_data(np.roll(x_path_1, -idx_apogee_path_1), np.roll(y_path_1, -idx_apogee_path_1))
        
        perigee_marker_1.set_data(r_p_interp_1, 0)
        apogee_marker_1.set_data(-r_a_interp_1, 0)
        

    x_ellip_1, y_ellip_1, theta_1, T_1 = compute_position(a_1_calc, e_1_calc, state['time_in_current_orbit_1'])
    satellite_ellip_1.set_data(x_ellip_1, y_ellip_1)
    
    if not state['transitioning_1']: # Only trigger if not already in a transition
        if np.abs(theta_1 - np.pi) < APOGEE_TOLERANCE_ANGLE:
            if not state['apogee_passage_detected_1']:
                state['apogee_passage_detected_1'] = True
                
                state['transitioning_1'] = True
                state['transition_progress_1'] = 0.0 # Reset progress for new transition
                
                state['r_p_start_1'] = state['current_perigee_1']
                state['r_a_start_1'] = state['current_apogee_1']
                
                next_orbit_idx_1 = (state['current_orbit_idx_1'] + 1) % len(perigees_1)
                state['r_p_end_1'] = perigees_1[next_orbit_idx_1]
                state['r_a_end_1'] = apogees_1[next_orbit_idx_1]
                
                print(f"Sat 1 Apogee DETECTED. Starting transition from (Rp={state['r_p_start_1']:.2f}, Ra={state['r_a_start_1']:.2f}) to (Rp={state['r_p_end_1']:.2f}, Ra={state['r_a_end_1']:.2f})")

        else: # Reset apogee_passage_detected_1 once moved past apogee
            state['apogee_passage_detected_1'] = False
            
    if not state['transitioning_1']:
        # If not transitioning, update the orbit line and markers to the current *fixed* orbit
        a_1_current_display = (state['current_perigee_1'] + state['current_apogee_1']) / 2
        e_1_current_display = (state['current_apogee_1'] - state['current_perigee_1']) / (state['current_apogee_1'] + state['current_perigee_1'])
        theta_path_1 = np.linspace(0, 2*np.pi, num_points_per_orbit_path)
        r_path_1 = a_1_current_display * (1 - e_1_current_display**2) / (1 + e_1_current_display * np.cos(theta_path_1))
        x_path_1 = r_path_1 * np.cos(theta_path_1)
        y_path_1 = r_path_1 * np.sin(theta_path_1)
        idx_apogee_path_1 = np.argmin(np.abs(theta_path_1 - np.pi))
        orbit_line_1.set_data(np.roll(x_path_1, -idx_apogee_path_1), np.roll(y_path_1, -idx_apogee_path_1))
        
        perigee_marker_1.set_data(state['current_perigee_1'], 0) # Use the stored current_perigee_1
        apogee_marker_1.set_data(-state['current_apogee_1'], 0)   # Use the stored current_apogee_1
    
    state['time_in_current_orbit_1'] += dt


    # --- Update Elliptical Satellite 2 (Mirroring Satellite 1's logic) ---
    a_2_calc = (state['current_perigee_2'] + state['current_apogee_2']) / 2
    e_2_calc = (state['current_apogee_2'] - state['current_perigee_2']) / (state['current_apogee_2'] + state['current_perigee_2'])

    if state['transitioning_2']:
        state['transition_progress_2'] += dt / (transition_frames * (frame_interval / 1000.0))
        alpha = np.clip(state['transition_progress_2'], 0, 1)

        r_p_interp_2 = (1 - alpha) * state['r_p_start_2'] + alpha * state['r_p_end_2']
        r_a_interp_2 = (1 - alpha) * state['r_a_start_2'] + alpha * state['r_a_end_2']
            
        a_2_calc = (r_p_interp_2 + r_a_interp_2) / 2
        e_2_calc = (r_a_interp_2 - r_p_interp_2) / (r_a_interp_2 + r_p_interp_2)

        if alpha >= 1.0:
            state['transitioning_2'] = False
            state['current_orbit_idx_2'] = (state['current_orbit_idx_2'] + 1) % len(perigees_2)
            state['current_perigee_2'] = perigees_2[state['current_orbit_idx_2']]
            state['current_apogee_2'] = apogees_2[state['current_orbit_idx_2']]

            a_new_2 = (state['current_perigee_2'] + state['current_apogee_2']) / 2
            # FIX: Corrected the denominator here
            e_new_2 = (state['current_apogee_2'] - state['current_perigee_2']) / (state['current_apogee_2'] + state['current_perigee_2'])
            T_new_2 = 2 * np.pi * np.sqrt(a_new_2**3 / mu)
            state['time_in_current_orbit_2'] = T_new_2 / 2.0 
            
            print(f"Sat 2 Transition COMPLETE. New Perigee: {state['current_perigee_2']:.2f} km, New Apogee: {state['current_apogee_2']:.2f} km")

        theta_path_2 = np.linspace(0, 2*np.pi, num_points_per_orbit_path)
        r_path_2 = a_2_calc * (1 - e_2_calc**2) / (1 + e_2_calc * np.cos(theta_path_2))
        x_path_2 = r_path_2 * np.cos(theta_path_2)
        y_path_2 = r_path_2 * np.sin(theta_path_2)
        idx_apogee_path_2 = np.argmin(np.abs(theta_path_2 - np.pi))
        orbit_line_2.set_data(np.roll(x_path_2, -idx_apogee_path_2), np.roll(y_path_2, -idx_apogee_path_2))

        perigee_marker_2.set_data(r_p_interp_2, 0)
        apogee_marker_2.set_data(-r_a_interp_2, 0)

    x_ellip_2, y_ellip_2, theta_2, T_2 = compute_position(a_2_calc, e_2_calc, state['time_in_current_orbit_2'])
    satellite_ellip_2.set_data(x_ellip_2, y_ellip_2)

    if not state['transitioning_2']:
        if np.abs(theta_2 - np.pi) < APOGEE_TOLERANCE_ANGLE:
            if not state['apogee_passage_detected_2']:
                state['apogee_passage_detected_2'] = True
                
                state['transitioning_2'] = True
                state['transition_progress_2'] = 0.0
                
                state['r_p_start_2'] = state['current_perigee_2']
                state['r_a_start_2'] = state['current_apogee_2']
                
                next_orbit_idx_2 = (state['current_orbit_idx_2'] + 1) % len(perigees_2)
                state['r_p_end_2'] = perigees_2[next_orbit_idx_2]
                state['r_a_end_2'] = apogees_2[next_orbit_idx_2]
                
                print(f"Sat 2 Apogee DETECTED. Starting transition from (Rp={state['r_p_start_2']:.2f}, Ra={state['r_a_start_2']:.2f}) to (Rp={state['r_p_end_2']:.2f}, Ra={state['r_a_end_2']:.2f})")

        else:
            state['apogee_passage_detected_2'] = False

    if not state['transitioning_2']:
        a_2_current_display = (state['current_perigee_2'] + state['current_apogee_2']) / 2
        # FIX: Corrected the denominator here
        e_2_current_display = (state['current_apogee_2'] - state['current_perigee_2']) / (state['current_apogee_2'] + state['current_perigee_2'])
        theta_path_2 = np.linspace(0, 2*np.pi, num_points_per_orbit_path)
        r_path_2 = a_2_current_display * (1 - e_2_current_display**2) / (1 + e_2_current_display * np.cos(theta_path_2))
        x_path_2 = r_path_2 * np.cos(theta_path_2)
        y_path_2 = r_path_2 * np.sin(theta_path_2)
        idx_apogee_path_2 = np.argmin(np.abs(theta_path_2 - np.pi))
        orbit_line_2.set_data(np.roll(x_path_2, -idx_apogee_path_2), np.roll(y_path_2, -idx_apogee_path_2))

        perigee_marker_2.set_data(state['current_perigee_2'], 0)
        apogee_marker_2.set_data(-state['current_apogee_2'], 0)
    
    state['time_in_current_orbit_2'] += dt

    # --- Update Circular Satellite ---
    elapsed_t_circ = frame / fps 
    theta_circ = 2 * np.pi * (elapsed_t_circ % T_circ) / T_circ
    x_circ = r_circ * np.cos(theta_circ)
    y_circ = r_circ * np.sin(theta_circ)
    satellite_circ.set_data(x_circ, y_circ)

    return (satellite_ellip_1, orbit_line_1, perigee_marker_1, apogee_marker_1,
            satellite_ellip_2, orbit_line_2, perigee_marker_2, apogee_marker_2,
            satellite_circ, circular_orbit_line)

ani = FuncAnimation(fig, update, frames=total_frames, interval=frame_interval, blit=True)
plt.show()