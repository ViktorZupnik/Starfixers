import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# === Constants ===
mu = 398600.4418  # Earth's gravitational parameter [km^3/s^2]
R_EARTH = 6371.0 # Earth's mean radius in km

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

# === Visual Scaling Function ===
# Define the physical altitude range of LEO orbits that we want to exaggerate
# These are altitudes ABOVE Earth's surface
PHYSICAL_LEO_MIN_ALTITUDE = 200.0  # Example LEO minimum altitude (km)
PHYSICAL_LEO_MAX_ALTITUDE = 2000.0 # Example LEO maximum altitude (km)

# Define the corresponding visual radius range (from Earth center)
# This range will be much larger to exaggerate the differences
VISUAL_MIN_RADIUS = R_EARTH + 2000 # Visual radius for PHYSICAL_LEO_MIN_ALTITUDE (e.g., 2000km above surface)
VISUAL_MAX_RADIUS = R_EARTH + 15000 # Visual radius for PHYSICAL_LEO_MAX_ALTITUDE (e.g., 15000km above surface)

def get_display_radius(physical_radius):
    """
    Maps a physical radius to a visually exaggerated display radius.
    Designed to highlight small differences in LEO altitudes.
    """
    physical_altitude = physical_radius - R_EARTH

    # Convert physical altitude to a scaled altitude within the defined display range
    if physical_altitude <= PHYSICAL_LEO_MIN_ALTITUDE:
        # For altitudes at or below the min LEO, scale proportionally to the visual min radius
        # Ensure it doesn't go below Earth's visual surface
        if physical_altitude < 0: # Cap at 0 altitude
            scaled_altitude = 0.0
        else:
            scaled_altitude = physical_altitude * ((VISUAL_MIN_RADIUS - R_EARTH) / PHYSICAL_LEO_MIN_ALTITUDE)
    elif physical_altitude >= PHYSICAL_LEO_MAX_ALTITUDE:
        # For altitudes at or above the max LEO, extend linearly from the max display radius
        scaled_altitude = (VISUAL_MAX_RADIUS - R_EARTH) + \
                          (physical_altitude - PHYSICAL_LEO_MAX_ALTITUDE) * \
                          ((VISUAL_MAX_RADIUS - VISUAL_MIN_RADIUS) / (PHYSICAL_LEO_MAX_ALTITUDE - PHYSICAL_LEO_MIN_ALTITUDE))
    else:
        # Perform linear scaling within the target LEO altitude range for exaggeration
        fraction = (physical_altitude - PHYSICAL_LEO_MIN_ALTITUDE) / \
                   (PHYSICAL_LEO_MAX_ALTITUDE - PHYSICAL_LEO_MIN_ALTITUDE)
        scaled_altitude = (VISUAL_MIN_RADIUS - R_EARTH) + fraction * \
                          ((VISUAL_MAX_RADIUS - R_EARTH) - (VISUAL_MIN_RADIUS - R_EARTH))
    
    return scaled_altitude + R_EARTH # Convert back to radius from Earth's center for display


# === Settings ===
num_points_per_orbit_path = 500  # Points to draw the full orbit path
frame_interval = 10          # ms between frames
speed_factor = 20000         # Speed multiplier (higher for faster animation)
transition_frames = 100      # Frames for smooth visual transition

# Define an apogee tolerance for triggering transitions (in radians)
APOGEE_TOLERANCE_ANGLE = np.deg2rad(5) # within 5 degrees of apogee (180 deg true anomaly)
ORBIT_PARAM_TOLERANCE = 1e-6 # Tolerance for checking if a/e have changed enough to redraw orbit path


# === Inputs for Elliptical Satellite 1 ===
print("--- Elliptical Satellite 1 ---")
# Now asking for altitude above surface
perigees_1_input = input("Enter perigee altitudes (comma-separated, km above surface) for Satellite 1: ")
# Convert altitudes to radii by adding R_EARTH
perigees_1_physical = np.array([float(r.strip()) + R_EARTH for r in perigees_1_input.split(",")])

# For simplicity, keeping apogees fixed at a certain altitude above surface
apogees_1_alt = 600.0 # Example fixed apogee altitude above surface
apogees_1_physical = np.array([apogees_1_alt + R_EARTH] * len(perigees_1_physical))

if len(perigees_1_physical) != len(apogees_1_physical):
    raise ValueError("Number of perigees and apogees for Satellite 1 must be the same.")


# === Inputs for Elliptical Satellite 2 ===
print("\n--- Elliptical Satellite 2 ---")
# Now asking for altitude above surface
perigees_2_input = input("Enter perigee altitudes (comma-separated, km above surface) for Satellite 2: ")
perigees_2_physical = np.array([float(r.strip()) + R_EARTH for r in perigees_2_input.split(",")])

apogees_2_input = input("Enter apogee altitudes (comma-separated, km above surface) for Satellite 2: ")
apogees_2_physical = np.array([float(r.strip()) + R_EARTH for r in apogees_2_input.split(",")])

if len(perigees_2_physical) != len(apogees_2_physical):
    raise ValueError("Number of perigees and apogees for Satellite 2 must be the same.")


# --- Initial states for dynamic orbital parameters ---
# Satellite 1
current_perigee_1_physical = perigees_1_physical[0]
current_apogee_1_physical = apogees_1_physical[0]
current_orbit_idx_1 = 0
time_in_current_orbit_1 = 0.0
transitioning_1 = False
transition_progress_1 = 0.0
apogee_passage_detected_1 = False

# Satellite 2
current_perigee_2_physical = perigees_2_physical[0]
current_apogee_2_physical = apogees_2_physical[0]
current_orbit_idx_2 = 0
time_in_current_orbit_2 = 0.0
transitioning_2 = False
transition_progress_2 = 0.0
apogee_passage_detected_2 = False

# --- Circular Satellite ---
r_circ_alt = 600.0 # Example fixed altitude for circular satellite in km
r_circ_physical = r_circ_alt # Physical radius from Earth's center
T_circ = 2 * np.pi * np.sqrt(r_circ_physical**3 / mu)

# --- Determine initial semi-major axes and eccentricities for state tracking ---
a_init_1_physical = (current_perigee_1_physical + current_apogee_1_physical) / 2
e_init_1_physical = (current_apogee_1_physical - current_perigee_1_physical) / (current_perigee_1_physical + current_apogee_1_physical)

a_init_2_physical = (current_perigee_2_physical + current_apogee_2_physical) / 2
e_init_2_physical = (current_apogee_2_physical - current_perigee_2_physical) / (current_perigee_2_physical + current_perigee_2_physical)

# --- Determine the overall maximum extent for plotting based on DISPLAY radii ---
# Calculate display radii for all possible physical radii
all_possible_display_radii = np.array([get_display_radius(r) for r in np.concatenate([perigees_1_physical, apogees_1_physical, perigees_2_physical, apogees_2_physical, [r_circ_physical]])])
max_r_display = np.max(all_possible_display_radii)
margin = 0.1 * max_r_display

# --- Precompute circular orbit path (fixed and visually scaled) ---
r_circ_display = get_display_radius(r_circ_physical)
theta_circ_path = np.linspace(0, 2*np.pi, num_points_per_orbit_path)
x_circ_path = r_circ_display * np.cos(theta_circ_path)
y_circ_path = r_circ_display * np.sin(theta_circ_path)


# --- Define fps and total_frames BEFORE update function ---
fps = 1000 / frame_interval
total_frames = 30000

# === Plot Setup ===
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect('equal')
central_body = plt.Circle((0, 0), R_EARTH, color='gold', zorder=5, label='Earth') # Use R_EARTH for central body
ax.add_patch(central_body)

# Elliptical Satellite 1 elements
orbit_line_1, = ax.plot([], [], 'b--', alpha=0.5, label='Elliptical Orbit 1')
satellite_ellip_1, = ax.plot([], [], 'ro', markersize=8, label='Elliptical Satellite 1')
perigee_marker_1, = ax.plot([], [], 'rx', markersize=10, mew=2, label='Perigee 1')
apogee_marker_1, = ax.plot([], [], 'r+', markersize=10, mew=2, label='Apogee 1')

# Elliptical Satellite 2 elements
orbit_line_2, = ax.plot([], [], 'm--', alpha=0.5, label='Elliptical Orbit 2')
satellite_ellip_2, = ax.plot([], [], 'co', markersize=8, label='Elliptical Satellite 2')
perigee_marker_2, = ax.plot([], [], 'cx', markersize=10, mew=2, label='Perigee 2')
apogee_marker_2, = ax.plot([], [], 'c+', markersize=10, mew=2, label='Apogee 2')

# Circular Satellite element
circular_orbit_line, = ax.plot(x_circ_path, y_circ_path, 'g--', alpha=0.3, label='Circular Orbit')
satellite_circ, = ax.plot([], [], 'go', markersize=8, label='Circular Satellite')

ax.set_xlim(-max_r_display - margin, max_r_display + margin)
ax.set_ylim(-max_r_display - margin, max_r_display + margin)
ax.set_title("Multiple Satellite Orbits (Perigee changes at Apogee with Visual Exaggeration)")
ax.legend(loc='upper right')
ax.grid(True, linestyle=':', alpha=0.7)

# Store state variables for easy access in update
state = {
    'current_perigee_1_physical': current_perigee_1_physical,
    'current_apogee_1_physical': current_apogee_1_physical,
    'current_orbit_idx_1': current_orbit_idx_1,
    'time_in_current_orbit_1': time_in_current_orbit_1,
    'transitioning_1': transitioning_1,
    'transition_progress_1': 0.0,
    'apogee_passage_detected_1': apogee_passage_detected_1,
    'r_p_start_1_physical': 0.0, 'r_a_start_1_physical': 0.0,
    'r_p_end_1_physical': 0.0, 'r_a_end_1_physical': 0.0,
    'last_a_1_display': get_display_radius(a_init_1_physical), # Initialize with first orbit's display a
    'last_e_1_display': e_init_1_physical, # E doesn't change from scaling

    'current_perigee_2_physical': current_perigee_2_physical,
    'current_apogee_2_physical': current_apogee_2_physical,
    'current_orbit_idx_2': current_orbit_idx_2,
    'time_in_current_orbit_2': time_in_current_orbit_2,
    'transitioning_2': transitioning_2,
    'transition_progress_2': 0.0,
    'apogee_passage_detected_2': apogee_passage_detected_2,
    'r_p_start_2_physical': 0.0, 'r_a_start_2_physical': 0.0,
    'r_p_end_2_physical': 0.0, 'r_a_end_2_physical': 0.0,
    'last_a_2_display': get_display_radius(a_init_2_physical), # Initialize with first orbit's display a
    'last_e_2_display': e_init_2_physical, # E doesn't change from scaling
}

def update(frame):
    dt = (frame_interval / 1000.0) * speed_factor

    # --- Update Elliptical Satellite 1 ---
    # Determine physical 'a' and 'e' for calculations
    if state['transitioning_1']:
        state['transition_progress_1'] += dt / (transition_frames * (frame_interval / 1000.0))
        alpha = np.clip(state['transition_progress_1'], 0, 1)

        r_p_interp_1_physical = (1 - alpha) * state['r_p_start_1_physical'] + alpha * state['r_p_end_1_physical']
        r_a_interp_1_physical = (1 - alpha) * state['r_a_start_1_physical'] + alpha * state['r_a_end_1_physical']
        
        a_1_physical = (r_p_interp_1_physical + r_a_interp_1_physical) / 2
        e_1_physical = (r_a_interp_1_physical - r_p_interp_1_physical) / (r_a_interp_1_physical + r_p_interp_1_physical)

        if alpha >= 1.0: # Transition complete
            state['transitioning_1'] = False
            state['current_orbit_idx_1'] = (state['current_orbit_idx_1'] + 1) % len(perigees_1_physical)
            state['current_perigee_1_physical'] = perigees_1_physical[state['current_orbit_idx_1']]
            state['current_apogee_1_physical'] = apogees_1_physical[state['current_orbit_idx_1']]
            
            a_new_1_physical = (state['current_perigee_1_physical'] + state['current_apogee_1_physical']) / 2
            e_new_1_physical = (state['current_apogee_1_physical'] - state['current_perigee_1_physical']) / (state['current_apogee_1_physical'] + state['current_perigee_1_physical'])
            T_new_1 = 2 * np.pi * np.sqrt(a_new_1_physical**3 / mu)
            
            state['time_in_current_orbit_1'] = T_new_1 / 2.0 
            
            print(f"Sat 1 Transition COMPLETE. New Perigee: {state['current_perigee_1_physical']-R_EARTH:.2f} km (alt), New Apogee: {state['current_apogee_1_physical']-R_EARTH:.2f} km (alt)")
            
            # Invalidate display parameters to force redraw of orbit path on next frame
            state['last_a_1_display'] = -1.0 
            state['last_e_1_display'] = -1.0 
    else: # Not transitioning, use current stable physical orbit parameters
        a_1_physical = (state['current_perigee_1_physical'] + state['current_apogee_1_physical']) / 2
        e_1_physical = (state['current_apogee_1_physical'] - state['current_perigee_1_physical']) / (state['current_apogee_1_physical'] + state['current_perigee_1_physical'])

    # Get display parameters from physical ones
    r_p_display_1 = get_display_radius(a_1_physical * (1 - e_1_physical))
    r_a_display_1 = get_display_radius(a_1_physical * (1 + e_1_physical))
    a_1_display = (r_p_display_1 + r_a_display_1) / 2
    e_1_display = (r_a_display_1 - r_p_display_1) / (r_a_display_1 + r_p_display_1) # Eccentricity is scale-independent if distances scale linearly

    # --- Conditional Redraw of Orbit Path 1 ---
    redraw_orbit_1 = False
    if state['transitioning_1'] or \
       abs(a_1_display - state['last_a_1_display']) > ORBIT_PARAM_TOLERANCE or \
       abs(e_1_display - state['last_e_1_display']) > ORBIT_PARAM_TOLERANCE:
        redraw_orbit_1 = True

    if redraw_orbit_1:
        theta_path_1 = np.linspace(0, 2*np.pi, num_points_per_orbit_path)
        r_path_1_display = a_1_display * (1 - e_1_display**2) / (1 + e_1_display * np.cos(theta_path_1))
        
        x_path_1_display = r_path_1_display * np.cos(theta_path_1)
        y_path_1_display = r_path_1_display * np.sin(theta_path_1)
        idx_apogee_path_1 = np.argmin(np.abs(theta_path_1 - np.pi))
        orbit_line_1.set_data(np.roll(x_path_1_display, -idx_apogee_path_1), np.roll(y_path_1_display, -idx_apogee_path_1))
        
        state['last_a_1_display'] = a_1_display
        state['last_e_1_display'] = e_1_display
    
    # Always update perigee/apogee markers based on current calculated orbit display parameters
    perigee_marker_1.set_data(r_p_display_1, 0)
    apogee_marker_1.set_data(-r_a_display_1, 0)

    # Update satellite position using physical parameters for timing, then scale for display
    x_ellip_1_physical, y_ellip_1_physical, theta_1, T_1 = compute_position(a_1_physical, e_1_physical, state['time_in_current_orbit_1'])
    
    # Apply display scaling to the satellite's current position
    # The current position (x_ellip, y_ellip) is based on 'r', which is related to 'a' and 'e'.
    # We need to scale this single point from its physical distance to display distance.
    r_ellip_1_physical = np.sqrt(x_ellip_1_physical**2 + y_ellip_1_physical**2)[0]
    r_ellip_1_display = get_display_radius(r_ellip_1_physical)
    
    # Re-project the position based on the scaled radius and original angle (theta_1)
    x_ellip_1_display = r_ellip_1_display * np.cos(theta_1)
    y_ellip_1_display = r_ellip_1_display * np.sin(theta_1)
    
    satellite_ellip_1.set_data(x_ellip_1_display, y_ellip_1_display)
    
    # Apogee detection logic to trigger the *next* transition (based on physical theta)
    if not state['transitioning_1']:
        if np.abs(theta_1 - np.pi) < APOGEE_TOLERANCE_ANGLE:
            if not state['apogee_passage_detected_1']:
                state['apogee_passage_detected_1'] = True
                
                state['transitioning_1'] = True
                state['transition_progress_1'] = 0.0
                
                state['r_p_start_1_physical'] = state['current_perigee_1_physical']
                state['r_a_start_1_physical'] = state['current_apogee_1_physical']
                
                next_orbit_idx_1 = (state['current_orbit_idx_1'] + 1) % len(perigees_1_physical)
                state['r_p_end_1_physical'] = perigees_1_physical[next_orbit_idx_1]
                state['r_a_end_1_physical'] = apogees_1_physical[next_orbit_idx_1]
                
                print(f"Sat 1 Apogee DETECTED. Starting transition from (Rp={state['r_p_start_1_physical']-R_EARTH:.2f} alt, Ra={state['r_a_start_1_physical']-R_EARTH:.2f} alt) to (Rp={state['r_p_end_1_physical']-R_EARTH:.2f} alt, Ra={state['r_a_end_1_physical']-R_EARTH:.2f} alt)")

        else:
            state['apogee_passage_detected_1'] = False
            
    state['time_in_current_orbit_1'] += dt


    # --- Update Elliptical Satellite 2 ---
    # Determine physical 'a' and 'e' for calculations
    if state['transitioning_2']:
        state['transition_progress_2'] += dt / (transition_frames * (frame_interval / 1000.0))
        alpha = np.clip(state['transition_progress_2'], 0, 1)

        r_p_interp_2_physical = (1 - alpha) * state['r_p_start_2_physical'] + alpha * state['r_p_end_2_physical']
        r_a_interp_2_physical = (1 - alpha) * state['r_a_start_2_physical'] + alpha * state['r_a_end_2_physical']
            
        a_2_physical = (r_p_interp_2_physical + r_a_interp_2_physical) / 2
        e_2_physical = (r_a_interp_2_physical - r_p_interp_2_physical) / (r_a_interp_2_physical + r_p_interp_2_physical)

        if alpha >= 1.0:
            state['transitioning_2'] = False
            state['current_orbit_idx_2'] = (state['current_orbit_idx_2'] + 1) % len(perigees_2_physical)
            state['current_perigee_2_physical'] = perigees_2_physical[state['current_orbit_idx_2']]
            state['current_apogee_2_physical'] = apogees_2_physical[state['current_orbit_idx_2']]

            a_new_2_physical = (state['current_perigee_2_physical'] + state['current_apogee_2_physical']) / 2
            e_new_2_physical = (state['current_apogee_2_physical'] - state['current_perigee_2_physical']) / (state['current_apogee_2_physical'] + state['current_perigee_2_physical'])
            T_new_2 = 2 * np.pi * np.sqrt(a_new_2_physical**3 / mu)
            state['time_in_current_orbit_2'] = T_new_2 / 2.0 
            
            print(f"Sat 2 Transition COMPLETE. New Perigee: {state['current_perigee_2_physical']-R_EARTH:.2f} km (alt), New Apogee: {state['current_apogee_2_physical']-R_EARTH:.2f} km (alt)")

            # Invalidate display parameters to force redraw of orbit path on next frame
            state['last_a_2_display'] = -1.0
            state['last_e_2_display'] = -1.0

    else: # Not transitioning, use current stable physical orbit parameters
        a_2_physical = (state['current_perigee_2_physical'] + state['current_apogee_2_physical']) / 2
        e_2_physical = (state['current_apogee_2_physical'] - state['current_perigee_2_physical']) / (state['current_apogee_2_physical'] + state['current_perigee_2_physical'])

    # Get display parameters from physical ones
    r_p_display_2 = get_display_radius(a_2_physical * (1 - e_2_physical))
    r_a_display_2 = get_display_radius(a_2_physical * (1 + e_2_physical))
    a_2_display = (r_p_display_2 + r_a_display_2) / 2
    e_2_display = (r_a_display_2 - r_p_display_2) / (r_a_display_2 + r_p_display_2)

    # --- Conditional Redraw of Orbit Path 2 ---
    redraw_orbit_2 = False
    if state['transitioning_2'] or \
       abs(a_2_display - state['last_a_2_display']) > ORBIT_PARAM_TOLERANCE or \
       abs(e_2_display - state['last_e_2_display']) > ORBIT_PARAM_TOLERANCE:
        redraw_orbit_2 = True

    if redraw_orbit_2:
        theta_path_2 = np.linspace(0, 2*np.pi, num_points_per_orbit_path)
        r_path_2_display = a_2_display * (1 - e_2_display**2) / (1 + e_2_display * np.cos(theta_path_2))
        x_path_2_display = r_path_2_display * np.cos(theta_path_2)
        y_path_2_display = r_path_2_display * np.sin(theta_path_2)
        idx_apogee_path_2 = np.argmin(np.abs(theta_path_2 - np.pi))
        orbit_line_2.set_data(np.roll(x_path_2_display, -idx_apogee_path_2), np.roll(y_path_2_display, -idx_apogee_path_2))

        state['last_a_2_display'] = a_2_display
        state['last_e_2_display'] = e_2_display

    # Always update perigee/apogee markers based on current calculated orbit display parameters
    perigee_marker_2.set_data(r_p_display_2, 0)
    apogee_marker_2.set_data(-r_a_display_2, 0)

    # Update satellite position using physical parameters for timing, then scale for display
    x_ellip_2_physical, y_ellip_2_physical, theta_2, T_2 = compute_position(a_2_physical, e_2_physical, state['time_in_current_orbit_2'])
    
    # Apply display scaling to the satellite's current position
    r_ellip_2_physical = np.sqrt(x_ellip_2_physical**2 + y_ellip_2_physical**2)[0]
    r_ellip_2_display = get_display_radius(r_ellip_2_physical)
    
    x_ellip_2_display = r_ellip_2_display * np.cos(theta_2)
    y_ellip_2_display = r_ellip_2_display * np.sin(theta_2)

    satellite_ellip_2.set_data(x_ellip_2_display, y_ellip_2_display)
    
    # Apogee detection logic to trigger the *next* transition (based on physical theta)
    if not state['transitioning_2']:
        if np.abs(theta_2 - np.pi) < APOGEE_TOLERANCE_ANGLE:
            if not state['apogee_passage_detected_2']:
                state['apogee_passage_detected_2'] = True
                
                state['transitioning_2'] = True
                state['transition_progress_2'] = 0.0
                
                state['r_p_start_2_physical'] = state['current_perigee_2_physical']
                state['r_a_start_2_physical'] = state['current_apogee_2_physical']
                
                next_orbit_idx_2 = (state['current_orbit_idx_2'] + 1) % len(perigees_2_physical)
                state['r_p_end_2_physical'] = perigees_2_physical[next_orbit_idx_2]
                state['r_a_end_2_physical'] = apogees_2_physical[next_orbit_idx_2]
                
                print(f"Sat 2 Apogee DETECTED. Starting transition from (Rp={state['r_p_start_2_physical']-R_EARTH:.2f} alt, Ra={state['r_a_start_2_physical']-R_EARTH:.2f} alt) to (Rp={state['r_a_end_2_physical']-R_EARTH:.2f} alt, Ra={state['r_a_end_2_physical']-R_EARTH:.2f} alt)")

        else:
            state['apogee_passage_detected_2'] = False
    
    state['time_in_current_orbit_2'] += dt

    # --- Update Circular Satellite ---
    elapsed_t_circ = frame / fps 
    theta_circ = 2 * np.pi * (elapsed_t_circ % T_circ) / T_circ
    
    # Calculate position based on the scaled circular radius for display
    x_circ_display = r_circ_display * np.cos(theta_circ)
    y_circ_display = r_circ_display * np.sin(theta_circ)
    satellite_circ.set_data(x_circ_display, y_circ_display)

    return (satellite_ellip_1, orbit_line_1, perigee_marker_1, apogee_marker_1,
            satellite_ellip_2, orbit_line_2, perigee_marker_2, apogee_marker_2,
            satellite_circ, circular_orbit_line)

ani = FuncAnimation(fig, update, frames=total_frames, interval=frame_interval, blit=True)
plt.show()