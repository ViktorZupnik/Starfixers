import numpy as np
import matplotlib.pyplot as plt


# Constants
R_e = 6371  # Earth radius in kilometers
h = 600  # Altitude in kilometers of start of transfer
R = R_e + h  # Radius of the orbit in kilometers
mu = 398600  # Gravitational parameter of Earth in km^3/s^2
g0 = 9.81  # Gravitational acceleration at sea level in m/s^2

# Input parameters
DV = -0.002  # Delta V in km/s applied to the spacecraft
theta_iterate = np.arange(-1, 1, 0.0001) # Range of angles in degrees to optimize contact time for
results = []

for theta_deg in theta_iterate:
    theta_0 = np.radians(theta_deg)  # Convert to radians
    v0 = np.sqrt(mu / R)  # Orbital velocity in km/s
    v1 = v0 + DV  # Velocity after delta V in km/s

    a = 1 / ((2 / R) - (v1**2 / mu))
    perigee = 2*a - R  # Perigee in kilometers
    apogee = R  # Apogee in kilometers
    v_perigee = np.sqrt(mu * (2 / perigee - 1 / a))  # Velocity at perigee in km/s

    theta = np.linspace(0, 2 * np.pi, 7200)  # Angle in radians
    theta_s = np.linspace(0, 2 * np.pi, 7200)  # Angle in radians

    # Earth coordinates
    x_e = R_e * np.cos(theta)
    y_e = R_e * np.sin(theta)

    # Target coordinates
    x_t = R * np.cos(theta)
    y_t = R * np.sin(theta)
    T_t = 2*np.pi * np.sqrt(R**3 / mu)  # Orbital period of target in seconds

    # Spacecraft coordinates
    T_s = 2*np.pi * np.sqrt(a**3 / mu)  # Orbital period of sc in seconds
    e = (apogee - perigee) / (apogee + perigee)
    b = a * np.sqrt(1 - e**2)
    x_s = a * np.cos(theta) - a * e  # shift so that focus is at origin
    y_s = b * np.sin(theta)



    # Time simulation setup
    dt = 1  # time step in seconds
    t_max = int(T_t)  # simulate up to 2 target orbits
    t = np.arange(0, t_max, dt)
    # Mean motions (angular rates)
    omega_t = 2 * np.pi / T_t  # rad/s for target
    omega_s = 2 * np.pi / T_s  # rad/s for spacecraft

    # Position arrays
    x_target = R * np.cos(omega_t * t)
    y_target = R * np.sin(omega_t * t)

    x_spacecraft = a * np.cos(omega_s * t + theta_0) + a * e
    y_spacecraft = b * np.sin(omega_s * t + theta_0)

    # Distance between spacecraft and target
    dist = np.sqrt((x_target - x_spacecraft)**2 + (y_target - y_spacecraft)**2)

    # Find when they're within 2 km
    close_indices = np.where(dist <= 2)[0]
    if len(close_indices) == 0:
        continue  # No contact time found
    time_of_crossing = t[close_indices]
    contact_time = time_of_crossing[-1]-time_of_crossing[0]
    if dist[-1] < 0.01:
        results.append((theta_deg, contact_time, dist[-1]))


results = sorted(results, key=lambda x: x[1], reverse=True)
print("Best contact time found:")
for theta, contact_time, final_dist in results:
    print(f"Theta: {theta} rad, Contact Time: {contact_time:} s, Final Distance: {final_dist} km")


def plot_distance(t, dist):
    plt.figure(figsize=(10, 8))
    plt.plot(t, dist, label='Distance between spacecraft and target')
    plt.axhline(y=2, color='r', linestyle='--', label='2 km threshold')
    plt.title('Distance between Spacecraft and Target Over Time')
    plt.show()


def plot_orbits(x_t, y_t, x_s, y_s, x_e, y_e):
    plt.figure(figsize=(8,6))
    ax = plt.gca()
    plt.plot(x_e, y_e, label='Earth', color='blue')
    earth = plt.Circle((0, 0), R_e, color='blue', zorder=0)  # zorder=0 to be behind
    ax.add_patch(earth)
    plt.plot(x_t, y_t, label='Target', color='red')
    plt.plot(x_s, y_s, label='Spacecraft', color='green')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel("X (km)")
    plt.ylabel("Y (km)")
    plt.legend()
    plt.show()

plot_orbits(x_t, y_t, x_s, y_s, x_e, y_e)