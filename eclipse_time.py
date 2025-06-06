import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Input
h = 381
inclination = np.radians(53)
solar_longitude_plot = np.radians(126.80)

# Constants
R_e = 6371
obliquity_angle = np.radians(23.439292)
mu = 398600.4
R = R_e + h

def plot_earth(ax, R):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = R * np.outer(np.cos(u), np.sin(v))
    y = R * np.outer(np.sin(u), np.sin(v))
    z = R * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(x, y, z, color='blue', alpha=0.5, edgecolor='none')

def plot_orbit(ax, R, alpha, color='red'):
    theta = np.linspace(0, 2 * np.pi, 100)
    x = R * np.cos(theta)
    y = R * np.sin(theta)
    z = np.zeros_like(theta)

    x_rot = x * np.cos(alpha) + z * np.sin(alpha)
    y_rot = y
    z_rot = -x * np.sin(alpha) + z * np.cos(alpha)

    ax.plot(x_rot, y_rot, z_rot, color=color, linewidth=2)

def plot_rotating_axis(ax, R, alpha):
    z = np.linspace(-R, R, 100)
    x = z * np.tan(alpha)
    y = np.zeros_like(z)
    ax.plot(x, y, z, color='black', linewidth=2)

def plot_shadow(ax, radius, length, angle_z):
    # Create a cylinder along X-axis
    theta = np.linspace(0, 2 * np.pi, 50)
    x = np.linspace(0, length, 2)
    Theta, X = np.meshgrid(theta, x)

    Y = radius * np.cos(Theta)
    Z = radius * np.sin(Theta)

    # Flatten for rotation
    Xf = X.flatten()
    Yf = Y.flatten()
    Zf = Z.flatten()

    # Rotate around Z-axis
    Xr = Xf * np.cos(angle_z) - Yf * np.sin(angle_z)
    Yr = Xf * np.sin(angle_z) + Yf * np.cos(angle_z)
    Zr = Zf

    # Reshape back
    Xr = Xr.reshape(X.shape)
    Yr = Yr.reshape(Y.shape)
    Zr = Zr.reshape(Z.shape)

    ax.plot_surface(Xr, Yr, Zr, color='gray', alpha=0.3, edgecolor='none')

def compute_shadow_time(R, R_e, mu, inclination, obliquity, shadow_length, solar_longitude):
    # Generate orbit points in XY plane
    theta = np.linspace(0, 2 * np.pi, 1000)
    x = R * np.cos(theta)
    y = R * np.sin(theta)
    z = np.zeros_like(theta)

    # Rotate orbit around Y-axis by (obliquity + inclination)
    alpha = obliquity + inclination
    x_rot = x * np.cos(alpha) + z * np.sin(alpha)
    y_rot = y
    z_rot = -x * np.sin(alpha) + z * np.cos(alpha)

    orbit_points = np.vstack((x_rot, y_rot, z_rot)).T

    # Compute Sun direction vector:
    # 1. Start with X-axis
    sun_dir = np.array([1, 0, 0])

    # 2. Rotate around Z (solar longitude)
    Rz_sun = np.array([
        [np.cos(solar_longitude), -np.sin(solar_longitude), 0],
        [np.sin(solar_longitude),  np.cos(solar_longitude), 0],
        [0, 0, 1]
    ])

    # 3. Rotate around Y (obliquity)
    Ry_obliq = np.array([
        [ np.cos(obliquity), 0, np.sin(obliquity)],
        [0, 1, 0],
        [-np.sin(obliquity), 0, np.cos(obliquity)]
    ])

    sun_dir = Rz_sun @ Ry_obliq @ sun_dir
    sun_dir = sun_dir / np.linalg.norm(sun_dir)

    # Project orbit points onto shadow axis
    rel_coords = orbit_points
    proj_lengths = np.dot(rel_coords, sun_dir)
    closest_points_on_axis = np.outer(proj_lengths, sun_dir)
    perp_vectors = rel_coords - closest_points_on_axis
    perp_distances = np.linalg.norm(perp_vectors, axis=1)

    # Check if within shadow cylinder (radius and length)
    inside_shadow = (perp_distances < R_e) & (np.abs(proj_lengths) < shadow_length / 2)

    # Fraction of orbit in shadow
    shadow_fraction = np.sum(inside_shadow) / len(inside_shadow)

    # Orbital period
    T = 2 * np.pi * np.sqrt(R**3 / mu)

    t_shadow = shadow_fraction * T / 2
    return t_shadow, T

# Solar longitudes from 0 to 2π
longitudes = np.linspace(0, 2 * np.pi, 3600)

shadow_times = []

for phi in longitudes:
    t_shadow, T = compute_shadow_time(R, R_e, mu, inclination, obliquity_angle, 2*R, solar_longitude=phi)
    shadow_times.append(t_shadow)
max_shadow_time_longitude = longitudes[np.argmax(shadow_times)]
max_shadow_time = np.max(shadow_times)

print(f"Max shadow time: {max_shadow_time:.2f} seconds")
print(f"Max shadow time occurs at solar longitude: {max_shadow_time_longitude*90:.2f} degrees")
print(T)
plt.plot(longitudes, shadow_times)
plt.xlabel("Solar Longitude (radians)")
plt.ylabel("Eclipse Time (seconds)")
plt.title("Orbit Time in Eclipse vs Solar Longitude")
plt.grid(True)
plt.show()

# Create 3D figure and axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot Earth, orbit, axis, and shadow
plot_earth(ax, R_e)
plot_orbit(ax, R, obliquity_angle + inclination)
plot_orbit(ax, R_e, obliquity_angle, color = "blue")
plot_rotating_axis(ax, R, obliquity_angle)
plot_shadow(ax, R_e, 2 * R, angle_z=solar_longitude_plot)

# Set same limits for all axes
ax.set_xlim([-R, R])
ax.set_ylim([-R, R])
ax.set_zlim([-R, R])
ax.set_box_aspect([1, 1, 1])

# Set plot properties
ax.set_title('Orbit and Earth Shadow')
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.view_init(elev=20, azim=-120)

plt.show()
