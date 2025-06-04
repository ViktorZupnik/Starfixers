import matplotlib.pyplot as plt
import numpy as np

# User-defined parameters
half_cone_angle_deg = 15  # in degrees
max_distance = 20         # in meters

# Convert angle to radians
half_cone_angle_rad = np.radians(half_cone_angle_deg)

# Calculate plume radius at max distance
plume_radius = max_distance * np.tan(half_cone_angle_rad)

# Triangle coordinates (plume edges)
x = [0, max_distance, max_distance]
y = [0, -plume_radius, plume_radius]

# Plot
plt.figure(figsize=(8, 6))
plt.fill(x, y, color='orange', alpha=0.7, label=f'{half_cone_angle_deg}° half cone angle')
plt.plot([0, max_distance], [0, -plume_radius], 'k--')
plt.plot([0, max_distance], [0, plume_radius], 'k--')

# Labels and styling
plt.title("Exhaust Plume Profile")
plt.xlabel("Distance from Nozzle (m)")
plt.ylabel("Plume Radius (m)")
plt.grid(True)
plt.legend()

# Set equal aspect ratio
plt.gca().set_aspect('equal', adjustable='box')

plt.show()