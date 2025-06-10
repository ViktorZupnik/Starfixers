import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Parameters
half_cone_angle_deg = 15
half_cone_angle_rad = np.radians(half_cone_angle_deg)
max_distance = 20

# Spacecraft dimensions
sc_width = 3
sc_height = 11
x_positions = [2, 10]

# Create 3D figure
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Create cone surface (plume)
theta = np.linspace(0, 2 * np.pi, 30)
x = np.linspace(0, max_distance, 50)
X, Theta = np.meshgrid(x, theta)
R = X * np.tan(half_cone_angle_rad)
Y = R * np.cos(Theta)
Z = R * np.sin(Theta)
ax.plot_surface(X, Y, Z, color='orange', alpha=0.5, edgecolor='none')

# Draw spacecraft boxes
def draw_spacecraft(ax, x_center, width, height):
    y0 = -height / 2
    y1 = height / 2
    z0 = -width / 2
    z1 = width / 2
    verts = [
        [[x_center, y0, z0], [x_center, y1, z0], [x_center, y1, z1], [x_center, y0, z1]],
        [[x_center, y0, z0], [x_center, y1, z0], [x_center+0.01, y1, z0], [x_center+0.01, y0, z0]],
        [[x_center, y1, z0], [x_center, y1, z1], [x_center+0.01, y1, z1], [x_center+0.01, y1, z0]],
        [[x_center, y0, z1], [x_center, y1, z1], [x_center+0.01, y1, z1], [x_center+0.01, y0, z1]],
        [[x_center, y0, z0], [x_center, y0, z1], [x_center+0.01, y0, z1], [x_center+0.01, y0, z0]],
        [[x_center+0.01, y0, z0], [x_center+0.01, y1, z0], [x_center+0.01, y1, z1], [x_center+0.01, y0, z1]],
    ]
    box = Poly3DCollection(verts, facecolors='gray', edgecolors='black', alpha=1.0)
    ax.add_collection3d(box)

    # Label
    label_y = y0 - 1.0
    ax.text(x_center, label_y, 0,
            'Closest Position' if x_center == x_positions[0] else 'Furthest Position',
            color='black', ha='center')

for x_pos in x_positions:
    draw_spacecraft(ax, x_pos, sc_width, sc_height)

# Title and labels
ax.set_title("3D Exhaust Plume and Spacecraft Positions")
ax.set_xlabel("X (Distance from Nozzle, m)")
ax.set_ylabel("Y (Height, m)")
ax.set_zlabel("Z (Width, m)")

# -------------------------
# Force equal aspect ratio
# -------------------------
def set_axes_equal(ax):
    """Make axes of 3D plot have equal scale so the cone base appears circular."""
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

set_axes_equal(ax)

plt.tight_layout()
plt.show()