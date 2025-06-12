import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, RegularGridInterpolator

# Known data
md_vals = np.array([250, 312.5, 375, 437.5, 500])
h_vals = np.array([550, 560, 570, 580, 590, 600, 610, 620, 630])

# hmin data for each md (each row corresponds to a mass, columns to altitudes)
hmin_matrix = np.array([
    [402, 398, 394, 390, 386, 383, 380, 377, 374],  # 250 kg
    [388, 383, 379, 376, 373, 370, 367, 364, 361],  # 312.5 kg
    [376, 372, 369, 365, 362, 359, 356, 353, 351],  # 375 kg
    [367, 363, 360, 356, 353, 350, 348, 345, 343],  # 437.5 kg
    [359, 355, 352, 349, 346, 343, 341, 338, 336]   # 500 kg
])

# Interpolator function
interp_func = RegularGridInterpolator(
    (md_vals, h_vals),
    hmin_matrix,
    method='linear',
    bounds_error=True
)

# Function to evaluate hmin
def get_hmin(md_query, h_query):
    if not (250 <= md_query <= 500):
        raise ValueError("md must be between 250 and 500 kg")
    if not (550 <= h_query <= 630):
        raise ValueError("h must be between 550 and 630 km")
    
    return float(interp_func([[md_query, h_query]]))

# Grid for plotting
md_range = np.linspace(250, 500, 50)
h_range = np.linspace(550, 630, 50)
MD, H = np.meshgrid(md_range, h_range)

Z = np.vectorize(get_hmin)(MD, H)

# Plotting
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(MD, H, Z, cmap='viridis')

ax.set_xlabel('md [kg]')
ax.set_ylabel('h [km]')
ax.set_zlabel('hmin [km]')
ax.set_title('Interpolated hmin(md, h)')
fig.colorbar(surf, shrink=0.5, aspect=10)

plt.tight_layout()
plt.show()

# ✅ Example use
print(get_hmin(260, 600))
