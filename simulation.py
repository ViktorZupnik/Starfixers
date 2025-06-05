import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Known data
md_vals = np.array([250, 500])
h_vals = np.array([550, 560, 570, 580, 590, 600, 610, 620, 630])

hmin_250 = np.array([402, 398, 394, 390, 386, 383, 380, 377, 374])
hmin_500 = np.array([359, 355, 352, 349, 346, 343, 341, 338, 336])

def get_hmin(md_query, h_query):
    if not (250 <= md_query <= 500):
        raise ValueError("md must be between 250 and 500 kg")
    if not (550 <= h_query <= 630):
        raise ValueError("h must be between 550 and 630 km")

    # Interpolation functions in h for both md values
    f250 = interp1d(h_vals, hmin_250, kind='linear')
    f500 = interp1d(h_vals, hmin_500, kind='linear')

    hmin250 = f250(h_query)
    hmin500 = f500(h_query)

    # Linear interpolation over md
    alpha = (md_query - 250) / (500 - 250)
    hmin = hmin250 + alpha * (hmin500 - hmin250)
    return float(hmin)

# 🔁 Grid for plotting
md_range = np.linspace(250, 500, 50)
h_range = np.linspace(550, 630, 50)
MD, H = np.meshgrid(md_range, h_range)

Z = np.vectorize(get_hmin)(MD, H)

# 🎨 Plotting
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

print (get_hmin(500, 570), get_hmin(250, 620))