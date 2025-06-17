import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


import math

def optimize_tank_thickness(P, R, sigma_allowable, corrosion_allowance=0.0, safety_factor=1.0):
    """
    Calculate minimum thicknesses for cylinder wall and hemispherical domed end caps
    for a pressure vessel storing gas nitrogen (GN2).
    
    Parameters:
    - P: Internal pressure in Pascals (Pa)
    - R: Inner radius of the tank in meters (m)
    - sigma_allowable: Allowable stress of the material in Pascals (Pa)
    - corrosion_allowance: Additional thickness for corrosion (m)
    - safety_factor: Safety factor (typically >1)
    
    Returns:
    - dict with thickness for cylinder and domed end caps (meters)
    """

    # Adjust allowable stress by safety factor
    sigma_adj = sigma_allowable / safety_factor

    # Cylinder thickness (thin-wall approx): t = P*R / sigma_allowable
    t_cyl = (P * R) / sigma_adj + corrosion_allowance

    # Hemispherical dome thickness: t = P*R / (2 * sigma_allowable)
    t_dome = (P * R) / (2 * sigma_adj) + corrosion_allowance

    return {
        "Cylinder Wall Thickness (m)": t_cyl,
        "Domed End Cap Thickness (m)": t_dome
    }

if __name__ == "__main__":
    # Example parameters for GN2 tank
    P_bar = 90           # pressure in bar
    P = P_bar * 1e5      # convert bar to Pa
    R = 0.05             # tank radius in meters (5 cm)
    sigma_allowable = 150e6  # allowable stress in Pa (e.g. Aluminum alloy)
    corrosion_allowance = 0.0005  # 0.5 mm corrosion allowance
    safety_factor = 2.0  # typical safety factor

    thicknesses = optimize_tank_thickness(P, R, sigma_allowable, corrosion_allowance, safety_factor)

    print("Optimized tank thicknesses for GN2 at 90 bar:")
    for part, thickness in thicknesses.items():
        print(f"{part}: {thickness*1000:.2f} mm")





import math

def optimize_tank_from_volume(P, V_internal, R, sigma_allowable, density, corrosion_allowance=0.0, safety_factor=1.0):
    """
    Calculate the cylinder length, thicknesses for cylinder and domed ends, and mass of tank walls
    based on total internal volume and radius.

    Parameters:
    - P: Internal pressure in Pascals (Pa)
    - V_internal: Required internal volume in cubic meters (m³)
    - R: Inner radius of the tank in meters (m)
    - sigma_allowable: Allowable stress of the material in Pascals (Pa)
    - density: Material density in kg/m³
    - corrosion_allowance: Additional thickness for corrosion (m)
    - safety_factor: Safety factor (typically >1)

    Returns:
    - dict with cylinder length, thicknesses (meters), tank wall volume (m³), and mass (kg)
    """
    # Volume of the two hemispherical end caps (a full sphere)
    V_ends = (4/3) * math.pi * R**3

    # Volume of cylindrical section needed
    V_cyl = V_internal - V_ends
    if V_cyl < 0:
        raise ValueError("Total volume too small for given radius! Increase volume or reduce radius.")

    # Cylinder length from volume
    L = V_cyl / (math.pi * R**2)

    sigma_adj = sigma_allowable / safety_factor

    # Thickness calculations (thin wall pressure vessel formulas)
    t_cyl = (P * R) / sigma_adj + corrosion_allowance
    t_dome = (P * R) / (2 * sigma_adj) + corrosion_allowance

    # Outer radii
    R_outer_cyl = R + t_cyl
    R_outer_dome = R + t_dome

    # Wall volumes
    V_cyl_shell = math.pi * (R_outer_cyl**2 - R**2) * L
    V_dome_shell = (4/3) * math.pi * (R_outer_dome**3 - R**3)  # volume of sphere shell (2 hemispheres combined)

    V_total_shell = V_cyl_shell + V_dome_shell
    mass = V_total_shell * density

    return {
        "Cylinder Length (m)": L,
        "Cylinder Wall Thickness (m)": t_cyl,
        "Domed End Cap Thickness (m)": t_dome,
        "Tank Wall Volume (m³)": V_total_shell,
        "Tank Wall Mass (kg)": mass
    }

if __name__ == "__main__":
    # Input parameters
    P_bar = 90                  # pressure in bar
    P = P_bar * 1e5             # convert bar to Pa
    V_internal = 0.0025         # desired internal volume in m³ (e.g. 2 liters)
    R = 0.08412                   # inner radius in meters (5 cm)
    sigma_allowable = 395e6     # allowable stress in Pa (e.g., Aluminum alloy)
    density = 2840             # aluminum density in kg/m³
    corrosion_allowance = 0.0005 # 0.5 mm corrosion allowance
    safety_factor = 2        # safety factor

    results = optimize_tank_from_volume(P, V_internal, R, sigma_allowable, density, corrosion_allowance, safety_factor)

    print("Optimized tank design based on total internal volume:")
    print(f"Cylinder Length: {results['Cylinder Length (m)']*1000:.1f} mm")
    print(f"Cylinder Wall Thickness: {results['Cylinder Wall Thickness (m)']*1000:.2f} mm")
    print(f"Domed End Cap Thickness: {results['Domed End Cap Thickness (m)']*1000:.2f} mm")
    print(f"Tank Wall Volume: {results['Tank Wall Volume (m³)']:.6f} m³")
    print(f"Tank Wall Mass: {results['Tank Wall Mass (kg)']:.2f} kg")



# Create figure and axes
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect([1, 1, 1])  # Keep cube aspect ratio

# Shift cube to be centered at (0, 0, 0)
def draw_cube(ax):
    offset = -0.5
    cube_verts = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]
    ]) + offset

    faces = [
        [0, 1, 2, 3], [4, 5, 6, 7],
        [0, 1, 5, 4], [2, 3, 7, 6],
        [1, 2, 6, 5], [0, 3, 7, 4]
    ]
    for face in faces:
        face_coords = cube_verts[face + [face[0]]]  # close the face
        ax.plot(face_coords[:, 0], face_coords[:, 1], face_coords[:, 2], color='gray')

def draw_thruster(ax, label, position, direction, color='purple'):
    direction = np.array(direction)
    direction = direction / np.linalg.norm(direction)
    ax.quiver(*position, *direction, length=0.15, color=color, linewidth=2)
    ax.text(*(np.array(position) + 0.02), label, fontsize=8, color='black')

# Thrusters now based on centered cube corners (±0.5)
thruster_config = {
    "T1": ([0.5, 0.5, 0.5], [1, 1, 1]),
    "T2": ([-0.5, -0.5, -0.5], [-1, -1, -1]),
    "T3": ([0.5, -0.5, 0.5], [1, -1, 1]),
    "T4": ([-0.5, 0.5, -0.5], [-1, 1, -1]),
    "T5": ([0.5, 0.5, -0.5], [1, 1, -1]),
    "T6": ([-0.5, -0.5, 0.5], [-1, -1, 1]),
    "T7": ([0.5, -0.5, -0.5], [1, -1, -1]),
    "T8": ([-0.5, 0.5, 0.5], [-1, 1, 1]),
}

# Plot the spacecraft
draw_cube(ax)
for label, (pos, vec) in thruster_config.items():
    draw_thruster(ax, label, pos, vec)

# Centered axis limits and labels
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
#ax.set_title('8-Thruster Skewed Configuration (Centered at Origin)', pad=20)

plt.tight_layout()
plt.show()




"""

distance = [9,10,11,12,13,14]
fuel = [350.9, 356.0, 363.2, 367.4, 375.3, 381.1]
fuel_onboard = [382, 382, 382, 382, 382, 382]

plt.scatter(distance, fuel, label='Fuel required')
plt.plot(distance, fuel_onboard, color='red', label='Fuel on board')
plt.xlabel('Distance to debris [m]')
plt.ylabel('Total fuel mass [kg]')
plt.grid(True)
plt.legend()
plt.show()



degree = np.arange(0,11,1)
fuel_deg = [350.9, 351.8, 354.4, 357.9, 363.5, 367.5, 372.5, 381.5, 388,388.2, 390.7]
fuel_on_board=[382, 382, 382, 382, 382, 382, 382, 382, 382, 382, 382]
plt.scatter(degree, fuel_deg, label='Fuel required')
plt.plot(degree, fuel_on_board, color='red', label='Fuel on board')
plt.xlabel('Thruster offset [°]')
plt.ylabel('Total fuel mass [kg]')
plt.grid(True)
plt.legend()
plt.show()
"""