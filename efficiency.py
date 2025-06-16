import numpy as np
import matplotlib.pyplot as plt

def equal_radius_rings(n, radius):
    # Radii for outer edge of each ring (equal radial spacing)
    ring_edges = [radius * i / n for i in range(1, n + 1)]
    ring_areas = []

    prev_r = 0.0
    for r in ring_edges:
        area = np.pi * (r**2 - prev_r**2)
        ring_areas.append(area)
        prev_r = r

    return ring_edges, ring_areas

def gaussian_weight(r, width): # https://descanso.jpl.nasa.gov/SciTechBook/series1/Goebel_08_Chap8_plumes.pdf?utm_source=chatgpt.com p403
    """Unnormalized Gaussian weight based on radial distance."""
    return np.exp(-(r / width)**2)

def calculate_efficiency_circle(half_cone_angle, object_radius, distance, samples, diffuse):
    """
    Calculate the momentum transfer efficiency at a given distance and half cone angle.
    """
    cone_base_radius = distance * np.tan(half_cone_angle)

    ring_edges, ring_areas = equal_radius_rings(samples, cone_base_radius)

    gauss_weights = np.array([gaussian_weight(r, cone_base_radius / 2) for r in ring_edges])
    ring_masses = gauss_weights * ring_areas
    total_gauss_mass = np.sum(ring_masses)
    
    efficiency = 0
    for i, r in enumerate(ring_edges):
        if r <= object_radius:
            factor1 = np.cos(np.arctan(r / distance))  # Angular efficiency
            factor2 = ring_masses[i] / total_gauss_mass  # Factor to account for the middle of the exhaust having more mass density, not uniform
            factor3 = diffuse # Factor to account for diffusivity of the exhaust plume
            factor = factor1 * factor2 * factor3
            efficiency += factor
    return efficiency

def calculate_efficiency_rect(half_cone_angle, object_width, object_length, distance, resolution, diffuse):
    res = resolution  # grid resolution
    cone_base_radius = distance * np.tan(half_cone_angle)
    extent = cone_base_radius
    x = np.linspace(-extent, extent, res)
    y = np.linspace(-extent, extent, res)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)  # radial distance from plume center
    object_mask = (np.abs(X) < object_width/2) & (np.abs(Y) < object_length/2)
    gauss = gaussian_weight(R, cone_base_radius/2)  # Plume width standard deviation
    angles = np.arctan2(R, distance)
    cosines = np.cos(angles)
    dA = (2 * extent / res) ** 2  # Area of each grid cell
    efficiency_map = gauss * cosines * object_mask * dA
    efficiency = np.sum(efficiency_map) / np.sum(gauss * dA)
    efficiency *= diffuse  # Apply diffusivity
    return efficiency

def calculate_efficiency(half_cone_angle, distance, resolution, diffuse, pointing_error_angle, object_radius=None, object_width=None, object_length=None):
    """
    Calculate the momentum transfer efficiency at a given distance and half cone angle.
    Supports both circular and rectangular objects.
    """
    distance = abs(distance)
    if object_radius is not None:
        return calculate_efficiency_circle(half_cone_angle, object_radius, distance, resolution, diffuse)
    elif object_width is not None and object_length is not None:
        return pointing_efficiency(half_cone_angle, object_width, object_length, distance, resolution, diffuse, pointing_error_angle)
        #return calculate_efficiency_rect(half_cone_angle, object_width, object_length, distance, resolution, diffuse)
    else:
        raise ValueError("Either object_radius or (object_width and object_length) must be provided.")


def pointing_efficiency(half_cone_angle, object_width, object_length, distance, resolution, diffuse, pointing_error_angle):
    """
    Calculate efficiency using overlap area between the offset plume and a rectangular target.
    Plume is offset by pointing_error_angle (in radians) in the x-direction.
    """
    cone_base_radius = distance * np.tan(half_cone_angle)
    extent = cone_base_radius * 1.5  # domain size, large enough to contain shifted plume
    x = np.linspace(-extent, extent, resolution)
    y = np.linspace(-extent, extent, resolution)
    X, Y = np.meshgrid(x, y)

    # Calculate plume center offset due to pointing error
    offset_x = distance * np.tan(pointing_error_angle)
    plume_center_x = 0 + offset_x
    plume_center_y = 0  # could also offset in Y if desired

    # Radial distance from *shifted* plume center
    R_shifted = np.sqrt((X - plume_center_x)**2 + (Y - plume_center_y)**2)

    # Gaussian distribution and cosine angle correction
    gauss = gaussian_weight(R_shifted, cone_base_radius / 2)
    angles = np.arctan2(R_shifted, distance)
    cosines = np.cos(angles)
    dA = (2 * extent / resolution) ** 2

    # Target mask: centered at origin
    target_mask = (np.abs(X) < object_width / 2) & (np.abs(Y) < object_length / 2)

    # Momentum efficiency map
    efficiency_map = gauss * cosines * target_mask * dA

    # Normalize by total available plume momentum
    total_plume_mass = np.sum(gauss * dA)
    efficiency = np.sum(efficiency_map) / total_plume_mass

    return efficiency * diffuse



# CALCULATE USING THRUST PLUME OFFSET
diffuse = 1
half_cone_angle = np.radians(15)
distances = np.arange(0.05, 10, 0.05)

# Rectangular object
object_width = 3.1  # Width of the object in meters
object_length = 10.9  # Length of the object in meters
resolution = 1000

eff0 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(0)) for d in distances]
eff01 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(0.1)) for d in distances]
eff05 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(0.5)) for d in distances]
eff1 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(1)) for d in distances]
eff2 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(2)) for d in distances]
eff3 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(3)) for d in distances]
eff4 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(4)) for d in distances]
eff5 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(5)) for d in distances]
eff6 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(6)) for d in distances]
eff7 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(7)) for d in distances]
eff8 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(8)) for d in distances]
eff9 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(9)) for d in distances]
eff10 = [pointing_efficiency(half_cone_angle, object_width, object_length, d, resolution, diffuse, np.radians(10)) for d in distances]

plt.plot(distances, eff0, label='PO 0°')
#plt.plot(distances, eff01, label='Offset 0.1°')
#plt.plot(distances, eff05, label='PO 0.5°')
plt.plot(distances, eff1, label='PO 1°')
plt.plot(distances, eff2, label='PO 2°')
plt.plot(distances, eff3, label='PO 3°')
plt.plot(distances, eff4, label='PO 4°')
plt.plot(distances, eff5, label='PO 5°')
plt.plot(distances, eff6, label='PO 6°')
plt.plot(distances, eff7, label='PO 7°')
plt.plot(distances, eff8, label='PO 8°')
plt.plot(distances, eff9, label='PO 9°')
plt.plot(distances, eff10, label='PO 10°')
plt.xlabel('Distance (m)')
plt.ylabel('Momentum Transfer Efficiency')
#plt.title('Momentum Transfer Efficiency vs Distance')
plt.grid(True)
plt.legend()
plt.show()



"""
if __name__ == "__main__":
    diffuse = 0.5  # Diffusivity factor, can be adjusted based on exhaust characteristics, 1 for specular reflection, 0 for diffuse
    half_cone_angle = np.radians(15)  # Half cone angle in radians
    distances = np.arange(0.05, 10, 0.05)

    # Circular object with radius 0.5 m
    object_radius = 0.5  # Radius of the object in meters
    efficiencies1 = [calculate_efficiency_circle(half_cone_angle, object_radius, d, 1000, diffuse) for d in distances]
    efficiencies2 = [calculate_efficiency_circle(np.radians(40), object_radius, d, 1000, diffuse) for d in distances]
    efficiencies3 = [calculate_efficiency_circle(np.radians(40), object_radius, d, 1000, 1) for d in distances]
    efficiencies4 = [calculate_efficiency_circle(half_cone_angle, object_radius, d, 1000, 1) for d in distances]

    plt.plot(distances, efficiencies1, label='AC 0.5; HCA 15°')
    plt.plot(distances, efficiencies2, label='AC 0.5; HCA 40°')
    plt.plot(distances, efficiencies3, label='AC 1; HCA 40°')
    plt.plot(distances, efficiencies4, label='AC 1; HCA 15°')
    plt.xlabel('Distance (m)')
    plt.ylabel('Momentum Transfer Efficiency')
    plt.title('Momentum Transfer Efficiency vs Distance')
    plt.grid(True)
    plt.legend()
    plt.show()


    # Rectangular object
    object_width = 3.1  # Width of the object in meters
    object_length = 10.9  # Length of the object in meters
    efficiencies5 = [calculate_efficiency_rect(half_cone_angle, object_width, object_length, d, 1000, diffuse) for d in distances]
    efficiencies6 = [calculate_efficiency_rect(np.radians(40), object_width, object_length, d, 1000, diffuse) for d in distances]
    efficiencies7 = [calculate_efficiency_rect(np.radians(40), object_width, object_length, d, 1000, 1) for d in distances]
    efficiencies8 = [calculate_efficiency_rect(half_cone_angle, object_width, object_length, d, 1000, 1) for d in distances]

    plt.plot(distances, efficiencies5, label='AC 0.5; HCA 15°')
    plt.plot(distances, efficiencies6, label='AC 0.5; HCA 40°')
    plt.plot(distances, efficiencies7, label='AC 1; HCA 40°')
    plt.plot(distances, efficiencies8, label='AC 1; HCA 15°')
    plt.xlabel('Distance (m)')
    plt.ylabel('Momentum Transfer Efficiency')
    plt.title('Momentum Transfer Efficiency vs Distance (Rectangular Object)')
    plt.grid(True)
    plt.legend()
    plt.show()

"""
