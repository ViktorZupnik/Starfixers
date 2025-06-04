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

def gaussian_weight(r, width):
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

def calculate_efficiency(half_cone_angle, distance, resolution, diffuse, object_radius=None, object_width=None, object_length=None):
    """
    Calculate the momentum transfer efficiency at a given distance and half cone angle.
    Supports both circular and rectangular objects.
    """
    if object_radius is not None:
        return calculate_efficiency_circle(half_cone_angle, object_radius, distance, resolution, diffuse)
    elif object_width is not None and object_length is not None:
        return calculate_efficiency_rect(half_cone_angle, object_width, object_length, distance, resolution, diffuse)
    else:
        raise ValueError("Either object_radius or (object_width and object_length) must be provided.")

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

    plt.plot(distances, efficiencies1, label='diffuse; 15 degrees')
    plt.plot(distances, efficiencies2, label='diffuse; 40 degrees')
    plt.plot(distances, efficiencies3, label='specular; 40 degrees')
    plt.plot(distances, efficiencies4, label='specular; 15 degrees')
    plt.xlabel('Distance (m)')
    plt.ylabel('Momentum Transfer Efficiency')
    plt.title('Momentum Transfer Efficiency vs Distance')
    plt.grid(True)
    plt.legend()
    plt.show()


    # Rectangular object with width 0.5 m and length 1 m
    object_width = 3  # Width of the object in meters
    object_length = 11  # Length of the object in meters
    efficiencies5 = [calculate_efficiency_rect(half_cone_angle, object_width, object_length, d, 1000, diffuse) for d in distances]
    efficiencies6 = [calculate_efficiency_rect(np.radians(40), object_width, object_length, d, 1000, diffuse) for d in distances]
    efficiencies7 = [calculate_efficiency_rect(np.radians(40), object_width, object_length, d, 1000, 1) for d in distances]
    efficiencies8 = [calculate_efficiency_rect(half_cone_angle, object_width, object_length, d, 1000, 1) for d in distances]

    plt.plot(distances, efficiencies5, label='diffuse; 15 degrees; rect')
    plt.plot(distances, efficiencies6, label='diffuse; 40 degrees; rect')
    plt.plot(distances, efficiencies7, label='specular; 40 degrees; rect')
    plt.plot(distances, efficiencies8, label='specular; 15 degrees; rect')
    plt.xlabel('Distance (m)')
    plt.ylabel('Momentum Transfer Efficiency')
    plt.title('Momentum Transfer Efficiency vs Distance (Rectangular Object)')
    plt.grid(True)
    plt.legend()
    plt.show()


