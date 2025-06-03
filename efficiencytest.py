import numpy as np
import matplotlib.pyplot as plt

def equal_area_rings(n, radius):
    # Radii for outer edge of each ring
    ring_edges = [radius * np.sqrt(i / n) for i in range(1, n + 1)]
    ring_areas = []

    prev_r = 0.0
    for r in ring_edges:
        area = np.pi * (r**2 - prev_r**2)
        ring_areas.append(area)
        prev_r = r

    return ring_edges, ring_areas

def calculate_efficiency(half_cone_angle, side, distance, samples):
    """
    Calculate the momentum transfer efficiency at a given distance and half cone angle.
    """

    object_area = side ** 2  # Area of the object (circular area)
    cone_base_radius = distance * np.tan(half_cone_angle)
    cone_base_area = np.pi * cone_base_radius ** 2

    ring_edges, ring_areas = equal_area_rings(samples, cone_base_radius)

    for ring_radius in ring_edges:
        if ring_radius > cone_base_radius:
            break

    print(f'Cone base area: {cone_base_area}')
    print(f'Object area: {object_area}')
    print(f'Cone base radius: {cone_base_radius}')
    print("Ring edges:", [round(float(r), 6) for r in ring_edges])
    print("Ring areas:", [round(float(a), 6) for a in ring_areas])

calculate_efficiency(np.radians(15), 0.4, 5, 30)

