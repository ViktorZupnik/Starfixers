import numpy as np
g = 9.80665
def section_moment_of_inertia(r, t):
    """Moment of inertia for a thin-walled hollow cylinder"""
    return np.pi * (r ** 4 - (r - t) ** 4) / 4


def bending_stress_at_x(x, total_mass, beam_length, r, t, g_lateral):
    """Max bending stress in a simply supported beam with uniform load"""
    w = (g_lateral * total_mass * g) / beam_length  # distributed load [N/m]
    I = section_moment_of_inertia(r, t)

    # Bending moment at position x (beam center: x = L/2 is max)
    Mx = (w * beam_length / 12) * (6 * x - beam_length) - (w * x ** 2) / 2
    stress = (Mx * r) / I
    return stress  # in Pascals


def max_deflection_simple(total_mass, beam_length, r, t, g_lateral, E):
    """Max deflection at center of simply supported beam under uniform load"""
    w = (g_lateral * total_mass * g) / beam_length
    I = section_moment_of_inertia(r, t)
    delta_max = (5 * w * beam_length ** 4) / (384 * E * I)
    return delta_max


def max_deflection_with_elastic_support(
    total_mass,
    beam_length,
    r,
    t,
    g_lateral,
    E_beam,
    E_support,
    A_support,
    L_support
):
    """Max deflection at mid-span with elastic support at center"""
    g = 9.80665
    w = (g_lateral * total_mass * g) / beam_length  # N/m
    I = np.pi * (r**4 - (r - t)**4) / 4  # Moment of inertia

    # Deflection without center support
    delta_0 = (5 * w * beam_length**4) / (384 * E_beam * I)

    # Stiffness of center support (spring)
    k_s = (E_support * A_support) / L_support  # N/m

    # Reduction factor due to spring
    reduction_factor = 1 + (384 * k_s * beam_length**3) / (np.pi**4 * E_beam * I)
    delta_adjusted = delta_0 / reduction_factor

    return delta_adjusted


# INPUT PARAMETERS
total_mass = 150  # kg
beam_length = 1.065  # m
r = 0.533 / 2  # m
t = 0.001  # m
g_lateral = 3  # g
E =71.7e9  # Pa (Aluminum)
A= np.pi*(0.02**2 - (0.018**2))
L= 1.065-2*5.33/3
# Evaluate center stress and deflections
x_center = beam_length / 2
stress_no_support = bending_stress_at_x(x_center, total_mass, beam_length, r, t, g_lateral)
deflection_no_support = max_deflection_simple(total_mass, beam_length, r, t, g_lateral, E)
deflection_with_support = max_deflection_with_elastic_support(total_mass, beam_length, r, t, g_lateral, E, E, A, L )

print(f"Bending stress at center (no support): {stress_no_support / 1e6} MPa")
print(f"Deflection at center (no support): {deflection_no_support * 1e3} mm")
print(f"Deflection at center (with central support): {deflection_with_support * 1e3} mm")