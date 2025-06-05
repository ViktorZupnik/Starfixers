import numpy as np
import matplotlib.pyplot as plt

E = 71.7e9              # Young's modulus in Pa
nu = 0.33             # Poisson's ratio
t_p = 0.002             # Panel thickness in meters
t_t = 0.002             # Tank thickness in meters
k = 6.0               # Buckling coefficient for clamped edges
SF = 1.5               # Safety Factor
Kd = 0.5                # Knockdown factor
r_outer = 0.533/2       #Tanks outer radius

walls_axial = {
    "Front": (1.13, 1.065),
    "Back": (1.13, 1.065),
    "Left": (1.13, 1.065),
    "Right": (1.13, 1.065)
}

walls_lateral = {
    "Top": (1.13, 1.13),
    "Bottom": (1.13, 1.13),
    "Front": (1.13, 1.065),
    "Back": (1.13, 1.065)
}

def sigma_cr_plate(E, nu, t_p, b, k, SF):

    return (1/SF)*(k * (np.pi**2) * E / (12 * (1 - nu**2))) * (t_p / b)**2

def cylinder_buckling(E, t_t, r_outer, Kd):

    sigma_cr_ideal = 0.605 * E * t_t / (r_outer - t_t/2)
    sigma_cr_real = Kd * sigma_cr_ideal
    return sigma_cr_ideal, sigma_cr_real


ideal, real = cylinder_buckling(E, t_t, r_outer, Kd)

print(f"Ideal buckling stress: {ideal/1e6:.2f} MPa")
print(f"Realistic (with knockdown): {real/1e6:.2f} MPa")


print("Critical Buckling Stress for Cube Walls in axial direction:\n")
for name, (a, b) in walls_axial.items():
    b_eff = min(a, b)  # Buckling depends on shorter dimension
    sigma_cr = sigma_cr_plate(E, nu, t_p, b_eff, k, SF)
    print(f"{name} wall: {sigma_cr/1e6:.2f} MPa")

print("\n\n\nCritical Buckling Stress for Cube Walls in lateral direction:\n")
for name, (a, b) in walls_lateral.items():
    b_eff = min(a, b)  # Buckling depends on shorter dimension
    sigma_cr = sigma_cr_plate(E, nu, t_p, b_eff, k, SF)
    print(f"{name} wall: {sigma_cr / 1e6:.2f} MPa")