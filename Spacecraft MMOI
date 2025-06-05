import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import matplotlib.patches as patches

def tank_volume(H, D, V):
    """
    Equation to solve for total tank height H:
    Volume = cylindrical part + 2 hemispherical ends (1 sphere)
    """
    r = D / 2
    h_cylinder = H - D  # cylinder height = total height - 2*r
    V_calc = np.pi * r**2 * h_cylinder + (4/3) * np.pi * r**3
    return V_calc - V

def calculate_mass(D, H, rho, t):
    """
    Calculate the mass of the tank.
    """
    A = (D/2)**2 * np.pi * 4 + D * np.pi * (H-D)
    return A * t * rho

# ---------------------------------- MMOI CALCULATION ---------------------------------- #

# ---------------------------------- BUS ---------------------------------- #

def calculate_mmoi_bus(w_bus, d_bus, h_bus, t_wall_bus, rho_bus_mat, xyz_cg_bus):

    # Approximate bus as hollow rectangular box
    outer_vol = w_bus * d_bus * h_bus
    inner_vol = (w_bus - 2*t_wall_bus) * (d_bus - 2*t_wall_bus) * (h_bus - 2*t_wall_bus)
    mass_bus = (outer_vol - inner_vol) * rho_bus_mat

    mmoi_bus_total = np.zeros(3)

    Ixx_bus = (1/12) * mass_bus * (d_bus**2 + h_bus**2)
    Iyy_bus = (1/12) * mass_bus * (w_bus**2 + h_bus**2)
    Izz_bus = (1/12) * mass_bus * (w_bus**2 + d_bus**2)

    mmoi_bus_total[0] = Ixx_bus
    mmoi_bus_total[1] = Iyy_bus
    mmoi_bus_total[2] = Izz_bus

    print("Bus Contribution:", mmoi_bus_total, "kgm^2")
    return mmoi_bus_total

# ---------------------------------- SOLAR PANELS ---------------------------------- #

def calculate_mmoi_sp(t_sp, w_sp, d_sp, rho_sp, xyz_cg_bus = np.array([0, 0, 0]), 
                      x_cg_sp = [0, 0], y_cg_sp = [1.345, -1.345], z_cg_sp = [0, 0]):

    vol_sp = t_sp * w_sp * d_sp
    m_sp = vol_sp * rho_sp
    sp_mmoi_total = np.zeros(3)

    for i in range(2):
        cg = np.array([x_cg_sp[i], y_cg_sp[i], z_cg_sp[i]])
        Ixx = (1/12) * m_sp * (d_sp**2 + t_sp**2)
        Iyy = (1/12) * m_sp * (w_sp**2 + t_sp**2)
        Izz = (1/12) * m_sp * (w_sp**2 + d_sp**2)
        
        I_local = np.array([Ixx, Iyy, Izz])
        d = cg - xyz_cg_bus
        d_squared = np.array([d[1]**2 + d[2]**2, d[0]**2 + d[2]**2, d[0]**2 + d[1]**2])
        sp_mmoi_total += I_local + m_sp * d_squared

    print("Solar Panel Contribution:", sp_mmoi_total, "kgm^2")
    return sp_mmoi_total

# ---------------------------------- FUEL TANKS ---------------------------------- #


def calculate_mmoi_tanks(d_bus, w_bus, d_cyl_tank, rho_tank, t_tank, total_fuel_mass, R=2, V=0.02, tank_mass_margin=0.2):      #tank_mass_margin in decimal (e.g. 0.2)
    """
    Calculate the MMOI of 4 symmetric tanks around the bus CG.
    d_bus, w_bus, h_bus: dimensions of bus
    d_cyl_tank: diameter of cylindrical tank
    rho_tank: tank material density
    t_tank: tank wall thickness
    R: length-to-diameter ratio
    V: internal volume of tank (used to solve for height)
    """

    # Solve for total tank height H using tank_volume equation
    H_initial_guess = d_cyl_tank * R
    H_solution = fsolve(lambda H: tank_volume(H, d_cyl_tank, V), H_initial_guess)[0]    

    # Mass of one tank (using total height H_solution)
    mass_empty_tank = calculate_mass(d_cyl_tank, H_solution, rho_tank, t_tank)
    mass_full_tank = mass_empty_tank + total_fuel_mass/4
    mass_full_tank_margin = mass_full_tank * (1 + tank_mass_margin)

    # Geometry of tank
    r = d_cyl_tank / 2
    L_cylinder = H_solution - d_cyl_tank  # exclude hemispherical ends

    # MMOI of one tank about its own CG
    # Cylinder component (assumed axis along z)
    Ixx_cyl = (1/12) * mass_full_tank_margin * (3*r**2 + L_cylinder**2)
    Iyy_cyl = Ixx_cyl
    Izz_cyl = (1/2) * mass_full_tank_margin * r**2

    # Approximate hemispheres as spheres for MMOI simplicity
    # This is conservative (slightly under/over estimates but fine for engineering estimate)

    I_local = np.array([Ixx_cyl, Iyy_cyl, Izz_cyl])

    # Tank CG positions (in bus CG coordinate system)
    x_cg_tanks = [ w_bus/2 - r,  w_bus/2 - r, -w_bus/2 + r, -w_bus/2 + r]
    y_cg_tanks = [ d_bus/2 - r, -d_bus/2 + r,  d_bus/2 - r, -d_bus/2 + r]
    z_cg_tanks = [0, 0, 0, 0]

    total_mmoi_tanks = np.zeros(3)

    for i in range(4):
        d_vec = np.array([x_cg_tanks[i], y_cg_tanks[i], z_cg_tanks[i]])
        d_sq = np.array([
            d_vec[1]**2 + d_vec[2]**2,  # Ixx
            d_vec[0]**2 + d_vec[2]**2,  # Iyy
            d_vec[0]**2 + d_vec[1]**2   # Izz
        ])
        total_mmoi_tanks += I_local + mass_full_tank_margin * d_sq

    print("Fuel Tanks Contribution:", total_mmoi_tanks, "kgm^2")
    return total_mmoi_tanks


def calculate_mmoi_sc_total(w_bus, d_bus, h_bus, t_wall_bus, rho_bus_mat, d_cyl_tank, rho_tank, t_tank, t_sp, w_sp, d_sp, rho_sp, total_fuel_mass,  xyz_cg_bus = np.array([0, 0, 0])):

    total_mmoi_sp = calculate_mmoi_sp(t_sp, w_sp, d_sp, rho_sp, xyz_cg_bus, 
                      x_cg_sp = [0, 0], y_cg_sp = [1.345, -1.345], z_cg_sp = [0, 0])
    
    total_mmoi_bus = calculate_mmoi_bus(w_bus, d_bus, h_bus, t_wall_bus, rho_bus_mat, xyz_cg_bus)

    total_mmoi_tanks = calculate_mmoi_tanks(d_bus, w_bus, d_cyl_tank, rho_tank, t_tank, total_fuel_mass, R=2, V=0.02, tank_mass_margin=0.2)

    mmoi_sc_total = total_mmoi_tanks + total_mmoi_sp + total_mmoi_bus

    print("Total S/C MMOI:", mmoi_sc_total, "kgm^2")

    return None


calculate_mmoi_sc_total(w_bus = 1.13, d_bus = 1.13, h_bus = 1.065, t_wall_bus = 0.004, rho_bus_mat = 4430, d_cyl_tank = 0.533, rho_tank = 4430, t_tank = 0.002, t_sp = 0.02, w_sp = 1.56, d_sp = 0.9, rho_sp = 750, total_fuel_mass = 700, xyz_cg_bus = np.array([0, 0, 0]))

