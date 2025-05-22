import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import matplotlib.patches as patches


def calculate_tank_volume(tank_capacity, plot=False):
    """
    Calculate the tank volume based on the tank capacity using linear regression in L.
    """
    # data points for tank capacity and corresponding tank volume
    x = np.array([218, 180, 235, 4.5, 15.4, 30, 75, 154])
    y = np.array([267.6935612, 219.0719684, 263.0456283, 6, 17.5, 37.3, 102.5, 204])
    

    # Linear regression
    coeffs = np.polyfit(x, y, 1)  # Degree 1 for linear
    linear_poly = np.poly1d(coeffs)

    # Generate fitted line
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = linear_poly(x_fit)

    # Format the equation string
    eq_text = f"y = {coeffs[0]:.2f}x + {coeffs[1]:.2f}"
    if plot:
        # Plot data and linear fit
        plt.scatter(x, y, color='blue', label='Data')
        plt.plot(x_fit, y_fit, color='red', label='Linear Fit')

        # Add equation to the plot
        plt.text(0.05, 0.95, eq_text, transform=plt.gca().transAxes,
                fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        plt.xlabel('Tank capacity')
        plt.ylabel('Tank volume')
        plt.title('Linear Regression')
        plt.legend()
        plt.grid(True)
        plt.show()
    return linear_poly(tank_capacity)

def calculate_tank_mass(tank_capacity, plot=False):
    x = np.array([218, 180, 235, 4.5, 15.4, 30, 75, 154])
    y = np.array([11, 21, 16, 1.3, 3.1, 3.9, 7.5, 17.5])
    # Linear regression
    coeffs = np.polyfit(x, y, 1)  # Degree 1 for linear
    linear_poly = np.poly1d(coeffs)

    # Generate fitted line
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = linear_poly(x_fit)

    # Format the equation string
    eq_text = f"y = {coeffs[0]:.2f}x + {coeffs[1]:.2f}"
    if plot:
        # Plot data and linear fit
        plt.scatter(x, y, color='blue', label='Data')
        plt.plot(x_fit, y_fit, color='red', label='Linear Fit')

        # Add equation to the plot
        plt.text(0.05, 0.95, eq_text, transform=plt.gca().transAxes,
                fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        plt.xlabel('Tank capacity')
        plt.ylabel('Tank Mass')
        plt.title('Linear Regression')
        plt.legend()
        plt.grid(True)
        plt.show()
    return linear_poly(tank_capacity)

def tank_volume(D, V, R):
    """
    Equation to solve for a given diameter D:
    Volume = cylinder + sphere
    """
    H = D / R
    r = D / 2
    h_cylinder = H - D  # height minus the two hemispherical ends
    V_calc = np.pi * r**2 * h_cylinder + (4/3) * np.pi * r**3
    return V_calc - V

def calculate_tank_dimensions(volume, ratio, D_guess=1.0):
    """
    Calculates diameter and height from volume and D/H ratio
    """
    D_solution, = fsolve(tank_volume, D_guess, args=(volume, ratio))
    H_solution = D_solution / ratio
    return D_solution, H_solution

def plot_tank(D_m, H_m):
    # Convert meters to millimeters
    D = D_m * 1000
    H = H_m * 1000

    r = D / 2
    h_cylinder = H - D  # height of the cylindrical section

    # Shift value to center tank vertically
    y_shift = - (H-D) / 2

    # Cylinder rectangle (centered on x=0, vertically shifted)
    cyl_x = [-r, r, r, -r, -r]
    cyl_y = [0, 0, h_cylinder, h_cylinder, 0]
    cyl_y = [y + y_shift for y in cyl_y]

    # Top hemisphere
    theta_top = np.linspace(0, np.pi, 100)
    top_x = r * np.cos(theta_top)
    top_y = h_cylinder + r * np.sin(theta_top) + y_shift

    # Bottom hemisphere
    theta_bot = np.linspace(np.pi, 2*np.pi, 100)
    bot_x = r * np.cos(theta_bot)
    bot_y = r * np.sin(theta_bot) + y_shift

    # Create figure
    plt.figure(figsize=(4, 8))

    # Draw tank
    plt.plot(cyl_x, cyl_y, 'b-', linewidth=2)
    plt.plot(top_x, top_y, 'b-', linewidth=2)
    plt.plot(bot_x, bot_y, 'b-', linewidth=2)

    # Fill with color
    plt.fill_between(cyl_x[:2], cyl_y[:2], cyl_y[2], color='skyblue', alpha=0.5)
    plt.fill_between(top_x, h_cylinder + y_shift, top_y, color='skyblue', alpha=0.5)
    plt.fill_between(bot_x, y_shift, bot_y, color='skyblue', alpha=0.5)

    # Set axis properties
    plt.gca().set_aspect('equal')
    plt.title('Propellant Tank Side View (Centered)')
    plt.xlabel('Width (mm)')
    plt.ylabel('Height (mm)')
    plt.grid(True)

    # Set limits to center the tank in view
    plt.xlim(-r * 1.5, r * 1.5)
    plt.ylim(-H/2 - r * 0.5, H/2 + r * 0.5)

    plt.show()

def calculate_mass(D, H, rho, t):
    """
    Calculate the mass of the tank.
    """
    A = (D/2)**2 * np.pi * 4 + D * np.pi * (H-D)
    return A * t * rho

def calculate_bottom_surface_side(D, L):
    return (D + L) * np.sqrt(2) / 2 + D

def plot_topview(deptwidth, D):
    half = deptwidth / 2
    square_x = [-half, half, half, -half, -half]
    square_y = [-half, -half, half, half, -half]

    # Define center points for the circles
    center_points = [
        (-half + D/2, -half + D/2),
        (half - D/2, -half + D/2),
        (half - D/2, half - D/2),
        (-half + D/2, half - D/2)
    ]

    # Create one figure and axes
    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot the square
    ax.plot(square_x, square_y, 'g-', linewidth=3, label='Satellite Boundary')

    # Plot the circles and their centers
    for i, (x, y) in enumerate(center_points):
        label = 'Propellant Tank' if i == 0 else None
        circle = plt.Circle((x, y), D/2, edgecolor='blue', facecolor='skyblue', alpha=0.5, linewidth=2, label=label)
        ax.add_patch(circle)
        ax.plot(x, y, 'ro') 

    # Axes settings
    ax.axhline(0, color='gray', linestyle='--')
    ax.axvline(0, color='gray', linestyle='--')
    ax.set_aspect('equal', adjustable='box')
    ax.set_title("Top View of Satellite Layout")
    ax.grid(True)

    # Set limits to fit all elements nicely
    buffer = D
    ax.set_xlim(-half - buffer, half + buffer)
    ax.set_ylim(-half - buffer, half + buffer)
    ax.legend(loc='upper right')
    plt.show()

def plot_side_view(D, H, depthwidth):
    # Convert meters to millimeters
    r = D / 2
    h_cylinder = H - D  # height of the cylindrical section

    # Centering shifts
    y_shift = - (H - D) / 2

    def draw_tank(x_shift):
        # Cylinder rectangle (shifted)
        cyl_x = [-r, r, r, -r, -r]
        cyl_x = [x + x_shift for x in cyl_x]
        cyl_y = [0, 0, h_cylinder, h_cylinder, 0]
        cyl_y = [y + y_shift for y in cyl_y]

        # Top hemisphere
        theta_top = np.linspace(0, np.pi, 100)
        top_x = r * np.cos(theta_top) + x_shift
        top_y = h_cylinder + r * np.sin(theta_top) + y_shift

        # Bottom hemisphere
        theta_bot = np.linspace(np.pi, 2*np.pi, 100)
        bot_x = r * np.cos(theta_bot) + x_shift
        bot_y = r * np.sin(theta_bot) + y_shift

        # Draw tank
        plt.plot(cyl_x, cyl_y, 'b-', linewidth=2)
        plt.plot(top_x, top_y, 'b-', linewidth=2)
        plt.plot(bot_x, bot_y, 'b-', linewidth=2)

        # Fill tank
        plt.fill_between(cyl_x[:2], cyl_y[:2], cyl_y[2], color='skyblue', alpha=0.5)
        plt.fill_between(top_x, h_cylinder + y_shift, top_y, color='skyblue', alpha=0.5)
        plt.fill_between(bot_x, y_shift, bot_y, color='skyblue', alpha=0.5)

    # Create figure
    plt.figure(figsize=(6, 8))

    # Draw satellite frame (rectangle)
    sat_x = [-depthwidth / 2, depthwidth / 2, depthwidth / 2, -depthwidth / 2, -depthwidth / 2]
    sat_y = [-H / 2, -H / 2, H / 2, H / 2, -H / 2]
    plt.plot(sat_x, sat_y, linewidth=2, label='Satellite Structure', color='green')

    # Draw both tanks
    x_offset = depthwidth / 2 - D / 2
    draw_tank(+x_offset)
    draw_tank(-x_offset)

    # Formatting
    plt.gca().set_aspect('equal')
    plt.title('Side View: Quad Propellant Tanks in Satellite')
    plt.xlabel('Width (m)')
    plt.ylabel('Height (m)')
    plt.grid(True)
    plt.legend()

    # Centered axis limits
    total_width = depthwidth * 1.2
    plt.xlim(-total_width / 2, total_width / 2)
    plt.ylim(-H / 2 - r * 0.5, H / 2 + r * 0.5)

    plt.show()

if __name__ == "__main__":
    # Things to input
    R = 0.5  # D/H ratio, choose yourself
    capacity = 129.920007  # in L, input the necessary tank capacity
    t = 0.002  # in m, input the necessary tank thickness
    rho = 4430  # in kg/m^3, input the necessary tank material density
    L = 0.312 # in m, input the free space length

    # Things that are calculated
    volume = calculate_tank_volume(capacity) # in L, input the necessary tank capacity
    mass_stat = calculate_tank_mass(capacity) #
    diameter, height = [x / 10 for x in calculate_tank_dimensions(volume, R)]
    mass = calculate_mass(diameter, height, rho, t)  # in kg
    depthwidth = calculate_bottom_surface_side(diameter, L)  # in m

    print(f"Diameter: {diameter:.3f} m")
    print(f"Height: {height:.3f} m")
    print(f"Mass: {mass:.3f} kg")
    print(f"Mass of the tank from statistics: {mass_stat:.3f} kg")
    print(f"Cube side: {depthwidth:.3f} m")

    plot_topview(depthwidth, diameter)
    plot_side_view(diameter, height, depthwidth)



    plot_tank(diameter, height)
