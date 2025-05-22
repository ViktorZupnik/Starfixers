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
    h_cylinder = H - D

    # Cylinder rectangle points
    cyl_x = [-r, r, r, -r, -r]
    cyl_y = [0, 0, h_cylinder, h_cylinder, 0]

    # Top hemisphere semicircle
    theta_top = np.linspace(0, np.pi, 100)
    top_x = r * np.cos(theta_top)
    top_y = h_cylinder + r * np.sin(theta_top)

    # Bottom hemisphere semicircle
    theta_bot = np.linspace(np.pi, 2*np.pi, 100)
    bot_x = r * np.cos(theta_bot)
    bot_y = r * np.sin(theta_bot)

    plt.figure(figsize=(4, 8))
    plt.plot(cyl_x, cyl_y, 'b-', linewidth=2)
    plt.plot(top_x, top_y, 'b-', linewidth=2)
    plt.plot(bot_x, bot_y, 'b-', linewidth=2)
    plt.fill_between(cyl_x[:2], cyl_y[:2], cyl_y[2], color='skyblue', alpha=0.5)
    plt.fill_between(top_x, h_cylinder, top_y, color='skyblue', alpha=0.5)
    plt.fill_between(bot_x, 0, bot_y, color='skyblue', alpha=0.5)

    plt.gca().set_aspect('equal')
    plt.title('Propellant tank Side View')
    plt.xlabel('Width (mm)')
    plt.ylabel('Height (mm)')
    plt.grid(True)
    plt.show()



def calculate_free_volume_and_area(D, H, A):
    """
    Calculate the free volume of the tank.
    """
    circle = np.pi * (D / 2)**2
    bigsquare = A**2
    return (bigsquare - circle * 4) / 4 * H, (bigsquare - circle * 4) / 4

def calculate_mass(D, H, rho, t):
    """
    Calculate the mass of the tank.
    """
    A = (D/2)**2 * np.pi * 4 + D * np.pi * (H-D)
    return A * t * rho


if __name__ == "__main__":
    R = 0.5  # D/H ratio, choose yourself
    capacity = 131.54  # in L, input the necessary tank capacity
    t = 0.001  # in m, input the necessary tank thickness

    volume = calculate_tank_volume(capacity) # in L, input the necessary tank capacity

    diameter, height = [x / 10 for x in calculate_tank_dimensions(volume, R)]
    mass = calculate_mass(diameter, height, 2700, 0.005)  # in kg
    print(f"Diameter: {diameter:.3f} m")
    print(f"Height: {height:.3f} m")
    print(f"Mass: {mass:.3f} kg")


    free_volume, free_area = calculate_free_volume_and_area(diameter, height, 2*diameter)




    plot_tank(diameter, height)
