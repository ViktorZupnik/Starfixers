import matplotlib.pyplot as plt
import numpy as np

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