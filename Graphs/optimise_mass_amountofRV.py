import numpy as np
import matplotlib.pyplot as plt

# Data
fuel_mass_2_3 = np.array([446, 444, 442.9, 442.9, 437, 443, 457, 477, 482, 501, 512, 521])
amount_of_RV_2_3 = np.array([172, 154, 143, 134, 125, 120, 113, 110, 104, 103, 101, 96])
distance_2_3 = ['4 m', '5 m', '6 m', '7 m', '8 m', '9 m', '10 m', '11 m', '12 m', '13 m', '14 m', '15 m']

fuel_mass_3_4 = np.array([456, 445, 446, 442, 439, 442, 462, 475, 496, 505, 513, 532])
amount_of_RV_3_4 = np.array([219, 176, 159, 145, 135, 128, 122, 116, 113, 109, 105, 102])
distance_3_4 = ['4 m', '5 m', '6 m', '7 m', '8 m', '9 m', '10 m', '11 m', '12 m', '13 m', '14 m', '15 m']

fuel_mass_4_5 = np.array([455, 448, 444, 443, 445, 465, 483, 500, 514, 526, 536])
amount_of_RV_4_5 = np.array([244, 188, 166, 152, 142, 134, 127, 121, 117, 114, 108])
distance_4_5 = ['5 m', '6 m', '7 m', '8 m', '9 m', '10 m', '11 m', '12 m', '13 m', '14 m', '15 m']

# Create scatter plot
plt.figure(figsize=(8, 6))

# Scatter plot for dataset 2_3
plt.scatter(amount_of_RV_2_3, fuel_mass_2_3, color='blue', edgecolors='black', label="2m - 3m closest distance")

# Annotate points for dataset 2_3
for i, label in enumerate(distance_2_3):
    plt.annotate(label, (amount_of_RV_2_3[i], fuel_mass_2_3[i]), textcoords="offset points", xytext=(5,5), ha='center', fontsize=6)

# Scatter plot for dataset 3_4
plt.scatter(amount_of_RV_3_4, fuel_mass_3_4, color='red', edgecolors='black', label="3m - 4m closest distance")

# Annotate points for dataset 3_4
for i, label in enumerate(distance_3_4):
    plt.annotate(label, (amount_of_RV_3_4[i], fuel_mass_3_4[i]), textcoords="offset points", xytext=(5,5), ha='center', fontsize=6)

# Scatter plot for dataset 4_5
plt.scatter(amount_of_RV_4_5, fuel_mass_4_5, color='green', edgecolors='black', label="4m - 5m closest distance")

# Annotate points for dataset 4_5
for i, label in enumerate(distance_4_5):
    plt.annotate(label, (amount_of_RV_4_5[i], fuel_mass_4_5[i]), textcoords="offset points", xytext=(5,5), ha='center', fontsize=6)

# Labels, title, and legend with reduced font size
plt.xlabel("Amount of RV", fontsize=10)
plt.ylabel("Fuel Mass", fontsize=10)
plt.title("Scatter Plot of Amount of RV vs. Fuel Mass for 1000 kg Wet Mass", fontsize=12)
plt.legend(fontsize=10)

# Show plot
plt.show()
