from Naca equation import *

x = np.linspace(0, 1, 300)

# Calculate coordinates for the NACA airfoil
x_U, x_L, y_U, y_L = NACA_equation(x, 2412)

# Plot the upper and lower airfoil coordinates
plt.plot(x_U, y_U)
plt.grid(True)
plt.axis("equal")
plt.hold(True)  # In Python, this command is not necessary in Matplotlib.
plt.plot(x_L, y_L)

# Use subprocess to call XFOIL
foil_name = "NACA2412"  # Replace with the desired NACA profile
alphas = list(range(-20, 21))  # List of angles of attack

for alpha in alphas:
    command = f"xfoil {foil_name}\n"
    command += f"oper\n"
    command += f"iter 400\n"
    command += f"alpha {alpha}\n"
    command += f"pacc output.txt\n"  # Output polar data to a file
    command += f"gdes\n"
    command += f"exec\n"
    command += "quit\n"

    # Run XFOIL using subprocess
    subprocess.run(command, shell=True)

# Plot the airfoil coordinates
plt.show()