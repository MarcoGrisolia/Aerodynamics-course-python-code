import numpy as np
from functions import *
from NacaEquation import Nacavectors_inputs

f = Functions()

accuracy = np.linspace(0, 1, 1000)
Nacaprofile = '0015'

x_U, x_L, y_U, y_L = Nacavectors_inputs(accuracy, Nacaprofile)

x = np.concatenate((x_U, x_L[::-1]))
y = np.concatenate((y_L[::-1], y_U))

x_control = np.full_like(x, 0)
y_control = np.full_like(y, 0)
theta = np.full_like(x, 0)

# Initialize Beta_ij and r_ij as empty lists
Beta_ij = []
r_ij = []

#FIRST CYCLE FOR THE NUMBER OF PANELS

for i in range(len(accuracy) - 1):
    x_control[i] = (x[i] + x[i+1]) / 2
    y_control[i] = (y[i] + y[i+1]) / 2
    theta[i] = np.arctan2(y[i] - y[i+1], x[i] - x[i+1])

    Beta = np.full_like(y, 0)
    r = np.full_like(y, 0)

    for j in range(len(accuracy) - 1):
        Beta[j] = np.arctan2(((x[j] - x_control[i]) * (y[j + 1] - y_control[i])) - ((y[j] - y_control[i]) * (x[j + 1] - x_control[i])), ((x[j] - x_control[i]) * (x[j + 1] - x_control[i])) - ((y[j] - y_control[i]) * (y[j + 1] - y_control[i])))
        r[j] = np.sqrt((y_control[i] - y[j])**2 + (x_control[i] - x[j])**2)

        if i == j:
            Beta = np.abs(Beta)

    # Append the calculated values to the lists
    Beta_ij.append(Beta)
    r_ij.append(r)

# Convert the lists to NumPy arrays
Beta_ij = np.array(Beta_ij)
r_ij = np.array(r_ij)

print(f"Beta_ij:\n{Beta_ij}")
print(f"r_ij:\n{r_ij}")





# normals = np.full_like(x, 0)
# tangents = np.full_like(x, 0)
# for i in range(len(x) - 1):
#     theta[i] = np.arctan2(y[i] - y[i+1], x[i] - x[i+1])
#     normals[i] = -np.sin(theta[i]) + np.cos(theta[i]) * 1j
#     tangents[i] = cos(theta[i]) + np.sin(theta[i]) * 1j










plt.plot(x_U, y_U)
plt.plot(x_L, y_L)
plt.legend()
plt.axis('equal')
plt.grid(True)


# Show the plot
plt.show()


