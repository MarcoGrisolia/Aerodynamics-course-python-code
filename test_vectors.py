import numpy as np
import matplotlib.pyplot as plt

# Creating different arrows
t = 3
steps = 100j
Y, X = np.mgrid[-t:t:steps, -t:t:steps]
U = np.full_like(X, 1.0)
V = np.full_like(X, 0.5)
speed = np.sqrt(U**2 + V**2)



fig, axs = plt.subplots()
# Using varying densities for streamplot(0)
axs.streamplot(X, Y, U, V, density=1)
axs.set_title('Density Variation')
plt.show()