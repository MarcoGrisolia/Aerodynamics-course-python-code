import numpy as np
import matplotlib.pyplot as plt
from functions import Functions

f = Functions()

Vinf = 1
C_x = 0.1
C_y = 0.1
radius = 1
AoA = np.pi / 18
Zc = C_x + 1j * C_y
accuracy = 500

thetaprime = np.linspace(0, 2 * np.pi, accuracy)

b, Beta, e, t_max = f.K_J_transform(C_x, C_y, radius)

Z = f.complex_grid()

X, Y = f.grid(field_extension=3)

Gamma = -(2 * Vinf * np.sin(np.pi + Beta + AoA) * (2 * np.pi * radius))

W = -Vinf * (((np.array(Z) - Zc) * np.exp(1j * AoA)) + ((radius ** 2) * np.exp(-1j * AoA)) / (np.array(Z) - Zc)) - (
        (1j * Gamma * (np.log((np.array(Z) - Zc) / radius))) / (2 * np.pi))

phi = np.real(W)
psi = np.imag(W)

# Calculate the velocity field components
U = np.real(W)
V = -np.imag(W)

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(f.Cylinder(radius, xC=C_x)[0], f.Cylinder(radius, yC=C_y)[1])
plt.streamplot(X, Y, U, V, density=1, color='r')
plt.axis('equal')
plt.title('Streamlines around the cylinder')

plt.show()