from functions import *
Vinf = 1
m = 0.9
alpha = 0
z = symbols("z")
f = Functions()

X, Y = f.grid()

# Define complex function
Z = X + 1j * Y
k = Vinf * np.exp(-1j * np.deg2rad(alpha))
W = k * Z**m

phi = np.real(W)
psi = np.imag(W)

# Define edge
beta1 = np.pi - np.pi / m
y1 = (-np.arctan(beta1) * X) * (X < 0) + 0 * (X > 0)

# Plot potential functions
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.contour(X, Y, phi * (Y > y1), 20, colors='r')
plt.plot(X, y1, linewidth=1)
plt.axis('equal')

plt.subplot(1, 2, 2)
plt.contour(X, Y, psi * (Y > y1), 20, colors='b')
plt.plot(X, y1,  linewidth=1)
plt.axis('equal')

plt.show()
