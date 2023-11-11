from functions import *

f = Functions()

Vinf = 1
C_x = .056
C_y = .056
radius = 1
AoA = np.pi/18
Zc = 0
accuracy = 500

thetaprime = np.linspace(0,2*np.pi, accuracy)


b, Beta, e, t_max = f.K_J_transform(C_x, C_y, radius)
print(b)

Z = f.complex_grid()

#X , Y , U , V , Gamma, L, D = f.flowFieldCylinder(circle_radius= radius, angle_of_attack=AoA, field_extension=3)

X, Y = f.grid(field_extension=3)

a = (b / np.cos(Beta))

Gamma = -(2* Vinf * np.sin(np.pi + Beta + AoA)*(2*np.pi*a))

print(Gamma)

W = - Vinf * (((Z - Zc) * np.exp(1j * AoA)) + ((a**2) *np.exp(-1j * AoA)) / (Z - Zc)) - (1j*Gamma * (np.log((Z- Zc)/a)))/(2*np.pi)

V = 2* Vinf * np.sin(thetaprime + AoA) + (Gamma/(2*np.pi*a))


phi = np.real(W)
psi = np.imag(W) 

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(f.Cylinder(radius)[0], f.Cylinder(radius)[1])
plt.contour(X, Y, psi, 70, color = "r")
plt.axis('equal')


plt.show()
