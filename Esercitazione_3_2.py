from functions import *

f = Functions()
x = symbols("x")
y = symbols("y")
theta = atan2(y,x)

# psi = -log(x**2 + y**2)
# phi = theta

X1 , Y1 , U1, V1 = f.simpleIrrotational_Vortex(1)
X2 , Y2 , U2, V2 = f.simpleIrrotational_Vortex(-1, x0 = 1e-6)

U = U1 + U2
V = V1 + V2

plt.streamplot(X1, Y1, U, V)


X1 , Y1 , U1, V1 = f.simpleIrrotational_Vortex(1)
X2 , Y2 , U2, V2 = f.simpleIrrotational_Vortex(-1, y0 = 1e-6)

U = U1 + U2
V = V1 + V2

plt.streamplot(X1, Y1, U, V)
plt.axis('equal')
plt.grid(True)


plt.show()
