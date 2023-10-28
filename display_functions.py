from functions import *

f = Functions()

x=symbols('x')
y=symbols('y')
r = symbols('r')
theta = symbols('theta')
circle_radius = 1
accuracy = 200j
field = 4



X_velocityField, Y_velocityField, U_velocityField, V_velocityField, V_inf = f.SimpleVelocityVectorField(intensity= -1 , field_extension= field, steps = accuracy)

X_doublet, Y_doublet, U_doublet, V_doublet = f.simpleDoublet(constant= ((V_inf)*(circle_radius**2)*2*np.pi) , field_extension= field, steps= accuracy) #don'r forget to put j after steps

print(V_inf)


Y, X = f.grid(4, accuracy) #don't forget to add the parameter j 
U = U_doublet + U_velocityField
V = V_doublet + V_velocityField

plt.plot(f.Cylinder(circle_radius)[0], f.Cylinder(circle_radius)[1])
plt.streamplot(X, Y, U , V, density= 2)

plt.axis('equal')
plt.grid(True)
plt.show()