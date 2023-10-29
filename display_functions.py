from functions import *

f = Functions()

x=symbols('x')
y=symbols('y')
r = symbols('r')
theta = symbols('theta')
circle_radius = 1
accuracy_j = 200j
accuracy = 200
field_width = 6
flow_Velocity = .2
a    = field_width / 2                                                                      # Horizontal axis half-length
b    = field_width / 2                                                                     # Vertical axis half-length
x0   = 0                                                                        # Ellipse center X coordinate
y0   = 0                                                                        # Ellipse center Y coordinate
numT = 400
angle_of_attack = 0

psi = 1 / (x**2 + y**2)


X_velocityField, Y_velocityField, U_velocityField, V_velocityField, V_inf = f.SimpleVelocityVectorField(AoA = angle_of_attack ,intensity= -flow_Velocity , field_extension= field_width, steps = accuracy_j)

X_doublet, Y_doublet, U_doublet, V_doublet = f.simpleDoublet(constant= ((V_inf)*(circle_radius**2)*2*np.pi) , field_extension = field_width, steps= accuracy_j) #don'r forget to put j after steps

#X_doublet, Y_doublet, U_doublet, V_doublet = f.flowLineVectorfield(psi, field_width, accuracy_j)

#print(f"V_inf = {V_inf}")


#SUPERPOSITION OF SOLUTIONS
Y, X = f.grid(field_width, accuracy_j)                                          #don't forget to add the parameter j 
U =  U_velocityField + U_doublet
V =  V_velocityField + V_doublet



# Create grid                                                                   # Number of Y points
X      = np.linspace(-field_width,field_width,accuracy)                         # X-point array
Y      = np.linspace(-field_width,field_width,accuracy)                         # Y-point array    



# %% CALCULATIONS

Gamma, xC, yC, UC, VC = f.compute_Circulation(a,b,x0,y0,accuracy,U,V,X,Y)       # Call circulation calculation
print("Circulation: ", Gamma)                                                   # Display circulation result




#PLOTTING THE STREAMLINES
plt.plot(f.Cylinder(circle_radius)[0], f.Cylinder(circle_radius)[1])
plt.streamplot(X, Y, U , V, density= 4)
plt.axis('equal')
plt.grid(True)
plt.show()
