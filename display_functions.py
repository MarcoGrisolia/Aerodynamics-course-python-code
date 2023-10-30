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
a_mathematical = np.sqrt((1)/(2 * np.pi* flow_Velocity))
# gammaTheorical = 4*np.pi * a_mathematical
print(a_mathematical)



# psi = 1 / (x**2 + y**2)


X_velocityField, Y_velocityField, U_velocityField, V_velocityField, V_inf = f.SimpleVelocityVectorField(AoA = angle_of_attack ,intensity= -flow_Velocity , field_extension= field_width, steps = accuracy_j)

X_doublet, Y_doublet, U_doublet, V_doublet = f.simpleDoublet(constant= ((V_inf)*(circle_radius**2)*2*np.pi) , field_extension = field_width, steps= accuracy_j)  #don'r forget to put j after steps

# X_doublet, Y_doublet, U_doublet, V_doublet = f.flowLineVectorfield(psi, field_width, accuracy_j)


#print(f"V_inf = {V_inf}")

# Create grid                                                                   # Number of Y points
X      = np.linspace(-field_width,field_width,accuracy)                         # X-point array
Y      = np.linspace(-field_width,field_width,accuracy)                         # Y-point array    

#SUPERPOSITION OF SOLUTIONS                                         #don't forget to add the parameter j 
U =  U_velocityField + U_doublet  
V =  V_velocityField + V_doublet

Gamma, xC, yC, UC, VC = f.compute_Circulation(circle_radius,circle_radius,x0,y0,accuracy,U,V,X,Y)       # Call circulation calculation
print("Circulation: ", Gamma)                                                                           # Display circulation result

X_vortex, Y_vortex, U_vortex, V_vortex = f.simpleIrrotational_Vortex(gamma_over2pi= 1/ (2 * np.pi) , field_extension = field_width, steps= accuracy_j)


U =  U + U_vortex  
V =  V + V_vortex

# %% CALCULATIONS                                               
c_p = f.compute_Cp_from_velocity(U, V, flow_Velocity)
L , D = f.compute_Lift_and_Drag_from_Cp(c_p, X, Y, accuracy, circle_radius, circle_radius,circle_radius)
print(f"Lift = {L}")
print(f"Drag = {D}")

                                                                                
#print(f"Theorical vorticity gamma = {4*np.pi * a_cylinder}")



#PLOTTING THE STREAMLINES
plt.plot(f.Cylinder(circle_radius)[0], f.Cylinder(circle_radius)[1])
plt.streamplot(X, Y, U , V, density= 4)
plt.axis('equal')
plt.grid(True)
plt.show()
