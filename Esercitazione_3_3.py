from functions import *

f = Functions()

x=symbols('x')
y=symbols('y')
r = symbols('r')
theta = symbols('theta')
circle_radius = 1
accuracy_j = 400j
accuracy = 400
field_width = 6
flow_Velocity = 1                                                                # Vertical axis half-length
x0   = 0                                                                        # Ellipse center X coordinate
y0   = 0                                                                        # Ellipse center Y coordinate
numT = 400
angle_of_attack = 0
a_mathematical = np.sqrt((1)/(2 * np.pi* flow_Velocity))
# gammaTheorical = 4*np.pi * a_mathematical



psi = x*y**3 / (x**2 + y**2)*(x + y)


X_velocityField, Y_velocityField, U_velocityField, V_velocityField, V_inf = f.SimpleVelocityVectorField(AoA = angle_of_attack ,intensity= -flow_Velocity , field_extension= field_width, steps = accuracy_j)

X_doublet, Y_doublet, U_doublet, V_doublet = f.simpleDoublet(constant= ((V_inf)*(circle_radius**2)*2*np.pi) , field_extension = field_width, steps= accuracy_j)  #don'r forget to put j after steps

# X_vortex, Y_doublet, U_doublet, V_doublet = f.flowLineVectorfield(psi, field_width, accuracy_j)


#print(f"V_inf = {V_inf}")

# Create grid                                                                   # Number of Y points
X      = np.linspace(-field_width,field_width,accuracy)                         # X-point array
Y      = np.linspace(-field_width,field_width,accuracy)                         # Y-point array    

#SUPERPOSITION OF SOLUTIONS                                         #don't forget to add the parameter j 
U =  U_velocityField + U_doublet  
V =  V_velocityField + V_doublet

                                                                          # Display circulation result

X_vortex, Y_vortex, U_vortex, V_vortex = f.simpleIrrotational_Vortex( gamma_over2pi= circle_radius*flow_Velocity , field_extension = field_width, steps= accuracy_j)


U =  U + U_vortex  
V =  V + V_vortex

Gamma, xC, yC, UC, VC = f.compute_Circulation(U,V,X,Y)       # Call circulation calculation
print("Circulation: ", Gamma) 
# %% CALCULATIONS                                               

c_p , cp_function = f.compute_Cp_from_velocity(X,Y,U,V, circle_radius, circle_radius, accuracy, flow_Velocity)
L, D = f.compute_Lift_and_Drag_from_cp(c_p, circle_radius, flow_Velocity, steps= accuracy)

print("Lift: ", L)
print("Drag: ", D)



                                                                                
#print(f"Theorical vorticity gamma = {4*np.pi * a_cylinder}")



#PLOTTING THE STREAMLINES
fig, axs = plt.subplots(1, 2)

axs[0].plot(f.Cylinder(circle_radius)[0], f.Cylinder(circle_radius)[1])
axs[0].streamplot(X, Y, U , V, density= 2)
axs[0].axis('equal')
axs[0].grid(True)

t = np.linspace( - np.pi / 2  , (3 * np.pi) / 2, len(X))
f_1 = np.cos(t)
f_3 = np.sin(t)

# Calculate F_y and F_x
cp_y = c_p * np.sin(t)
cp_x = c_p * np.cos(t)

# Calculate the magnitude of the force vector
F_tot = np.sqrt(cp_x**2 + cp_y**2)



axs[1].plot(f.Cylinder(circle_radius)[0], f.Cylinder(circle_radius)[1])
axs[1].quiver(f_1, f_3, cp_x, cp_y)
axs[1].axis('equal')
axs[1].grid(True)

plt.show()