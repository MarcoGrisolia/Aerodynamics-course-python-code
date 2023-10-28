import mpmath
import math
import decimal
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

class Functions():
    def Circle(self, radius = 1, steps = 500):
        theta = np.linspace(-np.pi, np.pi, steps)    
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        return x,y
    
    def SimpleVectorField(self, AoA = 0, intensity = 1, grid_points = 6, x_extension = [-3,3], y_extension = [-3,3]):
        U_constant = intensity * np.cos(AoA)
        V_constant = intensity * np.sin(AoA)
        x = np.linspace(x_extension[0],x_extension[1], grid_points)
        y = np.linspace(y_extension[0],y_extension[1], grid_points)
        X, Y = np.meshgrid(x, y)

        # Define vector components (U and V) at each grid point

        U = U_constant * np.ones(grid_points ** 2)
        V = V_constant * np.ones(grid_points ** 2)

        return X, Y, U, V
        

    def simpleFlowLinefunction(self, field_extension = 6, steps = 100j):
        # Creating different arrows

        Y, X = np.mgrid[-field_extension:field_extension:steps, -field_extension:field_extension:steps]
        U = np.full_like(X, 0)
        V = np.full_like(X, 0.5)
        module_speed = np.sqrt(U**2 + V**2)

        return X , Y , U , V, module_speed
    


    def FlowLinefunction(self, function, field_extension = 1, steps = 3j):
            # Creating different arrows
            X, Y = np.mgrid[-field_extension:field_extension:steps, -field_extension:field_extension:steps]
            
            #f_U, f_V = self.differenctiateVelocityfrom_Psi(function)

            

            # U = f_U[X]
            # V = np.full_like(Y, 0.5)


            return X, Y



    def differenctiateVelocityfrom_Psi(function = symbols('x') + symbols('y'), x = symbols('x'), y = symbols('y'), field_extension = 6, steps = 100):

        z = function
        f_U = - diff(z,y)
        f_V = diff(z,x)

        return f_U, f_V





    def compute_normal(self, y_of_x, x=symbols('x'), field_extension = 6, steps = 100):

        t_values = np.linspace(- field_extension ,field_extension ,steps)
        y = y_of_x
        y_prime = [diff(y[0],x), diff(y[1],x)]
        tanvector = [y_prime[0]/self.make_norm(y_prime), y_prime[1]/self.make_norm(y_prime)]
        tanprime = [diff(tanvector[0],x), diff(tanvector[1],x)]
        normalvector = [tanprime[0]/self.make_norm(tanprime), tanprime[1]/self.make_norm(tanprime)]
        normal_vector_functions = [lambdify(x, normalvector[0]),lambdify(x, normalvector[1])]
        value_functions = [lambdify(x, y[0]), lambdify(x, y[1])]
        X = []
        Y = []
        U = []
        V = []

        for i in range(1, len(t_values)):
            normal_location = t_values[i]
            X.append(value_functions[0](normal_location))
            Y.append(value_functions[1](normal_location))
            U.append(normal_vector_functions[0](normal_location))
            V.append(normal_vector_functions[1](normal_location))

        return X, Y, U, V


    def make_norm(vector_value_function):
        return sqrt(Pow(vector_value_function[0],2) + Pow(vector_value_function[1],2))



x=symbols('x')
y=symbols('y')

f = Functions()



#X, Y, U, V = f.SimpleVectorField(10, 0.01,x_extension= [-10, 10], y_extension= [-10, 10])
Psi = x**2 + y**2


X, Y, U, V , speed = f.simpleFlowLinefunction(steps= 3j, field_extension=1)



# X, Y = f.FlowLinefunction(function= Psi, steps= 3j, field_extension=1)


print(U)
print(V)


#plt.plot(f.Circle()[0], f.Circle()[1])
#plt.quiver(X, Y, U, V, scale=7, color='blue', width=0.005)
#plt.streamplot(X, Y, U, V, density = .5, color='blue', linewidth = None)


# # Add labels and a title
# plt.xlabel('X-Axis')
# plt.ylabel('Y-Axis')
# plt.title('Vector Field (Quiver Plot)')


# # Show the plot
# plt.axis('equal')
# plt.grid(True)
#plt.show()