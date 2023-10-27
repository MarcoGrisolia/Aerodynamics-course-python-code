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
        

    def FlowLineFunction(self):
        X, Y, U, V = self.SimpleVectorField()

    def make_norm(vector_value_function):
        return sqrt(Pow(vector_value_function[0],2) + Pow(vector_value_function[1],2))






x = symbols('x')


f = Functions()
X, Y, U, V = f.SimpleVectorField(0, 0.3)


plt.plot(f.Circle()[0], f.Circle()[1])


plt.quiver(X, Y, U, V, scale=7, color='blue', width=0.005)
f.FlowLineFunction()
# Add labels and a title
plt.xlabel('X-Axis')
plt.ylabel('Y-Axis')
plt.title('Vector Field (Quiver Plot)')


# Show the plot
plt.axis('equal')
plt.grid(True)
plt.show()