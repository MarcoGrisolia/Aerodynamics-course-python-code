import mpmath
import math
import decimal
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

x=symbols('x')


# toft = [1/sqrt(1 + 4*Pow(t,2)), 2*t/sqrt(1 + 4*Pow(t,2))]

# tprimet = [diff(toft[0],t), diff(toft[1],t)]

# noft = [tprimet[0]/sqrt(Pow(tprimet[0],2) + Pow(tprimet[1],2)), tprimet[1]/sqrt(Pow(tprimet[0],2) + Pow(tprimet[1],2))]

# import matplotlib.pyplot as plt
# #setup our t values array
# x_values = np.linspace(-4,4,100)

# #define the r(t) vector value function
# roft = [t, t**2]

# #turn our r(t) vector value symbolic math function into a python function
# r_value_functions = [lambdify(t, roft[0]), lambdify(t, roft[1])]

# #turn our T(t) vector value symbolic math function into a python function
# tangent_value_functions = [lambdify(t, toft[0]), lambdify(t, toft[1])]

# #turn our N(t) vector value symbolic math function into a python function
# normal_value_functions = [lambdify(t, noft[0]), lambdify(t, noft[1])]




def make_norm(vector_value_function):
  return sqrt(Pow(vector_value_function[0],2) + Pow(vector_value_function[1],2))




def compute_normal_graph(x, y, x_min = -5, x_max = 5):
  x_values = np.linspace(x_min, x_max ,100)

  y_prime = [diff(y[0],x), diff(y[1],x)]
  tanvector = [y_prime[0]/make_norm(y_prime), y_prime[1]/make_norm(y_prime)]
  tanprime = [diff(tanvector[0],x), diff(tanvector[1],x)]
  normalvector = [tanprime[0]/make_norm(tanprime), tanprime[1]/make_norm(tanprime)]
  normal_vector_functions = [lambdify(x, normalvector[0]),lambdify(x, normalvector[1])]
  value_functions = [lambdify(x, y[0]), lambdify(x, y[1])]
  plt.plot(value_functions[0](x_values), value_functions[1](x_values))

  for i in range(1,len(x_values)):
    if(i%5 == 0):
      normal_location = x_values[i]
      plt.quiver(value_functions[0](normal_location),value_functions[1](normal_location),normal_vector_functions[0](normal_location),normal_vector_functions[1](normal_location),color='r')









def compute_tangent_graph(x, y, x_min = -5, x_max = 5):
  x_values = np.linspace(x_min, x_max ,100)

  y_prime = [diff(y[0],x), diff(y[1],x)]
  tanvector = [y_prime[0]/make_norm(y_prime), y_prime[1]/make_norm(y_prime)]

  tan_vector_functions = [lambdify(x, tanvector[0]),lambdify(x, tanvector[1])]

  value_functions = [lambdify(x, y[0]), lambdify(x, y[1])]
  plt.plot(value_functions[0](x_values), value_functions[1](x_values))

  for i in range(1,len(x_values)):
    if(i%5 == 0):
      normal_location = x_values[i]
      plt.quiver(value_functions[0](normal_location),value_functions[1](normal_location),tan_vector_functions[0](normal_location),tan_vector_functions[1](normal_location),color='g')




compute_tangent_graph(symbols('x'), [x, x**2])






# plt.quiver(X, Y, U, V, scale=6, color='blue', width=0.005)

# # Add labels and a title
# plt.xlabel('X-Axis')
# plt.ylabel('Y-Axis')
# plt.title('Vector Field (Quiver Plot)')




# # Show the plot
# plt.axis('equal')
# plt.grid(True)
# plt.show()
