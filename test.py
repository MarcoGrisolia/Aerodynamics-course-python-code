import mpmath
import math
import decimal
import numpy as np
from sympy import *
import matplotlib.pyplot as plt




# y = [1/sqrt(1 + 2*Pow(x,2)), 2*x/sqrt(1 + 4*Pow(x,2))]

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



def computeVelocityfrom_Psi(function = symbols('x') + symbols('y'), x = symbols('x'), y = symbols('y'), field_extension = 6, steps = 100):

    z = function
    U = - diff(z,y)
    V = diff(z,x)

    return U, V







x=symbols('x')
y=symbols('y')


Psi = (atan(y/x) - y)

f_U , f_V = computeVelocityfrom_Psi(Psi)

f_U_function = lambdify([x,y], f_U)
f_V_function = lambdify([x,y], f_V)



field_extension = 10
steps = 300j


Y, X = np.mgrid[-field_extension:field_extension:steps, -field_extension:field_extension:steps]

# print(f_U)
# print(f_U_function(2,4))



#computes U vector 

f_U_vector = []
for i in range(0,len(X[0])):
    f_U_raw = []
    for j in range(0,len(X[0])):
      if not np.isnan(f_U_function(X[i][j],Y[i][j])):
        f_U_raw.append(f_U_function(X[i][j],Y[i][j]))      
      else:
         f_U_raw.append(0)   

    f_U_vector.append(f_U_raw)



#computes V vector
f_V_vector = []
for i in range(0,len(X[0])):
    f_V_raw = []
    for j in range(0,len(X[0])):
      if not np.isnan(f_V_function(X[i][j],Y[i][j])):
        f_V_raw.append(f_V_function(X[i][j],Y[i][j]))      

      else:
         f_V_raw.append(0)   

    f_V_vector.append(f_V_raw)






U = np.full_like(X, f_U_vector)
V = np.full_like(Y, f_V_vector)





print(U)
print(V)







plt.streamplot(X, Y, U, V, density = 1, color='blue', linewidth = None)


#plt.quiver(X, Y, U, V, scale=6, color='blue', width=0.005)

# # Add labels and a title
# plt.xlabel('X-Axis')
# plt.ylabel('Y-Axis')
# plt.title('Vector Field (Quiver Plot)')




# # Show the plot
# plt.axis('equal')
# plt.grid(True)
plt.show()
