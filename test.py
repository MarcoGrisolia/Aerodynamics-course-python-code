import mpmath
import math
import decimal
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

import mpmath
import math
import decimal
import numpy as np
from sympy import *


t=symbols('t')
toft = [1/sqrt(1 + 4*Pow(t,2)), 2*t/sqrt(1 + 4*Pow(t,2))]

tprimet = [diff(toft[0],t), diff(toft[1],t)]

noft = [tprimet[0]/sqrt(Pow(tprimet[0],2) + Pow(tprimet[1],2)), tprimet[1]/sqrt(Pow(tprimet[0],2) + Pow(tprimet[1],2))]

import matplotlib.pyplot as plt
#setup our t values array
t_values = np.linspace(-4,4,100)

#define the r(t) vector value function
roft = [t, t**2]

#turn our r(t) vector value symbolic math function into a python function
r_value_functions = [lambdify(t, roft[0]), lambdify(t, roft[1])]

#turn our T(t) vector value symbolic math function into a python function
tangent_value_functions = [lambdify(t, toft[0]), lambdify(t, toft[1])]

#turn our N(t) vector value symbolic math function into a python function
normal_value_functions = [lambdify(t, noft[0]), lambdify(t, noft[1])]

#plot r(t)
plt.plot(r_value_functions[0](t_values), r_value_functions[1](t_values))

#We are going to plot the vectors at the location of x=2 (which means t=2)
#plot the tangent unit vector in green
plt.quiver(r_value_functions[0](2),r_value_functions[1](2),tangent_value_functions[0](2),tangent_value_functions[1](2),color='g')

#plot the normal unit vector in red
plt.quiver(r_value_functions[0](2),r_value_functions[1](2),normal_value_functions[0](2),normal_value_functions[1](2),color='r')

#make sure that our axis grids are equal so that we can clearly see the perpendicular geometry of the Normal Vector and the tangent geometry of the Tangent Vector.
plt.axis('equal')
plt.grid()

#plot it
plt.show()







# plt.quiver(X, Y, U, V, scale=6, color='blue', width=0.005)

# # Add labels and a title
# plt.xlabel('X-Axis')
# plt.ylabel('Y-Axis')
# plt.title('Vector Field (Quiver Plot)')




# # Show the plot
# plt.axis('equal')
# plt.grid(True)
# plt.show()
