import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

eta = np.linspace(0, 7, 100)

# Define the Blasius ODE system
def blasius_eq(f, t):
    return (f[1], f[2], -0.5 * f[0] * f[2])

def blas0(x, ts):
    initial_conditions = np.array([0, 0, x[0]])  
    f = odeint(blasius_eq, initial_conditions, ts)

    #returs the value of the initial condition for f(0)
    return 1 - f[-1, 1]

value = fsolve(blas0, [0.1], args=(eta,))

print(f'The correct initial condition for f(7.0) = 0 is {value[0]}')

initial_conditions = [0, 0, value[0]]

f = odeint(blasius_eq, initial_conditions, eta)

# Plot the solution
plt.plot(f, eta)
plt.xlabel('eta (y/delta)')
plt.ylabel('f(eta)')
plt.title('Blasius Solution')
plt.grid(True)
plt.legend(["Psi", "u/U_inf", "Tau_wall"], loc="lower right")

plt.savefig("output.png")