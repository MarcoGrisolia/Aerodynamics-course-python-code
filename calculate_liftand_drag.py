import mpmath
import math
import decimal
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from functions import *


def c_p(x):
    return (1.814 * x - 0.271 * x**3 - 0.0471 * x**5)**2

# Generate an array of theta values
theta = np.linspace(-np.pi, np.pi, 300)


# Calculate Lift and Drag coefficients
Lift, _ = quad(lambda theta: c_p(theta) * np.sin(theta), -np.pi, np.pi)
Drag, _ = quad(lambda theta: c_p(theta) * np.cos(theta), -np.pi, np.pi)
Lift /= 2
Drag /= 2


print("Lift:", Lift)
print("Drag:", Drag)


