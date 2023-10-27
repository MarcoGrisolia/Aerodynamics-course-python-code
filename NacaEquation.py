import matplotlib.pyplot as plt
import numpy as np

def y_t_function(x, b):
        return 5 * b * ((0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4))
def function1(x, m, p):
    return (m / (p**2)) * (2 * p * x - x**2)
def function2(x, m, p):
    return (m / ((1 - p)**2)) * ((1 - 2 * p) + 2 * p * x - x**2)
def dy_over_dx_function1(x, p, m):
    return ((2 * m) / (p**2) * (p - x))
def dy_over_dx_function2(x, p, m):
    return ((2 * m) / (1 - p**2) * (p - x))

def Nacavectors_inputs(x_i,input):

    digits = str(input).zfill(4)
    d1, d2, d3, d4 = map(int, digits)
    b = (d3 * 10 + d4) / 100
    p = d2 / 10
    m = d1 / 100
    condition1 = (x_i >= 0) & (x_i < p)
    condition2 = (x_i >= p) & (x_i <= 1)
    y_t = y_t_function(x_i, b)

    # Create an array of boolean values indicating where x falls between two values          

    y_c = np.zeros(len(x_i))          # Initialize an array to store the y_c values
    dy_over_dx = np.zeros(len(x_i))
    # Calculate y_c values based on the condition
    if p != 0:
        y_c[condition1] = function1(x_i[condition1], m, p)
        dy_over_dx[condition1] = dy_over_dx_function1(x_i[condition1], p, m)
        y_c[condition2] = function2(x_i[condition2], m, p)
        dy_over_dx[condition2] = dy_over_dx_function2(x_i[condition2], p, m)
    else:
        y_c[condition2] = function2(x_i[condition2], m, p)
        dy_over_dx[condition2] = dy_over_dx_function2(x_i[condition2], p, m)


    theta = np.arctan(dy_over_dx)

    if d1 + d2 == 0:
        x_U = x_i
        x_L = x_i
        y_U = y_t
        y_L = -y_t
        
    else:
        x_U = x_i - y_t * np.sin(theta)
        x_L = x_i + y_t * np.sin(theta)
        y_U = y_c + y_t * np.cos(theta)
        y_L = y_c - y_t * np.cos(theta)

    return x_U, x_L, y_U, y_L





