import mpmath
import math
import decimal
import numpy as np
from sympy import *
import matplotlib.pyplot as plt



class Functions():

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

    def grid(self, field_extension = 10, steps = 500j):
        Y, X = np.mgrid[-field_extension:field_extension:steps, -field_extension:field_extension:steps]
        return Y, X

    def Cylinder(self, radius = 1, steps = 500):
        theta = np.linspace(-np.pi, np.pi, steps)    
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        return x,y
    
    def SimpleVelocityVectorField(self, AoA = 0, intensity = 1, field_extension = 10, steps = 500j):
        U_constant = intensity * np.cos(AoA)
        V_constant = intensity * np.sin(AoA)

        Y, X = self.grid(field_extension, steps)

        # Define vector components (U and V) at each grid point

        U = U_constant * np.full_like(X, 1)
        V = V_constant * np.full_like(Y, 0)

        V_inf = sqrt(U[0][0]**2 + V[0][0]**2)

        return X, Y, U, V, V_inf
    
    def differentiateVelocityfrom_Psi(self, function, x = symbols('x'), y = symbols('y')):

        z = function
        diff_U = diff(z,y)
        diff_V = - diff(z,x)

        return diff_U, diff_V


    def flowLineVectorfield(self, psiFunction, field_extension = 10, steps = 500j):
                    

        x=symbols('x')
        y=symbols('y')   
        diff_U , diff_V = self.differentiateVelocityfrom_Psi(psiFunction)

        f_U_function = lambdify([x,y], diff_U)
        f_V_function = lambdify([x,y], diff_V)

        Y, X = self.grid(field_extension, steps)


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


        return X, Y, U, V

    def simpleSource_Sink(self, constant = 1,  field_extension = 10, steps = 500j):
        
        x=symbols('x')
        y=symbols('y')
        psiFunction = constant * atan(y/x)
   
        diff_U , diff_V = self.differentiateVelocityfrom_Psi(function= psiFunction)

        f_U_function = lambdify([x,y], diff_U)
        f_V_function = lambdify([x,y], diff_V)

        Y, X = self.grid(field_extension, steps)


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


        return X, Y, U, V        

    def simpleDoublet(self, constant = 1,  field_extension = 10, steps = 500j):
        
        x=symbols('x')
        y=symbols('y')
        psiFunction = constant/(2*np.pi) * ( y / (x**2 + y**2))
   
        diff_U , diff_V = self.differentiateVelocityfrom_Psi(function= psiFunction)

        f_U_function = lambdify([x,y], diff_U)
        f_V_function = lambdify([x,y], diff_V)

        Y, X = self.grid(field_extension, steps)


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


        return X, Y, U, V









