from functions import *
f = Functions()
x=symbols('x')
y=symbols('y')  







psiFunction = atan(1/(x**2 + y**2))
field_extension = 10
steps = 500j

diff_U = - diff(psiFunction,y)
diff_V = diff(psiFunction,x)


f_U_function = lambdify([x,y], diff_U)
f_V_function = lambdify([x,y], diff_V)

Y, X = np.mgrid[-field_extension:field_extension:steps, -field_extension:field_extension:steps]


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