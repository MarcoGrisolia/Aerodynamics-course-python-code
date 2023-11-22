import numpy as np
from functions import *
from NacaEquation import Nacavectors_inputs

f = Functions()

accuracy = np.linspace(0, 1, 1000)
Nacaprofile = '0015'
Vinf = 1
AoA = np.deg2rad(0)


x_u, x_l, y_u, y_l = Nacavectors_inputs(accuracy, Nacaprofile)


# x_U = np.delete(x_l, -1)
# y_U = np.delete(y_u, -1)
x_L = np.delete(x_l, 0)
y_L = np.delete(y_l, 0)

# x_U = np.array(x_U)   
# y_L = np.array(y_L)   
x= np.append(x_L[::-1], x_u)
y= np.append(y_L[::-1], y_u)



x_control = []
y_control = []
theta = []

# Initialize Beta_ij and r_ij as empty lists
Beta_ij = []
r_ij = []

N = len(x) - 1

#FIRST CYCLE FOR THE NUMBER OF PANELS
for i in range(N):
    
    x_control.append((x[i + 1] + x[i]) / 2)
    y_control.append((y[i + 1] + y[i]) / 2)
    theta.append(np.arctan2(y[i + 1] - y[i], x[i + 1] - x[i]))

    Beta_row = []
    r = []
    #OPERATIONS TO DO WITH EVERY PANEL
    for j in range(N):
        dx1 = x_control[i] - x[j]
        dy1 = y_control[i] - y[j]
        dx2 = x_control[i] - x[j + 1]
        dy2 = y_control[i] - y[j + 1]

        Beta_row.append(np.arctan2(dy2 * dx1 - dx2 * dy1, dx1 * dx2 + dy1 * dy2))
        r.append(np.sqrt(dx1**2 + dy1**2))

        if i == j:
            Beta_row[j] = np.abs(Beta_row[j])

    r.append(np.sqrt(dx2**2 + dy2**2))

    # Append the calculated values to the lists
    Beta_ij.append(Beta_row)
    r_ij.append(r)







# Convert the lists to NumPy arrays
Beta_ij = np.array(Beta_ij)
r_ij = np.array(r_ij)
x_control = np.array(x_control)
y_control = np.array(y_control)
theta = np.array(theta)



A_i_j = []
b = []
A_i_j_lastRow = []
lastValue = []

for i in range(N):

    A_i_j_row = []
    summ_value = 0
    
    for j in range(N):

        A_i_j_row.append( (1/(2*np.pi)) * ((np.log((r_ij[i][j + 1])/(r_ij[i][j])) * np.sin(theta[i] - theta[j])) + (Beta_ij[i][j] * np.cos(theta[i] - theta[j]))))
        summ_value += ((np.log((r_ij[i][j + 1])/(r_ij[i][j])) * np.cos(theta[i] - theta[j])) - (Beta_ij[i][j] * np.sin(theta[i] - theta[j])))
    
        if i == N - 1:
            A_i_j_lastRow.append(1/(2*np.pi) * (((Beta_ij[0][j] * np.sin(theta[0] - theta[j]) - (np.log((r_ij[0][j+1])/(r_ij[0][j]))) * np.cos(theta[0] - theta[j]))) + 
                             ((Beta_ij[N - 1][j] * np.sin(theta[N - 1] - theta[j]) - (np.log((r_ij[N - 1][j+1])/(r_ij[N - 1][j]))) * np.cos(theta[N - 1] - theta[j])))))

            lastValue.append((1/(2*np.pi) * (((Beta_ij[0][j] * np.cos(theta[0] - theta[j]) + (np.log((r_ij[0][j+1])/(r_ij[0][j]))) * np.sin(theta[0] - theta[j]))) + 
                             ((Beta_ij[N - 1][j] * np.cos(theta[N - 1] - theta[j]) + (np.log((r_ij[N - 1][j+1])/(r_ij[N - 1][j]))) * np.sin(theta[N - 1] - theta[j]))))))
        
       

    A_i_j_row.append(summ_value /(2*np.pi))

    
    summ_value = 0


    if not sum(lastValue) == 0:
        A_i_j_lastRow.append(sum(lastValue))

    A_i_j.append(A_i_j_row)
    
    b.append(np.sin(theta[i] - AoA))
A_i_j.append(A_i_j_lastRow)





b.append(-(np.cos(theta[0] - AoA) + np.cos(theta[N - 1] - AoA)))


A_i_j = np.array(A_i_j)
b = np.array(b)


print(f"x: \n{x}")
print(f"y: \n{y}")
print(f"theta:\n{theta}")

print(f"x_control: \n {x_control}")
print(f"y_control: \n {y_control}")
print(f"Beta_ij:\n{Beta_ij}")
print(f"r_ij:\n{r_ij}")
print(f"A:\n{A_i_j}")
print(f"b: \n {b}")



solution = np.linalg.solve(A_i_j, b)
q = solution[:-1]  # All elements except the last one
Gamma = solution[-1]  # The last element



V_tang = Vinf *np.cos(theta - AoA) 


print(f"V_tang before= {V_tang}")


for i in range(0 , len(A_i_j) - 1):
    
    summ_value = 0

    for j in range(0 , len(A_i_j) - 1):
        summ_value += Gamma *((Beta_ij[i][j] * np.cos(theta[i] - theta[j])) +((np.sin(theta[i] - theta[j])) * (np.log((r_ij[i][j])/(r_ij[i][j])))))
        summ_value += q[j] * (Beta_ij[i][j] * np.sin(theta[i] - theta[j]) - ((np.cos(theta[i] - theta[j])) * (np.log((r_ij[i][j])/(r_ij[i][j])))))
    summ_value  *= 1 / (2 * np.pi)


    V_tang[i] += summ_value
    




print(f"q: \n {q}")
print(f"Gamma: \n {Gamma}")



print(f"V_tang:\n{V_tang}")



cp = 1 - (V_tang / Vinf)**2





fig, axs = plt.subplots(1,3)



axs[0].set_title("Airfoil")
axs[0].plot(x, y)
axs[0].plot(x_control, y_control, "o")
axs[0].axis('equal')
axs[0].grid(True)


axs[1].set_title("cp")
axs[1].plot(x_control, cp, "-o")
axs[1].axis('equal')
axs[1].grid(True)
axs[1].invert_yaxis()

axs[2].set_title("V/Vinf")
axs[2].plot(x_control, V_tang, "-o")
axs[2].axis('equal')
axs[2].grid(True)
axs[2].invert_yaxis()

 
#Show the plot
plt.show()


