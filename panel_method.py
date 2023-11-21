import numpy as np
from functions import *
from NacaEquation import Nacavectors_inputs

f = Functions()

accuracy = np.linspace(0, 1, 5)
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



#FIRST CYCLE FOR THE NUMBER OF PANELS
for i in range(1, len(x)):
    x_control.append((x[i] + x[i-1]) / 2)
    y_control.append((y[i] + y[i-1]) / 2)
    theta.append(np.arctan2(y[i] - y[i-1], x[i] - x[i-1]))

    Beta_row = []
    r = []
    #OPERATIONS TO DO WITH EVERY PANEL
    for j in range(len(x) - 1):
        Beta_row.append(np.arctan2(((x[j] - x_control[i - 1]) * (y[j + 1] - y_control[i - 1])) - ((y[j] - y_control[i - 1]) * (x[j + 1] - x_control[i - 1])), ((x[j] - x_control[i - 1]) * (x[j + 1] - x_control[i - 1])) - ((y[j] - y_control[i - 1]) * (y[j + 1] - y_control[i - 1]))))
        r.append(abs(y_control[i - 1] - y[j]) + abs((x_control[i - 1] - x[j])))

        if i == j + 1:
            Beta_row[j] = np.abs(Beta_row[j])



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

for i in range(0 , len(Beta_ij)):

    A_i_j_row = []

    for j in range(0 , len(Beta_ij)):

        A_i_j_row.append((np.log((r_ij[i - 1][j])/(r_ij[i- 1][j - 1])) * np.sin(theta[i- 1] - theta[j - 1])) + (Beta_ij[i- 1][j - 1] * np.cos(theta[i- 1] - theta[j - 1])))


    A_i_j_row.append((np.log((r_ij[i - 1][j])/(r_ij[i - 1][j - 1])) * np.cos(theta[i - 1] - theta[j - 1])) + (Beta_ij[i - 1 ][j - 1] * np.sin(theta[i - 1] - theta[j - 1])))
    A_i_j.append(A_i_j_row)
    b_value = 2 * np.pi * Vinf* np.sin(theta[i - 1] - AoA)

    A_i_j_row = []

    if i == len(Beta_ij) - 1:
        print(len(Beta_ij))
        for k in range(0, len(Beta_ij)):

            A_i_j_row.append((Beta_ij[i - 1][k - 1] * np.sin(theta[i - 1] - theta[k - 1]) - (np.log((r_ij[i - 1][k])/(r_ij[i- 1][k - 1]))) * np.cos(theta[i - 1] - theta[k - 1])))

            if k == len(Beta_ij) - 1:
                value = 0
                for l in range(0, len(Beta_ij) - 1):

                    value += (Beta_ij[i - 1][l - 1] * np.cos(theta[i - 1] - theta[l - 1]) + (np.log((r_ij[i - 1][l])/(r_ij[i- 1][l - 1]))) * np.sin(theta[i - 1] - theta[l - 1]))
                    
                b_value = - 2 * np.pi * Vinf* (np.cos(theta[i - 1] - AoA) + np.cos(theta[k - 1] - AoA))   
                A_i_j_row.append(value)
                
        A_i_j.append(A_i_j_row)
    b.append(b_value)

b.append(-2 * np.pi * Vinf* (np.cos(theta[len(Beta_ij) -  1] - AoA) + np.cos(theta[len(Beta_ij) - 1] - AoA)))




A_i_j = np.array(A_i_j)
# v_ss = np.array(v_ss)
# u_vortex = np.array(u_vortex)
# v_vortex = np.array(v_vortex)


solution = np.linalg.solve(A_i_j, b)
q = solution[:-1]  # All elements except the last one
Gamma = solution[-1]  # The last element



V_tang = Vinf *np.cos(theta - AoA) 


print(f"V_tang before= {V_tang}")


for i in range(0 , len(A_i_j) - 1):
    
    summ_value = 0

    for j in range(0 , len(A_i_j) - 1):
        summ_value += Gamma *((Beta_ij[i - 1][j - 1] * np.cos(theta[i - 1] - theta[j - 1])) +((np.sin(theta[i - 1] - theta[j - 1])) * (np.log((r_ij[i - 1][j])/(r_ij[i - 1][j - 1])))))
        summ_value += q[j] * (Beta_ij[i - 1][j - 1] * np.sin(theta[i - 1] - theta[j - 1]) - ((np.cos(theta[i - 1] - theta[j - 1])) * (np.log((r_ij[i - 1][j])/(r_ij[i - 1][j - 1])))))
    summ_value  *= 1 / (2 * np.pi)


    V_tang[i] += summ_value
    



# print(f"x: \n{x}")
# print(f"theta:\n{theta}")

print(f"x_control: \n {x_control}")
# print(f"Beta_ij:\n{Beta_ij}")
# print(f"r_ij:\n{r_ij}")
print(f"A:\n{A_i_j}")
print(f"b: \n {b}")

print(f"q: \n {q}")
print(f"Gamma: \n {Gamma}")
# print(f"v_ss:\n{v_ss}")


print(f"V_tang:\n{np.abs(V_tang)}")

# normals = np.full_lile(x, 0)
# tangents = np.full_like(x, 0)
# for i in range(len(x) - 1):
#     theta[i] = np.arctan2(y[i] - y[i+1], x[i] - x[i+1])
#     normals[i] = -np.sin(theta[i]) + np.cos(theta[i]) * 1j
#     tangents[i] = cos(theta[i]) + np.sin(theta[i]) * 1j








fig, axs = plt.subplots(2,2)



axs[0,0].set_title("Airfoil")
axs[0,0].plot(x, y)
axs[0,0].plot(x_control, y_control, "o")
axs[0,0].axis('equal')
axs[0,0].grid(True)


axs[0,1].set_title("Airfoil flow")
axs[0,1].plot(x_control, np.abs(V_tang / Vinf), "-o")
axs[0,1].axis('equal')
axs[0,1].grid(True)

# Show the plot
plt.show()


