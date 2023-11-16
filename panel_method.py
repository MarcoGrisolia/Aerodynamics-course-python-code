import numpy as np
from functions import *
from NacaEquation import Nacavectors_inputs

f = Functions()

accuracy = np.linspace(0, 1, 4)
Nacaprofile = '0015'

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




for i in range(len(x_control)):
    Beta_row = []
    r = []
    #OPERATIONS TO DO WITH EVERY PANEL
    for j in range(len(x) - 1):
        Beta_row.append(np.arctan2(((x[j] - x_control[i - 1]) * (y[j + 1] - y_control[i - 1])) - ((y[j] - y_control[i - 1]) * (x[j + 1] - x_control[i - 1])), ((x[j] - x_control[i - 1]) * (x[j + 1] - x_control[i - 1])) - ((y[j] - y_control[i - 1]) * (y[j + 1] - y_control[i - 1]))))
        r.append(abs(y_control[i - 1] - y[j]) + abs((x_control[i - 1] - x[j])))

        if i == j:
            Beta_row[j] = np.abs(Beta_row[j])



    # Append the calculated values to the lists
    Beta_ij.append(Beta_row)
    r_ij.append(r)







# Convert the lists to NumPy arrays
Beta_ij = np.array(Beta_ij)
r_ij = np.array(r_ij)
x_control = np.array(x_control)
y_control = np.array(r_ij)
theta = np.array(theta)




# u_ss = np.full_like(r_ij, 0)
# v_ss = np.full_like(Beta_ij, 0)
# u_vortex = np.full_like(Beta_ij, 0)
# v_vortex = np.full_like(r_ij, 0)
u_ss = []
v_ss = []
u_vortex = []
v_vortex = []

for i in range(0 , len(Beta_ij)):

    u_ss_row = []
    v_ss_row = []
    u_vortex_row = []
    v_vortex_row = []


    for j in range(0 , len(Beta_ij)):

        u_ss_row.append(- np.log((r_ij[i - 1][j])/(r_ij[i - 1][j - 1]))/ (2 * np.pi))
        v_ss_row.append(Beta_ij[i - 1][j - 1] / (2 * np.pi))
        u_vortex_row.append(Beta_ij[i - 1][j - 1])
        v_vortex_row.append(np.log((r_ij[i - 1][j])/(r_ij[i - 1][j - 1]))/ (2 * np.pi))


    u_ss.append(u_ss_row)
    v_ss.append(v_ss_row)
    u_vortex.append(u_vortex_row)
    v_vortex.append(v_vortex_row)



u_ss = np.array(u_ss)
v_ss = np.array(v_ss)
u_vortex = np.array(u_vortex)
v_vortex = np.array(v_vortex)





print(f"x= {x}")
print(f"x_control = {x_control}")
print(f"Beta_ij:\n{Beta_ij}")
print(f"r_ij:\n{r_ij}")

print(f"u_ss:\n{u_ss}")
print(f"v_ss:\n{v_ss}")




# normals = np.full_like(x, 0)
# tangents = np.full_like(x, 0)
# for i in range(len(x) - 1):
#     theta[i] = np.arctan2(y[i] - y[i+1], x[i] - x[i+1])
#     normals[i] = -np.sin(theta[i]) + np.cos(theta[i]) * 1j
#     tangents[i] = cos(theta[i]) + np.sin(theta[i]) * 1j










plt.plot(x, y)
plt.legend()
plt.axis('equal')
plt.grid(True)


# Show the plot
plt.show()


