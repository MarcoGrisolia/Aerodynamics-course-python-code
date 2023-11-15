from functions import *
from NacaEquation import *

accuracy = np.linspace(0,1,5)
Nacaprofile = '0015'

x_U , x_L , y_U, y_L = Nacavectors_inputs(accuracy, Nacaprofile)

x = np.concatenate((x_U, x_L[::-1]))

y = np.concatenate((y_L[::-1], y_U))

x_control = np.full_like(x, 0)
y_control = np.full_like(y, 0)
theta = np.full_like(x, 0)
Beta = np.full_like(x, 0)
r_ij = np.full_like(x, 0)

for i in range(len(x) - 1):
    x_control[i] = (x[i] + x[i+1])/2
    y_control[i] = (y[i] + y[i+1])/2
    theta[i] = np.arctan2(y[i] - y[i+1], x[i] - x[i+1])
    Beta[i] = np.arctan2(((x[i] - x_control[i - 1])*(y[i + 1] - y_control[i - 1])) - ((y[i] - y_control[i - 1])*(x[i + 1] - x_control[i -1])), 
                         ((x[i] - x_control[i - 1])*(x[i + 1] - x_control[i -1])) - ((y[i] - y_control[i -1])*(y[i + 1] - y_control[i - 1])))

    #devi calcolare rij e  poi fare un altro for per farlo per tutti i pannelli


print(f"Beta {Beta}")






# normals = np.full_like(x, 0)
# tangents = np.full_like(x, 0)
# for i in range(len(x) - 1):
#     theta[i] = np.arctan2(y[i] - y[i+1], x[i] - x[i+1])
#     normals[i] = -np.sin(theta[i]) + np.cos(theta[i]) * 1j
#     tangents[i] = cos(theta[i]) + np.sin(theta[i]) * 1j










plt.plot(x_U, y_U)
plt.plot(x_L, y_L)
plt.legend()
plt.axis('equal')
plt.grid(True)


# Show the plot
plt.show()


