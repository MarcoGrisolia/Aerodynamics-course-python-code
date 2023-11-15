import matplotlib.pyplot as plt
import numpy as np

from NacaEquation import Nacavectors_inputs

x = np.linspace(0,1,1000)
Nacaprofile = '6315'


x_U , x_L , y_U, y_L = Nacavectors_inputs(x, Nacaprofile)

plt.plot(x_U, y_U)
plt.plot(x_L, y_L)
plt.legend()
plt.axis('equal')
plt.grid(True)


# Show the plot
plt.show()

