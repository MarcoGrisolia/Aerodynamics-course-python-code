# import numpy as np
# import matplotlib.pyplot as plt
# from functions import *

# f = Functions()
# eps = 1e-6
# accuracy_j = 1000j
# field = 5

# X1, Y1, U1, V1 = f.simpleSource_Sink(1, 1, 0, field, accuracy_j)
# _, _, U2, V2 = f.simpleSource_Sink(-1, 1 + eps, 0, field, accuracy_j)

# U = U1 + U2
# V = V1 + V2

# speed = np.sqrt(U**2 + V**2)
# lw = 5 * speed / speed.max()

# plt.streamplot(X1, Y1, U, V, linewidth=lw, density=0.8, color='k', arrowsize=1.5)

# # Add labels and title
# plt.xlabel('X-axis')
# plt.ylabel('Y-axis')
# plt.axis("equal")
# plt.title('Stream Plot with Customized Arrows')
# plt.grid(True)
# # Show the plot
# plt.show()

import numpy as np
import matplotlib.pyplot as plt
from functions import *

f = Functions()
eps = 1e-6
accuracy_j = 1000j
field = 7
magnet_intensity = 1

X1, Y1, U1, V1 = f.simpleSource_Sink(magnet_intensity, 1, 0, field, accuracy_j)
_, _, U2, V2 = f.simpleSource_Sink(-magnet_intensity, 1 + eps, 0, field, accuracy_j)
_, _, U3, V3 = f.simpleSource_Sink(magnet_intensity, -1, 0, field, accuracy_j)
_, _, U4, V4 = f.simpleSource_Sink(-magnet_intensity, -1 + eps, 0, field, accuracy_j)

# Combine U and V components into arrays
U = U1 + U2 + U3 + U4
V = V1 + V2 + V3 + V4


plt.streamplot(X1, Y1, U, V, density=3, color='k', arrowstyle="fancy", linewidth=.6)



# Add labels and title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.axis("equal")
plt.title('Magnetic Field')
plt.grid(True)
# Show the plot
plt.tight_layout()
plt.show()