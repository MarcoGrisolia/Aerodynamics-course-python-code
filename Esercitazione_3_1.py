from functions import *

### ESRECIZIO 1 

f = Functions()
x = symbols("x")
y = symbols("y")
Q = -2
Gamma_over_2pi = 3



# computing graph 1

X_sink, Y_sink , U_sink, V_sink = f.simpleSource_Sink(Q)
X_vortex, Y_vortex , U_vortex, V_vortex = f.simpleIrrotational_Vortex(Gamma_over_2pi)

U = U_vortex + U_sink
V = V_vortex + V_sink


# graph 1

fig, axs = plt.subplots(1, 2)

axs[0].streamplot(X_sink, Y_sink, U, V , density= 2)
axs[1].streamplot(X_sink, Y_sink, U, V , density= 2)
axs[0].set_title("isopotential")
axs[0].axis('equal')
axs[0].grid(True)


X_sink, Y_sink , U_sink, V_sink = f.simpleSource_Sink(Q)
X_vortex, Y_vortex , U_vortex, V_vortex = f.simpleIrrotational_Vortex(- Gamma_over_2pi)

U = U_vortex + U_sink
V = V_vortex + V_sink


axs[1].streamplot(X_sink, Y_sink, U, V , density= 2)


axs[1].set_title("phi and psi functions")
axs[1].axis('equal')
axs[1].grid(True)


plt.show()