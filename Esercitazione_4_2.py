from functions import *
from mpl_toolkits.mplot3d import axes3d

f = Functions()

Vinf = 1
C_x = .15
C_y = .15
radius = 2
angle = 10
AoA = np.deg2rad(angle)
Zc = C_x + 1j * C_y
accuracy = 500
accuracyj = 500j
field_ext = 3 * radius
rho = 997

thetaprime = np.linspace(- np.pi ,np.pi, accuracy)


b, Beta, e, t_max = f.K_J_transform(C_x, C_y, radius)


Z = np.array(f.complex_grid(field_extension=field_ext, steps= accuracyj))

zitaC = Zc + (b**2)/Zc


X, Y = f.grid(field_extension=field_ext)


z_cilinder = f.z_cylinder(circle_radius= radius, zc = Zc, steps=accuracy)

zita_airfoil = z_cilinder + (b**2) / np.array(z_cilinder) 

zita = np.full_like(Z, 0)

for i in range(len(Z)):
    for j in range(len(Z)):
        if np.abs(Z[i][j] - Zc) >= radius:
            zita[i][j] = Z[i][j] + (b**2) / (Z[i][j])

        else :
            zita[i][j] = nan
        



Gamma = -(2* Vinf * np.sin(np.pi + Beta + AoA)*(2*np.pi*radius))



def W(t, s):
    return -Vinf * (((t - s) * np.exp(1j * AoA)) + ((radius ** 2) * np.exp(-1j * AoA)) / (t - s)) - (
            (1j * Gamma * (np.log((t - s) / radius))) / (2 * np.pi))
    

V = 2* Vinf * np.sin(thetaprime + AoA) + (Gamma/(2*np.pi*radius))


V_zita = V / (abs((1 - b**2)/(Z**2))) 

cp = 1 - (V_zita / Vinf)**2




phi = np.real(W(Z, Zc))
psi = np.imag(W(Z, Zc)) 


eps = np.real(zita)
eta = np.imag(zita)




alpha = np.linspace(-10, 10, 20)

Gamma_1 = -(2* Vinf * np.sin(np.pi + Beta + np.deg2rad(alpha))*(2*np.pi*radius))

Lift = rho * Vinf * Gamma_1

cl_zita = Lift / (0.5 * rho * Vinf**2 * 4*b)


if (eps.all() + eta.all() == zita_airfoil.all()):
    print(f"zita {zita}")

 



fig, axs = plt.subplots(2,2)

axs[0,0].set_title("Cylider flow")
axs[0,0].plot(np.real(z_cilinder), np.imag(z_cilinder))
axs[0,0].contour(X, Y, psi * ((X - C_x)**2 + (Y-C_y)**2 > radius**2),  np.linspace(-3,3,50) , colors = '#A2142F')
axs[0,0].axis("equal")
axs[0,0].grid(True)

axs[0,1].set_title("Airfoil flow")
axs[0,1].plot(np.real(zita_airfoil), np.imag(zita_airfoil))
axs[0,1].contour(eps, eta, psi, np.linspace(-field_ext,field_ext,50),  colors = '#A2142F')
axs[0,1].axis('equal')
axs[0,1].grid(True)


axs[1,0].set_title("Cl vs alpha")
axs[1,0].plot(np.deg2rad(alpha), cl_zita)
axs[1,0].axis('equal')
axs[1,0].grid(True)


axs[1,1].set_title("Cp vs eps")
axs[1,1].plot(eps, cp * (eps.all() + eta.all() == zita_airfoil.any()))
axs[1,1].axis('equal')
axs[1,1].grid(True)

plt.show()
