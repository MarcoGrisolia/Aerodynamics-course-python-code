from functions import *
from mpl_toolkits.mplot3d import axes3d

f = Functions()

Vinf = 1

C_x , C_y, radius , chord, e = f.conform_param(.15)



angle = 10
AoA = np.deg2rad(angle)
Zc = C_x + 1j * C_y
accuracy = 500
accuracyj = 500j
field_ext = 3 * radius
rho = 997

thetaprime = np.linspace(- np.pi ,np.pi, accuracy)

theta = np.linspace(0 , 2 * np.pi, accuracy)


b, Beta, e, t_max = f.conform_transformation(C_x, C_y, radius)


Z = np.array(f.complex_grid(field_extension=field_ext, steps= accuracyj))

zitaC = Zc + (b**2)/Zc


X, Y = f.grid(field_extension=field_ext, steps= accuracyj)


z_cilinder = f.z_cylinder(circle_radius= radius, zc = Zc, steps=accuracy)

zita_airfoil = z_cilinder + (b**2) / np.array(z_cilinder) 

min_zita_airfoil = np.min(np.real(zita_airfoil))
max_zita_airfoil = np.max(np.real(zita_airfoil))

chord_lenght = max_zita_airfoil - min_zita_airfoil



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


V_airfoil = np.abs(V / np.abs(1 - ((b**2)/(z_cilinder**2))))


cp = 1 - (V_airfoil / Vinf)**2


phi = np.real(W(Z, Zc))
psi = np.imag(W(Z, Zc)) 


eps = np.real(zita)
eta = np.imag(zita)




alpha = np.linspace(-10, 10, 20)

Gamma_1 = -(2* Vinf * np.sin(np.pi + Beta + np.deg2rad(alpha))*(2*np.pi*radius))

Lift = rho * Vinf * Gamma_1

cl_zita = Lift / (0.5 * rho * Vinf**2 * 4*b)



print(f"z_cilinder : \n {z_cilinder}")
print(f"V_airfoil : \n {V_airfoil}")



fig, axs = plt.subplots(2,3)

axs[0,0].set_title("Cylider flow")
axs[0,0].plot(np.real(z_cilinder), np.imag(z_cilinder))
axs[0,0].contour(X, Y, psi * ((X - C_x)**2 + (Y-C_y)**2 > radius**2),  np.linspace(-field_ext,field_ext,50) , colors = '#A2142F')
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


axs[1,1].set_title("V/Vinf vs Eps")
axs[1,1].plot(np.real(zita_airfoil) / (chord_lenght), V_airfoil  / Vinf,)
axs[1,1].axis('equal')
axs[1,1].grid(True)


axs[1,2].set_title("Cp vs Eps")
axs[1,2].plot(np.real(zita_airfoil) / (chord_lenght), cp)
axs[1,2].axis('equal')
axs[1,2].grid(True)
axs[1,2].invert_yaxis()



plt.show()
