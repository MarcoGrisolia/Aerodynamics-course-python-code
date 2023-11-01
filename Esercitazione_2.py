from functions import *

# ESERCIZIO NUMERO 1 

# initializing global variables
accuracy = 5000

Beta_1 = np.linspace(0, np.pi/6, 7)

Beta_2 = np.pi/4

b = 0.1

rho = 1000

V_1 = np.linspace(0, 15, accuracy)



for i in range(len(Beta_1)):

    #throug velocity triangle we find V_2
    V_2 = V_1 * (np.cos(Beta_1[i]) / np.cos(Beta_2))

    # throug bernoulli equations we find the pressure difference between the blades 
    delta_p =  0.5 * rho * (V_1**2 - V_2**2)

    # calculating the lift and drag from the excersice manually
    Drag = (delta_p * b)
    Lift = rho * b * V_1 * np.cos(Beta_1[i]) * (V_2 * np.sin(Beta_2) - V_1 * np.sin(Beta_1[i]))


    Force = np.sqrt(Drag**2 + Lift**2)


    plt.plot(V_1, Force, label = f"Beta = {round(Beta_1[i] * 180 / np.pi)}")


plt.legend(loc="upper left")
plt.xlabel("V_1")
plt.ylabel("Force")
plt.grid(True)
plt.show()
 



 












