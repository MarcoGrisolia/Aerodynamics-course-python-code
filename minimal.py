import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.special import erf

# Define the system of ODEs
def system(eta, y):
    f, f_prime = y
    dydt = [f_prime, -2 * eta * f_prime]
    return dydt

# Event function to stop integration when f approaches 1
def event(eta, y):
    return y[0] - 1

event.terminal = True  # Stop the integration when the event is triggered

# Desired boundary condition at infinity
f_inf_desired = 1

# Solve the system of ODEs with the shooting method
solution = solve_ivp(system, [0, np.inf], [0, 1], events=event, vectorized=True)

# Extract the solution at the event point
eta_event = solution.t_events[0][0]
f_event = solution.y_events[0][0, 0]

# Print the result
print(f"The solution reaches f({eta_event:.2f}) = {f_event:.2f}")

# Plot the solution
eta_values = np.linspace(0, eta_event, 1000)
f_values = solution.y[0, :len(eta_values)]

plt.plot(eta_values, f_values, label='f(eta)')
plt.xlabel('Eta')
plt.ylabel('f(eta)')
plt.legend()
plt.show()