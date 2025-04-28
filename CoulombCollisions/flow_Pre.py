import numpy as np
from matplotlib import pyplot as plt
from scipy.special import erf as erfc
def mu(x2):
    x = np.sqrt(x2)
    return erfc(x) - (2/np.sqrt(np.pi)) * x * np.exp(-x2)
T0 = 10*1.6e-19
epsilon_f = 9.11e-31 * (T0 / 9.11e-31 )/2

def dV_dt(V, t):
    factor = -(T0 / epsilon_f) ** 1.5 * mu(epsilon_f / T0)
    return factor * V

# RK4 method implementation
def rk4_step(func, V, t, dt):
    k1 = dt * func(V, t)
    k2 = dt * func(V + 0.5 * k1, t + 0.5 * dt)
    k3 = dt * func(V + 0.5 * k2, t + 0.5 * dt)
    k4 = dt * func(V + k3, t + dt)
    return V + (k1 + 2 * k2 + 2 * k3 + k4) / 6

# Solve the differential equation using RK4
def solve_rk4(func, V0, t0, t_end, dt):
    t_values = np.arange(t0, t_end, dt)
    V_values = np.zeros_like(t_values)
    V_values[0] = V0

    for i in range(1, len(t_values)):
        V_values[i] = rk4_step(func, V_values[i - 1], t_values[i - 1], dt)

    return t_values, V_values

V0 = (T0 / 9.11e-31 )  # Initial value of V
t0 = 0.0  # Start time
t_end = 11.0  # End time
dt = 0.01  # Time step

time1, V = solve_rk4(dV_dt, V0, t0, t_end, dt)


# plt.style.use('dark_background')
# Load data

data = np.loadtxt('Flow_Pre55.txt')

# Create subplots

time = data[:, 0]

plt.figure(figsize=(8,10))
plt.plot(time, data[:, 3]/data[0,3], color='black', linewidth=1.3, label='Experimental')
plt.plot(time1, V/V0, color='purple', linewidth=1, linestyle='dashed', label='Theoretical')
plt.ylabel('$V/V_0$', fontsize=12)
plt.xlabel('$\\hat{t}$', fontsize=12)
plt.tick_params(axis='both', direction='in', which='major', labelsize=11, top=True, right=True)
plt.tick_params(axis='both', direction='in', which='minor', top=True, right=True)
plt.tick_params(which='major', length=6)
plt.tick_params(which='minor', length=3)
plt.minorticks_on()
# Show the plot
plt.savefig(rf"./Figures/Figure7_Pre55.png")
