import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

def f(A, t):
    """ Function to solve: coth(A) - 1/A - e^t """
    return 1/np.tanh(A) - 1/A - np.exp(-t)

def df(A):
    """ Derivative of f(A): -csch^2(A) + 1/A^2 """
    return -1/np.sinh(A)**2 + 1/A**2

def newton_raphson(t, A0=.01, tol=1e-6, max_iter=100):
    """ Newton-Raphson method to find A(t) """
    A = A0
    for i in range(max_iter):
        f_val = f(A, t)
        df_val = df(A)
        
        if abs(f_val) < tol:
            return A  # Converged
        
        A -= f_val / df_val  # Newton step
        
        if A <= 0:  # Ensure A remains positive
            A = A0 / 2  # Reset to a safer value
    return A  # Return the last value if not converged

# Generate data
# t_values = np.logspace(-5, 1, 100000)  # Range of t values
# A_values = np.array([newton_raphson(t) for t in t_values])  # Compute A(t)
# A_analytical = np.zeros(len(t_values))
# for t in range(len(t_values)):
#     if t_values[t]<=0.01:
#         A_analytical[t] = 1./t_values[t]
#     elif t_values[t]>0.01 and t_values[t]<0.1:
#         A_analytical[t] = (84.708*t_values[t]+166.4288)/(166.4681*t_values[t]-0.0003)
#     elif t_values[t]>=0.1 and t_values[t]<0.6:
#         A_analytical[t] = 1.2528*t_values[t]**(-0.9208)+0.1115*t_values[t]
#     elif t_values[t]>=0.6 and t_values[t]<2:
#         A_analytical[t] = (-50.8179*t_values[t]+181.7747)/(88.0921*t_values[t]+20.5056)
#     else:
#         A_analytical[t] = 3*np.exp(-t_values[t])



# # Plot the numerical and fitted solutions
# plt.figure(figsize=(8, 5))

# plt.plot(t_values, A_values, label='Numerical Solution (Newton-Raphson)',linewidth=2, color='k')
# plt.plot(t_values, A_analytical, '--', linewidth=1.6, label='Fitted Analytical Solution', color='r')


t_values = np.logspace(-5, 1, 100000)  # Range of t values
A_values = np.array([newton_raphson(t) for t in t_values])  # Compute A(t)

# Split analytical solution into parts
A_analytical_1 = np.where(t_values <= 0.01, 1. / t_values, np.nan)
A_analytical_2 = np.where((t_values > 0.01) & (t_values < 0.1), 
                          (84.708 * t_values + 166.4288) / (166.4681 * t_values - 0.0003), np.nan)
A_analytical_3 = np.where((t_values >= 0.1) & (t_values < 0.6), 
                          1.2528 * t_values**(-0.9208) + 0.1115 * t_values, np.nan)
A_analytical_4 = np.where((t_values >= 0.6) & (t_values < 2), 
                          (-50.8179 * t_values + 181.7747) / (88.0921 * t_values + 20.5056), np.nan)
A_analytical_5 = np.where(t_values >= 2, 3 * np.exp(-t_values), np.nan)

# Plot the numerical and fitted solutions
plt.figure(figsize=(8, 5))

# Plot numerical solution
plt.plot(t_values, A_values, label='Numerical Solution (Newton-Raphson)', linewidth=2, color='k')

# Plot analytical parts with different colors
plt.plot(t_values, A_analytical_1, '--', label=r'$A(t) = \frac{1}{t}$ (t ≤ 0.01)', linewidth=1.6, color='#1E88E5')
plt.plot(t_values, A_analytical_2, '--', label=r'$A(t) = \frac{84.708t + 166.4288}{166.4681t - 0.0003}$ (0.01 < t < 0.1)', linewidth=1.6, color='#FFC107')
plt.plot(t_values, A_analytical_3, '--', label=r'$A(t) = 1.2528t^{-0.9208} + 0.1115t$ (0.1 ≤ t < 0.6)', linewidth=1.6, color='#D81B60')
plt.plot(t_values, A_analytical_4, '--', label=r'$A(t) = \frac{-50.8179t + 181.7747}{88.0921t + 20.5056}$ (0.6 ≤ t < 2)', linewidth=1.6, color='#004D40')
plt.plot(t_values, A_analytical_5, '--', label=r'$A(t) = 3e^{-t}$ (t ≥ 2)', linewidth=1.6, color='orange')




plt.yscale('log')
plt.xscale('log')
plt.xlabel('t')
plt.ylabel('A(t)')
plt.title('Numerical vs Fitted Analytical Solutions (coth(A) - 1/A = e^t)')
plt.legend()
plt.grid()
plt.show()
#plt.savefig('Figures/AnalyticalSolution.png', dpi=300, bbox_inches='tight')