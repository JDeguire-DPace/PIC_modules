import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

plt.style.use('dark_background')
# Load data

data2 = np.loadtxt('Wang2.txt')
data20 = np.loadtxt('Wang20.txt')
data200 = np.loadtxt('Wang200.txt')
data1 = np.loadtxt('Wang1.txt')
data_5 = np.loadtxt('Wang_5.txt')
data_50 = np.loadtxt('Wang_50.txt')

# time = data16[:, 0]  # Time in seconds
# tau = 1.31e-7
# nu = 1 # / tau  # Frequency in Hz

# normalized_diff = (data16[:, 2] - data16[:, 3]) / (data16[0, 2] - data16[0, 3])

# # Define the exponential function for fitting
# def exp_func(t, b):
#     return 1 * np.exp(-b * t)

# Fit the exponential function to the data
# popt, pcov = curve_fit(exp_func, time, normalized_diff, p0=(1/tau))
# b_fit = popt[0]
# print(1/b_fit)
# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(16, 6))  # Create a figure with 1 row and 2 columns

# First plot

axes[0].plot(data2[:,0], data2[:, 2], color='white', linewidth=1)
axes[0].plot(data2[:,0], data2[:, 3], color='white', linewidth=1)
axes[0].plot(data20[:,0], data20[:, 2], color='violet', linewidth=1)
axes[0].plot(data20[:,0], data20[:, 3], color='violet', linewidth=1)
axes[0].plot(data200[:,0], data200[:, 2], color='red', linewidth=1)
axes[0].plot(data200[:,0], data200[:, 3], color='red', linewidth=1)
axes[0].plot(data1[:,0], data1[:, 2], color='orange', linewidth=1)
axes[0].plot(data1[:,0], data1[:, 3], color='orange', linewidth=1)
axes[0].plot(data_5[:,0], data_5[:, 2], color='cyan', linewidth=1)
axes[0].plot(data_5[:,0], data_5[:, 3], color='cyan', linewidth=1)
axes[0].plot(data_50[:,0], data_50[:, 2], color='lawngreen', linewidth=1)
axes[0].plot(data_50[:,0], data_50[:, 3], color='lawngreen', linewidth=1)

#axes[0].set_xlim(0, 2e-6*nu/12.5)
#axes[0].set_ylim(7e-3, 10.5e-3)
axes[0].set_title('Temperature Evolution')
axes[0].set_xlabel('Time (s)')
axes[0].set_ylabel('$T_p,T_z (eV)$')
tau=1.75e-10
# Second plot
axes[1].plot(data2[:,0], (data2[:, 2] - data2[:, 3]) / (data2[0, 2] - data2[0, 3]), color='white', linewidth=1, label = fr'${tau/2}s$' )
axes[1].plot(data20[:,0], (data20[:, 2] - data20[:, 3]) / (data20[0, 2] - data20[0, 3]), color='violet', linewidth=1,label = fr'${tau/20}s$')
axes[1].plot(data200[:,0], (data200[:, 2] - data200[:, 3]) / (data200[0, 2] - data200[0, 3]), color='red', linewidth=1,label = fr'${tau/200}s$')
axes[1].plot(data1[:,0], (data1[:, 2] - data1[:, 3]) / (data1[0, 2] - data1[0, 3]), color='orange', linewidth=1,label = fr'${tau/2}s$')
axes[1].plot(data_5[:,0], (data_5[:, 2] - data_5[:, 3]) / (data_5[0, 2] - data_5[0, 3]), color='cyan', linewidth=1,label = fr'${5*tau}s$')
axes[1].plot(data_50[:,0], (data_50[:, 2] - data_50[:, 3]) / (data_50[0, 2] - data_50[0, 3]), color='lawngreen', linewidth=1,label = fr'${50*tau}s$')
axes[1].plot(data20[:,0], np.exp(-data20[:,0]/2.5e-8), color='yellow', linewidth=0.6, linestyle='dashed', label='theroretical')
# axes[1].plot(time * nu, exp_func(time, b_fit), color='black', linestyle='dashed', linewidth=1, label='Fit')
axes[1].set_title('Normalized Temperature Difference')
axes[1].set_xlabel('Time (s)')
axes[1].set_ylabel('Normalized Difference')
plt.legend()
# Adjust layout and show the plot
plt.tight_layout()
# Show the plot
plt.show()
