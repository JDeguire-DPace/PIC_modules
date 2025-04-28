import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

# Load data
dataRan2 = np.loadtxt('Pre55_A.txt')
data = np.loadtxt('Pre55_B.txt')
time = data[:, 0] 
timeRan2 = dataRan2[:,0]
# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(16, 9))  # Create a figure with 1 row and 2 columns

# First plot
axes[0].plot(time, data[:, 2], color='black', linewidth=0.9)
axes[0].plot(time, data[:, 3], color='black', linewidth=0.9)
axes[0].plot(time,1+(data[0,3]-1)*np.exp(-time),color='black',linestyle='dashed',linewidth=0.7)
axes[0].plot(time,1+(data[0,2]-1)*np.exp(-time),color='black',linestyle='dashed',linewidth=0.7)
axes[0].plot(time,1+(data[0,3]-1)*np.exp(-8 * time / (5 * np.sqrt(2 * np.pi))),linewidth=0.6,color='purple',linestyle='dashed')
axes[0].plot(time,1+(data[0,2]-1)*np.exp(-8 * time / (5 * np.sqrt(2 * np.pi))),linewidth=0.6,color='purple',linestyle='dashed')
axes[0].set_xlim(-0.3, 12.3)
axes[0].set_ylim(0.891, 1.209)
axes[0].set_yticks([0.9,1,1.1,1.2],['0.9','1.0','1.1','1.2'])
axes[0].set_ylabel('$T_x/T_e$ , $T_y/T_e$', fontsize=12)
axes[0].set_xlabel('$\\hat{t}$', fontsize=12)
axes[0].tick_params(axis='both', direction='in', which='major', labelsize=11, top=True, right=True)
axes[0].tick_params(axis='both', direction='in', which='minor', top=True, right=True)
axes[0].tick_params(which='major', length=6)
axes[0].tick_params(which='minor', length=3)
axes[0].minorticks_on()

# Second plot
axes[1].plot(time, (data[:, 2] - data[:, 3]) / (data[0, 2] - data[0, 3]), color='black', linewidth=0.9, label='Experimental')
axes[1].plot(time, np.exp(-time), linewidth=0.7, color='black', linestyle='dashed',
             label=r'$\frac{\Delta T}{\Delta T_0} = \exp(-\hat{t})$')
axes[1].plot(time, np.exp(-8 * time / (5 * np.sqrt(2 * np.pi))), linewidth=0.6, color='purple', linestyle='dashed',
             label=r'$\frac{\Delta T}{\Delta T_0} = \exp\left(-\frac{5}{8\sqrt{2\pi}}\hat{t}\right)$')
axes[1].set_ylabel('Normalized Temperature Difference', fontsize=12)
axes[1].set_xlabel('$\\hat{t}$', fontsize=12)
axes[1].tick_params(axis='both', direction='in', which='major', labelsize=11, top=True, right=True)
axes[1].tick_params(axis='both', direction='in', which='minor', top=True, right=True)
axes[1].tick_params(which='major', length=6)
axes[1].tick_params(which='minor', length=3)

axes[1].minorticks_on()

# Add legend
plt.legend(fontsize=15)

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
