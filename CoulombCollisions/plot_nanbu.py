import numpy as np
from matplotlib import pyplot as plt

plt.style.use('dark_background')
# Load data

data_6 = np.loadtxt('Coulomb_6.txt')
data_7 = np.loadtxt('Coulomb_7.txt')
data_8 = np.loadtxt('Coulomb_8.txt')
data_9 = np.loadtxt('Coulomb_9.txt')
data_10 = np.loadtxt('Coulomb_10.txt')
data_11 = np.loadtxt('Coulomb_11.txt')
# data_12 = np.loadtxt('Coulomb_12.txt')
# data_13 = np.loadtxt('Coulomb_13.txt')
# data_14 = np.loadtxt('Coulomb_14.txt')
# data_15 = np.loadtxt('Coulomb_15.txt')
# data_16 = np.loadtxt('Coulomb_16.txt')
# data_17 = np.loadtxt('Coulomb_17.txt')
# data_18 = np.loadtxt('Coulomb_18.txt')
# data_19 = np.loadtxt('Coulomb_19.txt')

# Create subplots
fig, ax = plt.subplots(1, 2, figsize=(17, 5))

# Plot for temperature

ax[0].plot(data_6[:, 0], data_6[:, 1], color='white', linestyle = 'dashed',  linewidth=1)
ax[0].plot(data_6[:, 0], data_6[:, 2], color='white', linestyle = 'dashed', linewidth=1)
ax[0].plot(data_7[:, 0], data_7[:, 1], color='orange', linewidth=1)
ax[0].plot(data_7[:, 0], data_7[:, 2], color='orange', linewidth=1)
ax[0].plot(data_8[:, 0], data_8[:, 1], color='skyblue', linewidth=1)
ax[0].plot(data_8[:, 0], data_8[:, 2], color='skyblue', linewidth=1)
ax[0].plot(data_9[:, 0], data_9[:, 1], color='fuchsia', linewidth=0.8)
ax[0].plot(data_9[:, 0], data_9[:, 2], color='fuchsia', linewidth=0.8)
ax[0].plot(data_10[:, 0], data_10[:, 1], color='lawngreen', linewidth=0.8)
ax[0].plot(data_10[:, 0], data_10[:, 2], color='lawngreen', linewidth=0.8)
ax[0].plot(data_11[:, 0], data_11[:, 1], color='red', linewidth=0.6)
ax[0].plot(data_11[:, 0], data_11[:, 2], color='red', linewidth=0.6)
# ax[0].plot(data_12[:, 0], data_12[:, 1], color='orange', linewidth=0.3)
# ax[0].plot(data_12[:, 0], data_12[:, 2], color='orange', linewidth=0.3)
# ax[0].plot(data_13[:, 0], data_13[:, 1], color='red', linewidth=0.3)
# ax[0].plot(data_13[:, 0], data_13[:, 2], color='red', linewidth=0.3)
# ax[0].plot(data_14[:, 0], data_14[:, 1], color='blue', linewidth=0.3)
# ax[0].plot(data_14[:, 0], data_14[:, 2], color='blue', linewidth=0.3)
# ax[0].plot(data_15[:, 0], data_15[:, 1], color='cyan', linewidth=0.3)
# ax[0].plot(data_15[:, 0], data_15[:, 2], color='cyan', linewidth=0.3)
# ax[0].plot(data_16[:, 0], data_16[:, 1], color='yellow')
# ax[0].plot(data_16[:, 0], data_16[:, 2], color='yellow')
# ax[0].plot(data_17[:, 0], data_17[:, 1], color='orange', linewidth=0.3)
# ax[0].plot(data_17[:, 0], data_17[:, 2], color='orange', linewidth=0.3)
# ax[0].plot(data_18[:, 0], data_18[:, 1], color='purple', linewidth=0.3)
# ax[0].plot(data_18[:, 0], data_18[:, 2], color='purple', linewidth=0.3)
# ax[0].plot(data_19[:, 0], data_19[:, 1], color='blue', linewidth=0.3)
# ax[0].plot(data_19[:, 0], data_19[:, 2], color='blue', linewidth=0.3)
ax[0].set_xlabel('time (s)')
ax[0].set_ylabel('$T_e, T_i$ (eV)')
ax[0].set_xlim(0, 1.5e-5)
ax[0].set_ylim(0, 1200)
#ax[0].set_xscale('log')





ax[1].plot(data_6[:, 0], (data_6[:, 1]-data_6[:, 2])/(data_6[0, 1]-data_6[0, 2]), color='white', linewidth=0.6,label='$\\Delta t = 10^{-6}$ s')
ax[1].plot(data_7[:, 0], (data_7[:, 1]-data_7[:, 2])/(data_7[0, 1]-data_7[0, 2]), color='orange', linewidth=1, label='$\\Delta t = 10^{-7}$ s')
ax[1].plot(data_8[:, 0], (data_8[:, 1]-data_8[:, 2])/(data_8[0, 1]-data_8[0, 2]), color='skyblue', linewidth=1,label='$\\Delta t = 10^{-8}$ s')
ax[1].plot(data_9[:, 0], (data_9[:, 1]-data_9[:, 2])/(data_9[0, 1]-data_9[0, 2]), color='fuchsia', linewidth=0.8, label='$\\Delta t = 10^{-9}$ s')
ax[1].plot(data_10[:, 0], (data_10[:, 1]-data_10[:, 2])/(data_10[0, 1]-data_10[0, 2]), color='lawngreen', linewidth=0.8, label='$\\Delta t = 10^{-10}$ s')
ax[1].plot(data_11[:, 0], (data_11[:, 1]-data_11[:, 2])/(data_11[0, 1]-data_11[0, 2]), color='red', linewidth=0.6, label='$\\Delta t = 10^{-11}$ s')
ax[1].plot(data_11[:,0], np.exp(-1.566642671644375*data_11[:,0]/1.7166122553395085e-06), color='yellow', linewidth=0.6, linestyle='dashed', label='theroretical')
# Plot for velocity ratio
# ax[1].plot(data_6[:, 0], data_6[:, 3], color='white', linestyle = 'dashed', linewidth=1, label='$\\Delta t = 10^{-6}$ s')
# ax[1].plot(data_6[:, 0], data_6[:, 4], color='white', linestyle = 'dashed', linewidth=1)
# ax[1].plot(data_7[:, 0], data_7[:, 3], color='orange', linewidth=1, label='$\\Delta t = 10^{-7}$ s')
# ax[1].plot(data_7[:, 0], data_7[:, 4], color='orange', linewidth=1)
# ax[1].plot(data_8[:, 0], data_8[:, 3], color='skyblue', linewidth=1, label='$\\Delta t = 10^{-8}$ s')
# ax[1].plot(data_8[:, 0], data_8[:, 4], color='skyblue', linewidth=1)
# ax[1].plot(data_9[:, 0], data_9[:, 3], color='fuchsia', linewidth=0.8, label='$\\Delta t = 10^{-9}$ s')
# ax[1].plot(data_9[:, 0], data_9[:, 4], color='fuchsia', linewidth=0.8)
# ax[1].plot(data_10[:, 0], data_10[:, 3], color='lawngreen', linewidth=0.8, label='$\\Delta t = 10^{-10}$ s')
# ax[1].plot(data_10[:, 0], data_10[:, 4], color='lawngreen', linewidth=0.8)
# ax[1].plot(data_11[:, 0], data_11[:, 3], color='red', linewidth=0.8, label='$\\Delta t = 10^{-11}$ s')
# ax[1].plot(data_11[:, 0], data_11[:, 4], color='red', linewidth=0.8)
# ax[1].plot(data_12[:, 0], data_12[:, 3], color='orange', linewidth=0.3, label='$\\Delta t = 10^{-12}$ s')
# ax[1].plot(data_12[:, 0], data_12[:, 4], color='orange', linewidth=0.3)
# ax[1].plot(data_13[:, 0], data_13[:, 3], color='red', linewidth=0.3, label='$\\Delta t = 10^{-13}$ s')
# ax[1].plot(data_13[:, 0], data_13[:, 4], color='red', linewidth=0.3)
# ax[1].plot(data_14[:, 0], data_14[:, 3], color='blue', linewidth=0.3, label='$\\Delta t = 10^{-14}$ s')
# ax[1].plot(data_14[:, 0], data_14[:, 4], color='blue', linewidth=0.3)
# ax[1].plot(data_15[:, 0], data_15[:, 3], color='cyan', linewidth=0.3, label='$\\Delta t = 10^{-15}$ s')
# ax[1].plot(data_15[:, 0], data_15[:, 4], color='cyan', linewidth=0.3)
# ax[1].plot(data_16[:, 0], data_16[:, 3], color='yellow', label='$\\Delta t = 10^{-16}$ s')
# ax[1].plot(data_16[:, 0], data_16[:, 4], color='yellow')
# ax[1].plot(data_17[:, 0], data_17[:, 3], color='orange', linewidth=0.3, label='$\\Delta t = 10^{-17}$ s')
# ax[1].plot(data_17[:, 0], data_17[:, 4], color='orange', linewidth=0.3)
# ax[1].plot(data_18[:, 0], data_18[:, 3], color='purple', linewidth=0.3, label='$\\Delta t = 10^{-18}$ s')
# ax[1].plot(data_18[:, 0], data_18[:, 4], color='purple', linewidth=0.3)
# ax[1].plot(data_19[:, 0], data_19[:, 3], color='blue', linewidth=0.3, label='$\\Delta t = 10^{-19}$ s')
# ax[1].plot(data_19[:, 0], data_19[:, 4], color='blue', linewidth=0.3)

ax[1].set_xlabel('time (s)')
ax[1].set_ylabel('$V/V_{0_{elec}}$')
# ax[1].set_xlim(0, 1.5e-5)
# ax[1].set_ylim(0, 1.2)
#ax[1].set_xscale('log')

# Collect labels from both subplots
handles, labels = [], []
for a in ax:
    h, l = a.get_legend_handles_labels()
    handles.extend(h)
    labels.extend(l)

# Adjust layout first
plt.tight_layout()

fig.legend(handles, labels, fontsize=15, loc='center right')

# Adjust margins to make space for the legend
plt.subplots_adjust(right=0.85)  # Increase if necessary
# Show the plot
plt.show()
