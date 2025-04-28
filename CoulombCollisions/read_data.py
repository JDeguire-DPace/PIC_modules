import numpy as np
from matplotlib import pyplot as plt

# Load data
data_5 = np.loadtxt('CeciEstUnTest5.txt')
data_6 = np.loadtxt('CeciEstUnTest6.txt')
data_7 = np.loadtxt('CeciEstUnTest7.txt')
data_8 = np.loadtxt('CeciEstUnTest8.txt')
data_9 = np.loadtxt('CeciEstUnTest9.txt')
data_10 = np.loadtxt('CeciEstUnTest10.txt')



# Create subplots
fig, ax = plt.subplots(1, 2, figsize=(17, 5))

# Plot for temperature
ax[0].plot(data_5[:, 0], data_5[:, 1], '--', color='purple', linewidth=0.8)
ax[0].plot(data_5[:, 0], data_5[:, 2], '--', color='purple', linewidth=0.8)
ax[0].plot(data_6[:, 0], data_6[:, 1], color='#D81B60', linewidth=0.8)
ax[0].plot(data_6[:, 0], data_6[:, 2], color='#D81B60', linewidth=0.8)
ax[0].plot(data_7[:, 0], data_7[:, 1], '--', color='#1E88E5', linewidth=0.8)
ax[0].plot(data_7[:, 0], data_7[:, 2], '--', color='#1E88E5', linewidth=0.8)
ax[0].plot(data_8[:, 0], data_8[:, 1], color='#FFC107', linewidth=0.8)
ax[0].plot(data_8[:, 0], data_8[:, 2], color='#FFC107', linewidth=0.8)
ax[0].plot(data_9[:, 0], data_9[:, 1], color='#004D40', linewidth=0.5)
ax[0].plot(data_9[:, 0], data_9[:, 2], color='#004D40', linewidth=0.5)
ax[0].plot(data_10[:, 0], data_10[:, 1], '--', color='lawngreen', linewidth=0.3)
ax[0].plot(data_10[:, 0], data_10[:, 2], '--', color='lawngreen', linewidth=0.3)
ax[0].set_xlabel('time (s)')
ax[0].set_ylabel('$k_B T$ (eV)')
ax[0].set_xlim(0, 5e-5)

# Plot for velocity ratio
ax[1].plot(data_5[:, 0], data_5[:, 3], '--', color='purple', linewidth=0.8, label='$\\Delta t = 10^{-5}$ s')
ax[1].plot(data_5[:, 0], data_5[:, 4], '--', color='purple', linewidth=0.8)
ax[1].plot(data_6[:, 0], data_6[:, 3], color='#D81B60', linewidth=0.8, label='$\\Delta t = 10^{-6}$ s')
ax[1].plot(data_6[:, 0], data_6[:, 4], color='#D81B60', linewidth=0.8)
ax[1].plot(data_7[:, 0], data_7[:, 3], '--', color='#1E88E5', linewidth=0.8, label='$\\Delta t = 10^{-7}$ s')
ax[1].plot(data_7[:, 0], data_7[:, 4], '--', color='#1E88E5', linewidth=0.8)
ax[1].plot(data_8[:, 0], data_8[:, 3], color='#FFC107', linewidth=0.8, label='$\\Delta t = 10^{-8}$ s')
ax[1].plot(data_8[:, 0], data_8[:, 4], color='#FFC107', linewidth=0.8)
ax[1].plot(data_9[:, 0], data_9[:, 3], color='#004D40', linewidth=0.5, label='$\\Delta t = 10^{-9}$ s')
ax[1].plot(data_9[:, 0], data_9[:, 4], color='#004D40', linewidth=0.5)
ax[1].plot(data_10[:, 0], data_10[:, 3], '--', color='lawngreen', linewidth=0.3, label='$\\Delta t = 10^{-10}$ s')
ax[1].plot(data_10[:, 0], data_10[:, 4], '--', color='lawngreen', linewidth=0.3)
ax[1].set_xlabel('time (s)')
ax[1].set_ylabel('$V/V_{0_{elec}}$')
ax[1].set_xlim(0, 5e-5)

# Collect labels from both subplots
handles, labels = [], []
for a in ax:
    h, l = a.get_legend_handles_labels()
    handles.extend(h)
    labels.extend(l)

# Adjust layout first
plt.tight_layout()

# Add the legend after adjusting layout
fig.legend(handles, labels, fontsize=15, loc='center right')

# Adjust margins to make space for the legend
plt.subplots_adjust(right=0.85)  # Increase if necessary

# Show the plot
plt.show()
