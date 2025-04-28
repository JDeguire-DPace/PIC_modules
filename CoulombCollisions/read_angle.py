import numpy as np
from matplotlib import pyplot as plt

# Load data
data_5 = np.loadtxt('Angle10_5.txt')
data_6 = np.loadtxt('Angle10_6.txt')
data_7 = np.loadtxt('Angle10_7.txt')
data_8 = np.loadtxt('Angle10_8.txt')
data_9 = np.loadtxt('Angle10_9.txt')
data_10 = np.loadtxt('Angle10_10.txt')

angle_5_ie = data_5[:,5]
angle_6_ie = data_6[:,5]
angle_7_ie = data_7[:,5]
angle_8_ie = data_8[:,5]
# angle_9_ie = data_9[:,5]
# angle_10_ie = data_10[:,5]

angle_5_ii = data_5[:,6]
angle_6_ii = data_6[:,6]
angle_7_ii = data_7[:,6]
angle_8_ii = data_8[:,6]
# angle_9_ii = data_9[:,6]
# angle_10_ii = data_10[:,6]

angle_5_ee = data_5[:,7]
angle_6_ee = data_6[:,7]
angle_7_ee = data_7[:,7]
angle_8_ee = data_8[:,7]
# angle_9_ee = data_9[:,7]
# angle_10_ee = data_10[:,7]


angles_ie=[np.mean(angle_5_ie),np.mean(angle_6_ie),np.mean(angle_7_ie),np.mean(angle_8_ie)]#,np.mean(angle_9_ie),np.mean(angle_10_ie)]
angles_ii=[np.mean(angle_5_ii),np.mean(angle_6_ii),np.mean(angle_7_ii),np.mean(angle_8_ii)]#,np.mean(angle_9_ii),np.mean(angle_10_ii)]
angles_ee=[np.mean(angle_5_ee),np.mean(angle_6_ee),np.mean(angle_7_ee),np.mean(angle_8_ee)]#,np.mean(angle_9_ee),np.mean(angle_10_ee)]
t=[10**-5,10**-6,10**-7,10**-8]#,10**-9,10**-10]


# plt.scatter(10**-6,np.mean(angle_6), color='red', s=14, label='$\\Delta t = 10^{-6}$ s')
# plt.scatter( 10**-7,np.mean(angle_7), color='blue', s=14, label='$\\Delta t = 10^{-7}$ s')
# plt.scatter(10**-8, np.mean(angle_8), color='green', s=14, label='$\\Delta t = 10^{-8}$ s')
# plt.scatter(10**-9, np.mean(data_9[:,5]), color='black', s=14, label='$\\Delta t = 10^{-9}$ s')
plt.figure(figsize=(8, 5))
plt.plot(t, angles_ie, '--', marker='X', linewidth=0.8, color='black', markersize=8, label='electrons-ions' , markerfacecolor='red')
plt.plot(t, angles_ii, '--', marker='X', linewidth=0.8, color='black', markersize=8, label='ions-ions' , markerfacecolor='blue')
plt.plot(t, angles_ee, '--', marker='X', linewidth=0.8, color='black', markersize=8, label='electrons-electrons' , markerfacecolor='green')
plt.xscale('log')
plt.ylabel('Average diffusion angle $\\chi$ (rad)')
plt.xlabel('Time step $\\Delta t$ (s)') 
plt.yticks(np.linspace(0, np.pi/2, 7), ['0', 'π/12', 'π/6', 'π/4', 'π/3', '5π/12', 'π/2'])
plt.legend()
plt.savefig('Figures/DiffusionAngleDT.png', dpi=300, bbox_inches='tight')
# Average6 = np.mean(data_6[:,-1], axis=1)
# print(Average6)