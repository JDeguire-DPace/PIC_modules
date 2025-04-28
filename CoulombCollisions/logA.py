import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
charge_electron = 1.602e-19
mass_electron = 9.11e-31
epsilon_0 = 8.854e-12

def logA(eV,density):
    return 4*pi*(epsilon_0*eV)**1.5/(density**.5*charge_electron**2)

color = ['red','blue','green','orange','purple','pink']
energy = []
for i in range (5,11):
    k = []
    for j in range(-3,3):
        k.append((logA(10**j,10**i)))
    plt.scatter(range(-3,3),k,color=color[i-5],label=f'$10^{i}$')

plt.legend()
plt.yscale('log')
plt.show()

