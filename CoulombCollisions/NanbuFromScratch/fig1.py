import numpy as np
from matplotlib import pyplot as plt

s=np.linspace(0,6,1000)
x=0.5*(1-np.exp(-s))

plt.plot(s,x,linestyle='dashed',linewidth=1,color='black')
plt.xlabel('s',fontsize=14)
plt.ylabel('$\\langle \\sin^2 \\frac{\\chi_N}{2} \\rangle$',fontsize=14)
plt.xlim(0,6)
plt.show()