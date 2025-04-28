import numpy as np
import matplotlib.pyplot as plt

Number_of_particles = 200000
Number_collisions = 10000

def f(x):
    if x==0:
        return 0
    else:
        premier = x**2*(np.log((x**2+4)/x**2))
        deuxieme = (x**2+4)*(np.atan(x/2))**2 - (x*np.pi/2)**2
        troiseme = 4*x*np.atan(x/2)
        return premier + deuxieme + troiseme


angle_min = np.linspace(0, 3*np.pi/180, 1000)
func = []
for i in angle_min:
    func.append(f(i))



angle_phi = 2*np.pi*np.random.rand(Number_collisions)
angle_theta = np.pi*np.random.rand(Number_collisions)
somme = 0
for k in range(Number_collisions):
    for l in range(0,k-1):
        phi_k = angle_phi[k]
        phi_l = angle_phi[l]
        theta_k = angle_theta[k]
        theta_l = angle_theta[l]
        somme += 1/2*theta_l*theta_k*np.cos(phi_k - phi_l)/Number_collisions
for k in range(Number_collisions):
    somme+= angle_theta[k]**2/4

print(somme)
#plt.plot(angle_min, func)
#plt.show()