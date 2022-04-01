import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special

def E(N,beta, omega, t, eps):
    E = 0
    for i in range(0,N+1):
        tau_i = (2*np.cosh(beta*omega/2))**N / (eps**2*N**2*2*np.cosh((i - N/2)*beta*omega))
        E += scipy.special.binom(N,i)*(2*np.sinh(beta*omega*(i - N/2))/(2*np.cosh(beta*omega/2))**N * (1 - np.exp(-t/tau_i)) * omega*((i - N/2)))
    return E 

t = np.linspace(0,5,10000)
N = 4
omega = 1.5
beta = np.linspace(1,10,10)
Energy = np.zeros([len(beta), len(t)])
for i in range(len(beta)):
    for j in range(len(t)):
        Energy[i,j] = E(N,beta[i], omega, t[j], np.sqrt(10))

plt.figure()
plt.plot(t, Energy[-1,:])
plt.show()
