import numpy as np
import matplotlib.pyplot as plt 
import scipy.special

beta_max = 10
P1_betas = []
P2_betas = []
P3_betas = []
P4_betas = []
for b in range(1,beta_max+1):
    data1 = np.loadtxt("Test_EQW_N = 1_beta = " + str(b) + ".txt")
    data2 = np.loadtxt("Test_EQW_N = 2_beta = " + str(b) + ".txt")
    data3 = np.loadtxt("Test_EQW_N = 3_beta = " + str(b) + ".txt")
    data4 = np.loadtxt("Test_EQW_N = 4_beta = " + str(b) + ".txt")



    t = data1[:,0]
    E1 = data1[:,1]
    E2 = data2[:,1]
    E3 = data3[:,1]
    E4 = data4[:,1]

    N = len(t)
    dt = t[1] - t[0]
    tcarga1 = 0
    tcarga2 = 0
    tcarga3 = 0
    tcarga4 = 0
    for i in range(N):
        if E1[i] < E1[-1]*0.99:
            tcarga1 = i 
        else:
            break 
    for i in range(N):
        if E2[i] < E2[-1]*0.99:
            tcarga2 = i 
        else:
            break 
    for i in range(N):
        if E3[i] < E3[-1]*0.99:
            tcarga3 = i 
        else:
            break 
    for i in range(N):
        if E4[i] < E4[-1]*0.99:
            tcarga4 = i 
        else:
            break 

    P1_betas.append(E1[-1]/(tcarga1*dt)) 
    P2_betas.append(E2[-1]/(tcarga2*dt))
    P3_betas.append(E3[-1]/(tcarga3*dt))
    P4_betas.append(E4[-1]/(tcarga4*dt))
    
P1_betas = np.array(P1_betas)
P2_betas = np.array(P2_betas)
P3_betas = np.array(P3_betas)
P4_betas = np.array(P4_betas)




def E(N,beta, omega, t):
    E = 0
    for i in range(0,N+1):
        tau_i = (2*np.cosh(beta*omega/2))**N / (10*N**2*2*np.cosh((i - N/2)*beta*omega))
        E += scipy.special.binom(N,i)*(2*np.sinh(beta*omega*(i - N/2))/(2*np.cosh(beta*omega/2))**N * (1 - np.exp(-t/tau_i)) * omega*((i - N/2)))
    return E 

data2 = np.loadtxt("Test_EQW_N = 4_beta = " + str(1) + ".txt")
E2 = data2[:,1]
data1 = np.loadtxt("Test_EQW_N = 1_beta = " + str(1) + ".txt")
E1 = data1[:,1]
En = E(4, 1, 1.5, t)
plt.figure()
plt.plot(t,En)
plt.plot(t,E2, '--')
plt.plot(t, 4*E1)
plt.show()


betas = np.arange(1,beta_max+1)
x = np.linspace(6,15,1000)
plt.figure()
plt.plot(betas*1.5, P2_betas/(2*P1_betas), 'o', label = r'$N= 2$')
plt.plot(x, np.ones(len(x))*2**2, '--C0')
plt.plot(betas*1.5, P3_betas/(3*P1_betas), 'o', label = r'$N= 3$')
plt.plot(x, np.ones(len(x))*3**2, '--C1')
#plt.plot(betas*1.5, P4_betas/(4*P1_betas), '.', label = r'$N= 4$')
#plt.plot(x, np.ones(len(x))*4**2, '--C2')
plt.legend()
plt.xlabel(r'$\beta\omega$')
plt.ylabel(r'$\Gamma$')
plt.show()


N = np.array([1,2,3,4])
P = np.array([P1_betas[-1], P2_betas[-1],P3_betas[-1],P4_betas[-1]])
x = np.linspace(0,np.log(4),1000)
plt.figure()
plt.plot(np.log(N), np.log(P/P[0]), '.')
plt.plot(x,x*3, 'r')
plt.show()
