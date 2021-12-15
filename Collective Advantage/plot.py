import numpy as np
import matplotlib.pyplot as plt 


beta_max = 10
P1_betas = []
P2_betas = []
P3_betas = []
P4_betas = []
for b in range(1,beta_max+1):
    data1 = np.loadtxt("Test_EQW_N = 1_beta = " + str(b) + ".txt")
    data2 = np.loadtxt("Test_EQW_N = 2_beta = " + str(b) + ".txt")
    data3 = np.loadtxt("Test_EQW_N = 5_beta = " + str(b) + ".txt")
    data4 = np.loadtxt("Test_EQW_N = 6_beta = " + str(b) + ".txt")



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
        if E1[i] < E1[-1]*0.95:
            tcarga1 = i 
        else:
            break 
    for i in range(N):
        if E2[i] < E2[-1]*0.95:
            tcarga2 = i 
        else:
            break 
    for i in range(N):
        if E3[i] < E3[-1]*0.95:
            tcarga3 = i 
        else:
            break 
    for i in range(N):
        if E4[i] < E4[-1]*0.95:
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

betas = np.arange(1,beta_max+1)
x = np.linspace(0,15,1000)
plt.figure()
plt.plot(betas*1.5, P2_betas/(2*P1_betas), '.', label = r'$N= 2$')
plt.plot(x, 2*(np.cosh(x)*4/(2*np.cosh(x/2))**2), '--C0')
plt.plot(betas*1.5, P3_betas/(5*P1_betas), '.', label = r'$N= 3$')
plt.plot(x, 2*(np.cosh(5*x/2)*25/(2*np.cosh(x/2))**5), '--C1')
plt.plot(betas*1.5, P4_betas/(6*P1_betas), '.', label = r'$N= 4$')
plt.plot(x, 2*(np.cosh(3*x)*36/(2*np.cosh(x/2))**6), '--C2')
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
