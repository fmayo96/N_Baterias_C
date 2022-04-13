from cProfile import label
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit


beta_cr = np.array([4.211055276381909, 4.572864321608041, 4.844221105527638, 5.07035175879397, 5.206030150753769, 5.341708542713568, 5.477386934673367, 5.613065326633166, 5.703517587939698, 5.793969849246231, 5.839195979899498, 5.9296482412060305])
Tc = 1/beta_cr
N = [3,4,5,6,7,8,9,10,11,12,13,14]
x = np.linspace(N[0], N[-1], 1000)
dbeta = 1/20
errTc = Tc**2 * dbeta


def f(x,a,b):
    return a*x + b

popt, pcov = curve_fit(f,np.log(N[1:]), np.log(Tc[1:]))

print(popt)
a = popt[0]
b = popt[1]

perr = np.sqrt(np.diag(pcov))
ymean = np.mean(np.log(Tc[1:]))




perr=np.sqrt(np.diag(pcov))


SSTOT=np.zeros(len(Tc[1:]))
SSRES=np.zeros(len(Tc[1:]))

SSTOT = (np.log(Tc[1:]) - ymean*np.ones(len(np.log(Tc[1:]))) ) ** 2
SSRES = (np.log(Tc[1:]) - f(np.log(N[1:]), a,b) ) ** 2

sstot=np.sum(SSTOT)
ssres=np.sum(SSRES)

Rcuadrado = 1 - (ssres/sstot)

print(r"$R^2 = $" + str(Rcuadrado))

plt.figure()
plt.errorbar(N,Tc, yerr=errTc, fmt = 'o', capsize=2, label = r'$T_c$')
#plt.plot(x,np.exp(popt[1])*np.exp(popt[0]*x))
plt.plot(x, np.exp(-1.25)*x**(-0.2), label = r'$N^{-0.2}$')
plt.xlabel('N')
plt.ylabel(r'$T_c$')
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize = 12)
plt.show()

