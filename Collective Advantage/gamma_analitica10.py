from cProfile import label
from telnetlib import GA
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special


def E(N,beta, omega, eps, t):
    E = 0
    for i in range(0,N+1):
        tau_i = (2*np.cosh(beta*omega/2))**N / (eps**2*N**2*2*np.cosh((i - N/2)*beta*omega))
        E += scipy.special.binom(N,i)*(2*np.sinh(beta*omega*(i - N/2))/(2*np.cosh(beta*omega/2))**N * (1 - np.exp(-t/tau_i)) * omega*((i - N/2)))
    return E 

N_betas = 500
Nt = 150000
omega = 1.5
eps = np.sqrt(10)
betas = np.linspace(1,50,N_betas)
t,dt = np.linspace(0,1,Nt,retstep=True)

Energy_N2 = np.zeros([N_betas,Nt])
Gamma_N2 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N2[i,j] = E(10,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga2 = 0
    for j in range(Nt):
        if Energy_N2[i,j] < Energy_N2[i,-1]*0.99:
            jcarga2 = j 
        else:
            break 
    tcarga2 = jcarga2*dt 
    tparall = 1/eps**2 * np.log((Energy_N2[i,-1] - Energy_N2[i,0])/(0.01*Energy_N2[i,-1]))
    Gamma_N2[i] = tparall/tcarga2 

Energy_N3 = np.zeros([N_betas,Nt])
Gamma_N3 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N3[i,j] = E(20,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga3 = 0
    for j in range(Nt):
        if Energy_N3[i,j] < Energy_N3[i,-1]*0.99:
            jcarga3 = j 
        else:
            break 
    tcarga3 = jcarga3*dt 
    tparall = 1/eps**2 * np.log((Energy_N3[i,-1] - Energy_N3[i,0])/(0.01*Energy_N3[i,-1]))
    Gamma_N3[i] = tparall/tcarga3 

Energy_N4 = np.zeros([N_betas,Nt])
Gamma_N4 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N4[i,j] = E(30,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga4 = 0
    for j in range(Nt):
        if Energy_N4[i,j] < Energy_N4[i,-1]*0.99:
            jcarga4 = j 
        else:
            break 
    tcarga4 = jcarga4*dt 
    tparall = 1/eps**2 * np.log((Energy_N4[i,-1] - Energy_N4[i,0])/(0.01*Energy_N4[i,-1]))
    Gamma_N4[i] = tparall/tcarga4 

Energy_N5 = np.zeros([N_betas,Nt])
Gamma_N5 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N5[i,j] = E(40,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga5 = 0
    for j in range(Nt):
        if Energy_N5[i,j] < Energy_N5[i,-1]*0.99:
            jcarga5 = j 
        else:
            break 
    tcarga5 = jcarga5*dt 
    tparall = 1/eps**2 * np.log((Energy_N5[i,-1] - Energy_N5[i,0])/(0.01*Energy_N5[i,-1]))
    Gamma_N5[i] = tparall/tcarga5

Energy_N6 = np.zeros([N_betas,Nt])
Gamma_N6 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N6[i,j] = E(50,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga6 = 0
    for j in range(Nt):
        if Energy_N6[i,j] < Energy_N6[i,-1]*0.99:
            jcarga6 = j 
        else:
            break 
    tcarga6 = jcarga6*dt 
    tparall = 1/eps**2 * np.log((Energy_N6[i,-1] - Energy_N6[i,0])/(0.01*Energy_N6[i,-1]))
    Gamma_N6[i] = tparall/tcarga6

Energy_N7 = np.zeros([N_betas,Nt])
Gamma_N7 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N7[i,j] = E(60,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga7 = 0
    for j in range(Nt):
        if Energy_N7[i,j] < Energy_N7[i,-1]*0.99:
            jcarga7 = j 
        else:
            break 
    tcarga7 = jcarga7*dt 
    tparall = 1/eps**2 * np.log((Energy_N7[i,-1] - Energy_N7[i,0])/(0.01*Energy_N7[i,-1]))
    Gamma_N7[i] = tparall/tcarga7

Energy_N8 = np.zeros([N_betas,Nt])
Gamma_N8 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N8[i,j] = E(70,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga8 = 0
    for j in range(Nt):
        if Energy_N8[i,j] < Energy_N8[i,-1]*0.99:
            jcarga8 = j 
        else:
            break 
    tcarga8 = jcarga8*dt 
    tparall = 1/eps**2 * np.log((Energy_N8[i,-1] - Energy_N8[i,0])/(0.01*Energy_N8[i,-1]))
    Gamma_N8[i] = tparall/tcarga8

Energy_N9 = np.zeros([N_betas,Nt])
Gamma_N9 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N9[i,j] = E(80,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga9 = 0
    for j in range(Nt):
        if Energy_N9[i,j] < Energy_N9[i,-1]*0.99:
            jcarga9 = j 
        else:
            break 
    tcarga9 = jcarga9*dt 
    tparall = 1/eps**2 * np.log((Energy_N9[i,-1] - Energy_N9[i,0])/(0.01*Energy_N9[i,-1]))
    Gamma_N9[i] = tparall/tcarga9


Energy_N10 = np.zeros([N_betas,Nt])
Gamma_N10 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N10[i,j] = E(90,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga10 = 0
    for j in range(Nt):
        if Energy_N10[i,j] < Energy_N10[i,-1]*0.99:
            jcarga10 = j 
        else:
            break 
    tcarga10 = jcarga10*dt 
    tparall = 1/eps**2 * np.log((Energy_N10[i,-1] - Energy_N10[i,0])/(0.01*Energy_N10[i,-1]))
    Gamma_N10[i] = tparall/tcarga10

Energy_N11 = np.zeros([N_betas,Nt])
Gamma_N11 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N11[i,j] = E(100,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga11 = 0
    for j in range(Nt):
        if Energy_N11[i,j] < Energy_N11[i,-1]*0.99:
            jcarga11 = j 
        else:
            break 
    tcarga11 = jcarga11*dt 
    tparall = 1/eps**2 * np.log((Energy_N11[i,-1] - Energy_N11[i,0])/(0.01*Energy_N11[i,-1]))
    Gamma_N11[i] = tparall/tcarga11


Energy_N12 = np.zeros([N_betas,Nt])
Gamma_N12 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N12[i,j] = E(110,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga12 = 0
    for j in range(Nt):
        if Energy_N12[i,j] < Energy_N12[i,-1]*0.99:
            jcarga12 = j 
        else:
            break 
    tcarga12 = jcarga12*dt 
    tparall = 1/eps**2 * np.log((Energy_N12[i,-1] - Energy_N12[i,0])/(0.01*Energy_N12[i,-1]))
    Gamma_N12[i] = tparall/tcarga12


Energy_N13 = np.zeros([N_betas,Nt])
Gamma_N13 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N13[i,j] = E(120,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga13 = 0
    for j in range(Nt):
        if Energy_N13[i,j] < Energy_N13[i,-1]*0.99:
            jcarga13 = j 
        else:
            break 
    tcarga13 = jcarga13*dt 
    tparall = 1/eps**2 * np.log((Energy_N13[i,-1] - Energy_N13[i,0])/(0.01*Energy_N13[i,-1]))
    Gamma_N13[i] = tparall/tcarga13


Energy_N14 = np.zeros([N_betas,Nt])
Gamma_N14 = np.zeros(N_betas)

for i in range(N_betas):
    for j in range(Nt):
        Energy_N14[i,j] = E(130,betas[i], omega, eps, t[j])

for i in range(N_betas):
    jcarga14 = 0
    for j in range(Nt):
        if Energy_N14[i,j] < Energy_N14[i,-1]*0.99:
            jcarga14 = j 
        else:
            break 
    tcarga14 = jcarga14*dt 
    tparall = 1/eps**2 * np.log((Energy_N14[i,-1] - Energy_N14[i,0])/(0.01*Energy_N14[i,-1]))
    Gamma_N14[i] = tparall/tcarga14


betas_cr = []
beta_cr = 0
for i in range(N_betas):
    if Gamma_N2[i] <0.99*Gamma_N2[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N3[i] <0.99*Gamma_N3[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N4[i] <0.99*Gamma_N4[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N5[i] <0.99*Gamma_N5[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N6[i] <0.99*Gamma_N6[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N7[i] <0.99*Gamma_N7[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N8[i] <0.99*Gamma_N8[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N9[i] <0.99*Gamma_N9[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N10[i] <0.99*Gamma_N10[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N11[i] <0.99*Gamma_N11[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N12[i] <0.99*Gamma_N12[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N13[i] <0.99*Gamma_N13[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)

beta_cr = 0
for i in range(N_betas):
    if Gamma_N14[i] <0.99*Gamma_N14[-1]:
        beta_cr = betas[i]
    else:
        break 
betas_cr.append(beta_cr)


print(betas_cr)