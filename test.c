#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrixoperations.h"
#include "propagators.h"
#include "RK4.h"
#include "time_evolution.h"
#include "energy.h"
#include "pauli.h"
#include "hamiltonian.h"
int main()
{
//-----------Declaración de variblaes--------------------------------
double tf = 2, dt = 0.0001, h = 1.5, beta = 50, eps, Z = 0;
int dim, N_qubits = 13, N = (int)tf/dt, i;
dim = pow(2, N_qubits);
eps = N_qubits*sqrt(10);
double complex *state, *hamiltonian, *bath_state, *interaction;
double *E, *t; 
state = (double complex *) calloc(dim, sizeof(double complex));
hamiltonian = (double complex *) calloc(dim, sizeof(double complex));
bath_state = (double complex *) calloc(dim, sizeof(double complex));
interaction = (double complex *) calloc(dim*dim, sizeof(double complex));
E = (double*) calloc(N, sizeof(double));
t = (double*) calloc(N, sizeof(double));
//----------Hamiltonianos y Condiciones iniciales-----------------------
Hamiltonian(hamiltonian, h, N_qubits);
Interaction(interaction, eps, N_qubits);

for(i = 0; i< N; i++)
{
    *(t + i) = i*dt;
}
for(i = 0; i < dim; i++)
{
    *(bath_state + i) = exp(-beta* *(hamiltonian + i));
}
for(i = 0; i < dim; i++)
{
    Z += *(bath_state + i);
}
for(i = 0; i < dim; i++)
{
    *(bath_state + i) = *(bath_state + i)/Z;
    *(state + i) = *(bath_state + i);
}
Open_evolution(state, hamiltonian, bath_state, interaction, dim, tf, dt, E);

char filename[255];
sprintf(filename,"Test_EQW_N = %d.txt", N_qubits);

FILE *fp=fopen(filename,"w");

for(i = 0; i < N; i++)
{
    fprintf(fp, "%lf %lf \n", *(t + i),*(E + i));
}

    fclose(fp);
    free(state);
    free(bath_state);
    free(hamiltonian);
    free(interaction);
    free(E);
    return 0;
}
