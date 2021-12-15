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
double tf = 2, dt = 0.0001, h = 1.5, eps, Z, beta;
int dim, N_qubits = 6, N = (int)tf/dt, i, b, N_betas = 10;
dim = pow(2, N_qubits);
eps = N_qubits*sqrt(10);
double complex *state, *hamiltonian, *bath_state, *interaction;
state = (double complex *) calloc(dim, sizeof(double complex));
hamiltonian = (double complex *) calloc(dim, sizeof(double complex));
bath_state = (double complex *) calloc(dim, sizeof(double complex));
interaction = (double complex *) calloc(dim*dim, sizeof(double complex));
//----------Hamiltonianos y Condiciones iniciales-----------------------
Hamiltonian(hamiltonian, h, N_qubits);
Interaction(interaction, eps, N_qubits);

for(b = 0; b < N_betas; b++)
{
    beta = (double)(b + 1);
    double *E, *t;
    E = (double*) calloc(N, sizeof(double));
    t = (double*) calloc(N, sizeof(double));
    Z = 0;
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
    sprintf(filename,"Test_EQW_N = %d_beta = %d.txt", N_qubits, (int)beta);

    FILE *fp=fopen(filename,"w");

    for(i = 0; i < N; i++)
    {
        fprintf(fp, "%lf %lf \n", *(t + i),*(E + i));
    }

        fclose(fp);
        free(E);
        free(t);
}
    free(state);
    free(bath_state);
    free(hamiltonian);
    free(interaction);
    
    return 0;
}
