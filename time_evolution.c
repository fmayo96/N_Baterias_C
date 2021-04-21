#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrixoperations.h"
#include "RK4.h"
#include "propagators.h"
#include "energy.h"
void Open_evolution(double complex *state, double complex *hamiltonian, double complex *bath_state, double complex *interaction, int dim, double tf, double dt, double *E)
{
    int N = (int)(tf/dt), step;
    double complex *dissipator;
    dissipator = (double complex*) calloc(dim, sizeof(double complex));
    *E  = calc_Energy(state, hamiltonian, dim);
    for(step = 1; step < N; step++)
    {
        RK4_open(dissipator, state, hamiltonian, bath_state, interaction, dt, dim);
        *(E + step) = calc_Energy(state, hamiltonian, dim) - *E;
        printf("\r Progress (time evolution): %.2lf", (double)(step + 1)/N *100);
    }
    *E = 0;
    free(dissipator);
}
