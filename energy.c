#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrixoperations.h"

double calc_Energy(double complex *state, double complex *hamiltonian, int dim)
{
    double E;
    double complex *state_hamiltonian;
    state_hamiltonian = (double complex *) calloc(dim, sizeof(double complex));
    dot_diag_diag(state_hamiltonian, state, hamiltonian, dim);
    E = trace_diag(state_hamiltonian, dim);
    free(state_hamiltonian);
    return E;
}
