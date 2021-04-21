#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrixoperations.h"
void Dissipator(double complex *dissipator, double complex *state, double complex *bath_state, double complex *interaction, int dim)
{   
    int i;
    double complex *rho, *comm_v_rho, *double_comm;
    rho = (double complex*)calloc(dim*dim, sizeof(double complex));
    comm_v_rho = (double complex*)calloc(dim*dim, sizeof(double complex));
    double_comm = (double complex*)calloc(dim*dim, sizeof(double complex));
    kron_diag_diag(rho, state, bath_state, dim);
    commutator_antidiag_diag(comm_v_rho, interaction, rho, dim*dim);
    commutator_antidiag_antidiag(double_comm, interaction, comm_v_rho, dim*dim);
    partial_trace_diag(dissipator, double_comm, dim);
    for(i = 0; i < dim; i++)
    { 
        *(dissipator + i) = -*(dissipator + i) * 0.5;
    }
    free(rho);
    free(comm_v_rho);
    free(double_comm);
}
