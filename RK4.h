#ifndef RK4_H
#define RK4_H
void RK4_open(double complex *dissipator, double complex *state, double complex *hamiltonian, double complex *bath_state, double complex *interaction, double dt, int dim);
#endif
