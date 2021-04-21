#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H
void Open_evolution(double complex *state, double complex *hamiltonian, double complex *bath_state, double complex *interaction, int dim, double tf, double dt, double *E);
#endif