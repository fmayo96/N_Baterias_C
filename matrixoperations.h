#ifndef MATRIXOPERATIONS_H
#define MATRIXOPERATIONS_H
void dot_diag_diag(double complex *out, double complex *A, double complex *B, int dim);
void dot_antidiag_antidiag(double complex *out, double complex *A, double complex *B, int dim);
void dot_diag_antidiag(double complex *out, double complex *A, double complex *B, int dim);
void diag_antidiag_diag(double complex *out, double complex *A, double complex *B, int dim);
void kron_diag_diag(double complex *out, double complex *A, double complex *B, int dim);
double trace_diag(double complex *A, int dim);
void partial_trace_diag(double complex *out, double complex *A, int dim);
void commutator_antidiag_diag(double complex *out, double complex *A, double complex *B, int dim);
void commutator_antidiag_antidiag(double complex *out, double complex *A, double complex *B, int dim);
#endif
