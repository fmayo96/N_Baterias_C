#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
void dot_diag_diag(double complex *out, double complex *A, double complex *B, int dim)
{
    int i;
    for(i = 0; i < dim; i++)
    {
        *(out + i) = *(A + i) * *(B + i);
    }
}
void dot_antidiag_antidiag(double complex *out, double complex *A, double complex *B, int dim)
{
    int i;
    for(i = 0; i < dim; i++)
    {
        *(out + i) = *(A + i) * *(B + dim - 1 - i);
    }
}
void dot_diag_antidiag(double complex *out, double complex *A, double complex *B, int dim)
{
    int i;
    for(i = 0; i < dim; i++)
    {
        *(out + i) = *(A + i) * *(B + i);
    }
}
void dot_antidiag_diag(double complex *out, double complex *A, double complex *B, int dim)
{
    int i;
    for(i = 0; i < dim; i++)
    {
        *(out + i) = *(A + i) * *(B + dim - 1 - i);
    }
}
void kron_diag_diag(double complex *out, double complex *A, double complex *B, int dim)
{
    int i, j;
    for(i = 0; i < dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            *(out + dim*i + j) = *(A + i) * *(B + j);
        }
    }
}
double trace_diag(double complex *A, int dim)
{
    int i;
    double trace = 0;
    for(i = 0; i < dim; i++)
        {
            trace += creal(*(A + i));
        }
    return trace;
}
void partial_trace_diag(double complex *out, double complex *A, int dim)
{
    int i, j;
    for(i = 0; i < dim; i++)
    {
        *(out + i) = 0;
    }
    for(i = 0; i< dim; i++)
    {
        for(j = 0; j < dim; j++)
        {
            *(out + i) += *(A + dim*i + j);
        }
    }
}
void commutator_antidiag_diag(double complex *out, double complex *A, double complex *B, int dim)// result is antidiag
{
    int i;
    double complex *C, *D;
    C = (double complex *) calloc(dim, sizeof(double complex));
    D = (double complex *) calloc(dim, sizeof(double complex));
    dot_antidiag_diag(C, A, B, dim);
    dot_diag_antidiag(D, B, A, dim);
    for(i = 0; i < dim; i++)
    {
        *(out + i) = *(C + i) - *(D + i);
    }
    free(C);
    free(D);
}
void commutator_antidiag_antidiag(double complex *out, double complex *A, double complex *B, int dim)// result is diag
{
    int i;
    double complex *C, *D;
    C = (double complex *) calloc(dim, sizeof(double complex));
    D = (double complex *) calloc(dim, sizeof(double complex));
    dot_antidiag_antidiag(C, A, B, dim);
    dot_antidiag_antidiag(D, B, A, dim);
    for(i = 0; i < dim; i++)
    {
        *(out + i) = *(C + i) - *(D + i);
    }
    free(C);
    free(D);
}