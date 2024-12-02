#include "glass_utils.h"

int *int_vector(int N)
{
    int *v = malloc( N * sizeof(int) );
    return v;
}

void free_int_vector(int *v)
{
    free(v);
}

int **int_matrix(int N, int M)
{
    int **m = malloc( N * sizeof(int *));
    for(int i=0; i<N; i++) m[i] = int_vector(M);
    return m;
}

void free_int_matrix(int **m, int N)
{
    for(int i=0; i<N; i++) free_int_vector(m[i]);
    free(m);
}

double *double_vector(int N)
{
    double *v = calloc( N , sizeof(double) );
    return v;
}

void free_double_vector(double *v)
{
    free(v);
}

double **double_matrix(int N, int M)
{
    double **m = malloc( N * sizeof(double *));
    for(int i=0; i<N; i++) m[i] = double_vector(M);
    return m;
}

void free_double_matrix(double **m, int N)
{
    for(int i=0; i<N; i++) free_double_vector(m[i]);
    free(m);
}

double ***double_tensor(int N, int M, int L)
{
    double ***t = malloc( N * sizeof(double **));
    for(int i=0; i<N; i++) t[i] = double_matrix(M,L);
    return t;
}

void free_double_tensor(double ***t, int N, int M)
{
    for(int i=0; i<N; i++) free_double_matrix(t[i],M);
    free(t);
}
