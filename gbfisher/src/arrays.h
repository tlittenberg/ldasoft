int *int_vector(int N);
void free_int_vector(int *v);

double *double_vector(int N);
void free_double_vector(double *v);

double **double_matrix(int N, int M);
void free_double_matrix(double **m, int N);

double ***double_tensor(int N, int M, int L);
void free_double_tensor(double ***t, int N, int M);
