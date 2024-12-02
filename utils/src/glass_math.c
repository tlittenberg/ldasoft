//
//  math.c
//  data
//
//  Created by Tyson Littenberg on 11/29/23.
//

#include "glass_utils.h"

void invert_noise_covariance_matrix(struct Noise *noise)
{
    int X,Y,Z,A,E;
    
    for(int n=0; n<noise->N; n++)
    {
        switch(noise->Nchannel)
        {
            case 1:
                X=0;
                noise->detC[n] = noise->C[X][X][n];
                noise->invC[X][X][n] = 1./noise->C[X][X][n];
                break;
            case 2:
                A=0, E=1;
                noise->detC[n] = noise->C[A][A][n]*noise->C[E][E][n];
                noise->invC[A][A][n] = 1./noise->C[A][A][n];
                noise->invC[E][E][n] = 1./noise->C[E][E][n];
                break;
            case 3:
                X=0, Y=1, Z=2;
                double cxx = noise->C[X][X][n];
                double cyy = noise->C[Y][Y][n];
                double czz = noise->C[Z][Z][n];
                double cxy = noise->C[X][Y][n];
                double cxz = noise->C[X][Z][n];
                double cyz = noise->C[Y][Z][n];
                noise->detC[n] = cxx*(czz*cyy - cyz*cyz) - cxy*(cxy*czz - cxz*cyz) + cxz*(cxy*cyz - cyy*cxz);
                double invdetC = 1./noise->detC[n];
                noise->invC[X][X][n] = (cyy*czz - cyz*cyz)*invdetC;
                noise->invC[Y][Y][n] = (czz*cxx - cxz*cxz)*invdetC;
                noise->invC[Z][Z][n] = (cxx*cyy - cxy*cxy)*invdetC;
                noise->invC[X][Y][n] = (cxz*cyz - czz*cxy)*invdetC;
                noise->invC[X][Z][n] = (cxy*cyz - cxz*cyy)*invdetC;
                noise->invC[Y][Z][n] = (cxy*cxz - cxx*cyz)*invdetC;
                noise->invC[Y][X][n] = noise->invC[X][Y][n];
                noise->invC[Z][X][n] = noise->invC[X][Z][n];
                noise->invC[Z][Y][n] = noise->invC[Y][Z][n];
                break;
        }
    }
}

double ipow(double x, int n)
{
    double xn = 1.0;
    for(int i=0; i<n; i++) xn*=x;
    return xn;
}

double chirpmass(double m1, double m2)
{
    return pow(m1*m2,3./5.)/pow(m1+m2,1./5.);
}

double amplitude(double Mc, double f0, double D)
{
    double f = f0;;
    double M = Mc*TSUN;
    double dL= D*PC/CLIGHT;
    
    return 2.*pow(pow(M,5)*pow(M_PI*f,2),1./3.)/dL;
}

double power_spectrum(double *data, int n)
{
    int i,j;
    double Re, Im;
    
    i = 2*n;
    j = i+1;
    Re = data[i];
    Im = data[j];
    
    return (Re*Re + Im*Im);
}

double fourier_nwip(double *a, double *b, double *invC, int N)
{
    int j,k;
    double arg, product;
    double ReA, ReB, ImA, ImB;
    
    arg = 0.0;
    for(int n=0; n<N; n++)
    {
        j = n * 2;
        k = j + 1;
        ReA = a[j]; ImA = a[k];
        ReB = b[j]; ImB = b[k];
        product = ReA*ReB + ImA*ImB;
        arg += product*invC[n];
    }
    
    return(2.0*arg);
}

double wavelet_nwip(double *a, double *b, double *invC, int *list, int N)
{
    double arg = 0.0;
    for(int n=0; n<N; n++)
    {
        if(list[n] > 0)
        {
            int k = list[n];
            arg += a[k]*b[k]*invC[k];
        }
        
    }
    return arg;
}


// Recursive binary search function.
// Return nearest smaller neighbor of x in array[nmin,nmax] is present,
// otherwise -1
int binary_search(double *array, int nmin, int nmax, double x)
{
    int next;
    if(nmax>nmin)
    {
        int mid = nmin + (nmax - nmin) / 2;
        
        //find next unique element of array
        next = mid;
        while(array[mid]==array[next]) next++;
        
        // If the element is present at the middle
        // itself
        if (x > array[mid])
        {
            if(x < array[next]) return mid;
        }
        
        // the element is in the lower half
        if (array[mid] > x)
            return binary_search(array, nmin, mid, x);
        
        // the element is in upper half
        return binary_search(array, next, nmax, x);
    }
    
    // We reach here when element is not
    // present in array
    return -1;
}

void matrix_eigenstuff(double **matrix, double **evector, double *evalue, int N)
{
    int i,j;
    
    // Don't let errors kill the program (yikes)
    gsl_set_error_handler_off ();
    int err=0;
    
    // Find eigenvectors and eigenvalues
    gsl_matrix *GSLfisher = gsl_matrix_alloc(N,N);
    gsl_matrix *GSLcovari = gsl_matrix_alloc(N,N);
    gsl_matrix *GSLevectr = gsl_matrix_alloc(N,N);
    gsl_vector *GSLevalue = gsl_vector_alloc(N);
    
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            if(matrix[i][j]!=matrix[i][j])fprintf(stderr,"nan matrix element at line %d in file %s\n", __LINE__, __FILE__);
            gsl_matrix_set(GSLfisher,i,j,matrix[i][j]);
        }
    }
    
    // sort and put them into evec
    gsl_eigen_symmv_workspace * workspace = gsl_eigen_symmv_alloc (N);
    gsl_permutation * permutation = gsl_permutation_alloc(N);
    err += gsl_eigen_symmv (GSLfisher, GSLevalue, GSLevectr, workspace);
    err += gsl_eigen_symmv_sort (GSLevalue, GSLevectr, GSL_EIGEN_SORT_ABS_ASC);
    
    // eigenvalues destroy matrix
    for(i=0; i<N; i++) for(j=0; j<N; j++) gsl_matrix_set(GSLfisher,i,j,matrix[i][j]);
    
    err += gsl_linalg_LU_decomp(GSLfisher, permutation, &i);
    err += gsl_linalg_LU_invert(GSLfisher, permutation, GSLcovari);
    
    if(err>0)
    {
        for(i=0; i<N; i++)for(j=0; j<N; j++)
        {
            evector[i][j] = 0.0;
            if(i==j)
            {
                evector[i][j]=1.0;
                evalue[i]=1./matrix[i][j];
            }
        }
        
    }
    else
    {
        
        //unpack arrays from gsl inversion
        for(i=0; i<N; i++)
        {
            evalue[i] = gsl_vector_get(GSLevalue,i);
            for(j=0; j<N; j++)
            {
                evector[i][j] = gsl_matrix_get(GSLevectr,i,j);
                if(evector[i][j] != evector[i][j]) evector[i][j] = 0.;
            }
        }
                
        //copy covariance matrix back into Fisher
        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                matrix[i][j] = gsl_matrix_get(GSLcovari,i,j);
            }
        }
        
        //cap minimum size eigenvalues
        for(i=0; i<N; i++)
        {
            if(evalue[i] != evalue[i] || evalue[i] <= 10.0) evalue[i] = 10.;
        }
    }
    
    gsl_vector_free (GSLevalue);
    gsl_matrix_free (GSLfisher);
    gsl_matrix_free (GSLcovari);
    gsl_matrix_free (GSLevectr);
    gsl_eigen_symmv_free (workspace);
    gsl_permutation_free (permutation);
}

void invert_matrix(double **matrix, int N)
{
    int i,j;
    
    // Don't let errors kill the program (yikes)
    gsl_set_error_handler_off ();
    int err=0;
    
    // Find eigenvectors and eigenvalues
    gsl_matrix *GSLmatrix = gsl_matrix_alloc(N,N);
    gsl_matrix *GSLinvrse = gsl_matrix_alloc(N,N);
    
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            if(matrix[i][j]!=matrix[i][j])fprintf(stderr,"nan matrix element at line %d in file %s\n", __LINE__, __FILE__);
            gsl_matrix_set(GSLmatrix,i,j,matrix[i][j]);
        }
    }
    
    gsl_permutation * permutation = gsl_permutation_alloc(N);
    
    err += gsl_linalg_LU_decomp(GSLmatrix, permutation, &i);
    err += gsl_linalg_LU_invert(GSLmatrix, permutation, GSLinvrse);
    
    if(err>0)
    {
        fprintf(stderr,"WARNING: singluar matrix at line %d in file %s\n", __LINE__, __FILE__);
        fflush(stderr);
        exit(1);
    }
    else
    {
        //copy inverse back into matrix
        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                matrix[i][j] = gsl_matrix_get(GSLinvrse,i,j);
            }
        }
    }
    
    gsl_matrix_free (GSLmatrix);
    gsl_matrix_free (GSLinvrse);
    gsl_permutation_free (permutation);
}

void matrix_multiply(double **A, double **B, double **AB, int N)
{
    //AB = A*B
    
    int i,j,k;
    
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            AB[i][j] = 0.0;
            for(k=0; k<N; k++)
            {
                AB[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    
}


void cholesky_decomp(double **A, double **L, int N)
{
    /*
     factorize matrix A into cholesky decomposition L.
     GSL overwrites the original matrix which we want
     to preserve
     */
    int i,j;
    gsl_matrix *GSLmatrix = gsl_matrix_alloc(N,N);
    
    //copy covariance matrix into workspace
    for(i=0; i<N; i++) for(j=0; j<N; j++) gsl_matrix_set(GSLmatrix,i,j,A[i][j]);
    
    //make the magic happen
    gsl_linalg_cholesky_decomp(GSLmatrix);
    
    //copy cholesky decomposition into output matrix
    for(i=0; i<N; i++) for(j=0; j<N; j++)  L[i][j] = gsl_matrix_get(GSLmatrix,i,j);
    
    //zero upper half of matrix (copy of A)
    for(i=0; i<N; i++) for(j=i+1; j<N; j++) L[i][j] = 0.0;
    
    gsl_matrix_free (GSLmatrix);
    
}

/*
static void window(double *data, int N)
{
    int i;
    double filter;
    double tau = LISA_CADENCE/FILTER_LENGTH;
    for(i=0; i<N; i++)
    {
        double x = i*LISA_CADENCE;
        filter = (0.5*(1+tanh(tau*(x-FILTER_LENGTH))))*(0.5*(1-tanh(tau*(x-N*LISA_CADENCE))));
        data[i]*=filter;
    }
}
*/
void tukey(double *data, double alpha, int N)
{
    int i, imin, imax;
    double filter;
    
    imin = (int)(alpha*(double)(N-1)/2.0);
    imax = (int)((double)(N-1)*(1.0-alpha/2.0));
    
    for(i=0; i< N; i++)
    {
        filter = 1.0;
        if(i < imin) filter = 0.5*(1.0+cos(M_PI*( (double)(i)/(double)(imin)-1.0 )));
        if(i>imax)   filter = 0.5*(1.0+cos(M_PI*( (double)(N-1-i)/(double)(imin)-1.0)));
        data[i] *= filter;
    }
}

/*
static double tukey_scale(double alpha, int N)
{
    int i, imin, imax;
    double scale = 0.0;
    double filter;
    
    imin = (int)(alpha*(double)(N-1)/2.0);
    imax = (int)((double)(N-1)*(1.0-alpha/2.0));
    
    int Nwin = N-imax;
    
    for(i=0; i< N; i++)
    {
        filter = 1.0;
        if(i<imin) filter = 0.5*(1.0+cos(M_PI*( (double)(i)/(double)(imin)-1.0 )));
        if(i>imax) filter = 0.5*(1.0+cos(M_PI*( (double)(i-imax)/(double)(Nwin))));
        scale += filter;
    }
    scale /= (double)(N);
    
    return scale;
}
*/

void unpack_gsl_rft_output(double *x, double *x_gsl, int N)
{
    x[0] = x_gsl[0];
    x[1] = 0.0;

    for(int n=1; n<N/2; n++)
    {
        x[2*n]   = x_gsl[2*n-1];
        x[2*n+1] = x_gsl[2*n];
    }
}

void unpack_gsl_fft_output(double *x, double *x_gsl, int N)
{
    x[0] = x_gsl[0];
    x[1] = 0.0;
    for(int n=1; n<N/2; n++)
    {
        x[2*n]   = x_gsl[n];
        x[2*n+1] = x_gsl[N-n];
    }
}

void pack_gsl_fft_input(double *x, double *x_gsl, int N)
{
    x_gsl[0]  = x[0];
    x_gsl[N/2] = 0.0;
    for(int n=1; n<N/2; n++)
    {
        x_gsl[n] = x[2*n];
        x_gsl[N-n] = x[2*n+1];
    }
}

void CubicSplineGSL(int N, double *x, double *y, int Nint, double *xint, double *yint)
{
    int n;
    
    /* do our own error catching from interpolator */
    gsl_set_error_handler_off();
    
    /* set up GSL spline */
    
    /* Standard cubic spline */
    //gsl_spline *cspline = gsl_spline_alloc(gsl_interp_cspline, N);
    
    /*
     Non-rounded Akima spline with natural boundary conditions.
     This method uses the non-rounded corner algorithm of Wodicka.
     Akima splines are ideal for fitting curves with rapidly
     changing second derivatives.  They are C1 differentiable.
     See
     https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_type.gsl_interp_akima
     */
    gsl_spline *cspline = gsl_spline_alloc(gsl_interp_akima, N);
    
    /*
     Steffen's splines are guaranteed to be monotonic between
     control points.  Local maxima and minima only occur at
     at control points. See
     https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_type.gsl_interp_steffen
     */
    //gsl_spline *cspline = gsl_spline_alloc(gsl_interp_steffen, N);
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    
    /* get derivatives */
    int status = gsl_spline_init(cspline,x,y,N);
    
    //if error, return values that will be rjected by sampler
    if(status) for(n=0; n<Nint; n++) yint[n]=1.0;
    
    //otherwise proceed w/ interpolation
    else for(n=0; n<Nint; n++) yint[n]=gsl_spline_eval(cspline,xint[n],acc);
    
    gsl_spline_free (cspline);
    gsl_interp_accel_free (acc);
    
}

static gsl_vector* gsl_vector_union(gsl_vector *N, gsl_vector *Np)
{
    // allocate memory for work space
    double *N_temp = malloc(N->size*sizeof(double));
    double *Np_temp = malloc(Np->size*sizeof(double));
    
    int Nsize = N->size;
    int Nmax = N->size+Np->size;
    
    //printf("N.size=%i, Np.size=%i, Nmax=%i\n",N->size, Np->size, Nmax);
    
    double *N_union = malloc(Nmax*sizeof(double));
    
    // store contents of vectors in work space
    for(int n=0; n<N->size; n++)
    {
        N_temp[n]  = gsl_vector_get(N,n);
        N_union[n] = N_temp[n];
    }
    for(int n=0; n<Np->size; n++) Np_temp[n] = gsl_vector_get(Np,n);
    
    // get union of work space arrays
    for(int i=0; i<Np->size; i++)
    {
        //traverse N to make sure it is a new point
        int flag=0;
        for(int l=0; l<N->size; l++)
        {
            if(N_temp[l]==Np_temp[i]) flag=1;
        }
        if(!flag)
        {
            N_union[Nsize]=Np_temp[i];
            Nsize++;
        }
    }
    

    // resize N to and copy work space union into GSL vector    
    gsl_vector_free(N);
    gsl_vector *Nnew = gsl_vector_alloc(Nsize);
    for(int n=0; n<Nsize; n++) gsl_vector_set(Nnew,n,N_union[n]);
    
    // clean up
    free(N_temp);
    free(Np_temp);
    free(N_union);
    
    return Nnew;
            
}

static gsl_vector * find_neighbors(gsl_vector *X, double P, double epsilon)
{
    /*
     return all points in X w/in epsilon of P, including P
     */
    int N[X->size];
    int Nsize=0;

    // check all points in X
    for(int n=0; n<X->size; n++)
    {
        // distance measure
        double D = fabs(gsl_vector_get(X,n) - P);
        
        // store points closer than epsilon
        if(D < epsilon)
        {
            N[Nsize] = n;
            Nsize++;
        }
    }
    
    //create GSL vector with list of points within epsilon of P
    gsl_vector *Neighbors = gsl_vector_alloc(Nsize);
    for(int n=0; n<Nsize; n++) gsl_vector_set(Neighbors,n,N[n]);
    
    return Neighbors;
}

void dbscan(gsl_vector *X, double eps, int min, int C[], int *K)
{
            
    //Step 1: initialize cluster index and mark all points as unvisited
    int visited[X->size];
    int cluster_index=0;
    for(int n=0; n<X->size; n++)
    {
        visited[n] = 0;
        C[n] = 0;
    }
    
    //loop through all data points
    for(int n=0; n<X->size; n++)
    {
        //skip visited points
        if(visited[n]) continue;

        //Step 2: Choose a random unvisited data point p and mark it as visited
        double P = gsl_vector_get(X,n);
        visited[n] = 1;
        
        //Step 3: Find neighbors of P within epsilon and store them in N
        gsl_vector *N = find_neighbors(X,P,eps);
                
        /*
         Step 4: If N has at least min_samples elements,
         then P is a core point and a new cluster C is formed with P and N.
         Otherwise, X[n] is a noise point and no cluster is formed.
         */
        if(N->size>=min)
        {
            
            //assign current data point to current cluster
            C[n] = cluster_index;
            
            //loop through all neighbors of current data point
            for(int m=0; m<N->size; m++)
            {
                int j = (int)gsl_vector_get(N,m);
                
                //skip visited neighbors
                if(visited[j]) continue;
                
                //mark neighbor as visited
                visited[j] = 1;
                
                //find neighbors of P's neighbors
                gsl_vector *Np = find_neighbors(X,gsl_vector_get(X,j),eps);

                //add neighbors' of P's neighbors to cluster
                if(Np->size>=min)
                {
                    // N = N U N'
                    N = gsl_vector_union(N,Np);
                }
                
                
                //double check that neighbor is not already assigned to any cluster
                if(C[j]==0) C[j] = cluster_index;
                
                //free memory holding Np for point P
                gsl_vector_free(Np);
            }

            //advance index to be ready for next cluster
            cluster_index++;

        }else C[n]=-1;
        
        //free memory holding list of neighbors for next point
        gsl_vector_free(N);

    }//end loop over data points
    
    //assign number of clusters
    *K = cluster_index;
}

void unwrap_phase(int N, double *phase)
{
    double u, v, q;
    int i;
    
    v = phase[0];
    for(i=0; i<N ;i++)
    {
        u = phase[i];
        q = rint(fabs(u-v)/(PI2));
        if(q > 0.0)
        {
           if(v > u) u += q*PI2;
           else      u -= q*PI2;
        }
        v = u;
        phase[i] = u;
    }
}

double simpson_integration_3(double f0, double f1, double f2, double h)
{
    return h*(f0 + 4.0*f1 + f2)/6.0;
}

double simpson_integration_5(double f0, double f1, double f2, double f3, double f4, double h)
{
    return h*(f0 + 4.0*f1 + 2.0*f2 + 4.0*f3 + f4)/12.0;
}

static int int_compare(const void* a, const void* b) 
{
   return (*(int*)a - *(int*)b);
}

void integer_sort(int *x, int N)
{
    qsort(x, N, sizeof(int), int_compare);
}

void list_union(int *A, int *B, int NA, int NB, int *AUB, int *NAUB)
{
    int i,j;
    int Ntemp = NA+NB;
    int *Utemp = int_vector(Ntemp);
    int *temp = int_vector(Ntemp);

    //first list
    j=0;
    for(i=0; i<NA; i++)
    {
        Utemp[j]=A[i];
        j++;
    }
    //second list
    for(i=0; i<NB; i++)
    {
        Utemp[j]=B[i];
        j++;
    }
    //sort master list
    integer_sort(Utemp,Ntemp);
    
    //find and remeove repeats
    j=0;
    for(i=0; i<Ntemp-1;i++)
    {
        if (Utemp[i] != Utemp[i + 1])
        {
            temp[j] = Utemp[i];
            j++;
        }
    }
    
    // Store the last element as whether 
    // it is unique or repeated, it hasn't 
    // stored previously
    temp[j] = Utemp[Ntemp - 1];
    j++;

    // Modify original array
    for (i=0; i<j; i++)
        AUB[i] = temp[i];

    *NAUB = j;
    free_int_vector(Utemp);
    free_int_vector(temp);
}
