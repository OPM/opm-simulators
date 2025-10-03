#include "bslv.h"
#include "vec.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#pragma GCC push_options
#pragma GCC target("avx2")

bslv_memory *bslv_new()
{
    bslv_memory *mem = malloc(sizeof(bslv_memory));
    return mem;
}

void bslv_info(bslv_memory *mem, int count)
{
    double * restrict e = mem->e;
    printf("bslv_info: iterations=%d reduction=%.2e\n",count,e[count]);
}


void bslv_init(bslv_memory *mem, double tol, int max_iter, bsr_matrix const *A, bool use_dilu)
{
    int n=(A->nrows * A->b);

    mem->use_dilu = use_dilu;

    mem->tol = tol;
    mem->max_iter = max_iter;
    mem->n = n;

    mem->e = (double*) malloc(max_iter*sizeof(double));

    int narrays=7;
    mem->dtmp = (double**) malloc(narrays*sizeof(double*));

    int np = 8*((n+7)/8); // padded to nearest multiple of 8
    for(int i=0;i<narrays;i++)
    {
        posix_memalign((void**)&(mem->dtmp[i]),64,np*sizeof(double));
        for(int k=8*(n/8);k<np;k++) mem->dtmp[i][k] = 0.0; //zeroing out padded section
    }

    mem->P = prec_new();
    prec_init(mem->P, A); // initialize structure of L,D,U components of P

    // random initialization of one-dimensinal shadow space
    //   - note we fix the random seed for reproducibility
    //   - also note that this done once at solver initialzation
    srand(0);
    double rmax=RAND_MAX;
    double *y = mem->dtmp[0];
    double norm=0;
    for(int i=0;i<n;i++)
    {
        double x = rand()/rmax;// - 0.5;
        y[i]=x;
        norm+=x*x;
    }
    norm=sqrt(norm);
    for(int i=0;i<n;i++) y[i]/=norm;
    //printf("RAND_MAX=%d\n",RAND_MAX);

}

double __attribute__((noinline)) vec_inner2(const double *a, const double *b, int n)
{
    const double *x = __builtin_assume_aligned(a,64);
    const double *y = __builtin_assume_aligned(b,64);

    int const N=8;
    double agg[N];
    for(int i=0;i<N;i++) agg[i]=0.0;
    for(int i=0;i<n;i+=N)
    {
        for(int j=0;j<N;j++) agg[j]+=x[i+j]*y[i+j];
    }
    //for(int j=0;j<8;j++) agg[j]+=agg[j+8];
    for(int j=0;j<4;j++) agg[j]+=agg[j+4];
    for(int j=0;j<2;j++) agg[j]+=agg[j+2];
    for(int j=0;j<1;j++) agg[j]+=agg[j+1];

    return agg[0];

}

int bslv_pbicgstab3m(bslv_memory *mem, bsr_matrix *A, const double *b, double *x)
{

    double tol = mem->tol;
    int max_iter = mem->max_iter;
    int n = mem->n;

    double * restrict e = mem->e;
    const double *r0 = b;
    //const double * restrict r0  = mem->dtmp[0]; //access randomly initialized one-dimensional shadow space
          double * restrict p_j = mem->dtmp[1];
          double * restrict q_j = mem->dtmp[2];
          double * restrict r_j = mem->dtmp[3];
          double * restrict s_j = mem->dtmp[4];
          double * restrict t_j = mem->dtmp[5];
          double * restrict v_j = mem->dtmp[6];
          double * restrict x_j = x;

    prec_t * restrict P = mem->P;
    mem->use_dilu ? prec_dilu_factorize(P,A) : prec_ilu0_factorize(P,A); // choose dilu or ilu0
    prec_downcast(P);

    vec_fill(x_j,0.0,n);
    vec_copy(r_j,b,n);
    vec_copy(p_j,b,n);

    vec_copy(q_j,p_j,n);
    //double norm_0 = sqrt(vec_inner(r_j,r_j,n));
    double norm_0 = sqrt(vec_inner2(r_j,r_j,n));

    //double rho_j = vec_inner(r0,r_j,n);
    double rho_j = vec_inner2(r0,r_j,n);
    int j;
    for(j=0;j<max_iter;j++)
    {
        //vec_copy(q_j,p_j,n);                                        //q_j=p_j
        prec_mapply3c(P,q_j);                                          //q_j=P.q_j;
        bsr_vmspmv3(A,q_j,v_j);                                        //v_j= A.q_j

        //double alpha_j = rho_j/vec_inner(r0,v_j,n);
        double alpha_j = rho_j/vec_inner2(r0,v_j,n);
        for (int k=0;k<n;k++) q_j[k] = s_j[k] = r_j[k]-alpha_j*v_j[k];       // r_j and s_j can overwrite each other
        //vec_axpy(-alpha_j,v_j,r_j,s_j,n);

        //vec_copy(q_j,s_j,n);                                        //q_j=s_j
        prec_mapply3c(P,q_j);                                          //q_j=P.q_j;
        bsr_vmspmv3(A,q_j,t_j);                                        //t_j= A.q_j

        //double w_j = vec_inner(s_j,t_j,n)/vec_inner(t_j,t_j,n);
        double w_j = vec_inner2(s_j,t_j,n)/vec_inner2(t_j,t_j,n);
        for (int k=0;k<n;k++) x_j[k] += alpha_j*p_j[k] + w_j*s_j[k];
        for (int k=0;k<n;k++) r_j[k]  =         s_j[k] - w_j*t_j[k];
        //vec_axpy(-w_j,t_j,s_j,r_j,n);
        //double norm_e = sqrt(vec_inner(r_j,r_j,n));
        double norm_e = sqrt(vec_inner2(r_j,r_j,n));
        e[j+1]=norm_e/norm_0;

        if (norm_e<tol*norm_0) break;                               //convergence check

        double rho_0 =rho_j;
        //rho_j=vec_inner(r0,r_j,n);
        rho_j=vec_inner2(r0,r_j,n);

        double beta_j = (alpha_j/w_j)*(rho_j/rho_0);
        for (int k=0;k<n;k++) q_j[k] = p_j[k] = r_j[k] + beta_j*(p_j[k] - w_j*v_j[k]);
    }
    prec_mapply3c(P,x_j);                                       //x_j=P.x_j;

    return j == max_iter ? j : ++j;
}


int bslv_pbicgstab3d(bslv_memory *mem, bsr_matrix *A, const double *b, double *x)
{

    double tol = mem->tol;
    int max_iter = mem->max_iter;
    int n = mem->n;

    double * restrict e = mem->e;
    const double *r0 = b;
    //const double * restrict r0  = mem->dtmp[0]; //access randomly initialized one-dimensional shadow space
          double * restrict p_j = mem->dtmp[1];
          double * restrict q_j = mem->dtmp[2];
          double * restrict r_j = mem->dtmp[3];
          double * restrict s_j = mem->dtmp[4];
          double * restrict t_j = mem->dtmp[5];
          double * restrict v_j = mem->dtmp[6];
          double * restrict x_j = x;

    prec_t * restrict P = mem->P;
    mem->use_dilu ? prec_dilu_factorize(P,A) : prec_ilu0_factorize(P,A); // choose dilu or ilu0
    prec_downcast(P);

    vec_fill(x_j,0.0,n);
    vec_copy(r_j,b,n);
    vec_copy(p_j,b,n);

    vec_copy(q_j,p_j,n);
    //double norm_0 = sqrt(vec_inner(r_j,r_j,n));
    double norm_0 = sqrt(vec_inner2(r_j,r_j,n));

    //double rho_j = vec_inner(r0,r_j,n);
    double rho_j = vec_inner2(r0,r_j,n);
    int j;
    for(j=0;j<max_iter;j++)
    {
        //vec_copy(q_j,p_j,n);                                        //q_j=p_j
        prec_dapply3c(P,q_j);                                          //q_j=P.q_j;
        bsr_vdspmv3(A,q_j,v_j);                                        //v_j= A.q_j

        //double alpha_j = rho_j/vec_inner(r0,v_j,n);
        double alpha_j = rho_j/vec_inner2(r0,v_j,n);
        for (int k=0;k<n;k++) q_j[k] = s_j[k] = r_j[k]-alpha_j*v_j[k];       // r_j and s_j can overwrite each other
        //vec_axpy(-alpha_j,v_j,r_j,s_j,n);

        //vec_copy(q_j,s_j,n);                                        //q_j=s_j
        prec_dapply3c(P,q_j);                                          //q_j=P.q_j;
        bsr_vdspmv3(A,q_j,t_j);                                        //t_j= A.q_j

        //double w_j = vec_inner(s_j,t_j,n)/vec_inner(t_j,t_j,n);
        double w_j = vec_inner2(s_j,t_j,n)/vec_inner2(t_j,t_j,n);
        for (int k=0;k<n;k++) x_j[k] += alpha_j*p_j[k] + w_j*s_j[k];
        for (int k=0;k<n;k++) r_j[k]  =         s_j[k] - w_j*t_j[k];
        //vec_axpy(-w_j,t_j,s_j,r_j,n);
        //double norm_e = sqrt(vec_inner(r_j,r_j,n));
        double norm_e = sqrt(vec_inner2(r_j,r_j,n));
        e[j+1]=norm_e/norm_0;

        if (norm_e<tol*norm_0) break;                               //convergence check

        double rho_0 =rho_j;
        //rho_j=vec_inner(r0,r_j,n);
        rho_j=vec_inner2(r0,r_j,n);

        double beta_j = (alpha_j/w_j)*(rho_j/rho_0);
        for (int k=0;k<n;k++) q_j[k] = p_j[k] = r_j[k] + beta_j*(p_j[k] - w_j*v_j[k]);
    }
    prec_dapply3c(P,x_j);                                       //x_j=P.x_j;

    return j == max_iter ? j : ++j;
}


#pragma GCC pop_options



