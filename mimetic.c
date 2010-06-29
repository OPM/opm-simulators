#include <stddef.h>

#ifndef MAT_SIZE_T
#define MAT_SIZE_T size_t
#endif

#include "blas_lapack.h"
#include "mimetic.h"


/* ------------------------------------------------------------------ */
static double
trace(size_t n, double *A);
/* ------------------------------------------------------------------ */


/* ------------------------------------------------------------------ */
static void
mim_ip_simple_fill_null(size_t nf, size_t d, double *C, double *A,
                        double *Binv, double *work, MAT_SIZE_T lwork);
/* ------------------------------------------------------------------ */


/* ------------------------------------------------------------------ */
void
mim_ip_simple(MAT_SIZE_T nf, MAT_SIZE_T d, double vol,
              double *A, double *K, double *C, double *N,
              double *Binv, double* work, MAT_SIZE_T lwork)
/* ------------------------------------------------------------------ */
{
    MAT_SIZE_T m, n, k, ld1, ld2, ldBinv;
    double     a1, a2;

    mim_ip_simple_fill_null(nf, d, C, A, Binv, work, lwork);

    /* Step 4) C <- N*K */
    m   = nf ;  n   = nf ;  k = d;
    ld1 = nf ;  ld2 = d  ;
    a1  = 1.0;  a2  = 0.0;
    dgemm_("No Transpose", "No Transpose", &m, &n, &k,
           &a1, N, &ld1, K, &ld1, &a2, C, &ld1);

    /*  Step 5) Binv <- (N*K*N' + t*A*(I-Q*Q')*A) / vol */
    a1 = 1.0                /      vol ;
    a2 = 6.0 * trace(nf, K) / (d * vol);
    ldBinv = nf;
    dgemm_("No Transpose", "Transpose", &m, &m, &n,
           &a1, C, &ld1, N, &ld1, &a2, Binv, &ldBinv);
}


/* ------------------------------------------------------------------ */
static void
mim_ip_simple_fill_null(size_t nf, size_t d, double *C, double *A,
                        double *Binv,
                        double *work, MAT_SIZE_T lwork)
/* ------------------------------------------------------------------ */
{
    MAT_SIZE_T m, n, k, ld, info;

    size_t i, j;
    double a1, a2;

    double tau[3] = { 0.0 };  /* No more than 3 spatial dimensions */

    /* Step 1) Binv <- I */
    for (i = 0; i < nf*nf; i++) { Binv[i           ] = 0.0; }
    for (i = 0; i < nf   ; i++) { Binv[i * (nf + 1)] = 1.0; }

    /* Step 2) C <- orth(A * C) */
    for (j = 0; j < d; j++) {
        for (i = 0; i < nf; i++) {
            C[i + j*nf] *= A[i];
        }
    }

    m = nf;  n = d;  ld = nf;  k = d;
    dgeqrf_(&m, &n,     C, &ld, tau, work, &lwork, &info);
    dorgqr_(&m, &n, &k, C, &ld, tau, work, &lwork, &info);

    /* Step 3) Binv <- A * (Binv - C*C') * A */
    a1 = -1.0;  a2 = 1.0;
    dsyrk_("Upper Triangular", "No Transpose",
           &m, &n, &a1, C, &ld, &a2, Binv, &ld);
    for (j = 0; j < nf; j++) {
        for (i = 0; i <= j; i++) {
            Binv[i + j*nf] *= A[i] * A[j];
        }
    }

    /* Account for DSYRK only assigning upper triangular part. */
    for (j = 0; j < nf; j++) {
        for (i = j + 1; i < nf; i++) {
            Binv[i + j*nf] = Binv[j + i*nf];
        }
    }
}


/* ------------------------------------------------------------------ */
static double
trace(size_t n, double *A)
/* ------------------------------------------------------------------ */
{
    size_t i;
    double t = 0.0;

    for (i = 0; i < n; i++) {
        t += A[i + i*n];
    }

    return t;
}
