#include <stddef.h>
#include <stdlib.h>

#ifndef MAT_SIZE_T
#define MAT_SIZE_T int
#endif

#include "blas_lapack.h"
#include "mimetic.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ------------------------------------------------------------------ */
void
mim_ip_simple_all(int ncells, int d, int max_ncf, int *ncf, int *cf,
                  int *fneighbour,  double *fcentroid, double *fnormal,
                  double *farea, double *ccentroid, double *cvol,
                  double *perm, double *Binv)
/* ------------------------------------------------------------------ */
{
    int i, j, c, f, nf, fpos, fpos2, lwork;

    double *C, *N, *A, *work, s;

    double cc[3] = { 0.0 };     /* No more than 3 space dimensions */

    lwork = 64 * (max_ncf * d);                   /* 64 from ILAENV() */
    C     = malloc((max_ncf * d) * sizeof *C);
    N     = malloc((max_ncf * d) * sizeof *N);
    A     = malloc(max_ncf       * sizeof *A);
    work  = malloc(lwork         * sizeof *work);

    if ((C != NULL) && (N != NULL) && (A != NULL) && (work != NULL)) {
        fpos = fpos2 = 0;

        for (c = 0; c < ncells; c++) {
            for (j = 0; j < d; j++) {
                cc[j] = ccentroid[j + c*d];
            }

            nf = ncf[c];

            for (i = 0; i < nf; i++) {
                f = cf[fpos + i];
                s = 2.0*(fneighbour[2 * f] == c) - 1.0;

                A[i] = farea[f];

                for (j = 0; j < d; j++) {
                    C[i + j*nf] = fcentroid  [j + f*d] - cc[j];
                    N[i + j*nf] = s * fnormal[j + f*d];
                }
            }

            mim_ip_simple(nf, d, cvol[c], A, &perm[c * d * d],
                          C, N, &Binv[fpos2], work, lwork);

            fpos  += nf;
            fpos2 += nf * nf;
        }
    }

    free(work);  free(A);  free(N);  free(C);
}


/* ------------------------------------------------------------------ */
void
mim_ip_simple(int nf, int d, double vol,
              double *A, double *K, double *C, double *N,
              double *Binv, double *work, int lwork)
/* ------------------------------------------------------------------ */
{
    mim_ip_span_nullspace(nf, d, C, A, Binv, work, lwork);
    mim_ip_linpress_exact(nf, d, vol, N, K, Binv, C);
}


/* ------------------------------------------------------------------ */
void
mim_ip_span_nullspace(int nf, int d, double *C, double *A,
                      double *X, double *work, int nwork)
/* ------------------------------------------------------------------ */
{
    MAT_SIZE_T m, n, k, ld, info, lwork;

    int    i, j;
    double a1, a2;

    double tau[3] = { 0.0 };  /* No more than 3 spatial dimensions */

    /* Step 1) X <- I */
    for (i = 0; i < nf*nf; i++) { X[i           ] = 0.0; }
    for (i = 0; i < nf   ; i++) { X[i * (nf + 1)] = 1.0; }

    /* Step 2) C <- orth(A * C) */
    for (j = 0; j < d; j++) {
        for (i = 0; i < nf; i++) {
            C[i + j*nf] *= A[i];
        }
    }

    m = nf;  n = d;  ld = nf;  k = d;  lwork = nwork;
    dgeqrf_(&m, &n,     C, &ld, tau, work, &lwork, &info);
    dorgqr_(&m, &n, &k, C, &ld, tau, work, &lwork, &info);

    /* Step 3) X <- A * (X - C*C') * A */
    a1 = -1.0;  a2 = 1.0;
    dsyrk_("Upper Triangular", "No Transpose",
           &m, &n, &a1, C, &ld, &a2, X, &ld);
    for (j = 0; j < nf; j++) {
        for (i = 0; i <= j; i++) {
            X[i + j*nf] *= A[i] * A[j];
        }
    }

    /* Account for DSYRK only assigning upper triangular part. */
    for (j = 0; j < nf; j++) {
        for (i = j + 1; i < nf; i++) {
            X[i + j*nf] = X[j + i*nf];
        }
    }
}


/* ------------------------------------------------------------------ */
void
mim_ip_linpress_exact(int nf, int d, double vol, double *N, double *K,
                      double *Binv, double *T)
/* ------------------------------------------------------------------ */
{
    MAT_SIZE_T m, n, k, ld1, ld2, ldBinv;
    int        i;
    double     a1, a2, t;

    t = 0.0;
    for (i = 0; i < d; i++) {
        t += K[i + i*d];
    }

    /* Step 4) T <- N*K */
    m   = nf ;  n   = d  ;  k = d;
    ld1 = nf ;  ld2 = d  ;
    a1  = 1.0;  a2  = 0.0;
    dgemm_("No Transpose", "No Transpose", &m, &n, &k,
           &a1, N, &ld1, K, &ld2, &a2, T, &ld1);

    /*  Step 5) Binv <- (N*K*N' + t*X) / vol */
    a1 = 1.0     /      vol ;
    a2 = 6.0 * t / (d * vol);
    ldBinv = nf;
    dgemm_("No Transpose", "Transpose", &m, &m, &n,
           &a1, T, &ld1, N, &ld1, &a2, Binv, &ldBinv);
}
