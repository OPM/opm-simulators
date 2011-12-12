/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>

#include <opm/core/linalg/blas_lapack.h>
#include <opm/core/pressure/mimetic/mimetic.h>

/* ------------------------------------------------------------------ */
void
mim_ip_simple_all(int ncells, int d, int max_nconn,
                  int *pconn, int *conn,
                  int *fneighbour, double *fcentroid, double *fnormal,
                  double *farea, double *ccentroid, double *cvol,
                  double *perm, double *Binv)
/* ------------------------------------------------------------------ */
{
    int i, j, c, f, nconn, fpos2, lwork;

    double *C, *N, *A, *work, s;

    double cc[3] = { 0.0 };     /* No more than 3 space dimensions */

    lwork = 64 * (max_nconn * d);                 /* 64 from ILAENV() */
    C     = malloc((max_nconn * d) * sizeof *C);
    N     = malloc((max_nconn * d) * sizeof *N);
    A     = malloc(max_nconn       * sizeof *A);
    work  = malloc(lwork           * sizeof *work);

    if ((C != NULL) && (N != NULL) && (A != NULL) && (work != NULL)) {
        fpos2 = 0;

        for (c = 0; c < ncells; c++) {
            for (j = 0; j < d; j++) {
                cc[j] = ccentroid[j + c*d];
            }

            nconn = pconn[c + 1] - pconn[c];

            for (i = 0; i < nconn; i++) {
                f = conn[pconn[c] + i];
                s = 2.0*(fneighbour[2 * f] == c) - 1.0;

                A[i] = farea[f];

                for (j = 0; j < d; j++) {
                    C[i + j*nconn] = fcentroid  [j + f*d] - cc[j];
                    N[i + j*nconn] = s * fnormal[j + f*d];
                }
            }

            mim_ip_simple(nconn, nconn, d, cvol[c], &perm[c * d * d],
                          C, A, N, &Binv[fpos2], work, lwork);

            fpos2 += nconn * nconn;
        }
    }

    free(work);  free(A);  free(N);  free(C);
}


/* ------------------------------------------------------------------ */
void
mim_ip_simple(int nf, int nconn, int d,
              double v, double *K, double *C,
              double *A, double *N,
              double *Binv,
              double *work, int lwork)
/* ------------------------------------------------------------------ */
{
    mim_ip_span_nullspace(nf, nconn, d,    C, A, Binv, work, lwork);
    mim_ip_linpress_exact(nf, nconn, d, v, K, N, Binv, work, lwork);
}


/* ------------------------------------------------------------------ */
void
mim_ip_span_nullspace(int nf, int nconn, int d,
                      double *C,
                      double *A,
                      double *X,
                      double *work, int nwork)
/* ------------------------------------------------------------------ */
{
    MAT_SIZE_T m, n, k, ldC, ldX, info, lwork;

    int    i, j;
    double a1, a2;

    double tau[3] = { 0.0 };  /* No more than 3 spatial dimensions */

    /* Step 1) X(1:nf, 1:nf) <- I_{nf} */
    for (j = 0; j < nf; j++) {
        for (i = 0; i < nf; i++) {
            X[i + j*nconn] = 0.0;
        }
        X[j * (nconn + 1)] = 1.0;
    }

    /* Step 2) C <- orth(A * C) */
    for (j = 0; j < d; j++) {
        for (i = 0; i < nf; i++) {
            C[i + j*nf] *= A[i];
        }
    }

    m = nf;  n = d;  ldC = nf;  k = d;  lwork = nwork;
    dgeqrf_(&m, &n,     C, &ldC, tau, work, &lwork, &info);
    dorgqr_(&m, &n, &k, C, &ldC, tau, work, &lwork, &info);

    /* Step 3) X <- A * (X - C*C') * A */
    ldX = nconn;
    a1 = -1.0;  a2 = 1.0;
    dsyrk_("Upper Triangular", "No Transpose",
           &m, &n, &a1, C, &ldC, &a2, X, &ldX);
    for (j = 0; j < nf; j++) {
        for (i = 0; i <= j; i++) {
            X[i + j*nconn] *= A[i] * A[j];
        }
    }

    /* Account for DSYRK only assigning upper triangular part. */
    for (j = 0; j < nf; j++) {
        for (i = j + 1; i < nf; i++) {
            X[i + j*nconn] = X[j + i*nconn];
        }
    }
}


/* ------------------------------------------------------------------ */
void
mim_ip_linpress_exact(int nf, int nconn, int d,
                      double vol, double *K,
                      double *N,
                      double *Binv,
                      double *work, int lwork)
/* ------------------------------------------------------------------ */
{
    MAT_SIZE_T m, n, k, ld1, ld2, ldBinv;
    int        i;
    double     a1, a2, t;

    assert (lwork >= d * nf);
#if defined(NDEBUG)
    /* Suppress warning about unused parameter. */
    (void) lwork;
#endif

    t = 0.0;
    for (i = 0; i < d; i++) {
        t += K[i + i*d];
    }

    /* Step 4) T <- N*K */
    m   = nf ;  n   = d  ;  k = d;
    ld1 = nf ;  ld2 = d  ;
    a1  = 1.0;  a2  = 0.0;
    dgemm_("No Transpose", "No Transpose", &m, &n, &k,
           &a1, N, &ld1, K, &ld2, &a2, work, &ld1);

    /*  Step 5) Binv <- (N*K*N' + t*X) / vol */
    a1 = 1.0     /      vol ;
    a2 = 6.0 * t / (d * vol);
    ldBinv = nconn;
    dgemm_("No Transpose", "Transpose", &m, &m, &n,
           &a1, work, &ld1, N, &ld1, &a2, Binv, &ldBinv);
}


/* ---------------------------------------------------------------------- */
void
mim_ip_compute_gpress(int nc, int d, const double *grav,
                      const int *pconn, const int *conn,
                      const double *fcentroid, const double *ccentroid,
                      double *gpress)
/* ---------------------------------------------------------------------- */
{
    int c, i, j;

    const double *cc, *fc;

    for (c = i = 0; c < nc; c++) {
        cc = ccentroid + (c * d);

        for (; i < pconn[c + 1]; i++) {
            fc = fcentroid + (conn[i] * d);

            gpress[i] = 0.0;
            for (j = 0; j < d; j++) {
                gpress[i] += grav[j] * (fc[j] - cc[j]);
            }
        }
    }
}


/* inv(B) <- \lambda_t(s)*inv(B)_0 */
/* ---------------------------------------------------------------------- */
void
mim_ip_mobility_update(int nc, const int *pconn, const double *totmob,
                       const double *Binv0, double *Binv)
/* ---------------------------------------------------------------------- */
{
    int c, i, n, p2;

    for (c = p2 = 0; c < nc; c++) {
        n = pconn[c + 1] - pconn[c];

        for (i = 0; i < n * n; i++) {
            Binv[p2 + i] = totmob[c] * Binv0[p2 + i];
        }

        p2 += n * n;
    }
}


/* G <- \sum_i \rho_i f_i(s) * G_0 */
/* ---------------------------------------------------------------------- */
void
mim_ip_density_update(int nc, const int *pconn, const double *omega,
                      const double *gpress0, double *gpress)
/* ---------------------------------------------------------------------- */
{
    int c, i;

    for (c = i = 0; c < nc; c++) {
        for (; i < pconn[c + 1]; i++) {
            gpress[i] = omega[c] * gpress0[i];
        }
    }
}
