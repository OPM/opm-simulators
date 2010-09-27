/*
 * Copyright (c) 2010 SINTEF ICT, Applied Mathematics
 */

#ifndef MIMETIC_H_INCLUDED
#define MIMETIC_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

void mim_ip_span_nullspace(int nf, int nconn, int d,
                           double *C,
                           double *A,
                           double *X,
                           double *work, int lwork);

void mim_ip_linpress_exact(int nf, int nconn, int d,
                           double vol, double *K,
                           double *N,
                           double *Binv,
                           double *work, int lwork);

void mim_ip_simple(int nf, int nconn, int d,
                   double v, double *K, double *C,
                   double *A, double *N,
                   double *Binv,
                   double *work, int lwork);

void mim_ip_simple_all(int ncells, int d, int max_ncf, int *ncf,
                       int *pconn, int *conn,
                       int *fneighbour, double *fcentroid, double *fnormal,
                       double *farea, double *ccentroid, double *cvol,
                       double *perm, double *Binv);

#ifdef __cplusplus
}
#endif

#endif /* MIMETIC_H_INCLUDED */
