#ifndef MIMETIC_H_INCLUDED
#define MIMETIC_H_INCLUDED


void mim_ip_span_nullspace(int nf, int d, double *C, double *A,
                           double *X, double *work, int lwork);

void mim_ip_linpress_exact(int nf, int d, double vol, double *N,
                           double *K, double *Binv, double *T);

void mim_ip_simple(int nf, int d, double v, double *A, double *K,
                   double *C, double *N, double *Binv,
                   double *work, int lwork);


#endif /* MIMETIC_H_INCLUDED */
