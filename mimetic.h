#ifndef MIMETIC_H_INCLUDED
#define MIMETIC_H_INCLUDED

void mim_ip_simple(MAT_SIZE_T nf, MAT_SIZE_T d, double v,
                   double *A, double *K, double *C, double *N,
                   double *Binv,
                   double *work, MAT_SIZE_T lwork);


#endif /* MIMETIC_H_INCLUDED */
