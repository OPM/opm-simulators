#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "bsr.h"

typedef
struct prec_t
{
    bsr_matrix *L;
    bsr_matrix *D;
    bsr_matrix *U;
    int noffsets;
    int(*offsets)[3];
}
prec_t;

prec_t *prec_new();
void prec_init(prec_t *P, bsr_matrix const *A);
int  prec_analyze(bsr_matrix *M, int (*offsets)[3]);
void prec_factorize(prec_t *P, bsr_matrix *A);
void prec_factorize2(prec_t *P, bsr_matrix *A);
void prec_apply3(prec_t *P, double *x);
void prec_apply3c(prec_t *P, double *x);
void prec_mapply3c(prec_t *P, double *x);
void prec_dapply3c(prec_t *P, double *x);
void prec_downcast(prec_t *P);
void prec_info(prec_t *P);

#ifdef __cplusplus
}
#endif
