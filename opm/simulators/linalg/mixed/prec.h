#pragma once

#ifdef __cplusplus
extern "C" {
#endif


#include "bsr.h"

#include <stdbool.h>

typedef
struct bildu_prec
{
    bsr_matrix *L;
    bsr_matrix *D;
    bsr_matrix *U;
    int noffsets;
    int(*offsets)[3];
}
bildu_prec;

bildu_prec *bildu_new();
void bildu_init(bildu_prec *P, bsr_matrix const *A);
int  bildu_analyze(bsr_matrix *M, int (*offsets)[3]);
void bildu_factorize(bildu_prec *P, bsr_matrix *A);
void bildu_factorize2(bildu_prec *P, bsr_matrix *A);
void bildu_apply3(bildu_prec *P, double *x);
void bildu_apply3c(bildu_prec *P, double *x);
void bildu_mapply3c(bildu_prec *P, double *x);
void bildu_dapply3c(bildu_prec *P, double *x);
void bildu_downcast(bildu_prec *P);
void bildu_info(bildu_prec *P);

#ifdef __cplusplus
}
#endif
