#pragma once

#ifdef __cplusplus
extern "C" {
#endif


#include "bsr.h"
#include "prec.h"

#include <stdbool.h>

typedef
struct bslv_memory
{
    bool use_dilu;

    double tol;
    int    max_iter;
    double *e;

    int n;
    double **dtmp;

    bildu_prec *P;
}
bslv_memory;

bslv_memory *bslv_new();
void bslv_info(bslv_memory *mem, int count);
void bslv_init(bslv_memory *mem, double tol, int max_iter, bsr_matrix const *A, bool use_dilu);
int  bslv_pbicgstab3m(bslv_memory *mem, bsr_matrix *A, const double *b, double *x);
int  bslv_pbicgstab3d(bslv_memory *mem, bsr_matrix *A, const double *b, double *x);

#ifdef __cplusplus
}
#endif
