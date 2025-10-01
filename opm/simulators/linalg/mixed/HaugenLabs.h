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
    float  **stmp;

    bildu_prec *P;
}
bslv_memory;

bslv_memory *bslv_new();
void bslv_init(bslv_memory *mem, double tol, int max_iter, bsr_matrix const *A, bool use_dilu);
int  bslv_pbicgstab3m(bslv_memory *mem, bsr_matrix *A, const double *b, double *x);
void bslv_info(bslv_memory *mem, int count);

double __attribute__((noinline)) vec_inner2(const double *a, const double *b, int n);
void vec_fill(double *y, double x, int n);
void vec_copy(double *y, double const * x, int n);
void vec_show(const double *x, int n, const char *name);
double vec_norm(double const *v, int n);

void headtail(double *x, int n, char const *name);

#ifdef __cplusplus
}
#endif
