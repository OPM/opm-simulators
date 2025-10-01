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
void bildu_downcast(bildu_prec *P);
void bildu_info(bildu_prec *P);

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
