#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef
struct bsr_matrix
{
    int nrows;
    int ncols;
    int nnz;
    int b;

    int *rowptr;
    int *colidx;
    double *dbl;
    float *flt;

} bsr_matrix;

bsr_matrix* bsr_new();
void bsr_init(bsr_matrix *A, int nrows, int nnz, int b);
void bsr_vmspmv3(bsr_matrix *A, const double *x, double *y);
void bsr_downcast(bsr_matrix *M);

void bsr_info(bsr_matrix *A);
void bsr_sparsity(const bsr_matrix *A, const char *name);
void bsr_nonzeros(bsr_matrix *A, const char *name);


typedef
struct bildu_prec
{
    bsr_matrix *L;
    bsr_matrix *D;
    bsr_matrix *U;
}
bildu_prec;

bildu_prec *bildu_new();
void bildu_init(bildu_prec *P, bsr_matrix const *A);
void bildu_factorize(bildu_prec *P, bsr_matrix *A);
void bildu_apply3(bildu_prec *P, double *x);
void bildu_apply3c(bildu_prec *P, double *x);
void bildu_downcast(bildu_prec *P);
void bildu_info(bildu_prec *P);

typedef
struct bslv_memory
{
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
void bslv_init(bslv_memory *mem, double tol, int max_iter, bsr_matrix const *A);
void bslv_pbicgstab3(bslv_memory *mem, bsr_matrix *A, const double *b, double *x);


double __attribute__((noinline)) vec_inner2(const double *a, const double *b, int n);
void vec_fill(double *y, double x, int n);
void vec_show(const double *x, int n, const char *name);

void headtail(double *x, int n, char const *name);

#ifdef __cplusplus
}
#endif
