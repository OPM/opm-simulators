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
    //float *flt;

} bsr_matrix;

bsr_matrix* bsr_new();
void bsr_init(bsr_matrix *A, int nrows, int nnz, int b);

void bsr_info(bsr_matrix *A);
void bsr_sparsity(const bsr_matrix *A, const char *name);
void bsr_nonzeros(bsr_matrix *A, const char *name);

void bsr_hello();


void vec_show(const double *x, int n, const char *name);

#ifdef __cplusplus
}
#endif
