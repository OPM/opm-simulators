#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * @brief Mixed-precision bsr matrix.
 */
typedef
struct bsr_matrix
{
    // number or rows
    int nrows;
    // number or columns
    int ncols;
    // number or non-zero blocks
    int nnz;
    // block-size
    int b;

    // pointer ot row offsets
    int *rowptr;
    // pointer to column indices
    int *colidx;
    // pointer to double-precision values
    double *dbl;
    // pointer to single-precision values
    float *flt;

} bsr_matrix;

/**
 * @brief Create empty bsr matrix.
 *
 * @return Pointer to bsr matrix.
 */
bsr_matrix* bsr_alloc();

/**
 * @brief Delete bsr matrix.
 *
 * @param A Pointer to bsr matrix.
 */
void bsr_free(bsr_matrix *A);

/**
 * @brief Initialize bsr matrix.
 *
 * @param A Pointer to bsr matrix.
 * @param nrows Number of rows.
 * @apram nnz Number of nonzero blocks.
 * @param b Block size.
 */
void bsr_init(bsr_matrix *A, int nrows, int nnz, int b);

/**
 * @brief Sparse matrix-vector multiplication in mixed precision.
 *
 * @note Function is specialized for 3x3 block-sparse matrices.
 * @note Function uses AVX2 intrinsics.
 *
 * @param A Pointer to bsr matrix.
 * @param x Pointer to input vector.
 * @param y Pointer to output vector.
 */
void bsr_vmspmv3(bsr_matrix *A, const double *x, double *y);

/**
 * @brief Sparse matrix-vector multiplication in double precision.
 *
 * @note Function is specialized for 3x3 block-sparse matrices.
 * @note Function uses AVX2 intrinsics.
 *
 * @param A Pointer to bsr matrix.
 * @param x Pointer to input vector.
 * @param y Pointer to output vector.
 */
void bsr_vdspmv3(bsr_matrix *A, const double *x, double *y);

/**
 * @brief Make single-precision copy of double-precision values.
 *
 * @param M Pointer to bsr matrix
 */
void bsr_downcast(bsr_matrix *M);

/**
 * @brief Display matrix statistics.
 *
 * @param A Pointer to bsr matrix
 */
void bsr_info(bsr_matrix *A);

/**
 * @brief Display spasity pattern of first few rows.
 *
 * @param A Pointer to bsr matrix
 * @param name String with desired name of bsr_matrix
 */
void bsr_sparsity(const bsr_matrix *A, const char *name);

/**
 * @brief Display nonzero blocks of first few rows.
 *
 * @param A Pointer to bsr matrix
 * @param name String with desired name of bsr_matrix
 */
void bsr_nonzeros(bsr_matrix *A, const char *name);

#ifdef __cplusplus
}
#endif
