#pragma once

#ifdef __cplusplus
extern "C" {
#endif


#include "bsr.h"
#include "prec.h"

#include <stdbool.h>

// Solver memory struct
typedef
struct bslv_memory
{
    bool use_dilu;

    double tol;
    int    max_iter;
    double *e;

    int n;
    double **dtmp;

    prec_t *P;
}
bslv_memory;

/**
 * @brief Create empty solver memory object.
 *
 * @return Pointer to solver memory object.
 */
bslv_memory *bslv_alloc();

/**
 * @brief Delete solver memroy object.
 *
 * @param mem Pointer to solver memory object.
 */
void bslv_free(bslv_memory *mem);

/**
 * @brief Display convergence information.
 *
 * @param mem Pointer to solver memory object.
 * @param count Number of linear iterations (returned by linear solver)
 */
void bslv_info(bslv_memory *mem, int count);

/**
 * @brief Initialize solver memory object.
 *
 * @param mem Pointer to solver memory object.
 * @param tol Linear solver tolerance.
 * @apram max_iter Maximum number of linear iterations.
 * @param A Pointer to coeffient matrix in bsr format.
 * @param use_dilu: Select DILU preconditioner if true. Otherwise use ILU9.
 */
void bslv_init(bslv_memory *mem, double tol, int max_iter, bsr_matrix const *A, bool use_dilu);

/**
 * @brief Preconditioned bicgstab in mixed-precision.
 *
 * @note Preconditioner is either ILU0 or DILU based on value of mem->use_dilu
 *
 * @param mem Pointer to solver memory object.
 * @param A Pointer to coeffient matrix in bsr format
 * @param b Pointer to right-hand side vector
 * @apram x Pointer to solution vector
 *
 * @return Number of linear iterations.
 */
int  bslv_pbicgstab3m(bslv_memory *mem, bsr_matrix *A, const double *b, double *x);

/**
 * @brief Preconditioned bicgstab in double-precision.
 *
 * @note Preconditioner is either ILU0 or DILU based on value of mem->use_dilu
 *
 * @param mem Pointer to solver memory object.
 * @param A Pointer to coeffient matrix in bsr format
 * @param b Pointer to right-hand side vector
 * @apram x Pointer to solution vector
 *
 * @return Number of linear iterations.
 */
int  bslv_pbicgstab3d(bslv_memory *mem, bsr_matrix *A, const double *b, double *x);

#ifdef __cplusplus
}
#endif
