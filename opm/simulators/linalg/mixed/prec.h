#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "bsr.h"

// preconditioner struct
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

/**
 * @brief Create empty preconditioner object.
 *
 * @return Pointer to preconditioner object.
 */
prec_t *prec_alloc();

/**
 * @brief Delete preconditioner object.
 *
 * @param A Pointer to preconditioner object.
 */
void prec_free(prec_t *P);

/**
 * @brief Initialize preconditioner object.
 *
 * @param P Pointer preconditioner object.
 * @param A Pointer to bsr matrix.
 */
void prec_init(prec_t *P, bsr_matrix const *A);

/**
 * @brief Identify off-diagonal ILU0 targets.
 *
 * @param M Pointer to bsr matrix.
 * @param offsets Pointer to offsets identifying off-diagonal targets.
 *
 * @return number of offdiagonal targets.
 */
int  prec_analyze(bsr_matrix *M, int (*offsets)[3]);

/**
 * @brief DILU factorization.
 *
 * @param P Pointer preconditioner object.
 * @param A Pointer to bsr matrix.
 */
void prec_dilu_factorize(prec_t *P, bsr_matrix *A);

/**
 * @brief ILU0 factorization.
 *
 * @param P Pointer preconditioner object.
 * @param A Pointer to bsr matrix.
 */
void prec_ilu0_factorize(prec_t *P, bsr_matrix *A);

/**
 * @brief Preconditioner application in mixed-precision.
 *
 * @note Algorithm onsists of lower and upper triangular solves
 *
 * @param P Pointer to preconditioner object.
 * @apram x Pointer to input/output vector
 */
void prec_mapply3c(prec_t *P, double *x);

/**
 * @brief Preconditioner applicationin double-precision.
 *
 * @note Algorithm onsists of lower and upper triangular solves
 *
 * @param P Pointer to preconditioner object.
 * @apram x Pointer to input/output vector
 */
void prec_dapply3c(prec_t *P, double *x);

/**
 * @brief Make single-precision copy of double-precision values.
 *
 * @param P Pointer to preconditioner object.
 */
void prec_downcast(prec_t *P);

/**
 * @brief Display preconditioner statistics.
 *
 * @param P Pointer to preconditioner object.
 */
void prec_info(prec_t *P);

#ifdef __cplusplus
}
#endif
