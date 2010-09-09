#ifndef SPARSE_SYS_H_INCLUDED
#define SPARSE_SYS_H_INCLUDED

#include <stddef.h>

#ifndef MAT_SIZE_T
#define MAT_SIZE_T int
#endif


struct CSRMatrix
{
    size_t      m;
    size_t      n;
    size_t      nnz;

    MAT_SIZE_T *ia;
    MAT_SIZE_T *ja;

    double     *sa;
};



struct CSRMatrix *
csrmatrix_new_count_nnz(size_t m);

size_t
csrmatrix_new_elms_pushback(struct CSRMatrix *A);

size_t
csrmatrix_elm_index(size_t i, MAT_SIZE_T j, const struct CSRMatrix *A);

void
csrmatrix_sortrows(struct CSRMatrix *A);

void
csrmatrix_delete(struct CSRMatrix *A);

void
csrmatrix_zero(struct CSRMatrix *A);

#endif  /* SPARSE_SYS_H_INCLUDED */
