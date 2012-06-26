/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_SPARSE_SYS_HEADER_INCLUDED
#define OPM_SPARSE_SYS_HEADER_INCLUDED

/**
 * \file
 * Data structure and operations to manage sparse matrices in CSR formats.
 */

#include <stddef.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Basic compressed-sparse row (CSR) matrix data structure.
 */
struct CSRMatrix
{
    size_t      m;    /**< Number of rows */
    size_t      nnz;  /**< Number of structurally non-zero elements */

    int        *ia;   /**< Row pointers */
    int        *ja;   /**< Column indices */

    double     *sa;   /**< Matrix elements */
};


/**
 * Allocate a matrix structure and corresponding row pointers, @c ia,
 * sufficiently initialised to support "count and push-back"
 * construction scheme.
 *
 * The matrix will be fully formed in csrmatrix_new_elms_pushback().
 *
 * \param[in] m Number of matrix rows.
 *
 * \return Allocated matrix structure with allocated row pointers and
 * valid @c m field.  The row pointer elements are initialised all
 * zero to simplify the non-zero element counting procedure.  The
 * @c ja and @c sa fields are @c NULL.  This function returns @c NULL
 * in case of allocation failure.
 */
struct CSRMatrix *
csrmatrix_new_count_nnz(size_t m);


/**
 * Allocate a matrix structure and all constituent fields to hold a
 * sparse matrix with a specified number of (structural) non-zero
 * elements.
 *
 * The contents of the individual matrix arrays is undefined.  In
 * particular, the sparsity pattern must be constructed through some
 * other, external, means prior to using the matrix in (e.g.,) a
 * global system assembly process.
 *
 * The memory resources should be released through the
 * csrmatrix_delete() function.
 *
 * \param[in] m   Number of matrix rows.
 * \param[in] nnz Number of structural non-zeros.
 *
 * \return Allocated matrix structure and constituent element arrays.
 * @c NULL in case of allocation failure.
 */
struct CSRMatrix *
csrmatrix_new_known_nnz(size_t m, size_t nnz);


/**
 * Set row pointers and allocate column index and matrix element
 * arrays of a matrix previous obtained from
 * csrmatrix_new_count_nnz().
 *
 * The memory resources should be released through the
 * csrmatrix_delete() function.
 *
 * This function assumes that, on input, the total number of
 * structurally non-zero elements of row @c i are stored in
 * <CODE>A->ia[i+1]</CODE> for all <CODE>i = 0, ..., A->m - 1</CODE>
 * and that <CODE>A->ia[0] == 0</CODE>.  If successful, then on output
 * the row \em end pointers <CODE>A->ia[i+1]</CODE> are positioned at
 * the \em start of the corresponding rows.  If not, then the
 * <CODE>A->ja</CODE> and <CODE>A->sa</CODE> arrays remain unallocated.
 *
 * \param[in,out] A Matrix.
 *
 * \return Total number of allocated non-zeros, <CODE>A->nnz ==
 * A->ia[A->m]</CODE> if successful and zero in case of allocation
 * failure.
 */
size_t
csrmatrix_new_elms_pushback(struct CSRMatrix *A);


/**
 * Compute non-zero index of specified matrix element.
 *
 * \param[in] i Row index.
 * \param[in] j Column index.  Must be in the structural non-zero
 *              element set of row @c i.
 * \param[in] A Matrix.
 *
 * \return Non-zero index, into @c A->ja and @c A->sa, of the
 * <CODE>(i,j)</CODE> matrix element.
 */
size_t
csrmatrix_elm_index(int i, int j, const struct CSRMatrix *A);


/**
 * Sort column indices within each matrix row in ascending order.
 *
 * The corresponding matrix elements (i.e., @c sa) are not referenced.
 * Consequently, following a call to csrmatrix_sortrows(), all
 * relations to any pre-existing matrix elements are lost and must be
 * rebuilt.
 *
 * After a call to csrmatrix_sortrows(), the following relation holds
 * <CODE>A->ja[k] < A->ja[k+1]</CODE> for all <CODE>k = A->ia[i], ...,
 * A->ia[i+1]-2</CODE> in each row <CODE>i = 0, ..., A->m - 1</CODE>.
 *
 * \param[in,out] A Matrix.
 */
void
csrmatrix_sortrows(struct CSRMatrix *A);


/**
 * Dispose of memory resources obtained through prior calls to
 * allocation routines.
 *
 * \param[in,out] A Matrix obtained from csrmatrix_new_count_nnz() +
 *                csrmatrix_new_elms_pushback() or
 *                csrmatrix_new_known_nnz().
 *
 * The pointer @c A is invalid following a call to csrmatrix_delete().
 */
void
csrmatrix_delete(struct CSRMatrix *A);


/**
 * Zero all matrix elements, typically in preparation of elemental
 * assembly.
 *
 * \param[in,out] A Matrix for which to zero the elements.
 */
void
csrmatrix_zero(struct CSRMatrix *A);


/**
 * Zero all vector elements.
 *
 * \param[in]  n Number of vector elements.
 * \param[out] v Vector for which to zero the elements.
 */
void
vector_zero(size_t n, double *v);


/**
 * Print matrix to file.
 *
 * The matrix content is printed in coordinate format with row and
 * column indices ranging from @c 1 to @c A->m.  This output format
 * facilitates simple processing through the @c spconvert function in
 * MATLAB© or Octave.
 *
 * This function is implemented in terms of csrmatrix_write_stream().
 *
 * \param[in] A  Matrix.
 * \param[in] fn Name of file to which matrix contents will be output.
 */
void
csrmatrix_write(const struct CSRMatrix *A, const char *fn);


/**
 * Print matrix to stream.
 *
 * The matrix content is printed in coordinate format with row and
 * column indices ranging from @c 1 to @c A->m.  This output format
 * facilitates simple processing through the @c spconvert function in
 * MATLAB© or Octave.
 *
 * \param[in]     A  Matrix.
 * \param[in,out] fp Open (text) stream to which matrix contents
 *                   will be output.
 */
void
csrmatrix_write_stream(const struct CSRMatrix *A, FILE *fp);


/**
 * Print vector to file.
 *
 * Elements are printed with one line (separated by <CODE>'\n'</CODE>)
 * per vector element.
 *
 * This function is implemented in terms of vector_write_stream().
 *
 * \param[in] n  Number of vector elements.
 * \param[in] v  Vector.
 * \param[in] fn Name of file to which vector contents will be output.
 */
void
vector_write(size_t n, const double *v, const char *fn);


/**
 * Print vector to stream.
 *
 * Elements are printed with one line (separated by <CODE>'\n'</CODE>)
 * per vector element.
 *
 * \param[in]     n  Number of vector elements.
 * \param[in]     v  Vector.
 * \param[in,out] fp Open (text) stream to which vector contents will be
 *                   output.
 */
void
vector_write_stream(size_t n, const double *v, FILE *fp);

#ifdef __cplusplus
}
#endif

#endif  /* OPM_SPARSE_SYS_HEADER_INCLUDED */
