/*
  Copyright 2019 Equinor ASA

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

#ifndef OPM_BLOCKED_MATRIX_HPP
#define OPM_BLOCKED_MATRIX_HPP

namespace Opm
{
namespace Accelerator
{

/// This struct resembles a blocked csr matrix, like Dune::BCRSMatrix.
/// The data is stored in contiguous memory, such that they can be copied to a device in one transfer.
class BlockedMatrix
{

public:

    /// Allocate BlockedMatrix and data arrays with given sizes
    /// \param[in] Nb               number of blockrows
    /// \param[in] nnzbs            number of nonzero blocks
    /// \param[in] block_size       the number of rows and columns for each block
    BlockedMatrix(int Nb_, int nnzbs_, unsigned int block_size_)
    : nnzValues(new double[nnzbs_*block_size_*block_size_]),
      colIndices(new int[nnzbs_*block_size_*block_size_]),
      rowPointers(new int[Nb_+1]),
      Nb(Nb_),
      nnzbs(nnzbs_),
      block_size(block_size_),
      deleteNnzs(true),
      deleteSparsity(true)
    {}

    /// Allocate BlockedMatrix, but copy sparsity pattern instead of allocating new memory
    /// \param[in] M              matrix to be copied
    BlockedMatrix(const BlockedMatrix& M)
    : nnzValues(new double[M.nnzbs*M.block_size*M.block_size]),
      colIndices(M.colIndices),
      rowPointers(M.rowPointers),
      Nb(M.Nb),
      nnzbs(M.nnzbs),
      block_size(M.block_size),
      deleteNnzs(true),
      deleteSparsity(false)
    {}

    /// Allocate BlockedMatrix, but let data arrays point to existing arrays
    /// \param[in] Nb             number of blockrows
    /// \param[in] nnzbs          number of nonzero blocks
    /// \param[in] block_size     the number of rows and columns for each block
    /// \param[in] nnzValues      array of nonzero values, contains nnzb*block_size*block_size scalars
    /// \param[in] colIndices     array of column indices, contains nnzb entries
    /// \param[in] rowPointers    array of row pointers, contains Nb+1 entries
    BlockedMatrix(int Nb_, int nnzbs_, unsigned int block_size_, double *nnzValues_, int *colIndices_, int *rowPointers_)
    : nnzValues(nnzValues_),
      colIndices(colIndices_),
      rowPointers(rowPointers_),
      Nb(Nb_),
      nnzbs(nnzbs_),
      block_size(block_size_),
      deleteNnzs(false),
      deleteSparsity(false)
    {}

    ~BlockedMatrix(){
        if (deleteNnzs) {
            delete[] nnzValues;
        }
        if (deleteSparsity) {
            delete[] colIndices;
            delete[] rowPointers;
        }
    }


    double *nnzValues;
    int *colIndices;
    int *rowPointers;
    int Nb;
    int nnzbs;
    unsigned int block_size;
    bool deleteNnzs;
    bool deleteSparsity;
};


/// Sort a row of matrix elements from a CSR-format, where the nonzeroes are ints
/// These ints aren't actually nonzeroes, but represent a mapping used later
/// \param[inout] colIndices     represent keys in sorting
/// \param[inout] data           sorted according to the colIndices
/// \param[in] left              lower index of data of row
/// \param[in] right             upper index of data of row
void sortRow(int *colIndices, int *data, int left, int right);

/// Multiply and subtract blocks
/// a = a - (b * c)
/// \param[inout] a              block to be subtracted from
/// \param[in] b                 input block
/// \param[in] c                 input block
/// \param[in] block_size        size of block
void blockMultSub(double *a, double *b, double *c, unsigned int block_size);

/// Perform a matrix-matrix multiplication on two blocks
/// resMat = mat1 * mat2
/// \param[in] mat1              input block 1
/// \param[in] mat2              input block 2
/// \param[out] resMat           output block
/// \param[in] block_size        size of block
void blockMult(double *mat1, double *mat2, double *resMat, unsigned int block_size);

} // namespace Accelerator
} // namespace Opm

#endif
