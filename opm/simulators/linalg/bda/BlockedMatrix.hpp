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

#ifndef BLOCKED_MATRIX_HPP
#define BLOCKED_MATRIX_HPP

namespace bda
{

/// This struct resembles a blocked csr matrix, like Dune::BCRSMatrix.
/// The data is stored in contiguous memory, such that they can be copied to a device in one transfer.
template<int BS>
struct BlockedMatrix{
    /// Allocate BlockedMatrix and data arrays with given sizes
    /// \param[in] Nb               number of blockrows
    /// \param[in] nnzbs            number of nonzero blocks
    BlockedMatrix(int Nb_, int nnzbs_)
    : nnzValues(new double[nnzbs_*BS*BS]),
      colIndices(new int[nnzbs_*BS*BS]),
      rowPointers(new int[Nb_+1]),
      Nb(Nb_),
      nnzbs(nnzbs_),
      deleteNnzs(true),
      deleteSparsity(true)
    {}
    /// Allocate BlockedMatrix, but copy sparsity pattern instead of allocating new memory
    /// \param[in] M              matrix to be copied
    BlockedMatrix(const BlockedMatrix& M)
    : nnzValues(new double[M.nnzbs*BS*BS]),
      colIndices(M.colIndices),
      rowPointers(M.rowPointers),
      Nb(M.Nb),
      nnzbs(M.nnzbs),
      deleteNnzs(true),
      deleteSparsity(false)
    {}
    /// Allocate BlockedMatrix, but let data arrays point to existing arrays
    /// \param[in] Nb             number of blockrows
    /// \param[in] nnzbs          number of nonzero blocks
    /// \param[in] nnzValues      array of nonzero values, contains nnzb*BS*BS scalars
    /// \param[in] colIndices     array of column indices, contains nnzb entries
    /// \param[in] rowPointers    array of row pointers, contains Nb+1 entries
    BlockedMatrix(int Nb_, int nnzbs_, double *nnzValues_, int *colIndices_, int *rowPointers_)
    : nnzValues(nnzValues_),
      colIndices(colIndices_),
      rowPointers(rowPointers_),
      Nb(Nb_),
      nnzbs(nnzbs_),
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
    bool deleteNnzs;
    bool deleteSparsity;
};


/// Sort a row of matrix elements from a blocked CSR-format
/// \param[inout] colIndices     
/// \param[inout] data           
/// \param[in] left              lower index of data of row
/// \param[in] right             upper index of data of row
template <unsigned int block_size>
void sortBlockedRow(int *colIndices, double *data, int left, int right);

/// Multiply and subtract blocks
/// a = a - (b * c)
/// \param[inout] a              block to be subtracted from
/// \param[in] b                 input block
/// \param[in] c                 input block
template <unsigned int block_size>
void blockMultSub(double *a, double *b, double *c);

/// Perform a 3x3 matrix-matrix multiplication on two blocks
/// \param[in] mat1              input block 1
/// \param[in] mat2              input block 2
/// \param[inout] resMat         output block
template <unsigned int block_size>
void blockMult(double *mat1, double *mat2, double *resMat);

/// Subtract blocks
/// a = b - c
/// \param[out] a                result block
/// \param[in] b                 input block
/// \param[in] c                 input block
template <unsigned int block_size>
void blockSub(double *a, double *b, double *c);

} // end namespace bda

#endif
