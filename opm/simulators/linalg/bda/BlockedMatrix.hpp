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

#if HAVE_FPGA
#include <vector>
namespace Opm
{
namespace Accelerator
{
    class Matrix;
}
}
#endif


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

#if HAVE_FPGA
    constexpr static double nnzThreshold = 1e-80;  // for unblocking, a nonzero must be bigger than this threshold to be considered a nonzero

    int countUnblockedNnzs();

    void unblock(Matrix *mat, bool isUMatrix);

    /// Converts this matrix to the dataformat used by the FPGA
    /// Is done every linear solve. The exact sparsity pattern can change every time since the zeros are removed during unblocking
    int toRDF(int numColors, int *nodesPerColor, bool isUMatrix,
        std::vector<std::vector<int> >& colIndicesInColor, int nnzsPerRowLimit, int *nnzValsSizes,
        std::vector<std::vector<double> >& nnzValues, short int *colIndices, unsigned char *NROffsets, int *colorSizes, int *valSize);

    /// Analyses the sparsity pattern and prepares for toRDF()
    /// Is only called once
    int findPartitionColumns(int numColors, int *nodesPerColor,
        int rowsPerColorLimit, int columnsPerColorLimit,
        std::vector<std::vector<int> >& colIndicesInColor, int *PIndicesAddr, int *colorSizes,
        std::vector<std::vector<int> >& LColIndicesInColor, int *LPIndicesAddr, int *LColorSizes,
        std::vector<std::vector<int> >& UColIndicesInColor, int *UPIndicesAddr, int *UColorSizes);
#endif

    double *nnzValues;
    int *colIndices;
    int *rowPointers;
    int Nb;
    int nnzbs;
    unsigned int block_size;
    bool deleteNnzs;
    bool deleteSparsity;
};


/// Sort a row of matrix elements from a blocked CSR-format
/// \param[inout] colIndices     
/// \param[inout] data           
/// \param[in] left              lower index of data of row
/// \param[in] right             upper index of data of row
/// \param[in] block_size        size of blocks in the row
void sortBlockedRow(int *colIndices, double *data, int left, int right, unsigned block_size);

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


#if HAVE_FPGA
template <unsigned int block_size>
void blockSub(double *mat1, double *mat2, double *resMat);

template <unsigned int block_size>
void blockVectMult(double *mat, double *vect, double scale, double *resVect, bool resetRes);

/// Convert a blocked inverse diagonal to the FPGA format.
/// This is the only blocked structure on the FPGA, since it needs blocked matrix-vector multiplication after the backwards substitution of U.
/// Since the rows of U are reversed, the rows of the diag are also reversed.
/// The cachelines can hold 8 doubles, a block has 9 doubles.
/// The format converts 3x3 blocks to 3x4 blocks, so 1 cacheline holds 2 unblocked rows.
/// Then 2 blocks (24 doubles) fit on 3 cachelines.
/// Example:
/// [1 2 3]    [1 2 3 0]              [1 2 3 0 4 5 6 0]
/// [4 5 6] -> [4 5 6 0] -> hardware: [7 8 9 0 block2 row1]
/// [7 8 9]    [7 8 9 0]              [block2 row2 block2 row3]
void blockedDiagtoRDF(double *blockedDiagVals, int rowSize, int numColors, std::vector<int>& rowsPerColor, double *RDFDiag);
#endif

} // namespace Accelerator
} // namespace Opm

#endif
