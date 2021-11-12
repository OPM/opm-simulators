/*
  Copyright 2020 Equinor ASA

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

#ifndef MATRIX_HEADER_INCLUDED
#define MATRIX_HEADER_INCLUDED

#include <vector>

#include <opm/simulators/linalg/bda/opencl.hpp>
#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

namespace Opm
{
namespace Accelerator
{

class Matrix;

/// This struct resembles a csr matrix, only doubles are supported
/// The matrix data is stored in OpenCL Buffers
template <unsigned int block_size>
class OpenclMatrix {
public:

    OpenclMatrix(cl::Context *context, int Nb_, int Mb_, int nnzbs_)
    : Nb(Nb_),
      Mb(Mb_),
      nnzbs(nnzbs_)
    {
        nnzValues = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * block_size * block_size * nnzbs);
        colIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * nnzbs);
        rowPointers = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (Nb + 1));
    }

    void upload(cl::CommandQueue *queue, double *vals, int *cols, int *rows);
    void upload(cl::CommandQueue *queue, Matrix *matrix);
    void upload(cl::CommandQueue *queue, BlockedMatrix<block_size> *matrix);

    cl::Buffer nnzValues;
    cl::Buffer colIndices;
    cl::Buffer rowPointers;
    int Nb, Mb;
    int nnzbs;
};


/// This struct resembles a csr matrix, only doubles are supported
/// The data is stored in contiguous memory, such that they can be copied to a device in one transfer.
class Matrix {

public:

    /// Allocate square Matrix and data arrays with given sizes
    /// \param[in] N               number of rows
    /// \param[in] nnzs            number of nonzeros
    Matrix(int N_, int nnzs_)
    : N(N_),
      nnzs(nnzs_)
    {
        nnzValues.resize(nnzs);
        colIndices.resize(nnzs);
        rowPointers.resize(N+1);
    }

    /// Allocate rectangular Matrix and data arrays with given sizes
    /// \param[in] N               number of rows
    /// \param[in] M               number of columns
    /// \param[in] nnzs            number of nonzeros
    Matrix(int N_, int M_, int nnzs_)
    : N(N_),
      M(M_),
      nnzs(nnzs_)
    {
        nnzValues.resize(nnzs);
        colIndices.resize(nnzs);
        rowPointers.resize(N+1);
    }

#if HAVE_FPGA
    /// Converts this matrix to the dataformat used by the FPGA.
    /// The FPGA uses a new data format called CSRO (Compressed Sparse Row Offset).
    /// The purpose of this format is to allow the data to be streamable.
    /// The rowPointers array has an unpredictable reading pattern/timing,
    /// it also needs a extra work if a row is shorter than a cacheline.
    /// The array of N+1 rowPointers is replaced by an array of nnz rowOffsets.
    /// The value of this offset is 0, unless the corresponding nnz is the first of a row,
    /// in that case it is 'the number of empty rows preceeding it + 1'.
    /// The FPGA can simply add the rowOffset to the current rowIdx to get the new rowIdx.
    /// Example:
    /// [1 0 0 3 0]    nnzValues   [1 3 2 2 1 4 3 4 1]
    /// [0 2 2 0 1]    colIndices  [0 3 1 2 4 0 1 2 4]
    /// [4 0 0 0 0] -> rowPointers [0 2 5 6 6 9]
    /// [0 0 0 0 0]    rowOffsets  [1 0 1 0 0 1 2 0 0]
    /// [0 3 4 0 1]
    /// The rowOffset is stored in 1 byte, meaning the maximum value is 255.
    int toRDF(int numColors, std::vector<int>& nodesPerColor,
        std::vector<std::vector<int> >& colIndicesInColor, int nnzsPerRowLimit, 
        std::vector<std::vector<double> >& ubNnzValues, short int *ubColIndices, int *nnzValsSizes, unsigned char *NROffsets, int *colorSizes);
#endif

    std::vector<double> nnzValues;
    std::vector<int> colIndices;
    std::vector<int> rowPointers;
    int N, M;
    int nnzs;
};

void sortRow(int *colIndices, double *data, int left, int right);

} // namespace Accelerator
} // namespace Opm

#endif // MATRIX_HEADER_INCLUDED
