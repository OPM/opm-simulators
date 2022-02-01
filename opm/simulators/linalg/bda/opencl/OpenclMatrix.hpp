/*
  Copyright 2021 Equinor ASA

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

#ifndef OPM_OPENCLMATRIX_HEADER_INCLUDED
#define OPM_OPENCLMATRIX_HEADER_INCLUDED

#include <vector>

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>

namespace Opm
{
namespace Accelerator
{

class Matrix;
class BlockedMatrix;

/// This struct resembles a csr matrix, only doubles are supported
/// The matrix data is stored in OpenCL Buffers
class OpenclMatrix {
public:

    OpenclMatrix(cl::Context *context, int Nb_, int Mb_, int nnzbs_, unsigned int block_size_)
    : Nb(Nb_),
      Mb(Mb_),
      nnzbs(nnzbs_),
      block_size(block_size_)
    {
        nnzValues = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * block_size * block_size * nnzbs);
        colIndices = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * nnzbs);
        rowPointers = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (Nb + 1));
    }

    void upload(cl::CommandQueue *queue, double *vals, int *cols, int *rows);
    void upload(cl::CommandQueue *queue, Matrix *matrix);
    void upload(cl::CommandQueue *queue, BlockedMatrix *matrix);

    cl::Buffer nnzValues;
    cl::Buffer colIndices;
    cl::Buffer rowPointers;
    int Nb, Mb;
    int nnzbs;
    unsigned int block_size;
};

} // namespace Accelerator
} // namespace Opm

#endif // OPM_OPENCLMATRIX_HEADER_INCLUDED
