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

#include <config.h>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/opencl/OpenclMatrix.hpp>
#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>
#include <opm/simulators/linalg/bda/Matrix.hpp>

namespace Opm
{
namespace Accelerator
{

void OpenclMatrix::upload(cl::CommandQueue *queue, double *vals, int *cols, int *rows) {
    std::vector<cl::Event> events(3);

    cl_int err = queue->enqueueWriteBuffer(nnzValues, CL_FALSE, 0, sizeof(double) * block_size * block_size * nnzbs, vals, nullptr, &events[0]);
    err |= queue->enqueueWriteBuffer(colIndices, CL_FALSE, 0, sizeof(int) * nnzbs, cols, nullptr, &events[1]);
    err |= queue->enqueueWriteBuffer(rowPointers, CL_FALSE, 0, sizeof(int) * (Nb + 1), rows, nullptr, &events[2]);

    cl::WaitForEvents(events);
    events.clear();
    if (err != CL_SUCCESS) {
        // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
        OPM_THROW(std::logic_error, "OpenclMatrix OpenCL enqueueWriteBuffer error");
    }
}

void OpenclMatrix::upload(cl::CommandQueue *queue, Matrix *matrix) {
    if (block_size != 1) {
        OPM_THROW(std::logic_error, "Error trying to upload a BlockedMatrix to OpenclMatrix with different block_size");
    }

    upload(queue, matrix->nnzValues.data(), matrix->colIndices.data(), matrix->rowPointers.data());
}

void OpenclMatrix::upload(cl::CommandQueue *queue, BlockedMatrix *matrix) {
    if (matrix->block_size != block_size) {
        OPM_THROW(std::logic_error, "Error trying to upload a BlockedMatrix to OpenclMatrix with different block_size");
    }

    upload(queue, matrix->nnzValues, matrix->colIndices, matrix->rowPointers);
}

} // namespace Accelerator
} // namespace Opm
