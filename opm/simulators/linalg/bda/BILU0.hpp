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

#ifndef BILU0_HPP
#define BILU0_HPP

#include <config.h> // CMake

#include <memory>

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

#if HAVE_OPENCL
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>            // up to OpenCL 1.2
#endif

namespace bda
{

    /// This class implementa a Blocked ILU0 preconditioner
    /// The decomposition is done on CPU, and reorders the rows of the matrix
    template <unsigned int block_size>
    class BILU0
    {

    private:
        int N;       // number of rows of the matrix
        int Nb;      // number of blockrows of the matrix
        int nnz;     // number of nonzeroes of the matrix (scalar)
        int nnzbs;   // number of blocks of the matrix
        BlockedMatrix *LMat, *UMat, *LUMat;
        BlockedMatrix *rMat = nullptr; // only used with PAR_SIM
        Block *invDiagVals;
        int *diagIndex, *rowsPerColor;
        int *toOrder, *fromOrder;
        int numColors;
        int verbosity;

        bool level_scheduling, graph_coloring;

        typedef struct {
            cl::Buffer Lvals, Uvals, invDiagVals;
            cl::Buffer Lcols, Lrows;
            cl::Buffer Ucols, Urows;
            cl::Buffer rowsPerColor;
        } GPU_storage;

        cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *ILU_apply1;
        cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *ILU_apply2;
        GPU_storage s;
        cl::Context *context;
        cl::CommandQueue *queue;
        int work_group_size = 0;
        int total_work_items = 0;
        int lmem_per_work_group = 0;
        bool pattern_uploaded = false;

    public:

        BILU0(bool level_scheduling, bool graph_coloring, int verbosity);

        ~BILU0();

        // analysis
        bool init(BlockedMatrix *mat);

        // ilu_decomposition
        bool create_preconditioner(BlockedMatrix *mat);

        // apply preconditioner, y = prec(x)
        void apply(cl::Buffer& x, cl::Buffer& y);

        void setOpenCLContext(cl::Context *context);
        void setOpenCLQueue(cl::CommandQueue *queue);
        void setKernelParameters(const unsigned int work_group_size, const unsigned int total_work_items, const unsigned int lmem_per_work_group);
        void setKernels(
            cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *ILU_apply1,
            cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *ILU_apply2
            );

        int* getToOrder()
        {
            return toOrder;
        }

        int* getFromOrder()
        {
            return fromOrder;
        }

        BlockedMatrix* getRMat()
        {
            return rMat;
        }

    };

} // end namespace bda

#endif

