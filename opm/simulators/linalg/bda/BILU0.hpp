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

#include <mutex>

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>
#include <opm/simulators/linalg/bda/ILUReorder.hpp>

#include <opm/simulators/linalg/bda/opencl.hpp>
#include <opm/simulators/linalg/bda/openclKernels.hpp>
#include <opm/simulators/linalg/bda/ChowPatelIlu.hpp>

// if CHOW_PATEL is 0, exact ILU decomposition is performed on CPU
// if CHOW_PATEL is 1, iterative ILU decomposition (FGPILU) is done, as described in:
//    FINE-GRAINED PARALLEL INCOMPLETE LU FACTORIZATION, E. Chow and A. Patel, SIAM 2015, https://doi.org/10.1137/140968896
// if CHOW_PATEL_GPU is 0, the decomposition is done on CPU
// if CHOW_PATEL_GPU is 1, the decomposition is done by bda::FGPILU::decomposition() on GPU
// the apply phase of the ChowPatelIlu uses two triangular matrices: L and U
// the exact decomposition uses a full matrix LU which is the superposition of L and U
// ChowPatelIlu could also operate on a full matrix LU when L and U are merged, but it is generally better to keep them split
#define CHOW_PATEL     0
#define CHOW_PATEL_GPU 1


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
        std::unique_ptr<BlockedMatrix<block_size> > LUmat = nullptr;
        std::shared_ptr<BlockedMatrix<block_size> > rmat = nullptr; // only used with PAR_SIM
#if CHOW_PATEL
        std::unique_ptr<BlockedMatrix<block_size> > Lmat = nullptr, Umat = nullptr;
#endif
        double *invDiagVals = nullptr;
        std::vector<int> diagIndex;
        std::vector<int> rowsPerColor;  // color i contains rowsPerColor[i] rows, which are processed in parallel
        std::vector<int> rowsPerColorPrefix;  // the prefix sum of rowsPerColor
        std::vector<int> toOrder, fromOrder;
        int numColors;
        int verbosity;
        std::once_flag pattern_uploaded;

        ILUReorder opencl_ilu_reorder;

        typedef struct {
            cl::Buffer invDiagVals;
            cl::Buffer diagIndex;
            cl::Buffer rowsPerColor;
#if CHOW_PATEL
            cl::Buffer Lvals, Lcols, Lrows;
            cl::Buffer Uvals, Ucols, Urows;
#else
            cl::Buffer LUvals, LUcols, LUrows;
#endif
        } GPU_storage;

        ilu_apply1_kernel_type *ILU_apply1;
        ilu_apply2_kernel_type *ILU_apply2;
        cl::make_kernel<const unsigned int, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                        cl::Buffer&, cl::Buffer&,
                                        const int, cl::LocalSpaceArg> *ilu_decomp_k;

        GPU_storage s;
        cl::Context *context;
        cl::CommandQueue *queue;
        std::vector<cl::Event> events;
        cl_int err;
        int work_group_size = 0;
        int total_work_items = 0;
        int lmem_per_work_group = 0;

        ChowPatelIlu chowPatelIlu;

        void chow_patel_decomposition();

    public:

        BILU0(ILUReorder opencl_ilu_reorder, int verbosity);

        ~BILU0();

        // analysis
        bool init(BlockedMatrix<block_size> *mat);

        // ilu_decomposition
        bool create_preconditioner(BlockedMatrix<block_size> *mat);

        // apply preconditioner, y = prec(x)
        void apply(cl::Buffer& x, cl::Buffer& y);

        void setOpenCLContext(cl::Context *context);
        void setOpenCLQueue(cl::CommandQueue *queue);
        void setKernelParameters(const unsigned int work_group_size, const unsigned int total_work_items, const unsigned int lmem_per_work_group);
        void setKernels(
            cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *ILU_apply1,
            cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg> *ILU_apply2,
            cl::make_kernel<const unsigned int, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const int, cl::LocalSpaceArg> *ilu_decomp_k
            );

        int* getToOrder()
        {
            return toOrder.data();
        }

        int* getFromOrder()
        {
            return fromOrder.data();
        }

        BlockedMatrix<block_size>* getRMat()
        {
            return rmat.get();
        }

    };

} // end namespace bda

#endif

