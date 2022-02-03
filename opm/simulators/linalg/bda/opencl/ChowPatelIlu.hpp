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

#ifndef CHOW_PATEL_ILU_HEADER_INCLUDED
#define CHOW_PATEL_ILU_HEADER_INCLUDED


#include <mutex>

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>

// Variables CHOW_PATEL, CHOW_PATEL_GPU and CHOW_PATEL_GPU_PARALLEL are set by CMake
// Pass -DUSE_CHOW_PATEL_ILU=1 to cmake to define CHOW_PATEL and use the iterative ILU decomposition
// Pass -DUSE_CHOW_PATEL_ILU_GPU=1 to run the ILU decomposition sweeps on the GPU
// Pass -DUSE_CHOW_PATEL_ILU_GPU_PARALLEL=1 to use more parallelisation in the GPU kernel, see ChowPatelIlu.cpp

// if CHOW_PATEL is 0, exact ILU decomposition is performed on CPU
// if CHOW_PATEL is 1, iterative ILU decomposition (FGPILU) is done, as described in:
//    FINE-GRAINED PARALLEL INCOMPLETE LU FACTORIZATION, E. Chow and A. Patel, SIAM 2015, https://doi.org/10.1137/140968896
// if CHOW_PATEL_GPU is 0, the decomposition is done on CPU
// if CHOW_PATEL_GPU is 1, the decomposition is done by gpu_decomposition() on GPU
// the apply phase of the ChowPatelIlu uses two triangular matrices: L and U
// the exact decomposition uses a full matrix LU which is the superposition of L and U
// ChowPatelIlu could also operate on a full matrix LU when L and U are merged, but it is generally better to keep them split

#if CHOW_PATEL

namespace Opm
{
namespace Accelerator
{

class BlockedMatrix;

// This class implements a blocked version on GPU of the Fine-Grained Parallel ILU (FGPILU) by Chow and Patel 2015:
//     FINE-GRAINED PARALLEL INCOMPLETE LU FACTORIZATION, E. Chow and A. Patel, SIAM 2015, https://doi.org/10.1137/140968896
// only blocksize == 3 is supported
// decomposition() allocates the cl::Buffers on the first call, these are C++ objects that deallocate automatically
template <unsigned int block_size>
class ChowPatelIlu
{
private:
    cl::Buffer d_Ut_vals, d_L_vals, d_LU_vals;
    cl::Buffer d_Ut_ptrs, d_Ut_idxs;
    cl::Buffer d_L_rows, d_L_cols;
    cl::Buffer d_LU_rows, d_LU_cols;
    cl::Buffer d_Ltmp, d_Utmp;

    cl::Event event;
    std::vector<cl::Event> events;
    cl_int err;
    std::once_flag initialize_flag;
    std::once_flag pattern_uploaded;
    int verbosity = 0;

    std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                    cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                    cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                    cl::Buffer&, cl::Buffer&,
                                    const int, cl::LocalSpaceArg, cl::LocalSpaceArg> > chow_patel_ilu_sweep_k;

public:
    /// Transposes the U matrix
    /// Executes the ChowPatelIlu sweeps for decomposition on CPU
    /// Also uploads the decomposition to the GPU (and sparsity pattern if needed)
    /// This function calls gpu_decomposition() if CHOW_PATEL_GPU is set
    void decomposition(
        cl::CommandQueue *queue, cl::Context *context,
        BlockedMatrix *LUmat, BlockedMatrix *Lmat, BlockedMatrix *Umat,
        double *invDiagVals, std::vector<int>& diagIndex,
        cl::Buffer& d_diagIndex, cl::Buffer& d_invDiagVals,
        cl::Buffer& d_Lvals, cl::Buffer& d_Lcols, cl::Buffer& d_Lrows,
        cl::Buffer& d_Uvals, cl::Buffer& d_Ucols, cl::Buffer& d_Urows);


    /// Executes the ChowPatelIlu sweeps for decomposition on GPU
    /// also copies data from CPU to GPU and GPU to CPU
    /// \param[in] queue         OpenCL commandqueue
    /// \param[in] context       OpenCL context
    /// \param[in] Ut_ptrs       BSC columnpointers
    /// \param[in] Ut_idxs       BSC rowindices
    /// \param[inout] Ut_vals    actual nonzeros for U
    /// \param[in] Ut_nnzbs      number of blocks in U
    /// \param[in] L_rows        BSR rowpointers
    /// \param[in] L_cols        BSR columnindices
    /// \param[inout] L_vals     actual nonzeroes for L
    /// \param[in] L_nnzbs       number of blocks in L
    /// \param[in] LU_rows       BSR rowpointers
    /// \param[in] LU_cols       BSR columnindices
    /// \param[in] LU_vals       actual nonzeroes for LU (original matrix)
    /// \param[in] LU_nnzbs      number of blocks in LU
    /// \param[in] Nb            number of blockrows
    /// \param[in] num_sweeps    number of sweeps to be done
    void gpu_decomposition(
        cl::CommandQueue *queue, cl::Context *context,
        int *Ut_ptrs, int *Ut_idxs, double *Ut_vals, int Ut_nnzbs,
        int *L_rows, int *L_cols, double *L_vals, int L_nnzbs,
        int *LU_rows, int *LU_cols, double *LU_vals, int LU_nnzbs,
        int Nb, int num_sweeps);

    /// Set the verbosity
    void setVerbosity(int verbosity_) {
        this->verbosity = verbosity_;
    }

};

} // namespace Accelerator
} // namespace Opm

#endif // CHOW_PATEL

#endif // CHOW_PATEL_ILU_HEADER_INCLUDED
