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

#include <opm/simulators/linalg/bda/opencl.hpp>


namespace bda
{

    // This class implements a blocked version on GPU of the Fine-Grained Parallel ILU (FGPILU) by Chow and Patel 2015:
    //     FINE-GRAINED PARALLEL INCOMPLETE LU FACTORIZATION, E. Chow and A. Patel, SIAM 2015, https://doi.org/10.1137/140968896
    // only blocksize == 3 is supported
    // decomposition() allocates the cl::Buffers on the first call, these are C++ objects that deallocate automatically
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

        std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                        cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                        cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                        cl::Buffer&, cl::Buffer&,
                                        const int, cl::LocalSpaceArg, cl::LocalSpaceArg> > chow_patel_ilu_sweep_k;

    public:
        /// Executes the ChowPatelIlu sweeps
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
        /// \param[in] verbosity     print verbosity
        void decomposition(
            cl::CommandQueue *queue, cl::Context *context,
            int *Ut_ptrs, int *Ut_idxs, double *Ut_vals, int Ut_nnzbs,
            int *L_rows, int *L_cols, double *L_vals, int L_nnzbs,
            int *LU_rows, int *LU_cols, double *LU_vals, int LU_nnzbs,
            int Nb, int num_sweeps, int verbosity);

    };

} // end namespace bda

#endif // CHOW_PATEL_ILU_HEADER_INCLUDED
