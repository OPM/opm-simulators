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

#ifndef OPM_OPENCLCPR_HPP
#define OPM_OPENCLCPR_HPP

#include <dune/istl/paamg/matrixhierarchy.hh>
#include <dune/istl/umfpack.hh>

#include <opm/simulators/linalg/gpubridge/opencl/opencl.hpp>
#include <opm/simulators/linalg/gpubridge/opencl/openclBILU0.hpp>
#include <opm/simulators/linalg/gpubridge/Matrix.hpp>
#include <opm/simulators/linalg/gpubridge/CprCreation.hpp>
#include <opm/simulators/linalg/gpubridge/opencl/OpenclMatrix.hpp>
#include <opm/simulators/linalg/gpubridge/opencl/openclPreconditioner.hpp>

#include <opm/simulators/linalg/gpubridge/opencl/openclSolverBackend.hpp>

#include <type_traits>

namespace Opm::Accelerator {

template<class Scalar> class BlockedMatrix;

/// This class implements a Constrained Pressure Residual (CPR) preconditioner
template<class Scalar, unsigned int block_size>
class openclCPR : public openclPreconditioner<Scalar,block_size>, public CprCreation<Scalar, block_size>
{
    using Base = openclPreconditioner<Scalar,block_size>;

    using Base::N;
    using Base::Nb;
    using Base::nnz;
    using Base::nnzb;
    using Base::verbosity;
    using Base::context;
    using Base::queue;
    using Base::events;
    using Base::err;

private:
    std::vector<OpenclMatrix<Scalar>> d_Amatrices, d_Rmatrices; // scalar matrices that represent the AMG hierarchy
    std::vector<cl::Buffer> d_PcolIndices;
    std::vector<cl::Buffer> d_invDiags;
    std::vector<cl::Buffer> d_t, d_f, d_u; // intermediate vectors used during amg cycle
    std::unique_ptr<cl::Buffer> d_rs;      // use before extracting the pressure
    std::unique_ptr<cl::Buffer> d_weights; // the quasiimpes weights, used to extract pressure
    std::unique_ptr<OpenclMatrix<Scalar>> d_mat;   // stores blocked matrix
    std::unique_ptr<cl::Buffer> d_coarse_y, d_coarse_x; // stores the scalar vectors
    std::once_flag opencl_buffers_allocated;  // only allocate OpenCL Buffers once

    std::unique_ptr<openclBILU0<Scalar,block_size>> bilu0;                    // Blocked ILU0 preconditioner

    std::unique_ptr<openclSolverBackend<Scalar,1>> coarse_solver; // coarse solver is scalar
    bool opencl_ilu_parallel;                   // whether ILU0 operation should be parallelized

    // Initialize and allocate matrices and vectors
    void init_opencl_buffers();

    // Copy matrices and vectors to GPU
    void opencl_upload();

    // apply pressure correction to vector
    void apply_amg(const cl::Buffer& y, cl::Buffer& x);

    void amg_cycle_gpu(const int level, cl::Buffer &y, cl::Buffer &x);

public:
    openclCPR(bool opencl_ilu_parallel, int verbosity);

    bool analyze_matrix(BlockedMatrix<Scalar>* mat) override;
    bool analyze_matrix(BlockedMatrix<Scalar>* mat,
                        BlockedMatrix<Scalar>* jacMat) override;

    // set own Opencl variables, but also that of the bilu0 preconditioner
    void setOpencl(std::shared_ptr<cl::Context>& context,
                   std::shared_ptr<cl::CommandQueue>& queue) override;

    // applies blocked ilu0
    // also applies amg for pressure component
    void apply(const cl::Buffer& y,
               cl::Buffer& x,
               WellContributions<Scalar>& wellContribs) override;

    bool create_preconditioner(BlockedMatrix<Scalar>* mat) override;
    bool create_preconditioner(BlockedMatrix<Scalar>* mat,
                               BlockedMatrix<Scalar>* jacMat) override;
};

} // namespace Opm::Accelerator

#endif
