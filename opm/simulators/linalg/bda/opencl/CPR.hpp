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

#ifndef OPM_CPR_HPP
#define OPM_CPR_HPP

#include <mutex>

#include <dune/istl/paamg/matrixhierarchy.hh>
#include <dune/istl/umfpack.hh>

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/opencl/BILU0.hpp>
#include <opm/simulators/linalg/bda/Matrix.hpp>
#include <opm/simulators/linalg/bda/opencl/OpenclMatrix.hpp>
#include <opm/simulators/linalg/bda/opencl/Preconditioner.hpp>

#include <opm/simulators/linalg/bda/opencl/openclSolverBackend.hpp>

namespace Opm
{
namespace Accelerator
{

class BlockedMatrix;

/// This class implements a Constrained Pressure Residual (CPR) preconditioner
template <unsigned int block_size>
class CPR : public Preconditioner<block_size>
{
    typedef Preconditioner<block_size> Base;

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
    int num_levels;
    std::vector<double> weights, coarse_vals, coarse_x, coarse_y;
    std::vector<Matrix> Amatrices, Rmatrices; // scalar matrices that represent the AMG hierarchy
    std::vector<OpenclMatrix> d_Amatrices, d_Rmatrices; // scalar matrices that represent the AMG hierarchy
    std::vector<std::vector<int> > PcolIndices; // prolongation does not need a full matrix, only store colIndices
    std::vector<cl::Buffer> d_PcolIndices;
    std::vector<std::vector<double> > invDiags; // inverse of diagonal of Amatrices
    std::vector<cl::Buffer> d_invDiags;
    std::vector<cl::Buffer> d_t, d_f, d_u; // intermediate vectors used during amg cycle
    std::unique_ptr<cl::Buffer> d_rs;      // use before extracting the pressure
    std::unique_ptr<cl::Buffer> d_weights; // the quasiimpes weights, used to extract pressure
    std::unique_ptr<OpenclMatrix> d_mat;   // stores blocked matrix
    std::unique_ptr<cl::Buffer> d_coarse_y, d_coarse_x; // stores the scalar vectors
    std::once_flag opencl_buffers_allocated;  // only allocate OpenCL Buffers once

    std::unique_ptr<BILU0<block_size> > bilu0;                    // Blocked ILU0 preconditioner
    BlockedMatrix *mat = nullptr;    // input matrix, blocked

    using DuneMat = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1> >;
    using DuneVec = Dune::BlockVector<Dune::FieldVector<double, 1> >;
    using MatrixOperator = Dune::MatrixAdapter<DuneMat, DuneVec, DuneVec>;
    using DuneAmg = Dune::Amg::MatrixHierarchy<MatrixOperator, Dune::Amg::SequentialInformation>;
    std::unique_ptr<DuneAmg> dune_amg;
    std::unique_ptr<DuneMat> dune_coarse;       // extracted pressure matrix, finest level in AMG hierarchy
    std::shared_ptr<MatrixOperator> dune_op;    // operator, input to Dune AMG
    std::vector<int> level_sizes;               // size of each level in the AMG hierarchy
    std::vector<std::vector<int> > diagIndices; // index of diagonal value for each level
    Dune::UMFPack<DuneMat> umfpack;             // dune/istl/umfpack object used to solve the coarsest level of AMG
    bool always_recalculate_aggregates = false; // OPM always reuses the aggregates by default
    bool recalculate_aggregates = true;         // only rerecalculate if true
    const int pressure_idx = 1;                 // hardcoded to mimic OPM
    unsigned num_pre_smooth_steps;              // number of Jacobi smooth steps before restriction
    unsigned num_post_smooth_steps;             // number of Jacobi smooth steps after prolongation

    std::unique_ptr<openclSolverBackend<1> > coarse_solver; // coarse solver is scalar
    bool opencl_ilu_parallel;                   // whether ILU0 operation should be parallelized

    // Analyze the AMG hierarchy build by Dune
    void analyzeHierarchy();

    // Analyze the aggregateMaps from the AMG hierarchy
    // These can be reused, so only use when recalculate_aggregates is true
    void analyzeAggregateMaps();

    // Initialize and allocate matrices and vectors
    void init_opencl_buffers();

    // Copy matrices and vectors to GPU
    void opencl_upload();

    // apply pressure correction to vector
    void apply_amg(const cl::Buffer& y, cl::Buffer& x);

    void amg_cycle_gpu(const int level, cl::Buffer &y, cl::Buffer &x);

    void create_preconditioner_amg(BlockedMatrix *mat);

public:

    CPR(int verbosity, bool opencl_ilu_parallel);

    bool analyze_matrix(BlockedMatrix *mat) override;
    bool analyze_matrix(BlockedMatrix *mat, BlockedMatrix *jacMat) override;

    // set own Opencl variables, but also that of the bilu0 preconditioner
    void setOpencl(std::shared_ptr<cl::Context>& context, std::shared_ptr<cl::CommandQueue>& queue) override;

    // applies blocked ilu0
    // also applies amg for pressure component
    void apply(const cl::Buffer& y, cl::Buffer& x) override;

    bool create_preconditioner(BlockedMatrix *mat) override;
    bool create_preconditioner(BlockedMatrix *mat, BlockedMatrix *jacMat) override;
};

// solve A^T * x = b
// A should represent a 3x3 matrix
// x and b are vectors with 3 elements
void solve_transposed_3x3(const double *A, const double *b, double *x);

} // namespace Accelerator
} // namespace Opm

#endif

