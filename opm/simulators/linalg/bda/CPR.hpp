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

#include <opm/simulators/linalg/bda/opencl.hpp>
#include <opm/simulators/linalg/bda/Matrix.hpp>
#include <opm/simulators/linalg/bda/OpenclMatrix.hpp>
#include <opm/simulators/linalg/bda/ILUReorder.hpp>

#include <opm/simulators/linalg/bda/openclKernels.hpp>
#include <opm/simulators/linalg/bda/ChowPatelIlu.hpp>
#include <opm/simulators/linalg/bda/openclSolverBackend.hpp>

namespace Opm
{
namespace Accelerator
{

template <unsigned int block_size>
class openclSolverBackend;

class BlockedMatrix;

/// This class implements a Constrained Pressure Residual (CPR) preconditioner
template <unsigned int block_size>
class CPR
{

private:
    int N;       // number of rows of the matrix
    int Nb;      // number of blockrows of the matrix
    int nnz;     // number of nonzeroes of the matrix (scalar)
    int nnzb;    // number of blocks of the matrix
    int verbosity;

    int num_levels;
    std::vector<double> weights, coarse_vals, coarse_x, coarse_y;
    std::vector<Matrix> Amatrices, Pmatrices, Rmatrices; // scalar matrices that represent the AMG hierarchy
    std::vector<OpenclMatrix> d_Amatrices, d_Pmatrices, d_Rmatrices; // scalar matrices that represent the AMG hierarchy
    std::vector<std::vector<double> > invDiags; // inverse of diagonal of Amatrices
    std::vector<cl::Buffer> d_invDiags;
    std::vector<cl::Buffer> d_t, d_f, d_u; // intermediate vectors used during amg cycle
    std::unique_ptr<cl::Buffer> d_rs;      // use before extracting the pressure
    std::unique_ptr<cl::Buffer> d_weights; // the quasiimpes weights, used to extract pressure
    std::unique_ptr<OpenclMatrix> d_mat;   // stores blocked matrix
    std::unique_ptr<cl::Buffer> d_coarse_y, d_coarse_x; // stores the scalar vectors
    std::once_flag opencl_buffers_allocated;  // only allocate OpenCL Buffers once

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

    std::unique_ptr<openclSolverBackend<1> > coarse_solver; // coarse solver is scalar
    ILUReorder opencl_ilu_reorder;                          // reordering strategy for ILU0 in coarse solver

    std::shared_ptr<cl::Context> context;
    std::shared_ptr<cl::CommandQueue> queue;
    std::vector<cl::Event> events;
    cl_int err;

    // Analyze the AMG hierarchy build by Dune
    void analyzeHierarchy();

    // Analyze the aggregateMaps from the AMG hierarchy
    // These can be reused, so only use when recalculate_aggregates is true
    void analyzeAggregateMaps();

    // Initialize and allocate matrices and vectors
    void init_opencl_buffers();

    // Copy matrices and vectors to GPU
    void opencl_upload();

    void amg_cycle_gpu(const int level, cl::Buffer &y, cl::Buffer &x);

public:

    CPR(int verbosity, ILUReorder opencl_ilu_reorder);

    void init(int Nb, int nnzb, std::shared_ptr<cl::Context>& context, std::shared_ptr<cl::CommandQueue>& queue);

    // apply preconditioner, x = prec(y)
    void apply(const cl::Buffer& y, cl::Buffer& x);

    void create_preconditioner(BlockedMatrix *mat);

};

} // namespace Accelerator
} // namespace Opm

#endif

