/*
  Copyright 2016 IRIS AS
  Copyright 2019, 2020 Equinor ASA
  Copyright 2020 SINTEF Digital, Mathematics and Cybernetics

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
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/linalg/ISTLSolverBda.hpp>

#include <dune/istl/schwarz.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <fmt/format.h>

#include <opm/simulators/linalg/bda/BdaBridge.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <ebos/alucartesianindexmapper.hh>
#endif // HAVE_DUNE_ALUGRID

#include <opm/grid/polyhedralgrid.hh>

namespace Opm {
namespace detail {

template<class Matrix, class Vector>
BdaSolverInfo<Matrix,Vector>::
BdaSolverInfo(const std::string& accelerator_mode,
              const int linear_solver_verbosity,
              const int maxit,
              const double tolerance,
              const int platformID,
              const int deviceID,
              const bool opencl_ilu_parallel,
              const std::string& linsolver)
    : bridge_(std::make_unique<Bridge>(accelerator_mode,
                                       linear_solver_verbosity, maxit,
                                       tolerance, platformID, deviceID,
                                       opencl_ilu_parallel, linsolver))
    , accelerator_mode_(accelerator_mode)
{}

template<class Matrix, class Vector>
BdaSolverInfo<Matrix,Vector>::~BdaSolverInfo() = default;

template<class Matrix, class Vector>
template<class Grid>
void BdaSolverInfo<Matrix,Vector>::
prepare(const Grid& grid,
        const Dune::CartesianIndexMapper<Grid>& cartMapper,
        const std::vector<Well>& wellsForConn,
        const std::vector<int>& cellPartition,
        const std::size_t nonzeroes,
        const bool useWellConn)
{
    if (numJacobiBlocks_ > 1) {
      detail::setWellConnections(grid, cartMapper, wellsForConn,
                                 useWellConn,
                                 wellConnectionsGraph_,
                                 numJacobiBlocks_);
      this->blockJacobiAdjacency(grid, cellPartition, nonzeroes);
    }
}

template<class Matrix, class Vector>
bool BdaSolverInfo<Matrix,Vector>::
apply(Vector& rhs,
      const bool useWellConn,
      [[maybe_unused]] WellContribFunc getContribs,
      const int rank,
      Matrix& matrix,
      Vector& x,
      Dune::InverseOperatorResult& result)
{
    bool use_gpu = bridge_->getUseGpu();
    if (use_gpu) {
        auto wellContribs = WellContributions::create(accelerator_mode_, useWellConn);
        bridge_->initWellContributions(*wellContribs, x.N() * x[0].N());

         // the WellContributions can only be applied separately with CUDA, OpenCL or rocsparse, not with amgcl or rocalution
#if HAVE_CUDA || HAVE_OPENCL || HAVE_ROCSPARSE
        if (!useWellConn) {
            getContribs(*wellContribs);
        }
#endif

        if (numJacobiBlocks_ > 1) {
            this->copyMatToBlockJac(matrix, *blockJacobiForGPUILU0_);
            // Const_cast needed since the CUDA stuff overwrites values for better matrix condition..
            bridge_->solve_system(&matrix, blockJacobiForGPUILU0_.get(),
                                  numJacobiBlocks_, rhs, *wellContribs, result);
        }
        else
            bridge_->solve_system(&matrix, &matrix,
                                  numJacobiBlocks_, rhs, *wellContribs, result);
        if (result.converged) {
            // get result vector x from non-Dune backend, iff solve was successful
            bridge_->get_result(x);
            return true;
        } else {
            // warn about CPU fallback
            // BdaBridge might have disabled its BdaSolver for this simulation due to some error
            // in that case the BdaBridge is disabled and flexibleSolver is always used
            // or maybe the BdaSolver did not converge in time, then it will be used next linear solve
            if (rank == 0) {
                OpmLog::warning(bridge_->getAccleratorName() + " did not converge, now trying Dune to solve current linear system...");
            }
        }
    }

    return false;
}

template<class Matrix, class Vector>
bool BdaSolverInfo<Matrix,Vector>::
gpuActive()
{
    return bridge_->getUseGpu();
}

template<class Matrix, class Vector>
template<class Grid>
void BdaSolverInfo<Matrix,Vector>::
blockJacobiAdjacency(const Grid& grid,
                     const std::vector<int>& cell_part,
                     std::size_t nonzeroes)
{
    using size_type = typename Matrix::size_type;
    using Iter = typename Matrix::CreateIterator;
    size_type numCells = grid.size(0);
    blockJacobiForGPUILU0_ = std::make_unique<Matrix>(numCells, numCells,
                                                      nonzeroes, Matrix::row_wise);

    const auto& lid = grid.localIdSet();
    const auto& gridView = grid.leafGridView();
    auto elemIt = gridView.template begin<0>(); // should never overrun, since blockJacobiForGPUILU0_ is initialized with numCells rows

    //Loop over cells
    for (Iter row = blockJacobiForGPUILU0_->createbegin(); row != blockJacobiForGPUILU0_->createend(); ++elemIt, ++row)
    {
        const auto& elem = *elemIt;
        size_type idx = lid.id(elem);
        row.insert(idx);

        // Add well non-zero connections
        for (const auto wc : wellConnectionsGraph_[idx]) {
            row.insert(wc);
        }

        int locPart = cell_part[idx];

        //Add neighbor if it is on the same part
        auto isend = gridView.iend(elem);
        for (auto is = gridView.ibegin(elem); is!=isend; ++is)
        {
            //check if face has neighbor
            if (is->neighbor())
            {
                size_type nid = lid.id(is->outside());
                int nabPart = cell_part[nid];
                if (locPart == nabPart) {
                    row.insert(nid);
                }
            }
        }
    }
}

template<class Matrix, class Vector>
void BdaSolverInfo<Matrix,Vector>::
copyMatToBlockJac(const Matrix& mat, Matrix& blockJac)
{
    auto rbegin = blockJac.begin();
    auto rend = blockJac.end();
    auto outerRow = mat.begin();
    for (auto row = rbegin; row != rend; ++row, ++outerRow) {
        auto outerCol = (*outerRow).begin();
        for (auto col = (*row).begin(); col != (*row).end(); ++col) {
            // outerRow is guaranteed to have all column entries that row has!
            while(outerCol.index() < col.index()) ++outerCol;
            assert(outerCol.index() == col.index());
            *col = *outerCol; // copy nonzero block
        }
    }
}

template<int Dim>
using BM = Dune::BCRSMatrix<MatrixBlock<double,Dim,Dim>>;
template<int Dim>
using BV = Dune::BlockVector<Dune::FieldVector<double,Dim>>;


#define INSTANCE_GRID(Dim, Grid)                   \
    template void BdaSolverInfo<BM<Dim>,BV<Dim>>:: \
    prepare(const Grid&, \
            const Dune::CartesianIndexMapper<Grid>&, \
            const std::vector<Well>&, \
            const std::vector<int>&, \
            const std::size_t, const bool);
using PolyHedralGrid3D = Dune::PolyhedralGrid<3, 3>;
#if HAVE_DUNE_ALUGRID
#if HAVE_MPI
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>;
#else
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI
#define INSTANCE(Dim) \
    template struct BdaSolverInfo<BM<Dim>,BV<Dim>>; \
    INSTANCE_GRID(Dim,Dune::CpGrid) \
    INSTANCE_GRID(Dim,ALUGrid3CN) \
    INSTANCE_GRID(Dim,PolyHedralGrid3D)
#else
#define INSTANCE(Dim) \
    template struct BdaSolverInfo<BM<Dim>,BV<Dim>>; \
    INSTANCE_GRID(Dim,Dune::CpGrid) \
    INSTANCE_GRID(Dim,PolyHedralGrid3D)
#endif
INSTANCE(1)
INSTANCE(2)
INSTANCE(3)
INSTANCE(4)
INSTANCE(5)
INSTANCE(6)

} // namespace detail
} // namespace Opm
