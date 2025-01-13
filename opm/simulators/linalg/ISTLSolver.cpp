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
#include <opm/simulators/linalg/ISTLSolver.hpp>

#include <dune/istl/schwarz.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <fmt/format.h>

#if COMPILE_GPU_BRIDGE
#include <opm/simulators/linalg/gpubridge/GpuBridge.hpp>
#include <opm/simulators/linalg/gpubridge/WellContributions.hpp>
#endif

namespace Opm {
namespace detail {

#ifdef HAVE_MPI
void copyParValues(std::any& parallelInformation, std::size_t size,
                   Dune::OwnerOverlapCopyCommunication<int,int>& comm)
{
  if (parallelInformation.type() == typeid(ParallelISTLInformation)) {
      const ParallelISTLInformation* parinfo = std::any_cast<ParallelISTLInformation>(&parallelInformation);
      assert(parinfo);
      parinfo->copyValuesTo(comm.indexSet(), comm.remoteIndices(), size, 1);
  }
}
#endif

template<class Matrix>
void makeOverlapRowsInvalid(Matrix& matrix,
                            const std::vector<int>& overlapRows)
{
    //value to set on diagonal
    const int numEq = Matrix::block_type::rows;
    typename Matrix::block_type diag_block(0.0);
    for (int eq = 0; eq < numEq; ++eq)
        diag_block[eq][eq] = 1.0;

    //loop over precalculated overlap rows and columns
    for (const auto row : overlapRows)
    {
        // Zero out row.
        matrix[row] = 0.0;

        //diagonal block set to diag(1.0).
        matrix[row][row] = diag_block;
    }
}

template<class Matrix, class Vector, class Comm>
void FlexibleSolverInfo<Matrix,Vector,Comm>::create(const Matrix& matrix,
                                                    bool parallel,
                                                    const PropertyTree& prm,
                                                    std::size_t pressureIndex,
                                                    std::function<Vector()> weightsCalculator,
                                                    const bool forceSerial,
                                                    [[maybe_unused]] Comm& comm)

{
    // Write sizes of linear systems on all ranks to debug log.
    if (!forceSerial) {
#if HAVE_MPI
        auto basic_comm = comm.communicator();
#else
        auto basic_comm = Dune::Communication<Dune::No_Comm>{};
#endif // HAVE_MPI
        std::ostringstream os;
        os << "Linear system ";
        if (basic_comm.size() > 1) {
            os << fmt::format("on MPI rank: {} ", basic_comm.rank());
        }
        // The static_cast of Matrix::block_type::rows is needed for fmt version 10.
        // TODO: Check if the cast is still needed in future versions.
        os << fmt::format("blocksize: {} size: {:7d} block nonzeroes: {:9d}",
                          static_cast<int>(Matrix::block_type::rows), matrix.N(), matrix.nonzeroes());
        DeferredLogger local_logger;
        local_logger.debug(os.str());
        auto global_logger = gatherDeferredLogger(local_logger, basic_comm);
        if (basic_comm.rank() == 0) {
            global_logger.logMessages();
        }
    }

    if (parallel) {
#if HAVE_MPI
        if (!wellOperator_) {
            using ParOperatorType = Opm::GhostLastMatrixAdapter<Matrix, Vector, Vector, Comm>;
            auto pop = std::make_unique<ParOperatorType>(matrix, comm);
            using FlexibleSolverType = Dune::FlexibleSolver<ParOperatorType>;
            auto sol = std::make_unique<FlexibleSolverType>(*pop, comm, prm,
                                                            weightsCalculator,
                                                            pressureIndex);
            this->pre_ = &sol->preconditioner();
            this->op_ = std::move(pop);
            this->solver_ = std::move(sol);
        } else {
            using ParOperatorType = WellModelGhostLastMatrixAdapter<Matrix, Vector, Vector, true>;
            auto pop = std::make_unique<ParOperatorType>(matrix, *wellOperator_,
                                                         interiorCellNum_);
            using FlexibleSolverType = Dune::FlexibleSolver<ParOperatorType>;
            auto sol = std::make_unique<FlexibleSolverType>(*pop, comm, prm,
                                                            weightsCalculator,
                                                            pressureIndex);
            this->pre_ = &sol->preconditioner();
            this->op_ = std::move(pop);
            this->solver_ = std::move(sol);
        }
#endif
    } else {
        if (!wellOperator_) {
            using SeqOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
            auto sop = std::make_unique<SeqOperatorType>(matrix);
            using FlexibleSolverType = Dune::FlexibleSolver<SeqOperatorType>;
            auto sol = std::make_unique<FlexibleSolverType>(*sop, prm,
                                                            weightsCalculator,
                                                            pressureIndex);
            this->pre_ = &sol->preconditioner();
            this->op_ = std::move(sop);
            this->solver_ = std::move(sol);
        } else {
            using SeqOperatorType = WellModelMatrixAdapter<Matrix, Vector, Vector, false>;
            auto sop = std::make_unique<SeqOperatorType>(matrix, *wellOperator_);
            using FlexibleSolverType = Dune::FlexibleSolver<SeqOperatorType>;
            auto sol = std::make_unique<FlexibleSolverType>(*sop, prm,
                                                            weightsCalculator,
                                                            pressureIndex);
            this->pre_ = &sol->preconditioner();
            this->op_ = std::move(sop);
            this->solver_ = std::move(sol);
        }
    }
}

template<class Scalar, int Dim>
using BM = Dune::BCRSMatrix<MatrixBlock<Scalar,Dim,Dim>>;
template<class Scalar, int Dim>
using BV = Dune::BlockVector<Dune::FieldVector<Scalar,Dim>>;

#if HAVE_MPI
using CommunicationType = Dune::OwnerOverlapCopyCommunication<int,int>;
#else
using CommunicationType = Dune::Communication<int>;
#endif

#define INSTANTIATE_FLEX(T,Dim)                                                           \
    template void makeOverlapRowsInvalid<BM<T,Dim>>(BM<T,Dim>&, const std::vector<int>&); \
    template struct FlexibleSolverInfo<BM<T,Dim>,BV<T,Dim>,CommunicationType>;

#define INSTANTIATE_TYPE(T) \
    INSTANTIATE_FLEX(T,1)   \
    INSTANTIATE_FLEX(T,2)   \
    INSTANTIATE_FLEX(T,3)   \
    INSTANTIATE_FLEX(T,4)   \
    INSTANTIATE_FLEX(T,5)   \
    INSTANTIATE_FLEX(T,6)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

}
}
