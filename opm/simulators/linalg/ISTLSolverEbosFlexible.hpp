/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_ISTLSOLVEREBOSFLEXIBLE_HEADER_INCLUDED
#define OPM_ISTLSOLVEREBOSFLEXIBLE_HEADER_INCLUDED

#include <ewoms/linear/matrixblock.hh>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>

#include <memory>
#include <utility>

BEGIN_PROPERTIES

NEW_TYPE_TAG(FlowIstlSolverFlexible, INHERITS_FROM(FlowIstlSolverParams));

NEW_PROP_TAG(LinearSolverConfiguration);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(SparseMatrixAdapter);
NEW_PROP_TAG(Simulator);

END_PROPERTIES


namespace Opm
{

//=====================================================================
// Implementation for ISTL-matrix based operator
//=====================================================================
/// This class solves the fully implicit black-oil system by
/// solving the reduced system (after eliminating well variables)
/// as a block-structured matrix (one block for all cell variables) for a fixed
/// number of cell variables.
///
/// The solvers and preconditioners used are run-time configurable.
template <class TypeTag>
class ISTLSolverEbosFlexible
{
    using SparseMatrixAdapter = typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter);
    using VectorType = typename GET_PROP_TYPE(TypeTag, GlobalEqVector);
    using Simulator = typename GET_PROP_TYPE(TypeTag, Simulator);
    using MatrixType = typename SparseMatrixAdapter::IstlMatrix;
#if HAVE_MPI
    using Communication = Dune::OwnerOverlapCopyCommunication<int, int>;
#endif
    using SolverType = Dune::FlexibleSolver<MatrixType, VectorType>;


public:
    static void registerParameters()
    {
        FlowLinearSolverParameters::registerParameters<TypeTag>();
    }

    explicit ISTLSolverEbosFlexible(const Simulator& simulator)
        : simulator_(simulator)
        , comm_(nullptr)
    {
        parameters_.template init<TypeTag>();
        prm_ = setupPropertyTree(parameters_);
        extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);
#if HAVE_MPI
        if (parallelInformation_.type() == typeid(ParallelISTLInformation)) {
            // Parallel case.
            const ParallelISTLInformation* parinfo = boost::any_cast<ParallelISTLInformation>(&parallelInformation_);
            assert(parinfo);
            comm_.reset(new Communication(parinfo->communicator()));
        }
#endif
    }

    void eraseMatrix()
    {
    }

    void prepare(const SparseMatrixAdapter& mat, VectorType& b)
    {
#if HAVE_MPI
        static bool firstcall = true;
        if (firstcall && parallelInformation_.type() == typeid(ParallelISTLInformation)) {
            // Parallel case.
            const ParallelISTLInformation* parinfo = boost::any_cast<ParallelISTLInformation>(&parallelInformation_);
            assert(parinfo);
            const size_t size = mat.istlMatrix().N();
            parinfo->copyValuesTo(comm_->indexSet(), comm_->remoteIndices(), size, 1);
            firstcall = false;
        }
#endif
        // Decide if we should recreate the solver or just do
        // a minimal preconditioner update.
        const int newton_iteration = this->simulator_.model().newtonMethod().numIterations();
        bool recreate_solver = false;
        if (this->parameters_.cpr_reuse_setup_ == 0) {
            recreate_solver = true;
        } else if (this->parameters_.cpr_reuse_setup_ == 1) {
            if (newton_iteration == 0) {
                recreate_solver = true;
            }
        } else if (this->parameters_.cpr_reuse_setup_ == 2) {
            if (this->iterations() > 10) {
                recreate_solver = true;
            }
        } else {
            assert(this->parameters_.cpr_reuse_setup_ == 3);
            assert(recreate_solver == false);
        }

        if (recreate_solver || !solver_) {
            if (isParallel()) {
                solver_.reset(new SolverType(prm_, mat.istlMatrix(), *comm_));
            } else {
                solver_.reset(new SolverType(prm_, mat.istlMatrix()));
            }
            rhs_ = b;
        } else {
            solver_->updatePreconditioner();
            rhs_ = b;
        }
    }

    bool solve(VectorType& x)
    {
        solver_->apply(x, rhs_, res_);
        return res_.converged;
    }

    bool isParallel() const
    {
#if HAVE_MPI
        return parallelInformation_.type() == typeid(ParallelISTLInformation);
#else
        return false;
#endif
    }

    int iterations() const
    {
        return res_.iterations;
    }

    void setResidual(VectorType& /* b */)
    {
        // rhs_ = &b; // Must be handled in prepare() instead.
    }

    void setMatrix(const SparseMatrixAdapter& /* M */)
    {
        // matrix_ = &M.istlMatrix(); // Must be handled in prepare() instead.
    }

protected:
    const Simulator& simulator_;

    std::unique_ptr<SolverType> solver_;
    FlowLinearSolverParameters parameters_;
    boost::property_tree::ptree prm_;
    VectorType rhs_;
    Dune::InverseOperatorResult res_;
    boost::any parallelInformation_;
    std::unique_ptr<Communication> comm_;
}; // end ISTLSolverEbosFlexible

} // namespace Opm

#endif // OPM_ISTLSOLVEREBOSFLEXIBLE_HEADER_INCLUDED
