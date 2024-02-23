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

#ifndef OPM_ISTLSOLVER_HEADER_INCLUDED
#define OPM_ISTLSOLVER_HEADER_INCLUDED

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solver.hh>

#include <ebos/eclbasevanguard.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/istlsparsematrixadapter.hh>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/WellOperators.hpp>
#include <opm/simulators/linalg/WriteSystemMatrixHelper.hpp>
#include <opm/simulators/linalg/findOverlapRowsAndColumns.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>

#include <any>
#include <cstddef>
#include <functional>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace Opm::Properties {

namespace TTag {
struct FlowIstlSolver {
    using InheritsFrom = std::tuple<FlowIstlSolverParams>;
};
}

template <class TypeTag, class MyTypeTag>
struct WellModel;

//! Set the type of a global jacobian matrix for linear solvers that are based on
//! dune-istl.
template<class TypeTag>
struct SparseMatrixAdapter<TypeTag, TTag::FlowIstlSolver>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    using Block = MatrixBlock<Scalar, numEq, numEq>;

public:
    using type = typename Linear::IstlSparseMatrixAdapter<Block>;
};

} // namespace Opm::Properties

namespace Opm
{


namespace detail
{

template<class Matrix, class Vector, class Comm>
struct FlexibleSolverInfo
{
    using AbstractSolverType = Dune::InverseOperator<Vector, Vector>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
    using AbstractPreconditionerType = Dune::PreconditionerWithUpdate<Vector, Vector>;

    void create(const Matrix& matrix,
                bool parallel,
                const PropertyTree& prm,
                std::size_t pressureIndex,
                std::function<Vector()> weightCalculator,
                const bool forceSerial,
                Comm& comm);

    std::unique_ptr<AbstractSolverType> solver_;
    std::unique_ptr<AbstractOperatorType> op_;
    std::unique_ptr<LinearOperatorExtra<Vector,Vector>> wellOperator_;
    AbstractPreconditionerType* pre_ = nullptr;
    std::size_t interiorCellNum_ = 0;
};


#ifdef HAVE_MPI
/// Copy values in parallel.
void copyParValues(std::any& parallelInformation, std::size_t size,
                   Dune::OwnerOverlapCopyCommunication<int,int>& comm);
#endif

/// Zero out off-diagonal blocks on rows corresponding to overlap cells
/// Diagonal blocks on ovelap rows are set to diag(1.0).
template<class Matrix>
void makeOverlapRowsInvalid(Matrix& matrix,
                            const std::vector<int>& overlapRows);

/// Create sparsity pattern for block-Jacobi matrix based on partitioning of grid.
/// Do not initialize the values, that is done in copyMatToBlockJac()
template<class Matrix, class Grid>
std::unique_ptr<Matrix> blockJacobiAdjacency(const Grid& grid,
                                             const std::vector<int>& cell_part,
                                             std::size_t nonzeroes,
                                             const std::vector<std::set<int>>& wellConnectionsGraph);
}

    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables) for a fixed
    /// number of cell variables np .
    template <class TypeTag>
    class ISTLSolver
    {
    protected:
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
        using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;
        using WellModel = GetPropType<TypeTag, Properties::WellModel>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Matrix = typename SparseMatrixAdapter::IstlMatrix;
        using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
        using AbstractSolverType = Dune::InverseOperator<Vector, Vector>;
        using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
        using AbstractPreconditionerType = Dune::PreconditionerWithUpdate<Vector, Vector>;
        using WellModelOperator = WellModelAsLinearOperator<WellModel, Vector, Vector>;
        using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
        constexpr static std::size_t pressureIndex = GetPropType<TypeTag, Properties::Indices>::pressureSwitchIdx;

#if HAVE_MPI
        using CommunicationType = Dune::OwnerOverlapCopyCommunication<int,int>;
#else
        using CommunicationType = Dune::CollectiveCommunication<int>;
#endif

    public:
        using AssembledLinearOperatorType = Dune::AssembledLinearOperator< Matrix, Vector, Vector >;

        static void registerParameters()
        {
            FlowLinearSolverParameters::registerParameters<TypeTag>();
        }

        /// Construct a system solver.
        /// \param[in] simulator   The opm-models simulator object
        /// \param[in] parameters  Explicit parameters for solver setup, do not
        ///                        read them from command line parameters.
        /// \param[in] forceSerial If true, will set up a serial linear solver only,
        ///                        local to the current rank, instead of creating a
        ///                        parallel (MPI distributed) linear solver.
        ISTLSolver(const Simulator& simulator,
                   const FlowLinearSolverParameters& parameters,
                   bool forceSerial = false)
            : simulator_(simulator),
              iterations_( 0 ),
              converged_(false),
              matrix_(nullptr),
              parameters_{parameters},
              forceSerial_(forceSerial)
        {
            initialize();
        }

        /// Construct a system solver.
        /// \param[in] simulator   The opm-models simulator object
        explicit ISTLSolver(const Simulator& simulator)
            : simulator_(simulator),
              iterations_( 0 ),
              solveCount_(0),
              converged_(false),
              matrix_(nullptr)
        {
            parameters_.resize(1);
            parameters_[0].template init<TypeTag>(simulator_.vanguard().eclState().getSimulationConfig().useCPR());
            initialize();
        }

        void initialize()
        {
            OPM_TIMEBLOCK(IstlSolver);

            if (parameters_[0].linsolver_ == "hybrid") {
                // Experimental hybrid configuration.
                // When chosen, will set up two solvers, one with CPRW
                // and the other with ILU0 preconditioner. More general
                // options may be added later.
                prm_.clear();
                parameters_.clear();
                {
                    FlowLinearSolverParameters para;
                    para.init<TypeTag>(false);
                    para.linsolver_ = "cprw";
                    parameters_.push_back(para);
                    prm_.push_back(setupPropertyTree(parameters_[0],
                                                     EWOMS_PARAM_IS_SET(TypeTag, int, LinearSolverMaxIter),
                                                     EWOMS_PARAM_IS_SET(TypeTag, double, LinearSolverReduction)));
                }
                {
                    FlowLinearSolverParameters para;
                    para.init<TypeTag>(false);
                    para.linsolver_ = "ilu0";
                    parameters_.push_back(para);
                    prm_.push_back(setupPropertyTree(parameters_[1],
                                                     EWOMS_PARAM_IS_SET(TypeTag, int, LinearSolverMaxIter),
                                                     EWOMS_PARAM_IS_SET(TypeTag, double, LinearSolverReduction)));
                }
                // ------------
            } else {
                // Do a normal linear solver setup.
                assert(parameters_.size() == 1);
                assert(prm_.empty());
                prm_.push_back(setupPropertyTree(parameters_[0],
                                                 EWOMS_PARAM_IS_SET(TypeTag, int, LinearSolverMaxIter),
                                                 EWOMS_PARAM_IS_SET(TypeTag, double, LinearSolverReduction)));
            }
            flexibleSolver_.resize(prm_.size());

            const bool on_io_rank = (simulator_.gridView().comm().rank() == 0);
#if HAVE_MPI
            comm_.reset( new CommunicationType( simulator_.vanguard().grid().comm() ) );
#endif
            extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);

            // For some reason simulator_.model().elementMapper() is not initialized at this stage
            //const auto& elemMapper = simulator_.model().elementMapper(); //does not work.
            // Set it up manually
            ElementMapper elemMapper(simulator_.vanguard().gridView(), Dune::mcmgElementLayout());
            detail::findOverlapAndInterior(simulator_.vanguard().grid(), elemMapper, overlapRows_, interiorRows_);
            useWellConn_ = EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions);
            const bool ownersFirst = EWOMS_GET_PARAM(TypeTag, bool, OwnerCellsFirst);
            if (!ownersFirst) {
                const std::string msg = "The linear solver no longer supports --owner-cells-first=false.";
                if (on_io_rank) {
                    OpmLog::error(msg);
                }
                OPM_THROW_NOLOG(std::runtime_error, msg);
            }

            const int interiorCellNum_ = detail::numMatrixRowsToUseInSolver(simulator_.vanguard().grid(), true);
            for (auto& f : flexibleSolver_) {
                f.interiorCellNum_ = interiorCellNum_;
            }

#if HAVE_MPI
            if (isParallel()) {
                const std::size_t size = simulator_.vanguard().grid().leafGridView().size(0);
                detail::copyParValues(parallelInformation_, size, *comm_);
            }
#endif

            // Print parameters to PRT/DBG logs.
            if (on_io_rank && parameters_[activeSolverNum_].linear_solver_print_json_definition_) {
                std::ostringstream os;
                os << "Property tree for linear solvers:\n";
                for (std::size_t i = 0; i<prm_.size(); i++) {
                    prm_[i].write_json(os, true);
                }
                OpmLog::note(os.str());
            }
        }

        // nothing to clean here
        void eraseMatrix()
        {
        }

        void setActiveSolver(const int num)
        {
            if (num > static_cast<int>(prm_.size()) - 1) {
                OPM_THROW(std::logic_error, "Solver number " + std::to_string(num) + " not available.");
            }
            activeSolverNum_ = num;
            if (simulator_.gridView().comm().rank() == 0) {
                OpmLog::debug("Active solver = " + std::to_string(activeSolverNum_)
                              + " (" + parameters_[activeSolverNum_].linsolver_ + ")");
            }
        }

        int numAvailableSolvers()
        {
            return flexibleSolver_.size();
        }

        void initPrepare(const Matrix& M, Vector& b)
        {
            const bool firstcall = (matrix_ == nullptr);

            // update matrix entries for solvers.
            if (firstcall) {
                // model will not change the matrix object. Hence simply store a pointer
                // to the original one with a deleter that does nothing.
                // Outch! We need to be able to scale the linear system! Hence const_cast
                matrix_ = const_cast<Matrix*>(&M);

                useWellConn_ = EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions);
                // setup sparsity pattern for jacobi matrix for preconditioner (only used for openclSolver)
            } else {
                // Pointers should not change
                if ( &M != matrix_ ) {
                        OPM_THROW(std::logic_error,
                                  "Matrix objects are expected to be reused when reassembling!");
                }
            }
            rhs_ = &b;

            // TODO: check all solvers, not just one.
            if (isParallel() && prm_[activeSolverNum_].template get<std::string>("preconditioner.type") != "ParOverILU0") {
                detail::makeOverlapRowsInvalid(getMatrix(), overlapRows_);
            }
        }

        void prepare(const SparseMatrixAdapter& M, Vector& b)
        {
            prepare(M.istlMatrix(), b);
        }

        void prepare(const Matrix& M, Vector& b)
        {
            OPM_TIMEBLOCK(istlSolverPrepare);

            initPrepare(M,b);

            prepareFlexibleSolver();
        }


        void setResidual(Vector& /* b */)
        {
            // rhs_ = &b; // Must be handled in prepare() instead.
        }

        void getResidual(Vector& b) const
        {
            b = *rhs_;
        }

        void setMatrix(const SparseMatrixAdapter& /* M */)
        {
            // matrix_ = &M.istlMatrix(); // Must be handled in prepare() instead.
        }

        int getSolveCount() const {
            return solveCount_;
        }

        void resetSolveCount() {
            solveCount_ = 0;
        }

        bool solve(Vector& x)
        {
            OPM_TIMEBLOCK(istlSolverSolve);
            ++solveCount_;
            // Write linear system if asked for.
            const int verbosity = prm_[activeSolverNum_].get("verbosity", 0);
            const bool write_matrix = verbosity > 10;
            if (write_matrix) {
                Helper::writeSystem(simulator_, //simulator is only used to get names
                                    getMatrix(),
                                    *rhs_,
                                    comm_.get());
            }

            // Solve system.
            Dune::InverseOperatorResult result;
            {
                OPM_TIMEBLOCK(flexibleSolverApply);
                assert(flexibleSolver_[activeSolverNum_].solver_);
                flexibleSolver_[activeSolverNum_].solver_->apply(x, *rhs_, result);
            }

            // Check convergence, iterations etc.
            checkConvergence(result);

            return converged_;
        }


        /// Solve the system of linear equations Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x

        /// \copydoc NewtonIterationBlackoilInterface::iterations
        int iterations () const { return iterations_; }

        /// \copydoc NewtonIterationBlackoilInterface::parallelInformation
        const std::any& parallelInformation() const { return parallelInformation_; }

        const CommunicationType* comm() const { return comm_.get(); }

    protected:
#if HAVE_MPI
        using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;
#endif

        void checkConvergence( const Dune::InverseOperatorResult& result ) const
        {
            // store number of iterations
            iterations_ = result.iterations;
            converged_ = result.converged;
            if(!converged_){
                if(result.reduction < parameters_[activeSolverNum_].relaxed_linear_solver_reduction_){
                    std::stringstream ss;
                    ss<< "Full linear solver tolerance not achieved. The reduction is:" << result.reduction
                      << " after " << result.iterations << " iterations ";
                    OpmLog::warning(ss.str());
                    converged_ = true;
                }
            }
            // Check for failure of linear solver.
            if (!parameters_[activeSolverNum_].ignoreConvergenceFailure_ && !converged_) {
                const std::string msg("Convergence failure for linear solver.");
                OPM_THROW_NOLOG(NumericalProblem, msg);
            }
        }
    protected:

        bool isParallel() const {
#if HAVE_MPI
            return !forceSerial_ && comm_->communicator().size() > 1;
#else
            return false;
#endif
        }

        void prepareFlexibleSolver()
        {
            OPM_TIMEBLOCK(flexibleSolverPrepare);
            if (shouldCreateSolver()) {
                if (!useWellConn_) {
                    auto wellOp = std::make_unique<WellModelOperator>(simulator_.problem().wellModel());
                    flexibleSolver_[activeSolverNum_].wellOperator_ = std::move(wellOp);
                }
                std::function<Vector()> weightCalculator = this->getWeightsCalculator(prm_[activeSolverNum_], getMatrix(), pressureIndex);
                OPM_TIMEBLOCK(flexibleSolverCreate);
                flexibleSolver_[activeSolverNum_].create(getMatrix(),
                                                         isParallel(),
                                                         prm_[activeSolverNum_],
                                                         pressureIndex,
                                                         weightCalculator,
                                                         forceSerial_,
                                                         *comm_);
            }
            else
            {
                OPM_TIMEBLOCK(flexibleSolverUpdate);
                flexibleSolver_[activeSolverNum_].pre_->update();
            }
        }


        /// Return true if we should (re)create the whole solver,
        /// instead of just calling update() on the preconditioner.
        bool shouldCreateSolver() const
        {
            // Decide if we should recreate the solver or just do
            // a minimal preconditioner update.
            if (flexibleSolver_.empty()) {
                return true;
            }
            if (!flexibleSolver_[activeSolverNum_].solver_) {
                return true;
            }
            if (this->parameters_[activeSolverNum_].cpr_reuse_setup_ == 0) {
                // Always recreate solver.
                return true;
            }
            if (this->parameters_[activeSolverNum_].cpr_reuse_setup_ == 1) {
                // Recreate solver on the first iteration of every timestep.
                const int newton_iteration = this->simulator_.model().newtonMethod().numIterations();
                return newton_iteration == 0;
            }
            if (this->parameters_[activeSolverNum_].cpr_reuse_setup_ == 2) {
                // Recreate solver if the last solve used more than 10 iterations.
                return this->iterations() > 10;
            }
            if (this->parameters_[activeSolverNum_].cpr_reuse_setup_ == 3) {
                // Recreate solver if the last solve used more than 10 iterations.
                return false;
            }
            if (this->parameters_[activeSolverNum_].cpr_reuse_setup_ == 4) {
                // Recreate solver every 'step' solve calls.
                const int step = this->parameters_[activeSolverNum_].cpr_reuse_interval_;
                const bool create = ((solveCount_ % step) == 0);
                return create;
            }

            // If here, we have an invalid parameter.
            const bool on_io_rank = (simulator_.gridView().comm().rank() == 0);
            std::string msg = "Invalid value: " + std::to_string(this->parameters_[activeSolverNum_].cpr_reuse_setup_)
                + " for --cpr-reuse-setup parameter, run with --help to see allowed values.";
            if (on_io_rank) {
                OpmLog::error(msg);
            }
            throw std::runtime_error(msg);

            // Never reached.
            return false;
        }


        // Weights to make approximate pressure equations.
        // Calculated from the storage terms (only) of the
        // conservation equations, ignoring all other terms.
        std::function<Vector()> getWeightsCalculator(const PropertyTree& prm,
                                                     const Matrix& matrix,
                                                     std::size_t pressIndex) const
        {
            std::function<Vector()> weightsCalculator;

            using namespace std::string_literals;

            auto preconditionerType = prm.get("preconditioner.type"s, "cpr"s);
            if (preconditionerType == "cpr" || preconditionerType == "cprt"
                || preconditionerType == "cprw" || preconditionerType == "cprwt") {
                const bool transpose = preconditionerType == "cprt" || preconditionerType == "cprwt";
                const auto weightsType = prm.get("preconditioner.weight_type"s, "quasiimpes"s);
                if (weightsType == "quasiimpes") {
                    // weights will be created as default in the solver
                    // assignment p = pressureIndex prevent compiler warning about
                    // capturing variable with non-automatic storage duration
                    weightsCalculator = [matrix, transpose, pressIndex]() {
                        return Amg::getQuasiImpesWeights<Matrix, Vector>(matrix,
                                                                         pressIndex,
                                                                         transpose);
                    };
                } else if ( weightsType == "trueimpes" ) {
                    weightsCalculator =
                        [this, pressIndex]
                        {
                            Vector weights(rhs_->size());
                            ElementContext elemCtx(simulator_);
                            Amg::getTrueImpesWeights(pressIndex, weights,
                                                             simulator_.vanguard().gridView(),
                                                             elemCtx, simulator_.model(),
                                                             ThreadManager::threadId());
                            return weights;
                        };
                } else if  (weightsType == "trueimpesanalytic" ) {
                    weightsCalculator =
                        [this, pressIndex]
                        {
                            Vector weights(rhs_->size());
                            ElementContext elemCtx(simulator_);
                            Amg::getTrueImpesWeightsAnalytic(pressIndex, weights,
                                                             simulator_.vanguard().gridView(),
                                                             elemCtx, simulator_.model(),
                                                             ThreadManager::threadId());
                            return weights;
                        };
                } else {
                    OPM_THROW(std::invalid_argument,
                              "Weights type " + weightsType +
                              "not implemented for cpr."
                              " Please use quasiimpes, trueimpes or trueimpesanalytic.");
                }
            }
            return weightsCalculator;
        }


        Matrix& getMatrix()
        {
            return *matrix_;
        }

        const Matrix& getMatrix() const
        {
            return *matrix_;
        }

        const Simulator& simulator_;
        mutable int iterations_;
        mutable int solveCount_;
        mutable bool converged_;
        std::any parallelInformation_;

        // non-const to be able to scale the linear system
        Matrix* matrix_;
        Vector *rhs_;

        int activeSolverNum_ = 0;
        std::vector<detail::FlexibleSolverInfo<Matrix,Vector,CommunicationType>> flexibleSolver_;
        std::vector<int> overlapRows_;
        std::vector<int> interiorRows_;

        bool useWellConn_;

        std::vector<FlowLinearSolverParameters> parameters_;
        bool forceSerial_ = false;
        std::vector<PropertyTree> prm_;

        std::shared_ptr< CommunicationType > comm_;
    }; // end ISTLSolver

} // namespace Opm

#endif // OPM_ISTLSOLVER_HEADER_INCLUDED
