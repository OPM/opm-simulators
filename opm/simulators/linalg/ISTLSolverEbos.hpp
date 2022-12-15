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

#ifndef OPM_ISTLSOLVER_EBOS_HEADER_INCLUDED
#define OPM_ISTLSOLVER_EBOS_HEADER_INCLUDED

#include <dune/istl/owneroverlapcopy.hh>

#include <ebos/eclbasevanguard.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>
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
struct EclWellModel;

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

#if COMPILE_BDA_BRIDGE
template<class Matrix, class Vector, int block_size> class BdaBridge;
class WellContributions;
#endif

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
                size_t pressureIndex,
                std::function<Vector()> trueFunc,
                Comm& comm);

    std::unique_ptr<AbstractSolverType> solver_;
    std::unique_ptr<AbstractOperatorType> op_;
    std::unique_ptr<LinearOperatorExtra<Vector,Vector>> wellOperator_;
    AbstractPreconditionerType* pre_ = nullptr;
    size_t interiorCellNum_ = 0;
};

#if COMPILE_BDA_BRIDGE
template<class Matrix, class Vector>
struct BdaSolverInfo
{
  using WellContribFunc = std::function<void(WellContributions&)>;
  using Bridge = BdaBridge<Matrix,Vector,Matrix::block_type::rows>;

  BdaSolverInfo(const std::string& accelerator_mode,
                const int linear_solver_verbosity,
                const int maxit,
                const double tolerance,
                const int platformID,
                const int deviceID,
                const bool opencl_ilu_parallel,
                const std::string& linsolver);

  ~BdaSolverInfo();

  template<class Grid>
  void prepare(const Grid& grid,
               const Dune::CartesianIndexMapper<Grid>& cartMapper,
               const std::vector<Well>& wellsForConn,
               const std::vector<int>& cellPartition,
               const size_t nonzeroes,
               const bool useWellConn);

  bool apply(Vector& rhs,
             const bool useWellConn,
             WellContribFunc getContribs,
             const int rank,
             Matrix& matrix,
             Vector& x,
             Dune::InverseOperatorResult& result);

  int numJacobiBlocks_ = 0;

private:
  /// Create sparsity pattern for block-Jacobi matrix based on partitioning of grid.
  /// Do not initialize the values, that is done in copyMatToBlockJac()
  template<class Grid>
  void blockJacobiAdjacency(const Grid& grid,
                            const std::vector<int>& cell_part,
                            size_t nonzeroes);

  void copyMatToBlockJac(const Matrix& mat, Matrix& blockJac);

  std::unique_ptr<Bridge> bridge_;
  std::string accelerator_mode_;
  std::unique_ptr<Matrix> blockJacobiForGPUILU0_;
  std::vector<std::set<int>> wellConnectionsGraph_;
};
#endif

#ifdef HAVE_MPI
/// Copy values in parallel.
void copyParValues(std::any& parallelInformation, size_t size,
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
                                             size_t nonzeroes,
                                             const std::vector<std::set<int>>& wellConnectionsGraph);
}

    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables) for a fixed
    /// number of cell variables np .
    template <class TypeTag>
    class ISTLSolverEbos
    {
    protected:
        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
        using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;
        using WellModel = GetPropType<TypeTag, Properties::EclWellModel>;
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
        /// \param[in] parallelInformation In the case of a parallel run
        ///                                with dune-istl the information about the parallelization.
        explicit ISTLSolverEbos(const Simulator& simulator)
            : simulator_(simulator),
              iterations_( 0 ),
              calls_( 0 ),
              converged_(false),
              matrix_()
        {
            const bool on_io_rank = (simulator.gridView().comm().rank() == 0);
#if HAVE_MPI
            comm_.reset( new CommunicationType( simulator_.vanguard().grid().comm() ) );
#endif
            parameters_.template init<TypeTag>();
            prm_ = setupPropertyTree(parameters_,
                                     EWOMS_PARAM_IS_SET(TypeTag, int, LinearSolverMaxIter),
                                     EWOMS_PARAM_IS_SET(TypeTag, int, CprMaxEllIter));

#if COMPILE_BDA_BRIDGE
            {
                std::string accelerator_mode = EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode);
                if ((simulator_.vanguard().grid().comm().size() > 1) && (accelerator_mode != "none")) {
                    if (on_io_rank) {
                        OpmLog::warning("Cannot use GPU with MPI, GPU are disabled");
                    }
                    accelerator_mode = "none";
                }
                const int platformID = EWOMS_GET_PARAM(TypeTag, int, OpenclPlatformId);
                const int deviceID = EWOMS_GET_PARAM(TypeTag, int, BdaDeviceId);
                const int maxit = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
                const double tolerance = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
                const bool opencl_ilu_parallel = EWOMS_GET_PARAM(TypeTag, bool, OpenclIluParallel);
                const int linear_solver_verbosity = parameters_.linear_solver_verbosity_;
                std::string linsolver = EWOMS_GET_PARAM(TypeTag, std::string, LinearSolver);
                bdaBridge = std::make_unique<detail::BdaSolverInfo<Matrix,Vector>>(accelerator_mode,
                                                                                   linear_solver_verbosity,
                                                                                   maxit,
                                                                                   tolerance,
                                                                                   platformID,
                                                                                   deviceID,
                                                                                   opencl_ilu_parallel,
                                                                                   linsolver);
            }
#else
            if (EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode) != "none") {
                OPM_THROW(std::logic_error,"Cannot use accelerated solver since CUDA, OpenCL and amgcl were not found by cmake");
            }
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

            flexibleSolver_.interiorCellNum_ = detail::numMatrixRowsToUseInSolver(simulator_.vanguard().grid(), true);

            // Print parameters to PRT/DBG logs.
            if (on_io_rank) {
                std::ostringstream os;
                os << "Property tree for linear solver:\n";
                prm_.write_json(os, true);
                OpmLog::note(os.str());
            }
        }

        // nothing to clean here
        void eraseMatrix()
        {
        }

        void prepare(const SparseMatrixAdapter& M, Vector& b)
        {
            static bool firstcall = true;
#if HAVE_MPI
            if (firstcall) {
                const size_t size = M.istlMatrix().N();
                detail::copyParValues(parallelInformation_, size, *comm_);
            }
#endif

            // update matrix entries for solvers.
            if (firstcall) {
                // ebos will not change the matrix object. Hence simply store a pointer
                // to the original one with a deleter that does nothing.
                // Outch! We need to be able to scale the linear system! Hence const_cast
                matrix_ = const_cast<Matrix*>(&M.istlMatrix());

                useWellConn_ = EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions);
                // setup sparsity pattern for jacobi matrix for preconditioner (only used for openclSolver)
#if HAVE_OPENCL
                bdaBridge->numJacobiBlocks_ = EWOMS_GET_PARAM(TypeTag, int, NumJacobiBlocks);
                bdaBridge->prepare(simulator_.vanguard().grid(),
                                   simulator_.vanguard().cartesianIndexMapper(),
                                   simulator_.vanguard().schedule().getWellsatEnd(),
                                   simulator_.vanguard().cellPartition(),
                                   getMatrix().nonzeroes(), useWellConn_);
#endif
            } else {
                // Pointers should not change
                if ( &(M.istlMatrix()) != matrix_ ) {
                        OPM_THROW(std::logic_error, "Matrix objects are expected to be reused when reassembling!"
                                  <<" old pointer was " << matrix_ << ", new one is " << (&M.istlMatrix()) );
                }
            }
            rhs_ = &b;

            if (isParallel() && prm_.get<std::string>("preconditioner.type") != "ParOverILU0") {
                detail::makeOverlapRowsInvalid(getMatrix(), overlapRows_);
            }
            prepareFlexibleSolver();
            firstcall = false;
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

        bool solve(Vector& x)
        {
            calls_ += 1;
            // Write linear system if asked for.
            const int verbosity = prm_.get<int>("verbosity", 0);
            const bool write_matrix = verbosity > 10;
            if (write_matrix) {
                Helper::writeSystem(simulator_, //simulator is only used to get names
                                    getMatrix(),
                                    *rhs_,
                                    comm_.get());
            }

            // Solve system.
            Dune::InverseOperatorResult result;

#if COMPILE_BDA_BRIDGE
            std::function<void(WellContributions&)> getContribs =
                [this](WellContributions& w)
                {
                    this->simulator_.problem().wellModel().getWellContributions(w);
                };
            if (!bdaBridge->apply(*rhs_, useWellConn_, getContribs,
                                  simulator_.gridView().comm().rank(),
                                  const_cast<Matrix&>(this->getMatrix()),
                                  x, result))
#endif
            {
                assert(flexibleSolver_.solver_);
                flexibleSolver_.solver_->apply(x, *rhs_, result);
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

    protected:
#if HAVE_MPI
        using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;
#endif

        void checkConvergence( const Dune::InverseOperatorResult& result ) const
        {
            // store number of iterations
            iterations_ = result.iterations;
            converged_ = result.converged;

            // Check for failure of linear solver.
            if (!parameters_.ignoreConvergenceFailure_ && !result.converged) {
                const std::string msg("Convergence failure for linear solver.");
                OPM_THROW_NOLOG(NumericalProblem, msg);
            }
        }
    protected:

        bool isParallel() const {
#if HAVE_MPI
            return comm_->communicator().size() > 1;
#else
            return false;
#endif
        }

        void prepareFlexibleSolver()
        {

            if (shouldCreateSolver()) {
                std::function<Vector()> trueFunc =
                    [this]
                    {
                        return this->getTrueImpesWeights(pressureIndex);
                    };

                if (!useWellConn_) {
                    auto wellOp = std::make_unique<WellModelOperator>(simulator_.problem().wellModel());
                    flexibleSolver_.wellOperator_ = std::move(wellOp);
                }

                flexibleSolver_.create(getMatrix(),
                                       isParallel(),
                                       prm_,
                                       pressureIndex,
                                       trueFunc,
                                       *comm_);
            }
            else
            {
                flexibleSolver_.pre_->update();
            }
        }


        /// Return true if we should (re)create the whole solver,
        /// instead of just calling update() on the preconditioner.
        bool shouldCreateSolver() const
        {
            // Decide if we should recreate the solver or just do
            // a minimal preconditioner update.
            if (!flexibleSolver_.solver_) {
                return true;
            }
            if (this->parameters_.cpr_reuse_setup_ == 0) {
                // Always recreate solver.
                return true;
            }
            if (this->parameters_.cpr_reuse_setup_ == 1) {
                // Recreate solver on the first iteration of every timestep.
                const int newton_iteration = this->simulator_.model().newtonMethod().numIterations();
                return newton_iteration == 0;
            }
            if (this->parameters_.cpr_reuse_setup_ == 2) {
                // Recreate solver if the last solve used more than 10 iterations.
                return this->iterations() > 10;
            }
            if (this->parameters_.cpr_reuse_setup_ == 3) {
                // Recreate solver if the last solve used more than 10 iterations.
                return false;
            }
            if (this->parameters_.cpr_reuse_setup_ == 4) {
                // Recreate solver every 'step' solve calls.
                const int step = this->parameters_.cpr_reuse_interval_;
                const bool create = ((calls_ % step) == 0);
                return create;
            }

            // If here, we have an invalid parameter.
            const bool on_io_rank = (simulator_.gridView().comm().rank() == 0);
            std::string msg = "Invalid value: " + std::to_string(this->parameters_.cpr_reuse_setup_)
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
        Vector getTrueImpesWeights(int pressureVarIndex) const
        {
            Vector weights(rhs_->size());
            ElementContext elemCtx(simulator_);
            Amg::getTrueImpesWeights(pressureVarIndex, weights,
                                     simulator_.vanguard().gridView(),
                                     elemCtx, simulator_.model(),
                                     ThreadManager::threadId());
            return weights;
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
        mutable int calls_;
        mutable bool converged_;
        std::any parallelInformation_;

        // non-const to be able to scale the linear system
        Matrix* matrix_;
        Vector *rhs_;

#if COMPILE_BDA_BRIDGE
        std::unique_ptr<detail::BdaSolverInfo<Matrix, Vector>> bdaBridge;
#endif

        detail::FlexibleSolverInfo<Matrix,Vector,CommunicationType> flexibleSolver_;
        std::vector<int> overlapRows_;
        std::vector<int> interiorRows_;

        bool useWellConn_;

        FlowLinearSolverParameters parameters_;
        PropertyTree prm_;

        std::shared_ptr< CommunicationType > comm_;
    }; // end ISTLSolver

} // namespace Opm
#endif
