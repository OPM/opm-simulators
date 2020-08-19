/*
  Copyright 2016 IRIS AS
  Copyright 2019 Equinor ASA

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

#include <opm/common/Exceptions.hpp>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/MatrixBlock.hpp>
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>
#include <opm/simulators/linalg/WellOperators.hpp>
#include <opm/simulators/linalg/WriteSystemMatrixHelper.hpp>
#include <opm/simulators/linalg/findOverlapRowsAndColumns.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>


#if HAVE_CUDA || HAVE_OPENCL
#include <opm/simulators/linalg/bda/BdaBridge.hpp>
#endif

BEGIN_PROPERTIES

NEW_TYPE_TAG(FlowIstlSolver, INHERITS_FROM(FlowIstlSolverParams));

template <class TypeTag, class MyTypeTag>
struct EclWellModel;

//! Set the type of a global jacobian matrix for linear solvers that are based on
//! dune-istl.
SET_PROP(FlowIstlSolver, SparseMatrixAdapter)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Opm::MatrixBlock<Scalar, numEq, numEq> Block;

public:
    typedef typename Opm::Linear::IstlSparseMatrixAdapter<Block> type;
};

END_PROPERTIES

namespace Opm
{

    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables) for a fixed
    /// number of cell variables np .
    template <class TypeTag>
    class ISTLSolverEbos
    {
    protected:
        using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
        using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
        using SparseMatrixAdapter = typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter);
        using Vector = typename GET_PROP_TYPE(TypeTag, GlobalEqVector);
        using WellModel = typename GET_PROP_TYPE(TypeTag, EclWellModel);
        using Simulator = typename GET_PROP_TYPE(TypeTag, Simulator);
        using Matrix = typename SparseMatrixAdapter::IstlMatrix;
        using ThreadManager = typename GET_PROP_TYPE(TypeTag, ThreadManager);
        using ElementContext = typename GET_PROP_TYPE(TypeTag, ElementContext);
        using FlexibleSolverType = Dune::FlexibleSolver<Matrix, Vector>;
        using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
        using WellModelOperator = WellModelAsLinearOperator<WellModel, Vector, Vector>;

#if HAVE_CUDA || HAVE_OPENCL
        static const unsigned int block_size = Matrix::block_type::rows;
        std::unique_ptr<BdaBridge<Matrix, Vector, block_size>> bdaBridge;
#endif

#if HAVE_MPI
        using CommunicationType = Dune::OwnerOverlapCopyCommunication<int,int>;
#else
        using CommunicationType = Dune::CollectiveCommunication< int >;
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
              converged_(false),
              matrix_()
        {
            const bool on_io_rank = (simulator.gridView().comm().rank() == 0);
#if HAVE_MPI
            comm_.reset( new CommunicationType( simulator_.vanguard().grid().comm() ) );
#endif
            parameters_.template init<TypeTag>();
            prm_ = setupPropertyTree<TypeTag>(parameters_);
            const auto& gridForConn = simulator_.vanguard().grid();
#if HAVE_CUDA || HAVE_OPENCL
            std::string gpu_mode = EWOMS_GET_PARAM(TypeTag, std::string, GpuMode);
            int platformID = EWOMS_GET_PARAM(TypeTag, int, OpenclPlatformId);
            int deviceID = EWOMS_GET_PARAM(TypeTag, int, BdaDeviceId);
            if (gridForConn.comm().size() > 1 && gpu_mode.compare("none") != 0) {
                if (on_io_rank) {
                    OpmLog::warning("Warning cannot use GPU with MPI, GPU is disabled");
                }
                gpu_mode = "none";
            }
            const int maxit = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
            const double tolerance = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
            const int linear_solver_verbosity = parameters_.linear_solver_verbosity_;
            bdaBridge.reset(new BdaBridge<Matrix, Vector, block_size>(gpu_mode, linear_solver_verbosity, maxit, tolerance, platformID, deviceID));
#else
            const std::string gpu_mode = EWOMS_GET_PARAM(TypeTag, std::string, GpuMode);
            if (gpu_mode.compare("none") != 0) {
                OPM_THROW(std::logic_error,"Error cannot use GPU solver since neither CUDA nor OpenCL was not found by cmake");
            }
#endif
            extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);

            // For some reason simulator_.model().elementMapper() is not initialized at this stage
            // Hence const auto& elemMapper = simulator_.model().elementMapper(); does not work.
            // Set it up manually
            using ElementMapper =
                Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
            ElementMapper elemMapper(simulator_.vanguard().grid().leafGridView(), Dune::mcmgElementLayout());
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

            interiorCellNum_ = detail::numMatrixRowsToUseInSolver(simulator_.vanguard().grid(), true);

            // Print parameters to PRT/DBG logs.
            if (on_io_rank) {
                std::ostringstream os;
                os << "Property tree for linear solver:\n";
                boost::property_tree::write_json(os, prm_, true);
                OpmLog::note(os.str());
            }
        }

        // nothing to clean here
        void eraseMatrix() {
        }

        void prepare(const SparseMatrixAdapter& M, Vector& b)
        {
            static bool firstcall = true;
#if HAVE_MPI
            if (firstcall && parallelInformation_.type() == typeid(ParallelISTLInformation)) {
                // Parallel case.
                const ParallelISTLInformation* parinfo = std::any_cast<ParallelISTLInformation>(&parallelInformation_);
                assert(parinfo);
                const size_t size = M.istlMatrix().N();
                parinfo->copyValuesTo(comm_->indexSet(), comm_->remoteIndices(), size, 1);
            }
#endif

            // update matrix entries for solvers.
            if (firstcall) {
                // ebos will not change the matrix object. Hence simply store a pointer
                // to the original one with a deleter that does nothing.
                // Outch! We need to be able to scale the linear system! Hence const_cast
                matrix_ = const_cast<Matrix*>(&M.istlMatrix());
            } else {
                // Pointers should not change
                if ( &(M.istlMatrix()) != matrix_ ) {
                        OPM_THROW(std::logic_error, "Matrix objects are expected to be reused when reassembling!"
                                  <<" old pointer was " << matrix_ << ", new one is " << (&M.istlMatrix()) );
                }
            }
            rhs_ = &b;

            if (isParallel()) {
                makeOverlapRowsInvalid(getMatrix());
            }
            prepareFlexibleSolver();
            firstcall = false;
        }


        void setResidual(Vector& /* b */) {
            // rhs_ = &b; // Must be handled in prepare() instead.
        }

        void getResidual(Vector& b) const {
            b = *rhs_;
        }

        void setMatrix(const SparseMatrixAdapter& /* M */) {
            // matrix_ = &M.istlMatrix(); // Must be handled in prepare() instead.
        }

        bool solve(Vector& x) {
            // Write linear system if asked for.
            const int verbosity = prm_.get<int>("verbosity", 0);
            const bool write_matrix = verbosity > 10;
            if (write_matrix) {
                Opm::Helper::writeSystem(simulator_, //simulator is only used to get names
                                         getMatrix(),
                                         *rhs_,
                                         comm_.get());
            }

            // Solve system.
            Dune::InverseOperatorResult result;
            bool gpu_was_used = false;

            // Use GPU if: available, chosen by user, and successful.
#if HAVE_CUDA || HAVE_OPENCL
            bool use_gpu = bdaBridge->getUseGpu();
            if (use_gpu) {
                WellContributions wellContribs;
                if (!useWellConn_) {
                    simulator_.problem().wellModel().getWellContributions(wellContribs);
                }
                // Const_cast needed since the CUDA stuff overwrites values for better matrix condition..
                bdaBridge->solve_system(const_cast<Matrix*>(&getMatrix()), *rhs_, wellContribs, result);
                if (result.converged) {
                    // get result vector x from non-Dune backend, iff solve was successful
                    bdaBridge->get_result(x);
                    gpu_was_used = true;
                } else {
                    // CPU fallback
                    use_gpu = bdaBridge->getUseGpu();  // update value, BdaBridge might have disabled cusparseSolver
                    if (use_gpu && simulator_.gridView().comm().rank() == 0) {
                        OpmLog::warning("cusparseSolver did not converge, now trying Dune to solve current linear system...");
                    }
                }
            }
#endif

            // Otherwise, use flexible istl solver.
            if (!gpu_was_used) {
                assert(flexibleSolver_);
                flexibleSolver_->apply(x, *rhs_, result);
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


        // 3x3 matrix block inversion was unstable at least 2.3 until and including
        // 2.5.0. There may still be some issue with the 4x4 matrix block inversion
        // we therefore still use the block inversion in OPM
        typedef ParallelOverlappingILU0<Dune::BCRSMatrix<Dune::MatrixBlock<typename Matrix::field_type,
                                                                            Matrix::block_type::rows,
                                                                            Matrix::block_type::cols> >,
                                                                            Vector, Vector> SeqPreconditioner;



#if HAVE_MPI
        typedef Dune::OwnerOverlapCopyCommunication<int, int> Comm;
        // 3x3 matrix block inversion was unstable from at least 2.3 until and
        // including 2.5.0
        typedef ParallelOverlappingILU0<Matrix,Vector,Vector,Comm> ParPreconditioner;
#endif

        void checkConvergence( const Dune::InverseOperatorResult& result ) const
        {
            // store number of iterations
            iterations_ = result.iterations;
            converged_ = result.converged;

            // Check for failure of linear solver.
            if (!parameters_.ignoreConvergenceFailure_ && !result.converged) {
                const std::string msg("Convergence failure for linear solver.");
                OPM_THROW_NOLOG(NumericalIssue, msg);
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

            std::function<Vector()> weightsCalculator = getWeightsCalculator();

            if (shouldCreateSolver()) {
                if (isParallel()) {
#if HAVE_MPI
                    if (useWellConn_) {
                        using ParOperatorType = Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector, Comm>;
                        linearOperatorForFlexibleSolver_ = std::make_unique<ParOperatorType>(getMatrix(), *comm_);
                        flexibleSolver_ = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_, *comm_, prm_, weightsCalculator);
                    } else {
                        using ParOperatorType = WellModelGhostLastMatrixAdapter<Matrix, Vector, Vector, true>;
                        wellOperator_ = std::make_unique<WellModelOperator>(simulator_.problem().wellModel());
                        linearOperatorForFlexibleSolver_ = std::make_unique<ParOperatorType>(getMatrix(), *wellOperator_, interiorCellNum_);
                        flexibleSolver_ = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_, *comm_, prm_, weightsCalculator);
                    }
#endif
                } else {
                    if (useWellConn_) {
                        using SeqOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
                        linearOperatorForFlexibleSolver_ = std::make_unique<SeqOperatorType>(getMatrix());
                        flexibleSolver_ = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_, prm_, weightsCalculator);
                    } else {
                        using SeqOperatorType = WellModelMatrixAdapter<Matrix, Vector, Vector, false>;
                        wellOperator_ = std::make_unique<WellModelOperator>(simulator_.problem().wellModel());
                        linearOperatorForFlexibleSolver_ = std::make_unique<SeqOperatorType>(getMatrix(), *wellOperator_);
                        flexibleSolver_ = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_, prm_, weightsCalculator);
                    }
                }
            }
            else
            {
                flexibleSolver_->preconditioner().update();
            }
        }


        /// Return true if we should (re)create the whole solver,
        /// instead of just calling update() on the preconditioner.
        bool shouldCreateSolver() const
        {
            // Decide if we should recreate the solver or just do
            // a minimal preconditioner update.
            if (!flexibleSolver_) {
                return true;
            }
            const int newton_iteration = this->simulator_.model().newtonMethod().numIterations();
            bool recreate_solver = false;
            if (this->parameters_.cpr_reuse_setup_ == 0) {
                // Always recreate solver.
                recreate_solver = true;
            } else if (this->parameters_.cpr_reuse_setup_ == 1) {
                // Recreate solver on the first iteration of every timestep.
                if (newton_iteration == 0) {
                    recreate_solver = true;
                }
            } else if (this->parameters_.cpr_reuse_setup_ == 2) {
                // Recreate solver if the last solve used more than 10 iterations.
                if (this->iterations() > 10) {
                    recreate_solver = true;
                }
            } else {
                assert(this->parameters_.cpr_reuse_setup_ == 3);
                assert(recreate_solver == false);
                // Never recreate solver.
            }
            return recreate_solver;
        }


        /// Return an appropriate weight function if a cpr preconditioner is asked for.
        std::function<Vector()> getWeightsCalculator() const
        {
            std::function<Vector()> weightsCalculator;

            auto preconditionerType = prm_.get("preconditioner.type", "cpr");
            if (preconditionerType == "cpr" || preconditionerType == "cprt") {
                const bool transpose = preconditionerType == "cprt";
                const auto weightsType = prm_.get("preconditioner.weight_type", "quasiimpes");
                const auto pressureIndex = this->prm_.get("preconditioner.pressure_var_index", 1);
                if (weightsType == "quasiimpes") {
                    // weighs will be created as default in the solver
                    weightsCalculator = [this, transpose, pressureIndex]() {
                        return Opm::Amg::getQuasiImpesWeights<Matrix, Vector>(this->getMatrix(), pressureIndex, transpose);
                    };
                } else if (weightsType == "trueimpes") {
                    weightsCalculator = [this, pressureIndex]() {
                        return this->getTrueImpesWeights(pressureIndex);
                    };
                } else {
                    OPM_THROW(std::invalid_argument,
                              "Weights type " << weightsType << "not implemented for cpr."
                              << " Please use quasiimpes or trueimpes.");
                }
            }
            return weightsCalculator;
        }


        // Weights to make approximate pressure equations.
        // Calculated from the storage terms (only) of the
        // conservation equations, ignoring all other terms.
        Vector getTrueImpesWeights(int pressureVarIndex) const
        {
            Vector weights(rhs_->size());
            ElementContext elemCtx(simulator_);
            Opm::Amg::getTrueImpesWeights(pressureVarIndex, weights, simulator_.vanguard().gridView(),
                                          elemCtx, simulator_.model(),
                                          ThreadManager::threadId());
            return weights;
        }


        /// Zero out off-diagonal blocks on rows corresponding to overlap cells
        /// Diagonal blocks on ovelap rows are set to diag(1.0).
        void makeOverlapRowsInvalid(Matrix& matrix) const
        {
            //value to set on diagonal
            const int numEq = Matrix::block_type::rows;
            typename Matrix::block_type diag_block(0.0);
            for (int eq = 0; eq < numEq; ++eq)
                diag_block[eq][eq] = 1.0;

            // Count overlap row elements.
            int overlap_diagonal = 0;
            int overlap_total = 0;
            for (auto row_index : overlapRows_) {
                const auto& row = matrix[row_index];
                for (auto it = row.begin(); it != row.end(); ++it) {
                    ++overlap_total;
                    if (it.index() == row_index) {
                        ++overlap_diagonal;
                    } else {
                        OpmLog::info("** " + std::to_string(row_index) + "  " + std::to_string(it.index()) + " **");
                    }
                }
            }
            std::ostringstream os;
            os << "makeOverlapRowsInvalid() found " << overlap_diagonal << " diagonal overlap entries and "
               << (overlap_total - overlap_diagonal) << " non-diagonal overlap entries.";
            OpmLog::info(os.str());

            //loop over precalculated overlap rows and columns
            for (auto row = overlapRows_.begin(); row != overlapRows_.end(); row++ )
                {
                    int lcell = *row;
                    // Zero out row.
                    matrix[lcell] = 0.0;

                    //diagonal block set to diag(1.0).
                    matrix[lcell][lcell] = diag_block;
                }
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
        mutable bool converged_;
        std::any parallelInformation_;

        // non-const to be able to scale the linear system
        Matrix* matrix_;
        Vector *rhs_;

        std::unique_ptr<FlexibleSolverType> flexibleSolver_;
        std::unique_ptr<AbstractOperatorType> linearOperatorForFlexibleSolver_;
        std::unique_ptr<WellModelAsLinearOperator<WellModel, Vector, Vector>> wellOperator_;
        std::vector<int> overlapRows_;
        std::vector<int> interiorRows_;
        std::vector<std::set<int>> wellConnectionsGraph_;

        bool useWellConn_;
        size_t interiorCellNum_;

        FlowLinearSolverParameters parameters_;
        boost::property_tree::ptree prm_;
        bool scale_variables_;

        std::shared_ptr< CommunicationType > comm_;
    }; // end ISTLSolver

} // namespace Opm
#endif
