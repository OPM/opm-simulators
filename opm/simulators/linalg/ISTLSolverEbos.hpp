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

#if HAVE_CUDA || HAVE_OPENCL || HAVE_FPGA || HAVE_AMGCL
#include <opm/simulators/linalg/bda/BdaBridge.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>
#endif

#include <dune/common/timer.hh>

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
    typedef MatrixBlock<Scalar, numEq, numEq> Block;

public:
    typedef typename Linear::IstlSparseMatrixAdapter<Block> type;
};

} // namespace Opm::Properties

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
        using FlexibleSolverType = Dune::FlexibleSolver<Matrix, Vector>;
        using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
        using WellModelOperator = WellModelAsLinearOperator<WellModel, Vector, Vector>;
        using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
        constexpr static std::size_t pressureIndex = GetPropType<TypeTag, Properties::Indices>::pressureSwitchIdx;

#if HAVE_CUDA || HAVE_OPENCL || HAVE_FPGA || HAVE_AMGCL
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
              matrix_(),
              solveCount_(0)
        {
            const bool on_io_rank = (simulator.gridView().comm().rank() == 0);
#if HAVE_MPI
            comm_.reset( new CommunicationType( simulator_.vanguard().grid().comm() ) );
#endif
            // please put me back!!! parameters_.template init<TypeTag>();
            // please put me back!!! prm_ = setupPropertyTree(parameters_,
            // please put me back!!!                          EWOMS_PARAM_IS_SET(TypeTag, int, LinearSolverMaxIter),
            // please put me back!!!                          EWOMS_PARAM_IS_SET(TypeTag, int, CprMaxEllIter));
            // ------------ the following is hard coded for testing purposes!!!
            prm_.clear();
            parameters_.clear();
            {
                FlowLinearSolverParameters para;
                para.init<TypeTag>();
                para.linsolver_ = "cpr";
                parameters_.push_back(para);
                prm_.push_back(setupPropertyTree(parameters_[0],
                            EWOMS_PARAM_IS_SET(TypeTag, int, LinearSolverMaxIter),
                            EWOMS_PARAM_IS_SET(TypeTag, int, CprMaxEllIter)
                            ));
            }
            {
                FlowLinearSolverParameters para;
                para.init<TypeTag>();
                para.linsolver_ = "ilu0";
                parameters_.push_back(para);
                prm_.push_back(setupPropertyTree(parameters_[1],
                            EWOMS_PARAM_IS_SET(TypeTag, int, LinearSolverMaxIter),
                            EWOMS_PARAM_IS_SET(TypeTag, int, CprMaxEllIter)
                            ));
            }

            flexibleSolver_.resize(prm_.size());
            linearOperatorForFlexibleSolver_.resize(prm_.size());
            wellOperator_.resize(prm_.size());
            // ------------
#if HAVE_CUDA || HAVE_OPENCL || HAVE_FPGA || HAVE_AMGCL
            {
                std::string accelerator_mode = EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode);
                if ((simulator_.vanguard().grid().comm().size() > 1) && (accelerator_mode != "none")) {
                    if (on_io_rank) {
                        OpmLog::warning("Cannot use GPU or FPGA with MPI, GPU/FPGA are disabled");
                    }
                    accelerator_mode = "none";
                }
                const int platformID = EWOMS_GET_PARAM(TypeTag, int, OpenclPlatformId);
                const int deviceID = EWOMS_GET_PARAM(TypeTag, int, BdaDeviceId);
                const int maxit = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
                const double tolerance = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
                const std::string opencl_ilu_reorder = EWOMS_GET_PARAM(TypeTag, std::string, OpenclIluReorder);
                const int linear_solver_verbosity = parameters_[activeSolverNum_].linear_solver_verbosity_;
                std::string fpga_bitstream = EWOMS_GET_PARAM(TypeTag, std::string, FpgaBitstream);
                std::string linsolver = EWOMS_GET_PARAM(TypeTag, std::string, Linsolver);
                bdaBridge.reset(new BdaBridge<Matrix, Vector, block_size>(accelerator_mode, fpga_bitstream, linear_solver_verbosity, maxit, tolerance, platformID, deviceID, opencl_ilu_reorder, linsolver));
            }
#else
            if (EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode) != "none") {
                OPM_THROW(std::logic_error,"Cannot use accelerated solver since CUDA, OpenCL and amgcl were not found by cmake and FPGA was not enabled");
            }
#endif
            extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);

            // For some reason simulator_.model().elementMapper() is not initialized at this stage
            // Hence const auto& elemMapper = simulator_.model().elementMapper(); does not work.
            // Set it up manually
            ElementMapper elemMapper(simulator_.vanguard().gridView(), Dune::mcmgElementLayout());
            detail::findOverlapAndInterior(simulator_.vanguard().grid(), elemMapper, overlapRows_, interiorRows_);

            useWellConn_ = EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions);
#if HAVE_FPGA
            // check usage of MatrixAddWellContributions: for FPGA they must be included
            if (EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode) == "fpga" && !useWellConn_) {
                OPM_THROW(std::logic_error,"fpgaSolver needs --matrix-add-well-contributions=true");
            }
#endif
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
                os << "Property tree for linear solvers:\n";
                for (std::size_t i = 0; i<prm_.size(); i++) {
                    prm_[i].write_json(os, true);
                    std::cerr<< "debug: ["<<i<<"] : " << os.str() <<std::endl;
                }
                OpmLog::note(os.str());
            }
        }

        // nothing to clean here
        void eraseMatrix() {
        }

        void setActiveSolver(const int num)
        {
            if (num>prm_.size()-1) {
                OPM_THROW(std::logic_error,"Solver number"+std::to_string(num)+"not available.");
            }
            activeSolverNum_ = num;
            auto cc = Dune::MPIHelper::getCollectiveCommunication();
            if (cc.rank() == 0) {
                OpmLog::note("Active solver = "+std::to_string(activeSolverNum_)+"\n");
                //std::ostringstream os;
                //os << "Setting active solver:\n";
                //prm_[activeSolverNum_].write_json(os, true);
                //std::cerr<< "debug: " << os.str() <<std::endl;
                //OpmLog::note(os.str());
            }
        }
        int numAvailableSolvers()
        {
            return activeSolverNum_;
        }

        //void chooseActiveSolver()
        //{
        //    activeSolverNum_ = (activeSolverNum_+1)%2;
        //    auto cc = Dune::MPIHelper::getCollectiveCommunication();
        //    if (cc.rank() == 0) {
        //        OpmLog::note("Active solver = "+std::to_string(activeSolverNum_)+"\n");
        //    }
        //}

        void prepare(const SparseMatrixAdapter& M, Vector& b)
        {
            //chooseActiveSolver();
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

            if (isParallel() && prm_[activeSolverNum_].template get<std::string>("preconditioner.type") != "ParOverILU0") {
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

        int getSolveCount() const {
            return solveCount_;
        }
        void resetSolveCount() {
            solveCount_ = 0;
        }

        bool solve(Vector& x) {
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
            bool accelerator_was_used = false;

            // Use GPU if: available, chosen by user, and successful.
            // Use FPGA if: support compiled, chosen by user, and successful.
#if HAVE_CUDA || HAVE_OPENCL || HAVE_FPGA || HAVE_AMGCL
            bool use_gpu = bdaBridge->getUseGpu();
            bool use_fpga = bdaBridge->getUseFpga();
            if (use_gpu || use_fpga) {
                const std::string accelerator_mode = EWOMS_GET_PARAM(TypeTag, std::string, AcceleratorMode);
                auto wellContribs = WellContributions::create(accelerator_mode, useWellConn_);
                bdaBridge->initWellContributions(*wellContribs);

                // the WellContributions can only be applied separately with CUDA or OpenCL, not with an FPGA or amgcl
#if HAVE_CUDA || HAVE_OPENCL
                if (!useWellConn_) {
                    simulator_.problem().wellModel().getWellContributions(*wellContribs);
                }
#endif

                // Const_cast needed since the CUDA stuff overwrites values for better matrix condition..
                bdaBridge->solve_system(const_cast<Matrix*>(&getMatrix()), *rhs_, *wellContribs, result);
                if (result.converged) {
                    // get result vector x from non-Dune backend, iff solve was successful
                    bdaBridge->get_result(x);
                    accelerator_was_used = true;
                } else {
                    // warn about CPU fallback
                    // BdaBridge might have disabled its BdaSolver for this simulation due to some error
                    // in that case the BdaBridge is disabled and flexibleSolver is always used
                    // or maybe the BdaSolver did not converge in time, then it will be used next linear solve
                    if (simulator_.gridView().comm().rank() == 0) {
                        OpmLog::warning(bdaBridge->getAccleratorName() + " did not converge, now trying Dune to solve current linear system...");
                    }
                }
            }
#endif

            // Otherwise, use flexible istl solver.
            if (!accelerator_was_used) {
                assert(flexibleSolver_[activeSolverNum_]);
                flexibleSolver_[activeSolverNum_]->apply(x, *rhs_, result);
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
            if (!parameters_[activeSolverNum_].ignoreConvergenceFailure_ && !result.converged) {
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
                        linearOperatorForFlexibleSolver_[activeSolverNum_] = std::make_unique<ParOperatorType>(getMatrix(), *comm_);
                        flexibleSolver_[activeSolverNum_] = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_[activeSolverNum_], *comm_, prm_[activeSolverNum_], weightsCalculator, pressureIndex);
                    } else {
                        using ParOperatorType = WellModelGhostLastMatrixAdapter<Matrix, Vector, Vector, true>;
                        wellOperator_[activeSolverNum_] = std::make_unique<WellModelOperator>(simulator_.problem().wellModel());
                        linearOperatorForFlexibleSolver_[activeSolverNum_] = std::make_unique<ParOperatorType>(getMatrix(), *wellOperator_[activeSolverNum_], interiorCellNum_);
                        flexibleSolver_[activeSolverNum_] = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_[activeSolverNum_], *comm_, prm_[activeSolverNum_], weightsCalculator, pressureIndex);
                    }
#endif
                } else {
                    if (useWellConn_) {
                        using SeqOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
                        linearOperatorForFlexibleSolver_[activeSolverNum_] = std::make_unique<SeqOperatorType>(getMatrix());
                        flexibleSolver_[activeSolverNum_] = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_[activeSolverNum_], prm_[activeSolverNum_], weightsCalculator, pressureIndex);
                    } else {
                        using SeqOperatorType = WellModelMatrixAdapter<Matrix, Vector, Vector, false>;
                        wellOperator_[activeSolverNum_] = std::make_unique<WellModelOperator>(simulator_.problem().wellModel());
                        linearOperatorForFlexibleSolver_[activeSolverNum_] = std::make_unique<SeqOperatorType>(getMatrix(), *wellOperator_[activeSolverNum_]);
                        flexibleSolver_[activeSolverNum_] = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_[activeSolverNum_], prm_[activeSolverNum_], weightsCalculator, pressureIndex);
                    }
                }
            }
            else
            {
                flexibleSolver_[activeSolverNum_]->preconditioner().update();
            }
        }


        /// Return true if we should (re)create the whole solver,
        /// instead of just calling update() on the preconditioner.
        bool shouldCreateSolver() const
        {
            // Decide if we should recreate the solver or just do
            // a minimal preconditioner update.
            if (flexibleSolver_.size()==0) {
                return true;
            }
            if (!flexibleSolver_[activeSolverNum_]) {
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

            // Otherwise, do not recreate solver.
            assert(this->parameters_[activeSolverNum_].cpr_reuse_setup_ == 3);

            return false;
        }


        /// Return an appropriate weight function if a cpr preconditioner is asked for.
        std::function<Vector()> getWeightsCalculator() const
        {
            std::function<Vector()> weightsCalculator;

            using namespace std::string_literals;

            auto preconditionerType = prm_[activeSolverNum_].template get<std::string>("preconditioner.type", "cpr");
            if (preconditionerType == "cpr" || preconditionerType == "cprt") {
                const bool transpose = preconditionerType == "cprt";
                const auto weightsType = prm_[activeSolverNum_].template get<std::string>("preconditioner.weight_type", "quasiimpes");
                if (weightsType == "quasiimpes") {
                    // weights will be created as default in the solver
                    // assignment p = pressureIndex prevent compiler warning about
                    // capturing variable with non-automatic storage duration
                    weightsCalculator = [this, transpose, p = pressureIndex]() {
                        return Amg::getQuasiImpesWeights<Matrix, Vector>(this->getMatrix(), p, transpose);
                    };
                } else if (weightsType == "trueimpes") {
                    // assignment p = pressureIndex prevent compiler warning about
                    // capturing variable with non-automatic storage duration
                    weightsCalculator = [this, p = pressureIndex]() {
                        return this->getTrueImpesWeights(p);
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
            Amg::getTrueImpesWeights(pressureVarIndex, weights, simulator_.vanguard().gridView(),
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

        int activeSolverNum_ = 0;

        std::vector<std::unique_ptr<FlexibleSolverType>> flexibleSolver_;
        std::vector<std::unique_ptr<AbstractOperatorType>> linearOperatorForFlexibleSolver_;
        std::vector<std::unique_ptr<WellModelAsLinearOperator<WellModel, Vector, Vector>>> wellOperator_;
        std::vector<int> overlapRows_;
        std::vector<int> interiorRows_;
        std::vector<std::set<int>> wellConnectionsGraph_;

        bool useWellConn_;
        size_t interiorCellNum_;

        std::vector<FlowLinearSolverParameters> parameters_;
        std::vector<PropertyTree> prm_;
        bool scale_variables_;

        std::shared_ptr< CommunicationType > comm_;
        int solveCount_;
    }; // end ISTLSolver

} // namespace Opm
#endif
