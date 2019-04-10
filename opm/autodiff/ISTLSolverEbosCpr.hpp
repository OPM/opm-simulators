/*
  Copyright 2016 IRIS AS

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

#ifndef OPM_ISTLSOLVERCPR_EBOS_HEADER_INCLUDED
#define OPM_ISTLSOLVERCPR_EBOS_HEADER_INCLUDED

#include <opm/autodiff/ISTLSolverEbos.hpp>
#include <opm/autodiff/BlackoilAmgCpr.hpp>
#include <utility>
#include <memory>

namespace Opm
{
//=====================================================================
// Implementation for ISTL-matrix based operator
//=====================================================================


    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables) for a fixed
    /// number of cell variables np .
    /// \tparam MatrixBlockType The type of the matrix block used.
    /// \tparam VectorBlockType The type of the vector block used.
    /// \tparam pressureIndex The index of the pressure component in the vector
    ///                       vector block. It is used to guide the AMG coarsening.
    ///                       Default is zero.
    template <class TypeTag>
    class ISTLSolverEbosCpr : public ISTLSolverEbos<TypeTag>
    {
    protected:
        // Types and indices from superclass.
        using SuperClass = ISTLSolverEbos<TypeTag>;
        using Matrix = typename SuperClass::Matrix;
        using Vector = typename SuperClass::Vector;
        using WellModel = typename SuperClass::WellModel;
        using Simulator = typename SuperClass::Simulator;
        using SparseMatrixAdapter = typename SuperClass::SparseMatrixAdapter;
        enum { pressureEqnIndex = SuperClass::pressureEqnIndex };
        enum { pressureVarIndex = SuperClass::pressureVarIndex };

        // New properties in this subclass.
        using Preconditioner            = Dune::Preconditioner<Vector, Vector>;
        using MatrixAdapter             = Dune::MatrixAdapter<Matrix,Vector, Vector>;

        using CouplingMetric            = Opm::Amg::Element<pressureEqnIndex,pressureVarIndex>;
        using CritBase                  = Dune::Amg::SymmetricCriterion<Matrix, CouplingMetric>;
        using Criterion                 = Dune::Amg::CoarsenCriterion<CritBase>;

        using ParallelMatrixAdapter     = Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector, Dune::OwnerOverlapCopyCommunication<int,int> >;
        using CprSmootherFine           = Opm::ParallelOverlappingILU0<Matrix, Vector, Vector, Dune::Amg::SequentialInformation>;
        using CprSmootherCoarse         = CprSmootherFine;
        using BlackoilAmgType           = BlackoilAmgCpr<MatrixAdapter,CprSmootherFine, CprSmootherCoarse, Criterion, Dune::Amg::SequentialInformation,
                                                 pressureEqnIndex, pressureVarIndex>;
        using ParallelCprSmootherFine   = Opm::ParallelOverlappingILU0<Matrix, Vector, Vector, Dune::OwnerOverlapCopyCommunication<int,int> >;
        using ParallelCprSmootherCoarse = ParallelCprSmootherFine;
        using ParallelBlackoilAmgType   = BlackoilAmgCpr<ParallelMatrixAdapter, ParallelCprSmootherFine, ParallelCprSmootherCoarse, Criterion,
                                                 Dune::OwnerOverlapCopyCommunication<int,int>, pressureEqnIndex, pressureVarIndex>;

        using OperatorSerial = WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, false>;
        using OperatorParallel = WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, true>;

    public:
        static void registerParameters()
        {
            FlowLinearSolverParameters::registerParameters<TypeTag>();
        }

        /// Construct a system solver.
        /// \param[in] parallelInformation In the case of a parallel run
        ///                                with dune-istl the information about the parallelization.
        explicit ISTLSolverEbosCpr(const Simulator& simulator)
            : SuperClass(simulator), oldMat()
        {
            extractParallelGridInformationToISTL(this->simulator_.vanguard().grid(), this->parallelInformation_);
            detail::findOverlapRowsAndColumns(this->simulator_.vanguard().grid(), this->overlapRowAndColumns_);
        }

        void prepare(const SparseMatrixAdapter& M, Vector& b)
        {
            if (oldMat != nullptr)
                std::cout << "old was "<<oldMat<<" new is "<<&M.istlMatrix()<<std::endl;
            oldMat = &M.istlMatrix();
            int newton_iteration = this->simulator_.model().newtonMethod().numIterations();
            if (newton_iteration < 1 or not(this->parameters_.cpr_reuse_setup_)) {
                SuperClass::matrix_.reset(new Matrix(M.istlMatrix()));
            } else {
                *SuperClass::matrix_ = M.istlMatrix();
            }
            SuperClass::rhs_ = &b;
            SuperClass::scaleSystem();
            const WellModel& wellModel = this->simulator_.problem().wellModel();

#if HAVE_MPI
            if( this->isParallel() ) {

                //remove ghost rows in local matrix without doing a copy.
                this->makeOverlapRowsInvalid(*(this->matrix_));

                if (newton_iteration < 1 or not(this->parameters_.cpr_reuse_setup_)) {
                    //Not sure what actual_mat_for_prec is, so put ebosJacIgnoreOverlap as both variables
                    //to be certain that correct matrix is used for preconditioning.
                    if( ! comm_ )
                    {
                        opAParallel_.reset(new OperatorParallel(*(this->matrix_), *(this->matrix_), wellModel,
                                                                this->parallelInformation_ ));
                        comm_ = opAParallel_->comm();
                        assert(comm_->indexSet().size()==0);
                        const size_t size = opAParallel_->getmat().N();

                        const ParallelISTLInformation& info =
                            boost::any_cast<const ParallelISTLInformation&>( this->parallelInformation_);

                        // As we use a dune-istl with block size np the number of components
                        // per parallel is only one.
                        info.copyValuesTo(comm_->indexSet(), comm_->remoteIndices(),
                                          size, 1);
                    }
                    else
                    {
                        opAParallel_.reset(new OperatorParallel(*(this->matrix_), *(this->matrix_), wellModel,
                                                                comm_ ));
                    }
                }

                using POrComm =  Dune::OwnerOverlapCopyCommunication<int,int>;

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
                constexpr Dune::SolverCategory::Category category=Dune::SolverCategory::overlapping;
                auto sp = Dune::createScalarProduct<Vector,POrComm>(*comm_, category);
                sp_ = std::move(sp);
#else
                constexpr int  category = Dune::SolverCategory::overlapping;
                typedef Dune::ScalarProductChooser<Vector, POrComm, category> ScalarProductChooser;
                typedef std::unique_ptr<typename ScalarProductChooser::ScalarProduct> SPPointer;
                SPPointer sp(ScalarProductChooser::construct(*comm_));
                sp_ = std::move(sp);
#endif

                using AMGOperator = Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector, POrComm>;
                // If clause is always execute as as Linearoperator is WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, false|true>;
                if( ! std::is_same< OperatorParallel, AMGOperator > :: value &&
                    ( newton_iteration < 1 or not(this->parameters_.cpr_reuse_setup_) ) ) {
                    // create new operator in case linear operator and matrix operator differ
                    opA_.reset( new AMGOperator( opAParallel_->getmat(), *comm_ ));
                }

                prepareSolver(*opAParallel_, *comm_);

            } else
#endif
            {

                if (newton_iteration < 1 or not(this->parameters_.cpr_reuse_setup_)) {
                    opASerial_.reset(new OperatorSerial(*(this->matrix_), *(this->matrix_), wellModel));
                }

                using POrComm = Dune::Amg::SequentialInformation;
                POrComm parallelInformation_arg;
                typedef  OperatorSerial LinearOperator;

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
                constexpr Dune::SolverCategory::Category category=Dune::SolverCategory::sequential;
                auto sp = Dune::createScalarProduct<Vector,POrComm>(parallelInformation_arg, category);
                sp_ = std::move(sp);
#else
                constexpr int  category = Dune::SolverCategory::sequential;
                typedef Dune::ScalarProductChooser<Vector, POrComm, category> ScalarProductChooser;
                typedef std::unique_ptr<typename ScalarProductChooser::ScalarProduct> SPPointer;
                SPPointer sp(ScalarProductChooser::construct(parallelInformation_arg));
                sp_ = std::move(sp);
#endif

                // If clause is always execute as as Linearoperator is WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, false|true>;
                if( ! std::is_same< LinearOperator, MatrixAdapter > :: value &&
                    ( newton_iteration < 1 or not(this->parameters_.cpr_reuse_setup_) ) ) {
                    // create new operator in case linear operator and matrix operator differ
                    opA_.reset( new MatrixAdapter( opASerial_->getmat()));//, parallelInformation_arg ) );
                }

                prepareSolver(*opASerial_, parallelInformation_arg);
            }
        }

        template<typename Operator, typename Comm>
        void prepareSolver(Operator& wellOpA, Comm& comm)
        {

            Vector& istlb = *(this->rhs_);
            comm.copyOwnerToAll(istlb, istlb);

            const double relax = this->parameters_.ilu_relaxation_;
            const MILU_VARIANT ilu_milu  = this->parameters_.ilu_milu_;
            using Matrix         = typename MatrixAdapter::matrix_type;
            const int verbosity    = ( this->parameters_.cpr_solver_verbose_ &&
                                       comm.communicator().rank()==0 ) ? 1 : 0;

            // TODO: revise choice of parameters
            // int coarsenTarget = 4000;
            int coarsenTarget = 1200;
            Criterion criterion(15, coarsenTarget);
            criterion.setDebugLevel( this->parameters_.cpr_solver_verbose_ ); // no debug information, 1 for printing hierarchy information
            criterion.setDefaultValuesIsotropic(2);
            criterion.setNoPostSmoothSteps( 1 );
            criterion.setNoPreSmoothSteps( 1 );
            //new guesses by hmbn
            //criterion.setAlpha(0.01); // criterion for connection strong 1/3 is default
            //criterion.setMaxLevel(2); //
            //criterion.setGamma(1); //  //1 V cycle 2 WW

            // Since DUNE 2.2 we also need to pass the smoother args instead of steps directly
            using AmgType           = typename std::conditional<std::is_same<Comm, Dune::Amg::SequentialInformation>::value,
                                                                BlackoilAmgType, ParallelBlackoilAmgType>::type;
            using SpType            = typename std::conditional<std::is_same<Comm, Dune::Amg::SequentialInformation>::value,
                                                                Dune::SeqScalarProduct<Vector>,
                                                                Dune::OverlappingSchwarzScalarProduct<Vector, Comm> >::type;
            using OperatorType      = typename std::conditional<std::is_same<Comm, Dune::Amg::SequentialInformation>::value,
                                                                MatrixAdapter, ParallelMatrixAdapter>::type;
            typedef typename AmgType::Smoother Smoother;
            typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments  SmootherArgs;
            SmootherArgs  smootherArgs;
            smootherArgs.iterations = 1;
            smootherArgs.relaxationFactor = relax;
            const Opm::CPRParameter& params(this->parameters_); // strange conversion
            ISTLUtility::setILUParameters(smootherArgs, ilu_milu);

            auto& opARef = reinterpret_cast<OperatorType&>(*opA_);
            int newton_iteration = this->simulator_.model().newtonMethod().numIterations();
            double dt = this->simulator_.timeStepSize();
            bool update_preconditioner = false;

            if (this->parameters_.cpr_reuse_setup_ < 1) {
                update_preconditioner = true;
            }
            if (this->parameters_.cpr_reuse_setup_ < 2) {
                if (newton_iteration < 1) {
                    update_preconditioner = true;
                }
            }
            if (this->parameters_.cpr_reuse_setup_ < 3) {
                if (this->iterations() > 10) {
                    update_preconditioner = true;
                }
            }

            if ( update_preconditioner or (amg_== 0) ) {
                amg_.reset( new AmgType( params, this->weights_, opARef, criterion, smootherArgs, comm ) );
            } else {
                if (this->parameters_.cpr_solver_verbose_) {
                    std::cout << " Only update amg solver " << std::endl;
                }
                reinterpret_cast<AmgType*>(amg_.get())->updatePreconditioner(opARef, smootherArgs, comm);
            }
            // Solve.
            //SuperClass::solve(linearOperator, x, istlb, *sp, *amg, result);
            //references seems to do something els than refering

            int verbosity_linsolve = 0;
            if (comm.communicator().rank() == 0) {
                verbosity_linsolve = this->parameters_.linear_solver_verbosity_;
            }

            linsolve_.reset(new Dune::BiCGSTABSolver<Vector>(wellOpA, reinterpret_cast<SpType&>(*sp_), reinterpret_cast<AmgType&>(*amg_),
                                                             this->parameters_.linear_solver_reduction_,
                                                             this->parameters_.linear_solver_maxiter_,
                                                             verbosity_linsolve));
        }

        bool solve(Vector& x)
        {
            // Solve system.
            Dune::InverseOperatorResult result;
            Vector& istlb = *(this->rhs_);
            linsolve_->apply(x, istlb, result);
            SuperClass::checkConvergence(result);
            if (this->parameters_.scale_linear_system_) {
                this->scaleSolution(x);
            }
	    return this->converged_;
        }


    protected:

      ///! \brief The dune-istl operator (either serial or parallel
      std::unique_ptr< Dune::LinearOperator<Vector, Vector> > opA_;
      ///! \brief Serial well matrix adapter
      std::unique_ptr< OperatorSerial > opASerial_;
      ///! \brief Parallel well matrix adapter
      std::unique_ptr< OperatorParallel > opAParallel_;
      ///! \brief The preconditoner to use (either serial or parallel CPR with AMG)
      std::unique_ptr< Preconditioner > amg_;
        
      using SPPointer = std::shared_ptr< Dune::ScalarProduct<Vector> >;
      SPPointer sp_;
      std::shared_ptr< Dune::BiCGSTABSolver<Vector> > linsolve_;
      const void* oldMat;
      using POrComm =  Dune::OwnerOverlapCopyCommunication<int,int>;
      std::shared_ptr<POrComm> comm_;
    }; // end ISTLSolver

} // namespace Opm
#endif
