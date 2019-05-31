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

#include <opm/simulators/linalg/MatrixBlock.hpp>
#include <opm/simulators/linalg/BlackoilAmg.hpp>
#include <opm/simulators/linalg/CPRPreconditioner.hpp>
#include <opm/autodiff/MPIUtilities.hpp>
#include <opm/simulators/linalg/ParallelRestrictedAdditiveSchwarz.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp>
#include <opm/autodiff/BlackoilDetails.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>
#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>

#include <ewoms/common/parametersystem.hh>
#include <ewoms/common/propertysystem.hh>

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

BEGIN_PROPERTIES

NEW_TYPE_TAG(FlowIstlSolver, INHERITS_FROM(FlowIstlSolverParams));

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(SparseMatrixAdapter);
NEW_PROP_TAG(Indices);
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(EclWellModel);

//! Set the type of a global jacobian matrix for linear solvers that are based on
//! dune-istl.
SET_PROP(FlowIstlSolver, SparseMatrixAdapter)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Ewoms::MatrixBlock<Scalar, numEq, numEq> Block;

public:
    typedef typename Ewoms::Linear::IstlSparseMatrixAdapter<Block> type;
};

END_PROPERTIES

namespace Opm
{
template <class DenseMatrix>
DenseMatrix transposeDenseMatrix(const DenseMatrix& M)
{
    DenseMatrix tmp;
    for (int i = 0; i < M.rows; ++i)
        for (int j = 0; j < M.cols; ++j)
            tmp[j][i] = M[i][j];

    return tmp;
}

//=====================================================================
// Implementation for ISTL-matrix based operator
//=====================================================================


/*!
   \brief Adapter to turn a matrix into a linear operator.

   Adapts a matrix to the assembled linear operator interface
 */
template<class M, class X, class Y, class WellModel, bool overlapping >
class WellModelMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
  typedef Dune::AssembledLinearOperator<M,X,Y> BaseType;

public:
  typedef M matrix_type;
  typedef X domain_type;
  typedef Y range_type;
  typedef typename X::field_type field_type;

#if HAVE_MPI
  typedef Dune::OwnerOverlapCopyCommunication<int,int> communication_type;
#else
  typedef Dune::CollectiveCommunication< int > communication_type;
#endif

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
  Dune::SolverCategory::Category category() const override
  {
    return overlapping ?
           Dune::SolverCategory::overlapping : Dune::SolverCategory::sequential;
  }
#else
  enum {
    //! \brief The solver category.
    category = overlapping ?
        Dune::SolverCategory::overlapping :
        Dune::SolverCategory::sequential
  };
#endif

  //! constructor: just store a reference to a matrix
  WellModelMatrixAdapter (const M& A,
                          const M& A_for_precond,
                          const WellModel& wellMod,
                          const boost::any& parallelInformation = boost::any() )
      : A_( A ), A_for_precond_(A_for_precond), wellMod_( wellMod ), comm_()
  {
#if HAVE_MPI
    if( parallelInformation.type() == typeid(ParallelISTLInformation) )
    {
      const ParallelISTLInformation& info =
          boost::any_cast<const ParallelISTLInformation&>( parallelInformation);
      comm_.reset( new communication_type( info.communicator() ) );
    }
#endif
  }
  WellModelMatrixAdapter (const M& A,
                          const M& A_for_precond,
                          const WellModel& wellMod,
                          std::shared_ptr<communication_type> comm )
      : A_( A ), A_for_precond_(A_for_precond), wellMod_( wellMod ), comm_(comm)
  {
  }


  virtual void apply( const X& x, Y& y ) const override
  {
    A_.mv( x, y );

    // add well model modification to y
    wellMod_.apply(x, y );

#if HAVE_MPI
    if( comm_ )
      comm_->project( y );
#endif
  }

  // y += \alpha * A * x
  virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const override
  {
    A_.usmv(alpha,x,y);

    // add scaled well model modification to y
    wellMod_.applyScaleAdd( alpha, x, y );

#if HAVE_MPI
    if( comm_ )
      comm_->project( y );
#endif
  }

  virtual const matrix_type& getmat() const override { return A_for_precond_; }

    std::shared_ptr<communication_type> comm()
  {
      return comm_;
  }

protected:
  const matrix_type& A_ ;
  const matrix_type& A_for_precond_ ;
  const WellModel& wellMod_;
  std::shared_ptr< communication_type > comm_;
};

    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables) for a fixed
    /// number of cell variables np .
    template <class TypeTag>
    class ISTLSolverEbos
    {
    protected:
        typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
        typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
        typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
        typedef typename GET_PROP_TYPE(TypeTag, EclWellModel) WellModel;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
        typedef typename SparseMatrixAdapter::IstlMatrix Matrix;

        typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlockType;
        typedef typename Vector::block_type BlockVector;
        typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
        typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;
        typedef typename GridView::template Codim<0>::Entity Element;
        typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
        // Due to miscibility oil <-> gas the water eqn is the one we can replace with a pressure equation.
        static const bool waterEnabled = Indices::waterEnabled;
        static const int pindex = (waterEnabled) ? BlackOilDefaultIndexTraits::waterCompIdx : BlackOilDefaultIndexTraits::oilCompIdx;
        enum { pressureEqnIndex = pindex };
        enum { pressureVarIndex = Indices::pressureSwitchIdx };
        static const int numEq = Indices::numEq;

    public:
        typedef Dune::AssembledLinearOperator< Matrix, Vector, Vector > AssembledLinearOperatorType;

        static void registerParameters()
        {
            FlowLinearSolverParameters::registerParameters<TypeTag>();
        }

        /// Construct a system solver.
        /// \param[in] parallelInformation In the case of a parallel run
        ///                                with dune-istl the information about the parallelization.
        ISTLSolverEbos(const Simulator& simulator)
            : simulator_(simulator),
              iterations_( 0 ),
              converged_(false)
        {
            parameters_.template init<TypeTag>();
            extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);
            detail::findOverlapRowsAndColumns(simulator_.vanguard().grid(),overlapRowAndColumns_);
        }

        // nothing to clean here
        void eraseMatrix() {
            matrix_for_preconditioner_.reset();
        }

        void prepare(const SparseMatrixAdapter& M, Vector& b)
        {
            matrix_.reset(new Matrix(M.istlMatrix()));
            rhs_ = &b;
            this->scaleSystem();
        }

        void scaleSystem()
        {
            const bool matrix_cont_added = EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions);

            if (matrix_cont_added) {
                bool form_cpr = true;
                if (parameters_.system_strategy_ == "quasiimpes") {
                    weights_ = getQuasiImpesWeights();
                } else if (parameters_.system_strategy_ == "trueimpes") {
                    weights_ = getStorageWeights();
                } else if (parameters_.system_strategy_ == "simple") {
                    BlockVector bvec(1.0);
                    weights_ = getSimpleWeights(bvec);
                } else if (parameters_.system_strategy_ == "original") {
                    BlockVector bvec(0.0);
                    bvec[pressureEqnIndex] = 1;
                    weights_ = getSimpleWeights(bvec);
                } else {
                    if (parameters_.system_strategy_ != "none") {
                        OpmLog::warning("unknown_system_strategy", "Unknown linear solver system strategy: '" + parameters_.system_strategy_ + "', applying 'none' strategy.");
                    }
                    form_cpr = false;
                }
                if (parameters_.scale_linear_system_) {
                    // also scale weights
                    this->scaleEquationsAndVariables(weights_);
                }
                if (form_cpr && !(parameters_.cpr_use_drs_)) {
                    scaleMatrixAndRhs(weights_);
                }
                if (weights_.size() == 0) {
                    // if weights are not set cpr_use_drs_=false;
                    parameters_.cpr_use_drs_ = false;
                }
            } else {
                if (parameters_.use_cpr_ && parameters_.cpr_use_drs_) {
                   OpmLog::warning("DRS_DISABLE", "Disabling DRS as matrix does not contain well contributions");
                }
                parameters_.cpr_use_drs_ = false;
                if (parameters_.scale_linear_system_) {
                    // also scale weights
                    this->scaleEquationsAndVariables(weights_);
                }
            }
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
            // Solve system.

            const WellModel& wellModel = simulator_.problem().wellModel();

            if( isParallel() )
            {
                typedef WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, true > Operator;

                auto ebosJacIgnoreOverlap = Matrix(*matrix_);
                //remove ghost rows in local matrix
                makeOverlapRowsInvalid(ebosJacIgnoreOverlap);

                //Not sure what actual_mat_for_prec is, so put ebosJacIgnoreOverlap as both variables
                //to be certain that correct matrix is used for preconditioning.
                Operator opA(ebosJacIgnoreOverlap, ebosJacIgnoreOverlap, wellModel,
                             parallelInformation_ );
                assert( opA.comm() );
                solve( opA, x, *rhs_, *(opA.comm()) );
            }
            else
            {
                const WellModel& wellModel = simulator_.problem().wellModel();
                typedef WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, false > Operator;
                Operator opA(*matrix_, *matrix_, wellModel);
                solve( opA, x, *rhs_ );
            }

            if (parameters_.scale_linear_system_) {
                scaleSolution(x);
            }

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
        const boost::any& parallelInformation() const { return parallelInformation_; }

    protected:
        /// \brief construct the CPR preconditioner and the solver.
        /// \tparam P The type of the parallel information.
        /// \param parallelInformation the information about the parallelization.
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
        template<Dune::SolverCategory::Category category=Dune::SolverCategory::sequential,
                 class LinearOperator, class POrComm>
#else
        template<int category=Dune::SolverCategory::sequential, class LinearOperator, class POrComm>
#endif
        void constructPreconditionerAndSolve(LinearOperator& linearOperator,
                                             Vector& x, Vector& istlb,
                                             const POrComm& parallelInformation_arg,
                                             Dune::InverseOperatorResult& result) const
        {
            // Construct scalar product.
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
            auto sp = Dune::createScalarProduct<Vector,POrComm>(parallelInformation_arg, category);
#else
            typedef Dune::ScalarProductChooser<Vector, POrComm, category> ScalarProductChooser;
            typedef std::unique_ptr<typename ScalarProductChooser::ScalarProduct> SPPointer;
            SPPointer sp(ScalarProductChooser::construct(parallelInformation_arg));
#endif

            // Communicate if parallel.
            parallelInformation_arg.copyOwnerToAll(istlb, istlb);

#if FLOW_SUPPORT_AMG // activate AMG if either flow_ebos is used or UMFPack is not available
            if( parameters_.linear_solver_use_amg_ || parameters_.use_cpr_)
            {
                typedef ISTLUtility::CPRSelector< Matrix, Vector, Vector, POrComm>  CPRSelectorType;
                typedef typename CPRSelectorType::Operator MatrixOperator;

                std::unique_ptr< MatrixOperator > opA;

                if( ! std::is_same< LinearOperator, MatrixOperator > :: value )
                {
                    // create new operator in case linear operator and matrix operator differ
                    opA.reset( CPRSelectorType::makeOperator( linearOperator.getmat(), parallelInformation_arg ) );
                }

                const double relax = parameters_.ilu_relaxation_;
                const MILU_VARIANT ilu_milu  = parameters_.ilu_milu_;
                if (  parameters_.use_cpr_ )
                {
                    using Matrix         = typename MatrixOperator::matrix_type;
                    using CouplingMetric = Opm::Amg::Element<pressureEqnIndex, pressureVarIndex>;
                    using CritBase       = Dune::Amg::SymmetricCriterion<Matrix, CouplingMetric>;
                    using Criterion      = Dune::Amg::CoarsenCriterion<CritBase>;
                    using AMG = typename ISTLUtility
                        ::BlackoilAmgSelector< Matrix, Vector, Vector,POrComm, Criterion, pressureEqnIndex, pressureVarIndex >::AMG;

                    std::unique_ptr< AMG > amg;
                    // Construct preconditioner.
                    Criterion crit(15, 2000);
                    constructAMGPrecond<Criterion>( linearOperator, parallelInformation_arg, amg, opA, relax, ilu_milu );

                    // Solve.
                    solve(linearOperator, x, istlb, *sp, *amg, result);
                }
                else
                {
                    typedef typename CPRSelectorType::AMG AMG;
                    std::unique_ptr< AMG > amg;

                    // Construct preconditioner.
                    constructAMGPrecond( linearOperator, parallelInformation_arg, amg, opA, relax, ilu_milu );

                    // Solve.
                    solve(linearOperator, x, istlb, *sp, *amg, result);
                }
            }
            else
#endif
            {
                // Construct preconditioner.
                auto precond = constructPrecond(linearOperator, parallelInformation_arg);

                // Solve.
                solve(linearOperator, x, istlb, *sp, *precond, result);
            }
        }


        // 3x3 matrix block inversion was unstable at least 2.3 until and including
        // 2.5.0. There may still be some issue with the 4x4 matrix block inversion
        // we therefore still use the block inversion in OPM
        typedef ParallelOverlappingILU0<Dune::BCRSMatrix<Dune::MatrixBlock<typename Matrix::field_type,
                                                                            Matrix::block_type::rows,
                                                                            Matrix::block_type::cols> >,
                                                                            Vector, Vector> SeqPreconditioner;


        template <class Operator>
        std::unique_ptr<SeqPreconditioner> constructPrecond(Operator& opA, const Dune::Amg::SequentialInformation&) const
        {
            const double relax   = parameters_.ilu_relaxation_;
            const int ilu_fillin = parameters_.ilu_fillin_level_;
            const MILU_VARIANT ilu_milu  = parameters_.ilu_milu_;
            const bool ilu_redblack = parameters_.ilu_redblack_;
            const bool ilu_reorder_spheres = parameters_.ilu_reorder_sphere_;
            std::unique_ptr<SeqPreconditioner> precond(new SeqPreconditioner(opA.getmat(), ilu_fillin, relax, ilu_milu, ilu_redblack, ilu_reorder_spheres));
            return precond;
        }

#if HAVE_MPI
        typedef Dune::OwnerOverlapCopyCommunication<int, int> Comm;
#if DUNE_VERSION_NEWER_REV(DUNE_ISTL, 2 , 5, 1)
        // 3x3 matrix block inversion was unstable from at least 2.3 until and
        // including 2.5.0
        typedef ParallelOverlappingILU0<Matrix,Vector,Vector,Comm> ParPreconditioner;
#else
        typedef ParallelOverlappingILU0<Dune::BCRSMatrix<Dune::MatrixBlock<typename Matrix::field_type,
                                                                           Matrix::block_type::rows,
                                                                           Matrix::block_type::cols> >,
                                        Vector, Vector, Comm> ParPreconditioner;
#endif
        template <class Operator>
        std::unique_ptr<ParPreconditioner>
        constructPrecond(Operator& opA, const Comm& comm) const
        {
            typedef std::unique_ptr<ParPreconditioner> Pointer;
            const double relax  = parameters_.ilu_relaxation_;
            const MILU_VARIANT ilu_milu  = parameters_.ilu_milu_;
            const bool ilu_redblack = parameters_.ilu_redblack_;
            const bool ilu_reorder_spheres = parameters_.ilu_reorder_sphere_;
            return Pointer(new ParPreconditioner(opA.getmat(), comm, relax, ilu_milu, ilu_redblack, ilu_reorder_spheres));
        }
#endif

        template <class LinearOperator, class MatrixOperator, class POrComm, class AMG >
        void
        constructAMGPrecond(LinearOperator& /* linearOperator */, const POrComm& comm, std::unique_ptr< AMG >& amg, std::unique_ptr< MatrixOperator >& opA, const double relax, const MILU_VARIANT milu) const
        {
            ISTLUtility::template createAMGPreconditionerPointer<pressureEqnIndex, pressureVarIndex>( *opA, relax, milu, comm, amg );
        }


        template <class C, class LinearOperator, class MatrixOperator, class POrComm, class AMG >
        void
        constructAMGPrecond(LinearOperator& /* linearOperator */, const POrComm& comm, std::unique_ptr< AMG >& amg, std::unique_ptr< MatrixOperator >& opA, const double relax,
                            const MILU_VARIANT /* milu */ ) const
        {
            ISTLUtility::template createAMGPreconditionerPointer<C>( *opA, relax,
                                                                     comm, amg, parameters_, weights_ );
        }


        /// \brief Solve the system using the given preconditioner and scalar product.
        template <class Operator, class ScalarProd, class Precond>
        void solve(Operator& opA, Vector& x, Vector& istlb, ScalarProd& sp, Precond& precond, Dune::InverseOperatorResult& result) const
        {
            // TODO: Revise when linear solvers interface opm-core is done
            // Construct linear solver.
            // GMRes solver
            int verbosity = 0;
            if (simulator_.gridView().comm().rank() == 0)
                verbosity = parameters_.linear_solver_verbosity_;

            if ( parameters_.newton_use_gmres_ ) {
                Dune::RestartedGMResSolver<Vector> linsolve(opA, sp, precond,
                          parameters_.linear_solver_reduction_,
                          parameters_.linear_solver_restart_,
                          parameters_.linear_solver_maxiter_,
                          verbosity);
                // Solve system.
                linsolve.apply(x, istlb, result);
            }
            else { // BiCGstab solver
                Dune::BiCGSTABSolver<Vector> linsolve(opA, sp, precond,
                          parameters_.linear_solver_reduction_,
                          parameters_.linear_solver_maxiter_,
                          verbosity);
                // Solve system.
                linsolve.apply(x, istlb, result);
            }
        }


        /// Solve the linear system Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] A   matrix A
        /// \param[inout] x  solution to be computed x
        /// \param[in] b   right hand side b
        void solve(Matrix& A, Vector& x, Vector& b ) const
        {
            // Parallel version is deactivated until we figure out how to do it properly.
#if HAVE_MPI
            if (parallelInformation_.type() == typeid(ParallelISTLInformation))
            {
                typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
                const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>( parallelInformation_);
                Comm istlComm(info.communicator());

                // Construct operator, scalar product and vectors needed.
                typedef Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector,Comm> Operator;
                Operator opA(A, istlComm);
                solve( opA, x, b, istlComm  );
            }
            else
#endif
            {
                // Construct operator, scalar product and vectors needed.
                Dune::MatrixAdapter< Matrix, Vector, Vector> opA( A );
                solve( opA, x, b );
            }
        }

        /// Solve the linear system Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] A   matrix A
        /// \param[inout] x  solution to be computed x
        /// \param[in] b   right hand side b
        template <class Operator, class Comm >
        void solve(Operator& opA, Vector& x, Vector& b, Comm& comm) const
        {
            Dune::InverseOperatorResult result;
            // Parallel version is deactivated until we figure out how to do it properly.
#if HAVE_MPI
            if (parallelInformation_.type() == typeid(ParallelISTLInformation))
            {
                const size_t size = opA.getmat().N();
                const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>( parallelInformation_);

                // As we use a dune-istl with block size np the number of components
                // per parallel is only one.
                info.copyValuesTo(comm.indexSet(), comm.remoteIndices(),
                                  size, 1);
                // Construct operator, scalar product and vectors needed.
                constructPreconditionerAndSolve<Dune::SolverCategory::overlapping>(opA, x, b, comm, result);
            }
            else
#endif
            {
                OPM_THROW(std::logic_error,"this method if for parallel solve only");
            }

            checkConvergence( result );
        }

        /// Solve the linear system Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] A   matrix A
        /// \param[inout] x  solution to be computed x
        /// \param[in] b   right hand side b
        template <class Operator>
        void solve(Operator& opA, Vector& x, Vector& b ) const
        {
            Dune::InverseOperatorResult result;
            // Construct operator, scalar product and vectors needed.
            Dune::Amg::SequentialInformation info;
            constructPreconditionerAndSolve(opA, x, b, info, result);
            checkConvergence( result );
        }

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
            return parallelInformation_.type() == typeid(ParallelISTLInformation);
#else
            return false;
#endif
        }

        /// Zero out off-diagonal blocks on rows corresponding to overlap cells
        /// Diagonal blocks on ovelap rows are set to diag(1e100).
        void makeOverlapRowsInvalid(Matrix& ebosJacIgnoreOverlap) const
        {
            //value to set on diagonal
            typedef Dune::FieldMatrix<Scalar, numEq, numEq >        MatrixBlockType;
            MatrixBlockType diag_block(0.0);
            for (int eq = 0; eq < numEq; ++eq)
                diag_block[eq][eq] = 1.0e100;

            //loop over precalculated overlap rows and columns
            for (auto row = overlapRowAndColumns_.begin(); row != overlapRowAndColumns_.end(); row++ )
            {
                int lcell = row->first;
                //diagonal block set to large value diagonal
                ebosJacIgnoreOverlap[lcell][lcell] = diag_block;

                //loop over off diagonal blocks in overlap row
                for (auto col = row->second.begin(); col != row->second.end(); ++col)
                {
                    int ncell = *col;
                    //zero out block
                    ebosJacIgnoreOverlap[lcell][ncell] = 0.0;
                }
            }
        }

        // Weights to make approximate pressure equations.
        // Calculated from the storage terms (only) of the
        // conservation equations, ignoring all other terms.
        Vector getStorageWeights() const
        {
            Vector weights(rhs_->size());
            BlockVector rhs(0.0);
            rhs[pressureVarIndex] = 1.0;
            int index = 0;
            ElementContext elemCtx(simulator_);
            const auto& vanguard = simulator_.vanguard();
            auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
            const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const Element& elem = *elemIt;
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                Dune::FieldVector<Evaluation, numEq> storage;
                unsigned threadId = ThreadManager::threadId();
                simulator_.model().localLinearizer(threadId).localResidual().computeStorage(storage,elemCtx,/*spaceIdx=*/0, /*timeIdx=*/0);
                Scalar extrusionFactor = elemCtx.intensiveQuantities(0, /*timeIdx=*/0).extrusionFactor();
                Scalar scvVolume = elemCtx.stencil(/*timeIdx=*/0).subControlVolume(0).volume() * extrusionFactor;
                Scalar storage_scale = scvVolume / elemCtx.simulator().timeStepSize();
                MatrixBlockType block;
                double pressure_scale = 50e5;
                for (int ii = 0; ii < numEq; ++ii) {
                    for (int jj = 0; jj < numEq; ++jj) {
                        block[ii][jj] = storage[ii].derivative(jj)/storage_scale;
                        if (jj == pressureVarIndex) {
                            block[ii][jj] *= pressure_scale;
                        }
                    }
                }
                BlockVector bweights;
                MatrixBlockType block_transpose = Opm::transposeDenseMatrix(block);
                block_transpose.solve(bweights, rhs);
                bweights /= 1000.0; // given normal densities this scales weights to about 1.
                weights[index] = bweights;
                ++index;
            }
            return weights;
        }

        // Interaction between the CPR weights (the function argument 'weights')
        // and the variable and equation weights from
        // simulator_.model().primaryVarWeight() and
        // simulator_.model().eqWeight() is nontrivial and does not work
        // at the moment. Possibly refactoring of ewoms weight treatment
        // is needed. In the meantime this function shows what needs to be
        // done to integrate the weights properly.
        void scaleEquationsAndVariables(Vector& weights)
        {
            // loop over primary variables
            const auto endi = matrix_->end();
            for (auto i = matrix_->begin(); i != endi; ++i) {
                const auto endj = (*i).end();
                BlockVector& brhs = (*rhs_)[i.index()];
                for (auto j = (*i).begin(); j != endj; ++j) {
                    MatrixBlockType& block = *j;
                    for (std::size_t ii = 0; ii < block.rows; ii++ ) {
                        for (std::size_t jj = 0; jj < block.cols; jj++) {
                            double var_scale = simulator_.model().primaryVarWeight(i.index(),jj);
                            block[ii][jj] /= var_scale;
                            block[ii][jj] *= simulator_.model().eqWeight(i.index(), ii);
                        }
                    }
                }
                for (std::size_t ii = 0; ii < brhs.size(); ii++) {
                    brhs[ii] *= simulator_.model().eqWeight(i.index(), ii);
                }
                if (weights.size() == matrix_->N()) {
                    BlockVector& bw = weights[i.index()];
                    for (std::size_t ii = 0; ii < brhs.size(); ii++) {
                        bw[ii] /= simulator_.model().eqWeight(i.index(), ii);
                    }
                    double abs_max =
                        *std::max_element(bw.begin(), bw.end(), [](double a, double b){ return std::abs(a) < std::abs(b); } );
                    bw /= abs_max;
                }
            }
        }

        void scaleSolution(Vector& x)
        {
            for (std::size_t i = 0; i < x.size(); ++i) {
                auto& bx = x[i];
                for (std::size_t jj = 0; jj < bx.size(); jj++) {
                    double var_scale = simulator_.model().primaryVarWeight(i,jj);
                    bx[jj] /= var_scale;
                }
            }
        }

        Vector getQuasiImpesWeights()
        {
            Matrix& A = *matrix_;
            Vector weights(rhs_->size());
            BlockVector rhs(0.0);
            rhs[pressureVarIndex] = 1;
            const auto endi = A.end();
            for (auto i = A.begin(); i!=endi; ++i) {
                const auto endj = (*i).end();
                MatrixBlockType diag_block(0.0);
                for (auto j=(*i).begin(); j!=endj; ++j) {
                    if (i.index() == j.index()) {
                        diag_block = (*j);
                        break;
                    }
                }
                BlockVector bweights;
                auto diag_block_transpose = Opm::transposeDenseMatrix(diag_block);
                diag_block_transpose.solve(bweights, rhs);
                double abs_max =
                    *std::max_element(bweights.begin(), bweights.end(), [](double a, double b){ return std::abs(a) < std::abs(b); } );
                bweights /= std::abs(abs_max);
                weights[i.index()] = bweights;
            }
            return weights;
        }

        Vector getSimpleWeights(const BlockVector& rhs)
        {
            Vector weights(rhs_->size(), 0);
            for (auto& bw : weights) {
                bw = rhs;
            }
            return weights;
        }

        void scaleMatrixAndRhs(const Vector& weights)
        {
            using Block = typename Matrix::block_type;
            const auto endi = matrix_->end();
            for (auto i = matrix_->begin(); i !=endi; ++i) {
                const BlockVector& bweights = weights[i.index()];
                BlockVector& brhs = (*rhs_)[i.index()];
                const auto endj = (*i).end();
                for (auto j = (*i).begin(); j != endj; ++j) {
                    // assume it is something on all rows
                    Block& block = (*j);
                    BlockVector neweq(0.0);
                    for (std::size_t ii = 0; ii < block.rows; ii++) {
                        for (std::size_t jj = 0; jj < block.cols; jj++) {
                            neweq[jj] += bweights[ii]*block[ii][jj];
                        }
                    }
                    block[pressureEqnIndex] = neweq;
                }
                Scalar newrhs(0.0);
                for (std::size_t ii = 0; ii < brhs.size(); ii++) {
                    newrhs += bweights[ii]*brhs[ii];
                }
                brhs[pressureEqnIndex] = newrhs;
            }
        }

        static void multBlocksInMatrix(Matrix& ebosJac, const MatrixBlockType& trans, const bool left = true)
        {
            const int n = ebosJac.N();
            for (int row_index = 0; row_index < n; ++row_index) {
                auto& row = ebosJac[row_index];
                auto* dataptr = row.getptr();
                for (int elem = 0; elem < row.N(); ++elem) {
                    auto& block = dataptr[elem];
                    if (left) {
                        block = block.leftmultiply(trans);
                    } else {
                        block = block.rightmultiply(trans);
                    }
                }
            }
        }

        static void multBlocksVector(Vector& ebosResid_cp, const MatrixBlockType& leftTrans)
        {
            for (auto& bvec : ebosResid_cp) {
                auto bvec_new = bvec;
                leftTrans.mv(bvec, bvec_new);
                bvec = bvec_new;
            }
        }

        static void scaleCPRSystem(Matrix& M_cp, Vector& b_cp, const MatrixBlockType& leftTrans)
        {
            multBlocksInMatrix(M_cp, leftTrans, true);
            multBlocksVector(b_cp, leftTrans);
        }

        const Simulator& simulator_;
        mutable int iterations_;
        mutable bool converged_;
        boost::any parallelInformation_;

        std::unique_ptr<Matrix> matrix_;
        Vector *rhs_;
        std::unique_ptr<Matrix> matrix_for_preconditioner_;

        std::vector<std::pair<int,std::vector<int>>> overlapRowAndColumns_;
        FlowLinearSolverParameters parameters_;
        Vector weights_;
        bool scale_variables_;
    }; // end ISTLSolver

} // namespace Opm
#endif
