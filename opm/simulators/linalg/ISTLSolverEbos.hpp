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
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
#include <opm/simulators/linalg/ParallelRestrictedAdditiveSchwarz.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/ExtractParallelGridInformationToISTL.hpp>
#include <opm/simulators/linalg/findOverlapRowsAndColumns.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/WriteSystemMatrixHelper.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>
#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>

#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#if HAVE_CUDA
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

  Dune::SolverCategory::Category category() const override
  {
    return overlapping ?
           Dune::SolverCategory::overlapping : Dune::SolverCategory::sequential;
  }

  //! constructor: just store a reference to a matrix
  WellModelMatrixAdapter (const M& A,
                          const WellModel& wellMod,
                          const std::shared_ptr< communication_type >& comm = std::shared_ptr< communication_type >())
      : A_( A ), wellMod_( wellMod ), comm_(comm)
  {}


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

  virtual const matrix_type& getmat() const override { return A_; }

protected:
  const matrix_type& A_ ;
  const WellModel& wellMod_;
  std::shared_ptr< communication_type > comm_;
};


/*!
   \brief Adapter to turn a matrix into a linear operator.
   Adapts a matrix to the assembled linear operator interface.
   We assume parallel ordering, where ghost rows are located after interior rows
 */
template<class M, class X, class Y, class WellModel, bool overlapping >
class WellModelGhostLastMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
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


    Dune::SolverCategory::Category category() const override
    {
        return overlapping ?
            Dune::SolverCategory::overlapping : Dune::SolverCategory::sequential;
    }

    //! constructor: just store a reference to a matrix
    WellModelGhostLastMatrixAdapter (const M& A,
                                     const WellModel& wellMod,
                                     const size_t interiorSize )
        : A_( A ), wellMod_( wellMod ), interiorSize_(interiorSize)
    {}

    virtual void apply( const X& x, Y& y ) const override
    {
        for (auto row = A_.begin(); row.index() < interiorSize_; ++row)
        {
            y[row.index()]=0;
            auto endc = (*row).end();
            for (auto col = (*row).begin(); col != endc; ++col)
                (*col).umv(x[col.index()], y[row.index()]);
        }

        // add well model modification to y
        wellMod_.apply(x, y );

        ghostLastProject( y );
    }

    // y += \alpha * A * x
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const override
    {
        for (auto row = A_.begin(); row.index() < interiorSize_; ++row)
        {
            auto endc = (*row).end();
            for (auto col = (*row).begin(); col != endc; ++col)
                (*col).usmv(alpha, x[col.index()], y[row.index()]);
        }
        // add scaled well model modification to y
        wellMod_.applyScaleAdd( alpha, x, y );

        ghostLastProject( y );
    }

    virtual const matrix_type& getmat() const override { return A_; }

protected:
    void ghostLastProject(Y& y) const
    {
        size_t end = y.size();
        for (size_t i = interiorSize_; i < end; ++i)
            y[i] = 0;
    }

    const matrix_type& A_ ;
    const WellModel& wellMod_;
    size_t interiorSize_;
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
        using FlexibleSolverType = Dune::FlexibleSolver<Matrix, Vector>;
        // Due to miscibility oil <-> gas the water eqn is the one we can replace with a pressure equation.
        static const bool waterEnabled = Indices::waterEnabled;
        static const int pindex = (waterEnabled) ? BlackOilDefaultIndexTraits::waterCompIdx : BlackOilDefaultIndexTraits::oilCompIdx;
        enum { pressureEqnIndex = pindex };
        enum { pressureVarIndex = Indices::pressureSwitchIdx };
        static const int numEq = Indices::numEq;

#if HAVE_CUDA
        std::unique_ptr<BdaBridge> bdaBridge;
#endif

#if HAVE_MPI
        typedef Dune::OwnerOverlapCopyCommunication<int,int> communication_type;
#else
        typedef Dune::CollectiveCommunication< int > communication_type;
#endif

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
              converged_(false),
              matrix_()
        {
#if HAVE_MPI
            comm_.reset( new communication_type( simulator_.vanguard().grid().comm() ) );
#endif
            parameters_.template init<TypeTag>();
            useFlexible_ = parameters_.use_cpr_ || EWOMS_PARAM_IS_SET(TypeTag, std::string, LinearSolverConfiguration);

            if (useFlexible_)
            {
                prm_ = setupPropertyTree<TypeTag>(parameters_);
            }
            const auto& gridForConn = simulator_.vanguard().grid();
#if HAVE_CUDA
            bool use_gpu = EWOMS_GET_PARAM(TypeTag, bool, UseGpu);
            if (gridForConn.comm().size() > 1 && use_gpu) {
                OpmLog::warning("Warning cannot use GPU with MPI, GPU is disabled");
                use_gpu = false;
            }
            const int maxit = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
            const double tolerance = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
            const int linear_solver_verbosity = parameters_.linear_solver_verbosity_;
            bdaBridge.reset(new BdaBridge(use_gpu, linear_solver_verbosity, maxit, tolerance));
#else
            const bool use_gpu = EWOMS_GET_PARAM(TypeTag, bool, UseGpu);
            if (use_gpu) {
                OPM_THROW(std::logic_error,"Error cannot use GPU solver since CUDA was not found during compilation");
            }
#endif
            extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);
            useWellConn_ = EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions);

            if (!useWellConn_ && useFlexible_)
            {
                OPM_THROW(std::logic_error, "Flexible solvers and CPR need the well contribution in the matrix. Please run with"
                          " --matrix-add-well-contributions=true");
            }

            ownersFirst_ = EWOMS_GET_PARAM(TypeTag, bool, OwnerCellsFirst);
            interiorCellNum_ = detail::numMatrixRowsToUseInSolver(simulator_.vanguard().grid(), ownersFirst_);

            if ( isParallel() && (!ownersFirst_ || parameters_.linear_solver_use_amg_  || useFlexible_ ) ) {
                detail::setWellConnections(gridForConn, simulator_.vanguard().schedule().getWellsatEnd(), useWellConn_, wellConnectionsGraph_);
                // For some reason simulator_.model().elementMapper() is not initialized at this stage
                // Hence const auto& elemMapper = simulator_.model().elementMapper(); does not work.
                // Set it up manually
                using ElementMapper =
                    Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
                ElementMapper elemMapper(simulator_.vanguard().gridView(), Dune::mcmgElementLayout());
                detail::findOverlapAndInterior(gridForConn, elemMapper, overlapRows_, interiorRows_);
                noGhostAdjacency();
                setGhostsInNoGhost(*noGhostMat_);
                if (ownersFirst_)
                    OpmLog::warning("OwnerCellsFirst option is true, but ignored.");
            }

            if (useFlexible_)
            {
                // Print parameters to PRT/DBG logs.
                if (simulator.gridView().comm().rank() == 0) {
                    std::ostringstream os;
                    os << "Property tree for linear solver:\n";
                    boost::property_tree::write_json(os, prm_, true);
                    OpmLog::note(os.str());
                }
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
            if (noGhostMat_)
            {
                copyJacToNoGhost(M.istlMatrix(), *noGhostMat_);
            }
            else
            {
                if (firstcall)
                {
                    // ebos will not change the matrix object. Hence simply store a pointer
                    // to the original one with a deleter that does nothing.
                    // Outch! We need to be able to scale the linear system! Hence const_cast
                    matrix_ = const_cast<Matrix*>(&M.istlMatrix());
                }
                else
                {
                    // Pointers should not change
                    if ( &(M.istlMatrix()) != matrix_ )
                    {
                        OPM_THROW(std::logic_error, "Matrix objects are expected to be reused when reassembling!"
                                  <<" old pointer was " << matrix_ << ", new one is " << (&M.istlMatrix()) );
                    }
                }
            }
            rhs_ = &b;

            if (useFlexible_)
            {
                prepareFlexibleSolver();
            }
            else
            {
                this->scaleSystem();
            }
            firstcall = false;
        }

        void scaleSystem()
        {
            if (useWellConn_) {
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
            const int verbosity = useFlexible_ ? prm_.get<int>("verbosity", 0) : parameters_.linear_solver_verbosity_;
            const bool write_matrix = verbosity > 10;
            // Solve system.

            if (useFlexible_)
            {
                Dune::InverseOperatorResult res;
                assert(flexibleSolver_);
                flexibleSolver_->apply(x, *rhs_, res);
                iterations_ = res.iterations;
                if (write_matrix) {
                    Opm::Helper::writeSystem(simulator_, //simulator is only used to get names
                                             getMatrix(),
                                             *rhs_,
                                             comm_.get());
                }

                return converged_ = res.converged;
            }

            const WellModel& wellModel = simulator_.problem().wellModel();

            if( isParallel() )
            {
                if ( ownersFirst_ && !parameters_.linear_solver_use_amg_ && !useFlexible_) {
                    typedef WellModelGhostLastMatrixAdapter< Matrix, Vector, Vector, WellModel, true > Operator;
                    assert(matrix_);
                    Operator opA(*matrix_, wellModel, interiorCellNum_);

                    solve( opA, x, *rhs_, *comm_ );
                }
                else {

                    typedef WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, true > Operator;
                    assert (noGhostMat_);
                    Operator opA(*noGhostMat_, wellModel,
                                 comm_ );

                    solve( opA, x, *rhs_, *comm_ );

                }
            }
            else
            {
                typedef WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, false > Operator;
                Operator opA(*matrix_, wellModel);
                solve( opA, x, *rhs_ );
            }

            if (parameters_.scale_linear_system_) {
                scaleSolution(x);
            }

            if (write_matrix) {
                Opm::Helper::writeSystem(simulator_, //simulator is only used to get names
                                         getMatrix(),
                                         *rhs_,
                                         comm_.get());
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
        const std::any& parallelInformation() const { return parallelInformation_; }

    protected:
        /// \brief construct the CPR preconditioner and the solver.
        /// \tparam P The type of the parallel information.
        /// \param parallelInformation the information about the parallelization.
        template<Dune::SolverCategory::Category category=Dune::SolverCategory::sequential,
                 class LinearOperator, class POrComm>
        void constructPreconditionerAndSolve(LinearOperator& linearOperator,
                                             Vector& x, Vector& istlb,
                                             const POrComm& parallelInformation_arg,
                                             Dune::InverseOperatorResult& result) const
        {
            // Construct scalar product.
            auto sp = Dune::createScalarProduct<Vector,POrComm>(parallelInformation_arg, category);

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
                    using MatrixType     = typename MatrixOperator::matrix_type;
                    using CouplingMetric = Opm::Amg::Element<pressureEqnIndex, pressureVarIndex>;
                    using CritBase       = Dune::Amg::SymmetricCriterion<MatrixType, CouplingMetric>;
                    using Criterion      = Dune::Amg::CoarsenCriterion<CritBase>;
                    using AMG = typename ISTLUtility
                        ::BlackoilAmgSelector< MatrixType, Vector, Vector,POrComm, Criterion, pressureEqnIndex, pressureVarIndex >::AMG;

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
                // tries to solve linear system
#if HAVE_CUDA
                bool use_gpu = bdaBridge->getUseGpu();
                if (use_gpu) {
                    WellContributions wellContribs;
                    if (!useWellConn_) {
                        simulator_.problem().wellModel().getWellContributions(wellContribs);
                    }
                    // Const_cast needed since the CUDA stuff overwrites values for better matrix condition..
                    bdaBridge->solve_system(const_cast<Matrix*>(&getMatrix()), istlb, wellContribs, result);
                    if (result.converged) {
                        // get result vector x from non-Dune backend, iff solve was successful
                        bdaBridge->get_result(x);
                    } else {
                        // CPU fallback
                        use_gpu = bdaBridge->getUseGpu();  // update value, BdaBridge might have disabled cusparseSolver
                        if (use_gpu) {
                            OpmLog::warning("cusparseSolver did not converge, now trying Dune to solve current linear system...");
                        }

                        // call Dune
                        auto precond = constructPrecond(linearOperator, parallelInformation_arg);
                        solve(linearOperator, x, istlb, *sp, *precond, result);
                    }
                } else { // gpu is not selected or disabled
                    auto precond = constructPrecond(linearOperator, parallelInformation_arg);
                    solve(linearOperator, x, istlb, *sp, *precond, result);
                }
#else
                // Construct preconditioner.
                auto precond = constructPrecond(linearOperator, parallelInformation_arg);

                // Solve.
                solve(linearOperator, x, istlb, *sp, *precond, result);
#endif
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
        // 3x3 matrix block inversion was unstable from at least 2.3 until and
        // including 2.5.0
        typedef ParallelOverlappingILU0<Matrix,Vector,Vector,Comm> ParPreconditioner;
        template <class Operator>
        std::unique_ptr<ParPreconditioner>
        constructPrecond(Operator& opA, const Comm& comm) const
        {
            typedef std::unique_ptr<ParPreconditioner> Pointer;
            const double relax  = parameters_.ilu_relaxation_;
            const MILU_VARIANT ilu_milu  = parameters_.ilu_milu_;
            const bool ilu_redblack = parameters_.ilu_redblack_;
            const bool ilu_reorder_spheres = parameters_.ilu_reorder_sphere_;
            return Pointer(new ParPreconditioner(opA.getmat(), comm, relax, ilu_milu, interiorCellNum_, ilu_redblack, ilu_reorder_spheres));
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
                // Construct operator, scalar product and vectors needed.
                typedef Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector,Comm> Operator;
                Operator opA(A, *comm_);
                solve( opA, x, b, *comm_  );
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
        void solve(Operator& opA OPM_UNUSED_NOMPI,
                   Vector& x OPM_UNUSED_NOMPI,
                   Vector& b OPM_UNUSED_NOMPI,
                   Comm& comm OPM_UNUSED_NOMPI) const
        {
            Dune::InverseOperatorResult result;
            // Parallel version is deactivated until we figure out how to do it properly.
#if HAVE_MPI
            if (parallelInformation_.type() == typeid(ParallelISTLInformation))
            {
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
            return comm_->communicator().size() > 1;
#else
            return false;
#endif
        }

        void prepareFlexibleSolver()
        {
            // Decide if we should recreate the solver or just do
            // a minimal preconditioner update.
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

            std::function<Vector()> weightsCalculator;

            auto preconditionerType = prm_.get("preconditioner.type", "cpr");
            if( preconditionerType  == "cpr" ||
                preconditionerType == "cprt"
                )
            {
                bool transpose = false;
                if(preconditionerType == "cprt"){
                    transpose = true;
                }

                auto weightsType = prm_.get("preconditioner.weight_type", "quasiimpes");
                auto pressureIndex = this->prm_.get("preconditioner.pressure_var_index", 1);
                if(weightsType == "quasiimpes") {
                    // weighs will be created as default in the solver
                    weightsCalculator =
                        [this, transpose, pressureIndex](){
                            return Opm::Amg::getQuasiImpesWeights<Matrix,
                                                                  Vector>(getMatrix(),
                                                                          pressureIndex,
                                                                          transpose);
                        };

                }else if(weightsType == "trueimpes"  ){
                    weightsCalculator =
                        [this](){
                            return this->getStorageWeights();
                        };
                }else{
                    OPM_THROW(std::invalid_argument, "Weights type " << weightsType << "not implemented for cpr."
                              << " Please use quasiimpes or trueimpes.");
                }
            }

            if (recreate_solver || !flexibleSolver_) {
                if (isParallel()) {
#if HAVE_MPI
                    assert(noGhostMat_);
                    flexibleSolver_.reset(new FlexibleSolverType(*noGhostMat_, *comm_, prm_, weightsCalculator));
#endif
                } else {
                    flexibleSolver_.reset(new FlexibleSolverType(*matrix_, prm_, weightsCalculator));
                }
            }
            else
            {
                flexibleSolver_->preconditioner().update();
            }
        }

        /// Create sparsity pattern of matrix without off-diagonal ghost entries.
        void noGhostAdjacency()
        {
            const auto& grid = simulator_.vanguard().grid();
            const auto& gridView = simulator_.vanguard().gridView();
            // For some reason simulator_.model().elementMapper() is not initialized at this stage.
            // Hence const auto& elemMapper = simulator_.model().elementMapper(); does not work.
            // Set it up manually
            using ElementMapper =
                Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
            ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());
            typedef typename Matrix::size_type size_type;
            size_type numCells = grid.size( 0 );
            noGhostMat_.reset(new Matrix(numCells, numCells, Matrix::random));

            std::vector<std::set<size_type>> pattern;
            pattern.resize(numCells);

            auto elemIt = gridView.template begin<0>();
            const auto& elemEndIt = gridView.template end<0>();

            //Loop over cells
            for (; elemIt != elemEndIt; ++elemIt)
            {
                const auto& elem = *elemIt;
                size_type idx = elemMapper.index(elem);
                pattern[idx].insert(idx);

                // Add well non-zero connections
                for (auto wc = wellConnectionsGraph_[idx].begin(); wc!=wellConnectionsGraph_[idx].end(); ++wc)
                    pattern[idx].insert(*wc);

                // Add just a single element to ghost rows
                if (elem.partitionType() != Dune::InteriorEntity)
                {
                    noGhostMat_->setrowsize(idx, pattern[idx].size());
                }
                else {
                    auto isend = gridView.iend(elem);
                    for (auto is = gridView.ibegin(elem); is!=isend; ++is)
                    {
                        //check if face has neighbor
                        if (is->neighbor())
                        {
                            size_type nid = elemMapper.index(is->outside());
                            pattern[idx].insert(nid);
                        }
                    }
                    noGhostMat_->setrowsize(idx, pattern[idx].size());
                }
            }
            noGhostMat_->endrowsizes();
            for (size_type dofId = 0; dofId < numCells; ++dofId)
            {
                auto nabIdx = pattern[dofId].begin();
                auto endNab = pattern[dofId].end();
                for (; nabIdx != endNab; ++nabIdx)
                {
                    noGhostMat_->addindex(dofId, *nabIdx);
                }
            }
            noGhostMat_->endindices();
        }

        /// Set the ghost diagonal in noGhost to diag(1.0)
        void setGhostsInNoGhost(Matrix& ng)
        {
            ng=0;
            typedef typename Matrix::block_type MatrixBlockTypeT;
            MatrixBlockTypeT diag_block(0.0);
            for (int eq = 0; eq < Matrix::block_type::rows; ++eq)
                diag_block[eq][eq] = 1.0;

            //loop over precalculated ghost rows and columns
            for (auto row = overlapRows_.begin(); row != overlapRows_.end(); row++ )
            {
                int lcell = *row;
                //diagonal block set to 1
                ng[lcell][lcell] = diag_block;
            }
        }

        /// Copy interior rows to noghost matrix
        void copyJacToNoGhost(const Matrix& jac, Matrix& ng)
        {
            //Loop over precalculated interior rows.
            for (auto row = interiorRows_.begin(); row != interiorRows_.end(); row++ )
            {
                //Copy row
                ng[*row] = jac[*row];
            }
        }

        // Weights to make approximate pressure equations.
        // Calculated from the storage terms (only) of the
        // conservation equations, ignoring all other terms.
        Vector getStorageWeights() const
        {
            Vector weights(rhs_->size());
            ElementContext elemCtx(simulator_);
            Opm::Amg::getTrueImpesWeights(pressureVarIndex, weights, simulator_.vanguard().gridView(),
                                          elemCtx, simulator_.model(),
                                          ThreadManager::threadId());
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
            const auto endi = getMatrix().end();
            for (auto i = getMatrix().begin(); i != endi; ++i) {
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
                if (weights.size() == getMatrix().N()) {
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
            return Amg::getQuasiImpesWeights<Matrix,Vector>(getMatrix(), pressureVarIndex, /* transpose=*/ true);
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
            const auto endi = getMatrix().end();
            for (auto i = getMatrix().begin(); i !=endi; ++i) {
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

        Matrix& getMatrix()
        {
            return noGhostMat_ ? *noGhostMat_ : *matrix_;
        }

        const Matrix& getMatrix() const
        {
            return noGhostMat_ ? *noGhostMat_ : *matrix_;
        }

        const Simulator& simulator_;
        mutable int iterations_;
        mutable bool converged_;
        std::any parallelInformation_;

        // non-const to be able to scale the linear system
        Matrix* matrix_;
        std::unique_ptr<Matrix> noGhostMat_;
        Vector *rhs_;

        std::unique_ptr<FlexibleSolverType> flexibleSolver_;
        std::vector<int> overlapRows_;
        std::vector<int> interiorRows_;
        std::vector<std::set<int>> wellConnectionsGraph_;

        bool ownersFirst_;
        bool useWellConn_;
        bool useFlexible_;
        size_t interiorCellNum_;

        FlowLinearSolverParameters parameters_;
        boost::property_tree::ptree prm_;
        Vector weights_;
        bool scale_variables_;

        std::shared_ptr< communication_type > comm_;
    }; // end ISTLSolver

} // namespace Opm
#endif
