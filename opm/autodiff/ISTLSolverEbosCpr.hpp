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
#include <utility>
#include <memory>
BEGIN_PROPERTIES

END_PROPERTIES

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
         typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
        typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
        typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
        typedef typename GET_PROP_TYPE(TypeTag, EclWellModel) WellModel;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
        typedef typename SparseMatrixAdapter::IstlMatrix Matrix;

      //typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlockType;
      //typedef typename Vector::block_type BlockVector;
      //typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
      //typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;
      //typedef typename GridView::template Codim<0>::Entity Element;
      //typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

      
        enum { pressureIndex = Indices::pressureSwitchIdx };
        static const int numEq = Indices::numEq;
      typedef WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, false, TypeTag > OperatorSerial;
      typedef ISTLSolverEbos<TypeTag> SuperClass;
      typedef Dune::Amg::SequentialInformation POrComm;
      typedef ISTLUtility::CPRSelector< Matrix, Vector, Vector, POrComm>  CPRSelectorType;
      typedef typename CPRSelectorType::Operator MatrixOperator;
      static constexpr int category = Dune::SolverCategory::sequential;
      typedef Dune::ScalarProductChooser<Vector, POrComm, category> ScalarProductChooser;
      


    public:
        typedef Dune::AssembledLinearOperator< Matrix, Vector, Vector > AssembledLinearOperatorType;

        static void registerParameters()
        {
            FlowLinearSolverParameters::registerParameters<TypeTag>();
        }

        /// Construct a system solver.
        /// \param[in] parallelInformation In the case of a parallel run
        ///                                with dune-istl the information about the parallelization.
        ISTLSolverEbosCpr(const Simulator& simulator)
	  : SuperClass(simulator)
         {
         }

        // nothing to clean here
        void eraseMatrix() {
            this->matrix_for_preconditioner_.reset();
        }

        void prepare(const Matrix& M, Vector& b) {
	  SuperClass::prepare(M,b);
	  const WellModel& wellModel = this->simulator_.problem().wellModel();
	  
	  
	  if( this->isParallel() )
	    {
#if HAVE_MPI			      
	      typedef WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, true ,TypeTag> Operator;
			
	      auto ebosJacIgnoreOverlap = Matrix(*(this->matrix_));
	      //remove ghost rows in local matrix
	      this->makeOverlapRowsInvalid(ebosJacIgnoreOverlap);
	      
	      //Not sure what actual_mat_for_prec is, so put ebosJacIgnoreOverlap as both variables
	      //to be certain that correct matrix is used for preconditioning.
	      Operator opA(ebosJacIgnoreOverlap, ebosJacIgnoreOverlap, wellModel,
			   this->parallelInformation_ );
	      assert( opA.comm() );
	      //SuperClass::solve( opA, x, *(this->rhs_), *(opA.comm()) );
	      typedef Dune::OwnerOverlapCopyCommunication<int,int>& comm = *(opA.comm());
	      const size_t size = opA.getmat().N();

	      const ParallelISTLInformation& info =
		boost::any_cast<const ParallelISTLInformation&>( this->parallelInformation_);

	      // As we use a dune-istl with block size np the number of components
	      // per parallel is only one.
	      info.copyValuesTo(comm.indexSet(), comm.remoteIndices(),
				size, 1);
	      // Construct operator, scalar product and vectors needed.
	      Dune::InverseOperatorResult result;
	      SuperClass::constructPreconditionerAndSolve<Dune::SolverCategory::overlapping>(opA, x, *(this->rhs_), comm, result);
	      SuperClass::checkConvergence(result);
#endif
            }
	  else
	    {
	      const WellModel& wellModel = this->simulator_.problem().wellModel();
	      //typedef WellModelMatrixAdapter< Matrix, Vector, Vector, WellModel, false ,TypeTag> OperatorSerial;
	      opASerial_.reset(new OperatorSerial(*(this->matrix_), *(this->matrix_), wellModel));
		
	      //Dune::Amg::SequentialInformation info;
	      typedef Dune::Amg::SequentialInformation POrComm;
	      POrComm parallelInformation_arg;
	      typedef  OperatorSerial LinearOperator;
	      
	      //SuperClass::constructPreconditionerAndSolve(opA, x, *(this->rhs_), info, result);
		
	      constexpr int  category = Dune::SolverCategory::sequential;
	      typedef Dune::ScalarProductChooser<Vector, POrComm, category> ScalarProductChooser;
	      typedef std::unique_ptr<typename ScalarProductChooser::ScalarProduct> SPPointer;
	      SPPointer sp(ScalarProductChooser::construct(parallelInformation_arg));
	      sp_ = std::move(sp);
	      Vector& istlb = *(this->rhs_);
	      parallelInformation_arg.copyOwnerToAll(istlb, istlb);
	      
	      
		
	      if( ! std::is_same< LinearOperator, MatrixOperator > :: value )
                {
		  // create new operator in case linear operator and matrix operator differ
		  opA_.reset( CPRSelectorType::makeOperator( opASerial_->getmat(), parallelInformation_arg ) );
                }

	      const double relax = this->parameters_.ilu_relaxation_;
	      //const MILU_VARIANT ilu_milu  = this->parameters_.ilu_milu_;
	      using Matrix         = typename MatrixOperator::matrix_type;
	      using CouplingMetric = Dune::Amg::Diagonal<pressureIndex>;
	      using CritBase       = Dune::Amg::SymmetricCriterion<Matrix, CouplingMetric>;
	      using Criterion      = Dune::Amg::CoarsenCriterion<CritBase>;
	      using AMG = typename ISTLUtility
		::BlackoilAmgSelector< Matrix, Vector, Vector,POrComm, Criterion, pressureIndex >::AMG;
		
	      //std::unique_ptr< AMG > amg;
	      // Construct preconditioner.
	      Criterion crit(15, 2000);
	      //SuperClass::constructAMGPrecond< Criterion >( linearOperator, parallelInformation_arg, amg, opA, relax, ilu_milu );
	      ISTLUtility::template createAMGPreconditionerPointer<Criterion>( *opA_,
									       relax,
									       parallelInformation_arg,
									       amg_,
									       this->parameters_,
									       this->weights_ );

	      
	      // Solve.
	      //SuperClass::solve(linearOperator, x, istlb, *sp, *amg, result);
	      int verbosity = ( this->isIORank_ ) ? this->parameters_.linear_solver_verbosity_ : 0;
	      linsolve_.reset(new Dune::BiCGSTABSolver<Vector>(*opASerial_, *sp_, *amg_,
						    this->parameters_.linear_solver_reduction_,
							       this->parameters_.linear_solver_maxiter_,
							       verbosity));

	    }	  
        }
      
        bool solve(Vector& x) {
	  //SuperClass::solve(x);

	  // Solve system.
	  Dune::InverseOperatorResult result;
	  Vector& istlb = *(this->rhs_);
	  linsolve_->apply(x, istlb, result);
	  SuperClass::checkConvergence(result);
	  
	  if(this->parameters_.scale_linear_system_){
	    this->scaleSolution(x);
	  }	    
	  return this->converged_;
        }
	

        /// Solve the system of linear equations Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x

        /// \copydoc NewtonIterationBlackoilInterface::iterations
        int iterations () const { return this->iterations_; }

        /// \copydoc NewtonIterationBlackoilInterface::parallelInformation
        const boost::any& parallelInformation() const { return this->parallelInformation_; }

    protected:
      
      //using Matrix         = typename MatrixOperator::matrix_type;
      using CouplingMetric = Dune::Amg::Diagonal<pressureIndex>;
      using CritBase       = Dune::Amg::SymmetricCriterion<Matrix, CouplingMetric>;
      using Criterion      = Dune::Amg::CoarsenCriterion<CritBase>;
      using AMG = typename ISTLUtility
	::BlackoilAmgSelector< Matrix, Vector, Vector,POrComm, Criterion, pressureIndex >::AMG;
      typedef std::unique_ptr<typename ScalarProductChooser::ScalarProduct> SPPointer;
      std::unique_ptr< MatrixOperator > opA_;
      std::unique_ptr< OperatorSerial > opASerial_;
      std::unique_ptr< AMG > amg_;
      SPPointer  sp_;
      std::shared_ptr< Dune::BiCGSTABSolver<Vector> > linsolve_;
      //std::shared_ptr< Dune::LinearOperator<Vector,Vector>  > op_;
      //std::shared_ptr< Dune::Preconditioner<Vector,Vector>  > prec_;
      //std::shared_ptr< Dune::ScalarProduct<Vector> > sp_;
      Vector solution_;
      
    }; // end ISTLSolver

} // namespace Opm
#endif
