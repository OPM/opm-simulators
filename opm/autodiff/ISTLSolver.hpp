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

#ifndef OPM_ISTLSOLVER_HEADER_INCLUDED
#define OPM_ISTLSOLVER_HEADER_INCLUDED

#include <opm/autodiff/MPIUtilities.hpp>
#include <opm/autodiff/BlackoilAmg.hpp>
#include <opm/autodiff/CPRPreconditioner.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterleaved.hpp>
#include <opm/autodiff/NewtonIterationUtilities.hpp>
#include <opm/autodiff/ParallelRestrictedAdditiveSchwarz.hpp>
#include <opm/autodiff/ParallelOverlappingILU0.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/MatrixBlock.hpp>

#include <opm/common/Exceptions.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

namespace Opm
{
    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables) for a fixed
    /// number of cell variables np .
    /// \tparam MatrixBlockType The type of the matrix block used.
    /// \tparam VectorBlockType The type of the vector block used.
    /// \tparam pressureIndex The index of the pressure component in the vector
    ///                       vector block. It is used to guide the AMG coarsening.
    ///                       Default is zero.
    template < class MatrixBlockType, class VectorBlockType, int pressureIndex=0 >
    class ISTLSolver : public NewtonIterationBlackoilInterface
    {
        typedef typename MatrixBlockType :: field_type  Scalar;

        typedef Dune::BCRSMatrix <MatrixBlockType>      Matrix;
        typedef Dune::BlockVector<VectorBlockType>      Vector;

    public:
        typedef Dune::AssembledLinearOperator< Matrix, Vector, Vector > AssembledLinearOperatorType;

        typedef NewtonIterationBlackoilInterface :: SolutionVector  SolutionVector;
        /// Construct a system solver.
        /// \param[in] param   parameters controlling the behaviour of the linear solvers
        /// \param[in] parallelInformation In the case of a parallel run
        ///                                with dune-istl the information about the parallelization.
        ISTLSolver(const NewtonIterationBlackoilInterleavedParameters& param,
                   const boost::any& parallelInformation_arg=boost::any())
        : iterations_( 0 ),
          parallelInformation_(parallelInformation_arg),
          isIORank_(isIORank(parallelInformation_arg)),
          parameters_( param )
        {
        }

        /// Construct a system solver.
        /// \param[in] param   ParameterGroup controlling the behaviour of the linear solvers
        /// \param[in] parallelInformation In the case of a parallel run
        ///                                with dune-istl the information about the parallelization.
        ISTLSolver(const ParameterGroup& param,
                   const boost::any& parallelInformation_arg=boost::any())
        : iterations_( 0 ),
          parallelInformation_(parallelInformation_arg),
          isIORank_(isIORank(parallelInformation_arg)),
          parameters_( param )
        {
        }

        // dummy method that is not implemented for this class
        SolutionVector computeNewtonIncrement(const LinearisedBlackoilResidual&) const
        {
            OPM_THROW(std::logic_error,"This method is not implemented");
            return SolutionVector();
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

    public:
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
                    // We should never end up here as this code is
                    // only part of flow_legacy  and if use_cpr_ is
                    // true for flow_legacy then solver_approach=cpr
                    // was specified and NewtonIterationBlackoilCPR
                    // is used as a solve and not ISTLSolver.
                    OPM_THROW(std::logic_error,
                              "This code path should bever be exectuded for parameters_.use_cpr_="
                              <<parameters_.use_cpr_<<" in flow_legacy.");
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


        template <class MatrixOperator, class POrComm, class AMG >
        void
        constructAMGPrecond(MatrixOperator& opA, const POrComm& comm, std::unique_ptr< AMG >& amg, std::unique_ptr< MatrixOperator >&, const double relax,
                            const MILU_VARIANT milu) const
        {
            ISTLUtility::template createAMGPreconditionerPointer<pressureIndex>( opA, relax,
                                                                                 milu, comm, amg );
        }

        /// \brief Solve the system using the given preconditioner and scalar product.
        template <class Operator, class ScalarProd, class Precond>
        void solve(Operator& opA, Vector& x, Vector& istlb, ScalarProd& sp, Precond& precond, Dune::InverseOperatorResult& result) const
        {
            // TODO: Revise when linear solvers interface opm-core is done
            // Construct linear solver.
            // GMRes solver
            int verbosity = ( isIORank_ ) ? parameters_.linear_solver_verbosity_ : 0;

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

            // Check for failure of linear solver.
            if (!parameters_.ignoreConvergenceFailure_ && !result.converged) {
                const std::string msg("Convergence failure for linear solver.");
                OPM_THROW_NOLOG(LinearSolverProblem, msg);
            }
        }
    protected:
        mutable int iterations_;
        boost::any parallelInformation_;
        bool isIORank_;

        NewtonIterationBlackoilInterleavedParameters parameters_;
    }; // end ISTLSolver

} // namespace Opm
#endif
