/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS
  Copyright 2015 IRIS AS

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

#include <config.h>

#include <opm/autodiff/DuneMatrix.hpp>

#include <opm/autodiff/NewtonIterationBlackoilInterleaved.hpp>
#include <opm/autodiff/NewtonIterationUtilities.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#if HAVE_UMFPACK
#include <Eigen/UmfPackSupport>
#else
#include <Eigen/SparseLU>
#endif
#include <opm/common/utility/platform_dependent/reenable_warnings.h>


namespace Opm
{


    typedef AutoDiffBlock<double> ADB;
    typedef ADB::V V;
    typedef ADB::M M;


    namespace detail {
        /**
         * Simple binary operator that always returns 0.1
         * It is used to get the sparsity pattern for our
         * interleaved system, and is marginally faster than using
         * operator+=.
         */
        template<typename Scalar> struct PointOneOp {
            EIGEN_EMPTY_STRUCT_CTOR(PointOneOp)
            Scalar operator()(const Scalar& a, const Scalar& b) const { return 0.1; }
        };
    }


    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables) for a fixed
    /// number of cell variables np .
    template <int np>
    class NewtonIterationBlackoilInterleavedImpl : public NewtonIterationBlackoilInterface
    {
        typedef Dune::FieldVector<double, np    >       VectorBlockType;
        typedef Dune::FieldMatrix<double, np, np>       MatrixBlockType;
        typedef Dune::BCRSMatrix <MatrixBlockType>      Mat;
        typedef Dune::BlockVector<VectorBlockType>      Vector;

    public:
        typedef NewtonIterationBlackoilInterface :: SolutionVector  SolutionVector;
        /// Construct a system solver.
        /// \param[in] param   parameters controlling the behaviour of
        ///                    the preconditioning and choice of
        ///                    linear solvers.
        ///                    Parameters:
        ///                        cpr_relax        (default 1.0) relaxation for the preconditioner
        ///                        cpr_ilu_n        (default 0) use ILU(n) for preconditioning of the linear system
        ///                        cpr_use_amg      (default false) if true, use AMG preconditioner for elliptic part
        ///                        cpr_use_bicgstab (default true)  if true, use BiCGStab (else use CG) for elliptic part
        /// \param[in] parallelInformation In the case of a parallel run
        ///                               with dune-istl the information about the parallelization.
        NewtonIterationBlackoilInterleavedImpl(const parameter::ParameterGroup& param,
                                               const boost::any& parallelInformation_arg=boost::any())
        : iterations_( 0 ),
          parallelInformation_(parallelInformation_arg),
          newton_use_gmres_( param.getDefault("newton_use_gmres", false ) ),
          linear_solver_reduction_( param.getDefault("linear_solver_reduction", 1e-2 ) ),
          linear_solver_maxiter_( param.getDefault("linear_solver_maxiter", 50 ) ),
          linear_solver_restart_( param.getDefault("linear_solver_restart", 40 ) ),
          linear_solver_verbosity_( param.getDefault("linear_solver_verbosity", 0 ))
        {
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
        template<int category=Dune::SolverCategory::sequential, class O, class POrComm>
        void constructPreconditionerAndSolve(O& opA,
                                             Vector& x, Vector& istlb,
                                             const POrComm& parallelInformation_arg,
                                             Dune::InverseOperatorResult& result) const
        {
            // Construct scalar product.
            typedef Dune::ScalarProductChooser<Vector, POrComm, category> ScalarProductChooser;
            typedef std::unique_ptr<typename ScalarProductChooser::ScalarProduct> SPPointer;
            SPPointer sp(ScalarProductChooser::construct(parallelInformation_arg));

            // Construct preconditioner.
            auto precond = constructPrecond(opA, parallelInformation_arg);

            // Communicate if parallel.
            parallelInformation_arg.copyOwnerToAll(istlb, istlb);

            // Solve.
            solve(opA, x, istlb, *sp, *precond, result);
        }


        typedef Dune::SeqILU0<Mat, Vector, Vector> SeqPreconditioner;

        template <class Operator>
        std::unique_ptr<SeqPreconditioner> constructPrecond(Operator& opA, const Dune::Amg::SequentialInformation&) const
        {
            const double relax = 1.0;
            std::unique_ptr<SeqPreconditioner> precond(new SeqPreconditioner(opA.getmat(), relax));
            return precond;
        }

#if HAVE_MPI
        typedef Dune::OwnerOverlapCopyCommunication<int, int> Comm;
        typedef Dune::BlockPreconditioner<Vector, Vector, Comm, SeqPreconditioner> ParPreconditioner;

        template <class Operator>
        std::unique_ptr<ParPreconditioner,
                        AdditionalObjectDeleter<SeqPreconditioner> >
        constructPrecond(Operator& opA, const Comm& comm) const
        {
            typedef AdditionalObjectDeleter<SeqPreconditioner> Deleter;
            typedef std::unique_ptr<ParPreconditioner, Deleter> Pointer;
            int ilu_setup_successful = 1;
            std::string message;
            const double relax = 1.0;
            SeqPreconditioner* seq_precond = nullptr;
            try {
                seq_precond = new SeqPreconditioner(opA.getmat(), relax);
            }
            catch ( Dune::MatrixBlockError error )
            {
                message = error.what();
                std::cerr<<"Exception occured on process " <<
                    comm.communicator().rank() << " during " <<
                    "setup of ILU0 preconditioner with message: " <<
                    message<<std::endl;
                ilu_setup_successful = 0;
            }
            // Check whether there was a problem on some process
            if ( comm.communicator().min(ilu_setup_successful) == 0 )
            {
                if ( seq_precond ) // not null if constructor succeeded
                {
                    // prevent memory leak
                    delete seq_precond;
                    throw Dune::MatrixBlockError();
                }
            }
            return Pointer(new ParPreconditioner(*seq_precond, comm),
                                  Deleter(*seq_precond));
        }
#endif

        /// \brief Solve the system using the given preconditioner and scalar product.
        template <class Operator, class ScalarProd, class Precond>
        void solve(Operator& opA, Vector& x, Vector& istlb, ScalarProd& sp, Precond& precond, Dune::InverseOperatorResult& result) const
        {
            // TODO: Revise when linear solvers interface opm-core is done
            // Construct linear solver.
            // GMRes solver
            if ( newton_use_gmres_ ) {
                Dune::RestartedGMResSolver<Vector> linsolve(opA, sp, precond,
                          linear_solver_reduction_, linear_solver_restart_, linear_solver_maxiter_, linear_solver_verbosity_);
                // Solve system.
                linsolve.apply(x, istlb, result);
            }
            else { // BiCGstab solver
                Dune::BiCGSTABSolver<Vector> linsolve(opA, sp, precond,
                          linear_solver_reduction_, linear_solver_maxiter_, linear_solver_verbosity_);
                // Solve system.
                linsolve.apply(x, istlb, result);
            }
        }

        void formInterleavedSystem(const std::vector<LinearisedBlackoilResidual::ADB>& eqs,
                                   Mat& istlA) const
        {
            assert( np == int(eqs.size()) );
            // Find sparsity structure as union of basic block sparsity structures,
            // corresponding to the jacobians with respect to pressure.
            // Use our custom PointOneOp to get to the union structure.
            // Note that we only iterate over the pressure derivatives on purpose.
            Eigen::SparseMatrix<double, Eigen::ColMajor> col_major = eqs[0].derivative()[0].getSparse();
            detail::PointOneOp<double> point_one;
            for (int phase = 1; phase < np; ++phase) {
                const AutoDiffMatrix::SparseRep& mat = eqs[phase].derivative()[0].getSparse();
                col_major = col_major.binaryExpr(mat, point_one);
            }

            // Automatically convert the column major structure to a row-major structure
            Eigen::SparseMatrix<double, Eigen::RowMajor> row_major = col_major;

            const int size = row_major.rows();
            assert(size == row_major.cols());

            // Create ISTL matrix with interleaved rows and columns (block structured).
            istlA.setSize(row_major.rows(), row_major.cols(), row_major.nonZeros());
            istlA.setBuildMode(Mat::row_wise);
            const int* ia = row_major.outerIndexPtr();
            const int* ja = row_major.innerIndexPtr();
            const typename Mat::CreateIterator endrow = istlA.createend();
            for (typename Mat::CreateIterator row = istlA.createbegin(); row != endrow; ++row) {
                const int ri = row.index();
                for (int i = ia[ri]; i < ia[ri + 1]; ++i) {
                    row.insert(ja[i]);
                }
            }

            // Set all blocks to zero.
            for (int row = 0; row < size; ++row) {
                for (int col_ix = ia[row]; col_ix < ia[row + 1]; ++col_ix) {
                    const int col = ja[col_ix];
                    istlA[row][col] = 0.0;
                }
            }


            /**
             * Go through all jacobians, and insert in correct spot
             *
             * The straight forward way to do this would be to run through each
             * element in the output matrix, and set all block entries by gathering
             * from all "input matrices" (derivatives).
             *
             * A faster alternative is to instead run through each "input matrix" and
             * insert its elements in the correct spot in the output matrix.
             *
             */
            for (int col = 0; col < size; ++col) {
                for (int p1 = 0; p1 < np; ++p1) {
                    for (int p2 = 0; p2 < np; ++p2) {
                        // Note that that since these are CSC and not CSR matrices,
                        // ja contains row numbers instead of column numbers.
                        const AutoDiffMatrix::SparseRep& s = eqs[p1].derivative()[p2].getSparse();
                        const int* ia = s.outerIndexPtr();
                        const int* ja = s.innerIndexPtr();
                        const double* sa = s.valuePtr();
                        for (int elem_ix = ia[col]; elem_ix < ia[col + 1]; ++elem_ix) {
                            const int row = ja[elem_ix];
                            istlA[row][col][p1][p2] = sa[elem_ix];
                        }
                    }
                }
            }
        }


        /// Solve the linear system Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x
        SolutionVector computeNewtonIncrement(const LinearisedBlackoilResidual& residual) const
        {
            // Build the vector of equations.
            //const int np = residual.material_balance_eq.size();
            assert( np == int(residual.material_balance_eq.size()) );
            std::vector<ADB> eqs;
            eqs.reserve(np + 2);
            for (int phase = 0; phase < np; ++phase) {
                eqs.push_back(residual.material_balance_eq[phase]);
            }

            // check if wells are present
            const bool hasWells = residual.well_flux_eq.size() > 0 ;
            std::vector<ADB> elim_eqs;
            if( hasWells )
            {
                eqs.push_back(residual.well_flux_eq);
                eqs.push_back(residual.well_eq);

                // Eliminate the well-related unknowns, and corresponding equations.
                elim_eqs.reserve(2);
                elim_eqs.push_back(eqs[np]);
                eqs = eliminateVariable(eqs, np); // Eliminate well flux unknowns.
                elim_eqs.push_back(eqs[np]);
                eqs = eliminateVariable(eqs, np); // Eliminate well bhp unknowns.
                assert(int(eqs.size()) == np);
            }

            // Scale material balance equations.
            for (int phase = 0; phase < np; ++phase) {
                eqs[phase] = eqs[phase] * residual.matbalscale[phase];
            }

            // calculating the size for b
            int size_b = 0;
            for (int elem = 0; elem < np; ++elem) {
                const int loc_size = eqs[elem].size();
                size_b += loc_size;
            }

            V b(size_b);

            int pos = 0;
            for (int elem = 0; elem < np; ++elem) {
                const int loc_size = eqs[elem].size();
                b.segment(pos, loc_size) = eqs[elem].value();
                pos += loc_size;
            }
            assert(pos == size_b);

            // Create ISTL matrix with interleaved rows and columns (block structured).
            Mat istlA;
            formInterleavedSystem(eqs, istlA);

            // Solve reduced system.
            SolutionVector dx(SolutionVector::Zero(b.size()));

            // Right hand side.
            const int size = istlA.N();
            Vector istlb(size);
            for (int i = 0; i < size; ++i) {
                for( int p = 0, idx = i; p<np; ++p, idx += size ) {
                    istlb[i][p] = b(idx);
                }
            }

            // System solution
            Vector x(istlA.M());
            x = 0.0;

            Dune::InverseOperatorResult result;
            // Parallel version is deactivated until we figure out how to do it properly.
#if HAVE_MPI
            if (parallelInformation_.type() == typeid(ParallelISTLInformation))
            {
                typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
                const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>( parallelInformation_);
                Comm istlComm(info.communicator());
                // As we use a dune-istl with block size np the number of components
                // per parallel is only one.
                info.copyValuesTo(istlComm.indexSet(), istlComm.remoteIndices(),
                                  size, 1);
                // Construct operator, scalar product and vectors needed.
                typedef Dune::OverlappingSchwarzOperator<Mat,Vector,Vector,Comm> Operator;
                Operator opA(istlA, istlComm);
                constructPreconditionerAndSolve<Dune::SolverCategory::overlapping>(opA, x, istlb, istlComm, result);
            }
            else
#endif
            {
                // Construct operator, scalar product and vectors needed.
                typedef Dune::MatrixAdapter<Mat,Vector,Vector> Operator;
                Operator opA(istlA);
                Dune::Amg::SequentialInformation info;
                constructPreconditionerAndSolve(opA, x, istlb, info, result);
            }

            // store number of iterations
            iterations_ = result.iterations;

            // Check for failure of linear solver.
            if (!result.converged) {
                OPM_THROW(LinearSolverProblem, "Convergence failure for linear solver.");
            }

            // Copy solver output to dx.
            for (int i = 0; i < size; ++i) {
                for( int p=0, idx = i; p<np; ++p, idx += size ) {
                    dx(idx) = x[i][p];
                }
            }

            if ( hasWells ) {
                // Compute full solution using the eliminated equations.
                // Recovery in inverse order of elimination.
                dx = recoverVariable(elim_eqs[1], dx, np);
                dx = recoverVariable(elim_eqs[0], dx, np);
            }
            return dx;
        }

    protected:
        mutable int iterations_;
        boost::any parallelInformation_;

        const bool newton_use_gmres_;
        const double linear_solver_reduction_;
        const int    linear_solver_maxiter_;
        const int    linear_solver_restart_;
        const int    linear_solver_verbosity_;

    }; // end NewtonIterationBlackoilInterleavedImpl



    /// Construct a system solver.
    NewtonIterationBlackoilInterleaved::NewtonIterationBlackoilInterleaved(const parameter::ParameterGroup& param,
                                                                           const boost::any& parallelInformation_arg)
      : newtonIncrement_( 6 ),
        param_( param ),
        parallelInformation_(parallelInformation_arg),
        iterations_( 0 )
    {
    }



    /// Solve the linear system Ax = b, with A being the
    /// combined derivative matrix of the residual and b
    /// being the residual itself.
    /// \param[in] residual   residual object containing A and b.
    /// \return               the solution x
    NewtonIterationBlackoilInterleaved::SolutionVector
    NewtonIterationBlackoilInterleaved::computeNewtonIncrement(const LinearisedBlackoilResidual& residual) const
    {
        // get np and call appropriate template method
        const int np = residual.material_balance_eq.size();
        switch( np )
        {
            case 2:
                return computeNewtonIncrementImpl< 2 >( residual );
            case 3:
                return computeNewtonIncrementImpl< 3 >( residual );
            case 4:
                return computeNewtonIncrementImpl< 4 >( residual );
            case 5:
                return computeNewtonIncrementImpl< 5 >( residual );
            case 6:
                return computeNewtonIncrementImpl< 6 >( residual );
            default:
                OPM_THROW(std::runtime_error,"computeNewtonIncrement: number of variables not supported yet. Consult the code to add a case for np = " << np);
                return SolutionVector();
        }
    }


    /// Solve the linear system Ax = b, with A being the
    /// combined derivative matrix of the residual and b
    /// being the residual itself.
    /// \param[in] residual   residual object containing A and b.
    /// \return               the solution x
    template <int np>
    NewtonIterationBlackoilInterleaved::SolutionVector
    NewtonIterationBlackoilInterleaved::computeNewtonIncrementImpl(const LinearisedBlackoilResidual& residual) const
    {
        assert( np < int(newtonIncrement_.size()) );
        // create NewtonIncrement with fixed np
        if( ! newtonIncrement_[ np ] )
            newtonIncrement_[ np ].reset( new NewtonIterationBlackoilInterleavedImpl< np >( param_, parallelInformation_ ) );

        // compute newton increment
        SolutionVector dx = newtonIncrement_[ np ]->computeNewtonIncrement( residual );
        // get number of linear iterations
        iterations_ = newtonIncrement_[ np ]->iterations();
        return dx;
    }


    const boost::any& NewtonIterationBlackoilInterleaved::parallelInformation() const
    {
        return parallelInformation_;
    }



} // namespace Opm

