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

// Define making clear that the simulator supports AMG
#define FLOW_SUPPORT_AMG !defined(HAVE_UMFPACK)

#include <opm/autodiff/CPRPreconditioner.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterleaved.hpp>
#include <opm/autodiff/NewtonIterationUtilities.hpp>
#include <opm/autodiff/ParallelRestrictedAdditiveSchwarz.hpp>
#include <opm/autodiff/ParallelOverlappingILU0.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/DuneMatrix.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>

#include <opm/autodiff/ISTLSolver.hpp>

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#if HAVE_UMFPACK
#include <Eigen/UmfPackSupport>
#else
#include <Eigen/SparseLU>
#endif
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

namespace Opm
{

    namespace detail {
        /**
         * Simple binary operator that always returns 0.1
         * It is used to get the sparsity pattern for our
         * interleaved system, and is marginally faster than using
         * operator+=.
         */
        template<typename Scalar> struct PointOneOp {
            EIGEN_EMPTY_STRUCT_CTOR(PointOneOp)
            Scalar operator()(const Scalar&, const Scalar&) const { return 0.1; }
        };
    }


    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables) for a fixed
    /// number of cell variables np .
    template <int np, class ScalarT = double >
    class NewtonIterationBlackoilInterleavedImpl : public NewtonIterationBlackoilInterface
    {
        typedef ScalarT                                 Scalar;
        typedef Dune::FieldVector<Scalar, np    >       VectorBlockType;

        typedef Dune::MatrixBlock<Scalar, np, np >      MatrixBlockType;

        typedef Dune::BCRSMatrix <MatrixBlockType>      Mat;
        typedef Dune::BlockVector<VectorBlockType>      Vector;

        typedef Opm::ISTLSolver< MatrixBlockType, VectorBlockType > ISTLSolverType;

    public:
        typedef NewtonIterationBlackoilInterface :: SolutionVector  SolutionVector;
        /// Construct a system solver.
        /// \param[in] param   parameters controlling the behaviour of the linear solvers
        /// \param[in] parallelInformation In the case of a parallel run
         ///                               with dune-istl the information about the parallelization.
        NewtonIterationBlackoilInterleavedImpl(const NewtonIterationBlackoilInterleavedParameters& param,
                                               const boost::any& parallelInformation_arg=boost::any())
        : istlSolver_( param, parallelInformation_arg ),
          parameters_( param )
        {
        }

        /// Solve the system of linear equations Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x

        /// \copydoc NewtonIterationBlackoilInterface::iterations
        int iterations () const { return istlSolver_.iterations(); }

        /// \copydoc NewtonIterationBlackoilInterface::parallelInformation
        const boost::any& parallelInformation() const { return istlSolver_.parallelInformation(); }

    public:
        void formInterleavedSystem(const std::vector<LinearisedBlackoilResidual::ADB>& eqs,
                                   Mat& istlA) const
        {
            assert( np == int(eqs.size()) );
            // Find sparsity structure as union of basic block sparsity structures,
            // corresponding to the jacobians with respect to pressure.
            // Use our custom PointOneOp to get to the union structure.
            // As default we only iterate over the pressure derivatives.
            Eigen::SparseMatrix<double, Eigen::ColMajor> col_major = eqs[0].derivative()[0].getSparse();
            detail::PointOneOp<double> point_one;
            for (int phase = 1; phase < np; ++phase) {
                const AutoDiffMatrix::SparseRep& mat = eqs[phase].derivative()[0].getSparse();
                col_major = col_major.binaryExpr(mat, point_one);
            }
            // For some cases (for instance involving Solvent flow) the reasoning for only adding
            // the pressure derivatives fails. As getting the sparsity pattern is non-trivial, in terms
            // of work, the full sparsity pattern is only added when required.
            if (parameters_.require_full_sparsity_pattern_) {
                for (int p1 = 0; p1 < np; ++p1) {
                    for (int p2 = 1; p2 < np; ++p2) { // pressure is already added
                        const AutoDiffMatrix::SparseRep& mat = eqs[p1].derivative()[p2].getSparse();
                        col_major = col_major.binaryExpr(mat, point_one);
                    }
                }
            }

            // Automatically convert the column major structure to a row-major structure
            Eigen::SparseMatrix<double, Eigen::RowMajor> row_major = col_major;

            const int size = row_major.rows();
            assert(size == row_major.cols());

            {
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
            }

            /*
            // not neeeded since MatrixBlock initially zeros all elements during construction
            // Set all blocks to zero.
            for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
                for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                    *col = 0.0;
                }
            }
            */

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
            for (int p1 = 0; p1 < np; ++p1) {
                for (int p2 = 0; p2 < np; ++p2) {
                    // Note that that since these are CSC and not CSR matrices,
                    // ja contains row numbers instead of column numbers.
                    const AutoDiffMatrix::SparseRep& s = eqs[p1].derivative()[p2].getSparse();
                    const int* ia = s.outerIndexPtr();
                    const int* ja = s.innerIndexPtr();
                    const double* sa = s.valuePtr();
                    for (int col = 0; col < size; ++col) {
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
            typedef LinearisedBlackoilResidual::ADB  ADB;
            typedef ADB::V   V;

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

            // solve linear system using ISTL methods
            istlSolver_.solve( istlA, x, istlb );

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
        ISTLSolverType istlSolver_;
        NewtonIterationBlackoilInterleavedParameters parameters_;
    }; // end NewtonIterationBlackoilInterleavedImpl



    /// Construct a system solver.
    NewtonIterationBlackoilInterleaved::NewtonIterationBlackoilInterleaved(const ParameterGroup& param,
                                                                           const boost::any& parallelInformation_arg)
      : newtonIncrementDoublePrecision_(),
        newtonIncrementSinglePrecision_(),
        parameters_( param ),
        parallelInformation_(parallelInformation_arg),
        iterations_( 0 )
    {
    }

    namespace detail {

        template< int NP, class Scalar >
        struct NewtonIncrement
        {
            template <class NewtonIncVector>
            static const NewtonIterationBlackoilInterface&
            get( NewtonIncVector& newtonIncrements,
                 const NewtonIterationBlackoilInterleavedParameters& param,
                 const boost::any& parallelInformation,
                 const int np )
            {
                if( np == NP )
                {
                    assert( np < int(newtonIncrements.size()) );
                    // create NewtonIncrement with fixed np
                    if( ! newtonIncrements[ NP ] )
                        newtonIncrements[ NP ].reset( new NewtonIterationBlackoilInterleavedImpl< NP, Scalar >( param, parallelInformation ) );
                    return *(newtonIncrements[ NP ]);
                }
                else
                {
                    return NewtonIncrement< NP-1, Scalar >::get(newtonIncrements, param, parallelInformation, np );
                }
            }
        };

        template<class Scalar>
        struct NewtonIncrement< 0, Scalar >
        {
            template <class NewtonIncVector>
            static const NewtonIterationBlackoilInterface&
            get( NewtonIncVector&,
                 const NewtonIterationBlackoilInterleavedParameters&,
                 const boost::any&,
                 const int np )
            {
                OPM_THROW(std::runtime_error,"NewtonIncrement::get: number of variables not supported yet. Adjust maxNumberEquations appropriately to cover np = " << np);
            }
        };





        std::pair<NewtonIterationBlackoilInterleaved::SolutionVector, Dune::InverseOperatorResult>
        computePressureIncrement(const LinearisedBlackoilResidual& residual)
        {
            typedef LinearisedBlackoilResidual::ADB ADB;

            // Build the vector of equations (should be just a single material balance equation
            // in which the pressure equation is stored).
            const int np = residual.material_balance_eq.size();
            assert(np == 1);
            std::vector<ADB> eqs;
            eqs.reserve(np + 2);
            for (int phase = 0; phase < np; ++phase) {
                eqs.push_back(residual.material_balance_eq[phase]);
            }

            // Check if wells are present.
            const bool hasWells = residual.well_flux_eq.size() > 0 ;
            std::vector<ADB> elim_eqs;
            if (hasWells) {
                // Eliminate the well-related unknowns, and corresponding equations.
                eqs.push_back(residual.well_flux_eq);
                eqs.push_back(residual.well_eq);
                elim_eqs.reserve(2);
                elim_eqs.push_back(eqs[np]);
                eqs = eliminateVariable(eqs, np); // Eliminate well flux unknowns.
                elim_eqs.push_back(eqs[np]);
                eqs = eliminateVariable(eqs, np); // Eliminate well bhp unknowns.
                assert(int(eqs.size()) == np);
            }

            // Solve the linearised oil equation.
            Eigen::SparseMatrix<double, Eigen::RowMajor> eigenA = eqs[0].derivative()[0].getSparse();
            DuneMatrix opA(eigenA);
            const int size = eqs[0].size();
            typedef Dune::BlockVector<Dune::FieldVector<double, 1> > Vector1;
            Vector1 x;
            x.resize(size);
            x = 0.0;
            Vector1 b;
            b.resize(size);
            b = 0.0;
            std::copy_n(eqs[0].value().data(), size, b.begin());

            // Solve with AMG solver.
            typedef Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1> > Mat;
            typedef Dune::MatrixAdapter<Mat, Vector1, Vector1> Operator;
            Operator sOpA(opA);

            typedef Dune::Amg::SequentialInformation ParallelInformation;
            typedef Dune::SeqILU0<Mat,Vector1,Vector1> EllipticPreconditioner;
            typedef EllipticPreconditioner Smoother;
            typedef Dune::Amg::AMG<Operator, Vector1, Smoother, ParallelInformation> AMG;
            typedef Dune::Amg::FirstDiagonal CouplingMetric;
            typedef Dune::Amg::SymmetricCriterion<Mat, CouplingMetric> CritBase;
            typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;

            // TODO: revise choice of parameters
            const int coarsenTarget = 1200;
            Criterion criterion(15, coarsenTarget);
            criterion.setDebugLevel(0); // no debug information, 1 for printing hierarchy information
            criterion.setDefaultValuesIsotropic(2);
            criterion.setNoPostSmoothSteps(1);
            criterion.setNoPreSmoothSteps(1);

            // for DUNE 2.2 we also need to pass the smoother args
            typedef typename AMG::Smoother Smoother;
            typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments  SmootherArgs;
            SmootherArgs  smootherArgs;
            smootherArgs.iterations = 1;
            smootherArgs.relaxationFactor = 1.0;

            AMG precond(sOpA, criterion, smootherArgs);

            const int verbosity = 0;
            const int maxit = 30;
            const double tolerance = 1e-5;

            // Construct linear solver.
            Dune::BiCGSTABSolver<Vector1> linsolve(sOpA, precond, tolerance, maxit, verbosity);

            // Solve system.
            Dune::InverseOperatorResult result;
            linsolve.apply(x, b, result);

            // Check for failure of linear solver.
            if (!result.converged) {
                const std::string msg("Convergence failure for linear solver in computePressureIncrement().");
                OpmLog::problem(msg);
                OPM_THROW_NOLOG(LinearSolverProblem, msg);
            }

            // Copy solver output to dx.
            NewtonIterationBlackoilInterleaved::SolutionVector dx(size);
            for (int i = 0; i < size; ++i) {
                dx(i)          = x[i];
            }

            if (hasWells) {
                // Compute full solution using the eliminated equations.
                // Recovery in inverse order of elimination.
                dx = recoverVariable(elim_eqs[1], dx, np);
                dx = recoverVariable(elim_eqs[0], dx, np);
            }
            return std::make_pair(dx, result);
        }


    } // end namespace detail


    NewtonIterationBlackoilInterleaved::SolutionVector
    NewtonIterationBlackoilInterleaved::computeNewtonIncrement(const LinearisedBlackoilResidual& residual) const
    {
        // get np and call appropriate template method
        const int np = residual.material_balance_eq.size();
        if (np == 1) {
            auto result = detail::computePressureIncrement(residual);
            iterations_ = result.second.iterations;
            return result.first;
        }

        const NewtonIterationBlackoilInterface& newtonIncrement = residual.singlePrecision ?
            detail::NewtonIncrement< maxNumberEquations_, float  > :: get( newtonIncrementSinglePrecision_, parameters_, parallelInformation_, np ) :
            detail::NewtonIncrement< maxNumberEquations_, double > :: get( newtonIncrementDoublePrecision_, parameters_, parallelInformation_, np );

        // compute newton increment
        SolutionVector dx = newtonIncrement.computeNewtonIncrement( residual );
        // get number of linear iterations
        iterations_ = newtonIncrement.iterations();
        return dx;
    }

    const boost::any& NewtonIterationBlackoilInterleaved::parallelInformation() const
    {
        return parallelInformation_;
    }



} // namespace Opm

