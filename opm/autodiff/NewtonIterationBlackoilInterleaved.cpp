/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

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





    /// Construct a system solver.
    NewtonIterationBlackoilInterleaved::NewtonIterationBlackoilInterleaved(const parameter::ParameterGroup& param,
                                                                           const boost::any& parallelInformation_arg)
      : iterations_( 0 ),
        parallelInformation_(parallelInformation_arg),
        newton_use_gmres_( param.getDefault("newton_use_gmres", false ) ),
        linear_solver_reduction_( param.getDefault("linear_solver_reduction", 1e-2 ) ),
        linear_solver_maxiter_( param.getDefault("linear_solver_maxiter", 50 ) ),
        linear_solver_restart_( param.getDefault("linear_solver_restart", 40 ) ),
        linear_solver_verbosity_( param.getDefault("linear_solver_verbosity", 0 ))
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
        // Build the vector of equations.
        const int np = residual.material_balance_eq.size();
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
            istlb[i][0] = b(i);
            istlb[i][1] = b(size + i);
            istlb[i][2] = b(2*size + i);
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
            dx(i)          = x[i][0];
            dx(size + i)   = x[i][1];
            dx(2*size + i) = x[i][2];
        }

        if ( hasWells ) {
            // Compute full solution using the eliminated equations.
            // Recovery in inverse order of elimination.
            dx = recoverVariable(elim_eqs[1], dx, np);
            dx = recoverVariable(elim_eqs[0], dx, np);
        }
        return dx;
    }




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


    void NewtonIterationBlackoilInterleaved::formInterleavedSystem(const std::vector<ADB>& eqs,
                                                                   Mat& istlA) const
    {
        const int np = eqs.size();
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
        {
            assert(np == 3);
            istlA.setSize(row_major.rows(), row_major.cols(), row_major.nonZeros());
            istlA.setBuildMode(Mat::row_wise);
            const int* ia = row_major.outerIndexPtr();
            const int* ja = row_major.innerIndexPtr();
            for (Mat::CreateIterator row = istlA.createbegin(); row != istlA.createend(); ++row) {
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





    const boost::any& NewtonIterationBlackoilInterleaved::parallelInformation() const
    {
        return parallelInformation_;
    }





} // namespace Opm

