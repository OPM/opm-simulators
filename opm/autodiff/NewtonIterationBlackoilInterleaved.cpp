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

#if HAVE_UMFPACK
#include <Eigen/UmfPackSupport>
#else
#include <Eigen/SparseLU>
#endif


namespace Opm
{


    typedef AutoDiffBlock<double> ADB;
    typedef ADB::V V;
    typedef ADB::M M;





    /// Construct a system solver.
    NewtonIterationBlackoilInterleaved::NewtonIterationBlackoilInterleaved(const parameter::ParameterGroup& param,
                                                                           const boost::any& parallelInformation)
      : iterations_( 0 ),
        parallelInformation_(parallelInformation),
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
        const double matbalscale[3] = { 1.1169, 1.0031, 0.0031 }; // HACK hardcoded instead of computed.
        for (int phase = 0; phase < np; ++phase) {
            eqs[phase] = eqs[phase] * matbalscale[phase];
        }

        // Form modified system.
        Eigen::SparseMatrix<double, Eigen::RowMajor> A;
        V b;
        formEllipticSystem(np, eqs, A, b);

        // Create ISTL matrix with interleaved rows and columns (block structured).
        Mat istlA;
        formInterleavedSystem(eqs, A, istlA);

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
            info.copyValuesTo(istlComm.indexSet(), istlComm.remoteIndices(),
                              size, np);
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





    void NewtonIterationBlackoilInterleaved::formInterleavedSystem(const std::vector<ADB>& eqs,
                                                                   const Eigen::SparseMatrix<double, Eigen::RowMajor>& A,
                                                                   Mat& istlA) const
    {
        const int np = eqs.size();

        // Find sparsity structure as union of basic block sparsity structures,
        // corresponding to the jacobians with respect to pressure.
        // Use addition to get to the union structure.
        Eigen::SparseMatrix<double> structure = eqs[0].derivative()[0];
        for (int phase = 0; phase < np; ++phase) {
            structure += eqs[phase].derivative()[0];
        }
        Eigen::SparseMatrix<double, Eigen::RowMajor> s = structure;

        // Create ISTL matrix with interleaved rows and columns (block structured).
        assert(np == 3);
        istlA.setSize(s.rows(), s.cols(), s.nonZeros());
        istlA.setBuildMode(Mat::row_wise);
        const int* ia = s.outerIndexPtr();
        const int* ja = s.innerIndexPtr();
        for (Mat::CreateIterator row = istlA.createbegin(); row != istlA.createend(); ++row) {
            int ri = row.index();
            for (int i = ia[ri]; i < ia[ri + 1]; ++i) {
                row.insert(ja[i]);
            }
        }
        const int size = s.rows();
        Span span[3] = { Span(size, 1, 0),
                         Span(size, 1, size),
                         Span(size, 1, 2*size) };
        for (int row = 0; row < size; ++row) {
            for (int col_ix = ia[row]; col_ix < ia[row + 1]; ++col_ix) {
                const int col = ja[col_ix];
                MatrixBlockType block;
                for (int p1 = 0; p1 < np; ++p1) {
                    for (int p2 = 0; p2 < np; ++p2) {
                        block[p1][p2] = A.coeff(span[p1][row], span[p2][col]);
                    }
                }
                istlA[row][col] = block;
            }
        }
    }





    const boost::any& NewtonIterationBlackoilInterleaved::parallelInformation() const
    {
        return parallelInformation_;
    }





} // namespace Opm

