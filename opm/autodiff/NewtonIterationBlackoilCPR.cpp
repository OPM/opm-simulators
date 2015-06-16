/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
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

#include <opm/autodiff/NewtonIterationBlackoilCPR.hpp>
#include <opm/autodiff/NewtonIterationUtilities.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>


#include <opm/core/utility/platform_dependent/disable_warnings.h>

// #include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/kamg.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <opm/core/utility/platform_dependent/reenable_warnings.h>

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
    typedef Dune::FieldMatrix<double, 1, 1> MatrixBlockType;
    typedef Dune::BCRSMatrix <MatrixBlockType>        Mat;





    /// Construct a system solver.
    NewtonIterationBlackoilCPR::NewtonIterationBlackoilCPR(const parameter::ParameterGroup& param,
                                                           const boost::any& parallelInformation)
      : cpr_param_( param ),
        iterations_( 0 ),
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
    NewtonIterationBlackoilCPR::SolutionVector
    NewtonIterationBlackoilCPR::computeNewtonIncrement(const LinearisedBlackoilResidual& residual) const
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

        // Add material balance equations (or other manipulations) to
        // form pressure equation in top left of full system.
        Eigen::SparseMatrix<double, Eigen::RowMajor> A;
        V b;
        formEllipticSystem(np, eqs, A, b);

        // Scale pressure equation.
        const double pscale = 200*unit::barsa;
        const int nc = residual.material_balance_eq[0].size();
        A.topRows(nc) *= pscale;
        b.topRows(nc) *= pscale;

        // Solve reduced system.
        SolutionVector dx(SolutionVector::Zero(b.size()));

        // Create ISTL matrix.
        DuneMatrix istlA( A );

        // Create ISTL matrix for elliptic part.
        DuneMatrix istlAe( A.topLeftCorner(nc, nc) );

        // Right hand side.
        Vector istlb(istlA.N());
        std::copy_n(b.data(), istlb.size(), istlb.begin());
        // System solution
        Vector x(istlA.M());
        x = 0.0;

        Dune::InverseOperatorResult result;
#if HAVE_MPI
        if(parallelInformation_.type()==typeid(ParallelISTLInformation))
        {
            typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
            const ParallelISTLInformation& info =
                boost::any_cast<const ParallelISTLInformation&>( parallelInformation_);
            Comm istlComm(info.communicator());
            Comm istlAeComm(info.communicator());
            info.copyValuesTo(istlAeComm.indexSet(), istlAeComm.remoteIndices());
            info.copyValuesTo(istlComm.indexSet(), istlComm.remoteIndices(),
                              istlAe.N(), istlA.N()/istlAe.N());
            // Construct operator, scalar product and vectors needed.
            typedef Dune::OverlappingSchwarzOperator<Mat,Vector,Vector,Comm> Operator;
            Operator opA(istlA, istlComm);
            constructPreconditionerAndSolve<Dune::SolverCategory::overlapping>(opA, istlAe, x, istlb, istlComm, istlAeComm, result);
        }
        else
#endif
        {
            // Construct operator, scalar product and vectors needed.
            typedef Dune::MatrixAdapter<Mat,Vector,Vector> Operator;
            Operator opA(istlA);
            Dune::Amg::SequentialInformation info;
            constructPreconditionerAndSolve(opA, istlAe, x, istlb, info, info, result);
        }

        // store number of iterations
        iterations_ = result.iterations;

        // Check for failure of linear solver.
        if (!result.converged) {
            OPM_THROW(LinearSolverProblem, "Convergence failure for linear solver.");
        }

        // Copy solver output to dx.
        std::copy(x.begin(), x.end(), dx.data());

        if( hasWells )
        {
            // Compute full solution using the eliminated equations.
            // Recovery in inverse order of elimination.
            dx = recoverVariable(elim_eqs[1], dx, np);
            dx = recoverVariable(elim_eqs[0], dx, np);
        }
        return dx;
    }





    const boost::any& NewtonIterationBlackoilCPR::parallelInformation() const
    {
        return parallelInformation_;
    }





} // namespace Opm

