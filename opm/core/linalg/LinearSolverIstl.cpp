/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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


#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/core/linalg/LinearSolverIstl.hpp>

#include <opm/core/utility/have_boost_redef.hpp>

// Silence compatibility warning from DUNE headers since we don't use
// the deprecated member anyway (in this compilation unit)
#define DUNE_COMMON_FIELDVECTOR_SIZE_IS_METHOD 1

// TODO: clean up includes.
#include <dune/common/deprecated.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/io.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>

#include <stdexcept>


namespace Opm
{

    using namespace Dune; // While not great, it's okay in a cpp file like this.

    namespace {
        typedef FieldVector<double, 1   > VectorBlockType;
        typedef FieldMatrix<double, 1, 1> MatrixBlockType;
        typedef BCRSMatrix <MatrixBlockType>        Mat;
        typedef BlockVector<VectorBlockType>        Vector;
        typedef MatrixAdapter<Mat,Vector,Vector> Operator;

        LinearSolverInterface::LinearSolverReport
        solveCG_ILU0(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity);

        LinearSolverInterface::LinearSolverReport
        solveCG_AMG(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity);

        LinearSolverInterface::LinearSolverReport
        solveBiCGStab_ILU0(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity);
    } // anonymous namespace




    LinearSolverIstl::LinearSolverIstl()
        : linsolver_residual_tolerance_(1e-8),
          linsolver_verbosity_(0),
          linsolver_type_(CG_AMG),
          linsolver_save_system_(false),
          linsolver_max_iterations_(0),
          linsolver_smooth_steps_(2),
          linsolver_prolongate_factor_(1.6)
    {
    }




    LinearSolverIstl::LinearSolverIstl(const parameter::ParameterGroup& param)
        : linsolver_residual_tolerance_(1e-8),
          linsolver_verbosity_(0),
          linsolver_type_(CG_AMG),
          linsolver_save_system_(false),
          linsolver_max_iterations_(0),
          linsolver_smooth_steps_(2),
          linsolver_prolongate_factor_(1.6)
    {
        linsolver_residual_tolerance_ = param.getDefault("linsolver_residual_tolerance", linsolver_residual_tolerance_);
        linsolver_verbosity_ = param.getDefault("linsolver_verbosity", linsolver_verbosity_);
        linsolver_type_ = LinsolverType(param.getDefault("linsolver_type", int(linsolver_type_)));
        linsolver_save_system_ = param.getDefault("linsolver_save_system", linsolver_save_system_);
        if (linsolver_save_system_) {
            linsolver_save_filename_ = param.getDefault("linsolver_save_filename", std::string("linsys"));
        }
        linsolver_max_iterations_ = param.getDefault("linsolver_max_iterations", linsolver_max_iterations_);
        linsolver_smooth_steps_ = param.getDefault("linsolver_smooth_steps", linsolver_smooth_steps_);
        linsolver_prolongate_factor_ = param.getDegfault("linsolver_prolongate_factor", linsolver_prolongate_factor_);
        
    }




    LinearSolverIstl::~LinearSolverIstl()
    {
    }




    LinearSolverInterface::LinearSolverReport
    LinearSolverIstl::solve(const int size,
                            const int nonzeros,
                            const int* ia,
                            const int* ja,
                            const double* sa,
                            const double* rhs,
                            double* solution) const
    {
        // Build Istl structures from input.
        // System matrix
        Mat A(size, size, nonzeros, Mat::row_wise);
        for (Mat::CreateIterator row = A.createbegin(); row != A.createend(); ++row) {
            int ri = row.index();
            for (int i = ia[ri]; i < ia[ri + 1]; ++i) {
                row.insert(ja[i]);
            }
        }
        for (int ri = 0; ri < size; ++ri) {
            for (int i = ia[ri]; i < ia[ri + 1]; ++i) {
                A[ri][ja[i]] = sa[i];
            }
        }
        // System RHS
        Vector b(size);
        std::copy(rhs, rhs + size, b.begin());
        // System solution
        Vector x(size);
        x = 0.0;

        if (linsolver_save_system_)
        {
            // Save system to files.
            writeMatrixToMatlab(A, linsolver_save_filename_ + "-mat");
            std::string rhsfile(linsolver_save_filename_ + "-rhs");
            std::ofstream rhsf(rhsfile.c_str());
            rhsf.precision(15);
            rhsf.setf(std::ios::scientific | std::ios::showpos);
            std::copy(b.begin(), b.end(),
                      std::ostream_iterator<VectorBlockType>(rhsf, "\n"));
        }
        
        int maxit = linsolver_max_iterations_;
        if (maxit == 0) {
            maxit = A.N();
        }

        LinearSolverReport res;
        switch (linsolver_type_) {
        case CG_ILU0:
            res = solveCG_ILU0(A, x, b, linsolver_residual_tolerance_, maxit, linsolver_verbosity_);
            break;
        case CG_AMG:
            res = solveCG_AMG(A, x, b, linsolver_residual_tolerance_, maxit, linsolver_verbosity_);
            break;
        case BiCGStab_ILU0:
            res = solveBiCGStab_ILU0(A, x, b, linsolver_residual_tolerance_, maxit, linsolver_verbosity_);
            break;
        default:
            std::cerr << "Unknown linsolver_type: " << int(linsolver_type_) << '\n';
            throw std::runtime_error("Unknown linsolver_type");
        }
        std::copy(x.begin(), x.end(), solution);
        return res;
    }

    void LinearSolverIstl::setTolerance(const double tol)
    {
        linsolver_residual_tolerance_ = tol;
    }

    double LinearSolverIstl::getTolerance() const
    {
        return linsolver_residual_tolerance_;
    }



    namespace
    {

    LinearSolverInterface::LinearSolverReport
    solveCG_ILU0(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity)
    {
        Operator opA(A);

        // Construct preconditioner.
        SeqILU0<Mat,Vector,Vector> precond(A, 1.0);

        // Construct linear solver.
        CGSolver<Vector> linsolve(opA, precond, tolerance, maxit, verbosity);

        // Solve system.
        InverseOperatorResult result;
        linsolve.apply(x, b, result);

        // Output results.
        LinearSolverInterface::LinearSolverReport res;
        res.converged = result.converged;
        res.iterations = result.iterations;
        res.residual_reduction = result.reduction;
        return res;
    }




    LinearSolverInterface::LinearSolverReport
    solveCG_AMG(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity)
    {
        // Solve with AMG solver.
#define FIRST_DIAGONAL 1
#define SYMMETRIC 1
#define SMOOTHER_ILU 1
#define ANISOTROPIC_3D 0

#if FIRST_DIAGONAL
        typedef Amg::FirstDiagonal CouplingMetric;
#else
        typedef Amg::RowSum        CouplingMetric;
#endif

#if SYMMETRIC
        typedef Amg::SymmetricCriterion<Mat,CouplingMetric>   CriterionBase;
#else
        typedef Amg::UnSymmetricCriterion<Mat,CouplingMetric> CriterionBase;
#endif

#if SMOOTHER_ILU
        typedef SeqILU0<Mat,Vector,Vector>        Smoother;
#else
        typedef SeqSSOR<Mat,Vector,Vector>        Smoother;
#endif
        typedef Amg::CoarsenCriterion<CriterionBase> Criterion;
        typedef Amg::AMG<Operator,Vector,Smoother>   Precond;

        Operator opA(A);

        // Construct preconditioner.
        double relax = 1;
        Precond::SmootherArgs smootherArgs;
        smootherArgs.relaxationFactor = relax;
        Criterion criterion;
        criterion.setDebugLevel(verbosity);
#if ANISOTROPIC_3D
        criterion.setDefaultValuesAnisotropic(3, 2);
#endif
        criterion.setProlongateDampingFactor(linsolve_prolongate_factor_);
        Precond precond(opA, criterion, smootherArgs, 1, linsolve_smooth_steps_,
                        linsolve_smooth_steps_);

        // Construct linear solver.
        CGSolver<Vector> linsolve(opA, precond, tolerance, maxit, verbosity);

        // Solve system.
        InverseOperatorResult result;
        linsolve.apply(x, b, result);

        // Output results.
        LinearSolverInterface::LinearSolverReport res;
        res.converged = result.converged;
        res.iterations = result.iterations;
        res.residual_reduction = result.reduction;
        return res;
    }



    LinearSolverInterface::LinearSolverReport
    solveBiCGStab_ILU0(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity)
    {
        Operator opA(A);

        // Construct preconditioner.
        SeqILU0<Mat,Vector,Vector> precond(A, 1.0);

        // Construct linear solver.
        BiCGSTABSolver<Vector> linsolve(opA, precond, tolerance, maxit, verbosity);

        // Solve system.
        InverseOperatorResult result;
        linsolve.apply(x, b, result);

        // Output results.
        LinearSolverInterface::LinearSolverReport res;
        res.converged = result.converged;
        res.iterations = result.iterations;
        res.residual_reduction = result.reduction;
        return res;
    }




    } // anonymous namespace


} // namespace Opm

