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
#include <dune/istl/paamg/kamg.hh>

#include <stdexcept>


namespace Opm
{


    namespace {
        typedef Dune::FieldVector<double, 1   > VectorBlockType;
        typedef Dune::FieldMatrix<double, 1, 1> MatrixBlockType;
        typedef Dune::BCRSMatrix <MatrixBlockType>        Mat;
        typedef Dune::BlockVector<VectorBlockType>        Vector;
        typedef Dune::MatrixAdapter<Mat,Vector,Vector> Operator;

        LinearSolverInterface::LinearSolverReport
        solveCG_ILU0(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity);

        LinearSolverInterface::LinearSolverReport
        solveCG_AMG(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity,
                    double prolongateFactor, int smoothsteps);

#ifdef HAS_DUNE_FAST_AMG
        LinearSolverInterface::LinearSolverReport
        solveKAMG(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity,
                  double prolongateFactor, int smoothsteps);

        LinearSolverInterface::LinearSolverReport
        solveFastAMG(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity,
                     double prolongateFactor, int smoothsteps);
#endif
    
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
        linsolver_prolongate_factor_ = param.getDefault("linsolver_prolongate_factor", linsolver_prolongate_factor_);
        
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
            maxit = 5000;
        }

        LinearSolverReport res;        
        switch (linsolver_type_) {
        case CG_ILU0:
            res = solveCG_ILU0(A, x, b, linsolver_residual_tolerance_, maxit, linsolver_verbosity_);
            break;
        case CG_AMG:
            res = solveCG_AMG(A, x, b, linsolver_residual_tolerance_, maxit, linsolver_verbosity_,
                              linsolver_prolongate_factor_, linsolver_smooth_steps_);
            break;
        case KAMG:
#ifdef HAS_DUNE_FAST_AMG
            res = solveKAMG(A, x, b, linsolver_residual_tolerance_, maxit, linsolver_verbosity_,
                            linsolver_prolongate_factor_, linsolver_smooth_steps_);
#else
            throw std::runtime_error("KAMG not supported with this version of DUNE");
#endif
            break;
        case FastAMG:
#ifdef HAS_DUNE_FAST_AMG
            res = solveFastAMG(A, x, b, linsolver_residual_tolerance_, maxit, linsolver_verbosity_,
                               linsolver_prolongate_factor_, linsolver_smooth_steps_);
#else
	    if(linsolver_verbosity_)
	      std::cerr<<"Fast AMG is not available; falling back to CG preconditioned with the normal one"<<std::endl;
	    res = solveCG_AMG(A, x, b, linsolver_residual_tolerance_, maxit, linsolver_verbosity_);
#endif
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
        Dune::SeqILU0<Mat,Vector,Vector> precond(A, 1.0);

        // Construct linear solver.
        Dune::CGSolver<Vector> linsolve(opA, precond, tolerance, maxit, verbosity);

        // Solve system.
        Dune::InverseOperatorResult result;
        linsolve.apply(x, b, result);

        // Output results.
        LinearSolverInterface::LinearSolverReport res;
        res.converged = result.converged;
        res.iterations = result.iterations;
        res.residual_reduction = result.reduction;
        return res;
    }



#define FIRST_DIAGONAL 1
#define SYMMETRIC 1
#define SMOOTHER_ILU 0
#define ANISOTROPIC_3D 0

    LinearSolverInterface::LinearSolverReport
    solveCG_AMG(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity,
                double linsolver_prolongate_factor, int linsolver_smooth_steps)
    {
        // Solve with AMG solver.

#if FIRST_DIAGONAL
        typedef Dune::Amg::FirstDiagonal CouplingMetric;
#else
        typedef Dune::Amg::RowSum        CouplingMetric;
#endif

#if SYMMETRIC
        typedef Dune::Amg::SymmetricCriterion<Mat,CouplingMetric>   CriterionBase;
#else
        typedef Dune::Amg::UnSymmetricCriterion<Mat,CouplingMetric> CriterionBase;
#endif

#if SMOOTHER_ILU
        typedef Dune::SeqILU0<Mat,Vector,Vector>        Smoother;
#else
        typedef Dune::SeqSOR<Mat,Vector,Vector>        Smoother;
#endif
        typedef Dune::Amg::CoarsenCriterion<CriterionBase> Criterion;
        typedef Dune::Amg::AMG<Operator,Vector,Smoother>   Precond;

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
        Precond precond(opA, criterion, smootherArgs, 1, linsolver_smooth_steps,
                        linsolver_smooth_steps);
        
        // Construct linear solver.
        Dune::CGSolver<Vector> linsolve(opA, precond, tolerance, maxit, verbosity);

        // Solve system.
        Dune::InverseOperatorResult result;
        linsolve.apply(x, b, result);

        // Output results.
        LinearSolverInterface::LinearSolverReport res;
        res.converged = result.converged;
        res.iterations = result.iterations;
        res.residual_reduction = result.reduction;
        return res;
    }


#ifdef HAS_DUNE_FAST_AMG
    LinearSolverInterface::LinearSolverReport
    solveKAMG(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity,
              double linsolver_prolongate_factor, int linsolver_smooth_steps)
    {
        // Solve with AMG solver.

#if FIRST_DIAGONAL
        typedef Dune::Amg::FirstDiagonal CouplingMetric;
#else
        typedef Dune::Amg::RowSum        CouplingMetric;
#endif

#if SYMMETRIC
        typedef Dune::Amg::SymmetricCriterion<Mat,CouplingMetric>   CriterionBase;
#else
        typedef Dune::Amg::UnSymmetricCriterion<Mat,CouplingMetric> CriterionBase;
#endif

#if SMOOTHER_ILU
        typedef Dune::SeqILU0<Mat,Vector,Vector>        Smoother;
#else
        typedef Dune::SeqSOR<Mat,Vector,Vector>        Smoother;
#endif
        typedef Dune::Amg::CoarsenCriterion<CriterionBase> Criterion;
        typedef Dune::Amg::KAMG<Operator,Vector,Smoother,Dune::Amg::SequentialInformation>   Precond;
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
        criterion.setProlongationDampingFactor(linsolver_prolongate_factor);
        Precond precond(opA, criterion, smootherArgs, 1, linsolver_smooth_steps,
                        linsolver_smooth_steps);

        // Construct linear solver.
        Dune::GeneralizedPCGSolver<Vector> linsolve(opA, precond, tolerance, maxit, verbosity);

        // Solve system.
        Dune::InverseOperatorResult result;
        linsolve.apply(x, b, result);

        // Output results.
        LinearSolverInterface::LinearSolverReport res;
        res.converged = result.converged;
        res.iterations = result.iterations;
        res.residual_reduction = result.reduction;
        return res;
    }

    LinearSolverInterface::LinearSolverReport
    solveFastAMG(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity,
                 double linsolver_prolongate_factor, int linsolver_smooth_steps)
    {
        // Solve with AMG solver.

#if FIRST_DIAGONAL
        typedef Dune::Amg::FirstDiagonal CouplingMetric;
#else
        typedef Dune::Amg::RowSum        CouplingMetric;
#endif

#if SYMMETRIC
        typedef Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricMatrixDependency<Mat,CouplingMetric> > CriterionBase;
#else
	typedef Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricMatrixDependency<Mat,CouplingMetric> > CriterionBase;
#endif

        typedef Dune::Amg::CoarsenCriterion<CriterionBase> Criterion;
        typedef Dune::Amg::FastAMG<Operator,Vector>   Precond;

        Operator opA(A);
        
        // Construct preconditioner.
	Criterion criterion;
        criterion.setDebugLevel(verbosity);
#if ANISOTROPIC_3D
        criterion.setDefaultValuesAnisotropic(3, 2);
#else 
        criterion.setDefaultValuesIsotropic(3,2);
#endif
        criterion.setProlongationDampingFactor(1.6);
	Dune::Amg::Parameters parms;
	parms.setDebugLevel(verbosity);
	parms.setNoPreSmoothSteps(1);
	parms.setNoPostSmoothSteps(1);
        Precond precond(opA, criterion, parms);

        // Construct linear solver.
        Dune::GeneralizedPCGSolver<Vector> linsolve(opA, precond, tolerance, maxit, verbosity);

        // Solve system.
        Dune::InverseOperatorResult result;
        linsolve.apply(x, b, result);

        // Output results.
        LinearSolverInterface::LinearSolverReport res;
        res.converged = result.converged;
        res.iterations = result.iterations;
        res.residual_reduction = result.reduction;
        return res;
    }
#endif


    LinearSolverInterface::LinearSolverReport
    solveBiCGStab_ILU0(const Mat& A, Vector& x, Vector& b, double tolerance, int maxit, int verbosity)
    {
        Operator opA(A);

        // Construct preconditioner.
        Dune::SeqILU0<Mat,Vector,Vector> precond(A, 1.0);

        // Construct linear solver.
        Dune::BiCGSTABSolver<Vector> linsolve(opA, precond, tolerance, maxit, verbosity);

        // Solve system.
        Dune::InverseOperatorResult result;
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

