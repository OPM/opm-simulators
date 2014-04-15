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
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

// Silence compatibility warning from DUNE headers since we don't use
// the deprecated member anyway (in this compilation unit)
#define DUNE_COMMON_FIELDVECTOR_SIZE_IS_METHOD 1

// TODO: clean up includes.
#include <dune/common/deprecated.hh>
#include <dune/common/version.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/io.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/kamg.hh>
#include <dune/istl/paamg/pinfo.hh>

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
#include <dune/istl/paamg/fastamg.hh>
#endif

#include <stdexcept>
#include <iostream>
#include <type_traits>

namespace Opm
{


    namespace {
        typedef Dune::FieldVector<double, 1   > VectorBlockType;
        typedef Dune::FieldMatrix<double, 1, 1> MatrixBlockType;
        typedef Dune::BCRSMatrix <MatrixBlockType>        Mat;
        typedef Dune::BlockVector<VectorBlockType>        Vector;
        typedef Dune::MatrixAdapter<Mat,Vector,Vector> Operator;

        template<class O, class S, class C>
        LinearSolverInterface::LinearSolverReport
        solveCG_ILU0(O& A, Vector& x, Vector& b, S& sp, const C& comm, double tolerance, int maxit, int verbosity);

        template<class O, class S, class C>
        LinearSolverInterface::LinearSolverReport
        solveCG_AMG(O& A, Vector& x, Vector& b, S& sp, const C& comm, double tolerance, int maxit, int verbosity,
                    double prolongateFactor, int smoothsteps);

#if defined(HAS_DUNE_FAST_AMG) || DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
       template<class O, class S, class C>
        LinearSolverInterface::LinearSolverReport
        solveKAMG(O& A, Vector& x, Vector& b, S& sp, const C& comm, double tolerance, int maxit, int verbosity,
                  double prolongateFactor, int smoothsteps);

       template<class O, class S, class C>
        LinearSolverInterface::LinearSolverReport
        solveFastAMG(O& A, Vector& x, Vector& b, S& sp, const C& comm, double tolerance, int maxit, int verbosity,
                     double prolongateFactor);
#endif

        template<class O, class S, class C>
        LinearSolverInterface::LinearSolverReport
        solveBiCGStab_ILU0(O& A, Vector& x, Vector& b, S& sp, const C& comm, double tolerance, int maxit, int verbosity);
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
    {}

    LinearSolverInterface::LinearSolverReport
    LinearSolverIstl::solve(const int size,
                            const int nonzeros,
                            const int* ia,
                            const int* ja,
                            const double* sa,
                            const double* rhs,
                            double* solution,
                            const boost::any& comm) const
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

        int maxit = linsolver_max_iterations_;
        if (maxit == 0) {
            maxit = 5000;
        }
#if HAVE_MPI
        if(comm.type()==typeid(ParallelISTLInformation))
        {
            typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
            const ParallelISTLInformation& info = boost::any_cast<const ParallelISTLInformation&>(comm);
            Comm istlComm(info.communicator());
            info.copyValuesTo(istlComm.indexSet(), istlComm.remoteIndices());
            Dune::OverlappingSchwarzOperator<Mat,Vector,Vector, Comm>
                opA(A, istlComm);
            Dune::OverlappingSchwarzScalarProduct<Vector,Comm> sp(istlComm);
            return solveSystem(opA, solution, rhs, sp, istlComm, maxit);
        }
        else
#endif
        {
            (void) comm; // Avoid warning for unused argument if no MPI.
            Dune::SeqScalarProduct<Vector> sp;
            Dune::Amg::SequentialInformation seq_comm;
            Operator opA(A);
            return solveSystem(opA, solution, rhs, sp, seq_comm, maxit);
        }
    }

    template<class O, class S, class C>
    LinearSolverInterface::LinearSolverReport
    LinearSolverIstl::solveSystem (O& opA, double* solution, const double* rhs,
                                   S& sp, const C& comm, int maxit) const
    {
                // System RHS
        Vector b(opA.getmat().N());
        std::copy(rhs, rhs+b.size(), b.begin());
        // System solution
        Vector x(opA.getmat().M());
        x = 0.0;

        if (linsolver_save_system_)
        {
            // Save system to files.
            writeMatrixToMatlab(opA.getmat(), linsolver_save_filename_ + "-mat");
            std::string rhsfile(linsolver_save_filename_ + "-rhs");
            std::ofstream rhsf(rhsfile.c_str());
            rhsf.precision(15);
            rhsf.setf(std::ios::scientific | std::ios::showpos);
            std::copy(b.begin(), b.end(),
                      std::ostream_iterator<VectorBlockType>(rhsf, "\n"));
        }

        LinearSolverReport res;
        switch (linsolver_type_) {
        case CG_ILU0:
            res = solveCG_ILU0(opA, x, b, sp, comm, linsolver_residual_tolerance_, maxit, linsolver_verbosity_);
            break;
        case CG_AMG:
            res = solveCG_AMG(opA, x, b, sp, comm, linsolver_residual_tolerance_, maxit, linsolver_verbosity_,
                              linsolver_prolongate_factor_, linsolver_smooth_steps_);
            break;
        case KAMG:
#if defined(HAS_DUNE_FAST_AMG) || DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
            res = solveKAMG(opA, x, b, sp, comm, linsolver_residual_tolerance_, maxit, linsolver_verbosity_,
                            linsolver_prolongate_factor_, linsolver_smooth_steps_);
#else
            throw std::runtime_error("KAMG not supported with this version of DUNE");
#endif
            break;
        case FastAMG:
#if defined(HAS_DUNE_FAST_AMG) || DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)

#if HAVE_MPI
            if(std::is_same<C,Dune::OwnerOverlapCopyCommunication<int,int> >::value)
            {
                OPM_THROW(std::runtime_error, "Trying to use sequential FastAMG solver for a parallel problem!");
            }
#endif // HAVE_MPI

            res = solveFastAMG(opA, x, b, sp, comm, linsolver_residual_tolerance_, maxit, linsolver_verbosity_,
                               linsolver_prolongate_factor_);
#else
            if(linsolver_verbosity_)
              std::cerr<<"Fast AMG is not available; falling back to CG preconditioned with the normal one"<<std::endl;
            res = solveCG_AMG(opA, x, b, sp, comm, linsolver_residual_tolerance_, maxit, linsolver_verbosity_,
                               linsolver_prolongate_factor_, linsolver_smooth_steps_);
#endif
            break;
        case BiCGStab_ILU0:
            res = solveBiCGStab_ILU0(opA, x, b, sp, comm, linsolver_residual_tolerance_, maxit, linsolver_verbosity_);
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
    template<class P, class O, class C>
    struct SmootherChooser
    {
        typedef P Type;
    };

#if HAVE_MPI
    template<class P, class O>
    struct SmootherChooser<P, O, Dune::OwnerOverlapCopyCommunication<int,int> >
    {
        typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
        typedef Dune::BlockPreconditioner<typename O::domain_type, typename O::range_type,
                                          Comm, P>
        Type;
    };
#endif


    template<class P, class O, class C>
    struct PreconditionerTraits
    {
        typedef typename SmootherChooser<P,O,C>::Type SmootherType;
        typedef std::shared_ptr<SmootherType> PointerType;
    };

    template<class P, class O, class C>
    typename PreconditionerTraits<P,O,C>::PointerType
    makePreconditioner(O& opA, double relax, const C& comm, int iterations=1)
    {
        typedef typename SmootherChooser<P,O,C>::Type SmootherType;
        typedef typename PreconditionerTraits<P,O,C>::PointerType PointerType;
        typename Dune::Amg::SmootherTraits<SmootherType>::Arguments args;
        typename Dune::Amg::ConstructionTraits<SmootherType>::Arguments cargs;
        cargs.setMatrix(opA.getmat());
        args.iterations=iterations;
        args.relaxationFactor=relax;
        cargs.setArgs(args);
        cargs.setComm(comm);
        return PointerType(Dune::Amg::ConstructionTraits<SmootherType>::construct(cargs));
    }

    template<class O, class S, class C>
    LinearSolverInterface::LinearSolverReport
    solveCG_ILU0(O& opA, Vector& x, Vector& b, S& sp, const C& comm, double tolerance, int maxit, int verbosity)
    {

        // Construct preconditioner.
        typedef Dune::SeqILU0<Mat,Vector,Vector> Preconditioner;
        auto precond = makePreconditioner<Preconditioner>(opA, 1.0, comm);

        // Construct linear solver.
        Dune::CGSolver<Vector> linsolve(opA, sp, *precond, tolerance, maxit, verbosity);

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

    template<typename C>
    void setUpCriterion(C& criterion, double linsolver_prolongate_factor,
                        int verbosity, std::size_t linsolver_smooth_steps)
    {
        criterion.setDebugLevel(verbosity);
#if ANISOTROPIC_3D
        criterion.setDefaultValuesAnisotropic(3, 2);
#endif
        criterion.setProlongationDampingFactor(linsolver_prolongate_factor);
        criterion.setNoPreSmoothSteps(linsolver_smooth_steps);
        criterion.setNoPostSmoothSteps(linsolver_smooth_steps);
        criterion.setGamma(1); // V-cycle; this is the default
    }

    template<class O, class S, class C>
    LinearSolverInterface::LinearSolverReport
    solveCG_AMG(O& opA, Vector& x, Vector& b, S& sp, const C& comm, double tolerance, int maxit, int verbosity,
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
        typedef Dune::SeqILU0<Mat,Vector,Vector>        SeqSmoother;
#else
        typedef Dune::SeqSOR<Mat,Vector,Vector>        SeqSmoother;
#endif
        typedef typename SmootherChooser<SeqSmoother, O, C>::Type Smoother;
        typedef Dune::Amg::CoarsenCriterion<CriterionBase> Criterion;
        typedef Dune::Amg::AMG<O,Vector,Smoother,C>   Precond;

        // Construct preconditioner.
        Criterion criterion;
        typename Precond::SmootherArgs smootherArgs;
        setUpCriterion(criterion, linsolver_prolongate_factor, verbosity,
                       linsolver_smooth_steps);
        Precond precond(opA, criterion, smootherArgs, comm);

        // Construct linear solver.
        Dune::CGSolver<Vector> linsolve(opA, sp, precond, tolerance, maxit, verbosity);

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


#if defined(HAS_DUNE_FAST_AMG) || DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
    template<class O, class S, class C>
    LinearSolverInterface::LinearSolverReport
    solveKAMG(O& opA, Vector& x, Vector& b, S& /* sp */, const C& /* comm */, double tolerance, int maxit, int verbosity,
              double linsolver_prolongate_factor, int linsolver_smooth_steps)
    {
        // Solve with AMG solver.
        Dune::MatrixAdapter<typename O::matrix_type,Vector,Vector> sOpA(opA.getmat());

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

        // Construct preconditioner.
        Precond::SmootherArgs smootherArgs;
        Criterion criterion;
        setUpCriterion(criterion, linsolver_prolongate_factor, verbosity,
                       linsolver_smooth_steps);
        Precond precond(sOpA, criterion, smootherArgs);

        // Construct linear solver.
        Dune::GeneralizedPCGSolver<Vector> linsolve(sOpA, precond, tolerance, maxit, verbosity);

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

    template<class O, class S, class C>
    LinearSolverInterface::LinearSolverReport
    solveFastAMG(O& opA, Vector& x, Vector& b, S& /* sp */, const C& /* comm */, double tolerance, int maxit, int verbosity,
                 double linsolver_prolongate_factor)
    {
        // Solve with AMG solver.
        typedef Dune::MatrixAdapter<typename O::matrix_type, Vector, Vector> Operator;
        Operator sOpA(opA.getmat());

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

        // Construct preconditioner.
        Criterion criterion;
        const int smooth_steps = 1;
        setUpCriterion(criterion, linsolver_prolongate_factor, verbosity, smooth_steps);
        Dune::Amg::Parameters parms;
        parms.setDebugLevel(verbosity);
        parms.setNoPreSmoothSteps(smooth_steps);
        parms.setNoPostSmoothSteps(smooth_steps);
        parms.setProlongationDampingFactor(linsolver_prolongate_factor);
        Precond precond(sOpA, criterion, parms);

        // Construct linear solver.
        Dune::GeneralizedPCGSolver<Vector> linsolve(sOpA, precond, tolerance, maxit, verbosity);

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

    template<class O, class S, class C>
    LinearSolverInterface::LinearSolverReport
    solveBiCGStab_ILU0(O& opA, Vector& x, Vector& b, S& sp, const C& comm, double tolerance, int maxit, int verbosity)
    {

        // Construct preconditioner.
        typedef Dune::SeqILU0<Mat,Vector,Vector> Preconditioner;
        auto precond = makePreconditioner<Preconditioner>(opA, 1.0, comm);

        // Construct linear solver.
        Dune::BiCGSTABSolver<Vector> linsolve(opA, sp, *precond, tolerance, maxit, verbosity);

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
