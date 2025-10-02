/*
  Copyright 2019, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2020 Equinor.

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

#ifndef OPM_FLEXIBLE_SOLVER_IMPL_HEADER_INCLUDED
#define OPM_FLEXIBLE_SOLVER_IMPL_HEADER_INCLUDED

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/ilufirstelement.hh>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/WellOperators.hpp>
#include <opm/simulators/linalg/PreconditionerFactoryGPUIncludeWrapper.hpp>
#include <opm/simulators/linalg/is_gpu_operator.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <type_traits>

#if HAVE_CUDA
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/SolverAdapter.hpp>
#else
#include <opm/simulators/linalg/gpuistl/SolverAdapter.hpp>
#endif
#endif

namespace Dune
{
    /// Create a sequential solver.
    template <class Operator>
    FlexibleSolver<Operator>::
    FlexibleSolver(Operator& op,
                   const Opm::PropertyTree& prm,
                   const std::function<VectorType()>& weightsCalculator,
                   std::size_t pressureIndex)
    {
        init(op, Dune::Amg::SequentialInformation(), prm, weightsCalculator,
             pressureIndex);
    }

    /// Create a parallel solver (if Comm is e.g. OwnerOverlapCommunication).
    template <class Operator>
    template <class Comm>
    FlexibleSolver<Operator>::
    FlexibleSolver(Operator& op,
                   const Comm& comm,
                   const Opm::PropertyTree& prm,
                   const std::function<VectorType()>& weightsCalculator,
                   std::size_t pressureIndex)
    {
        init(op, comm, prm, weightsCalculator, pressureIndex);
    }

    template <class Operator>
    void
    FlexibleSolver<Operator>::
    apply(VectorType& x, VectorType& rhs, Dune::InverseOperatorResult& res)
    {
        if (direct_solver_) {
            recreateDirectSolver();
        }
        linsolver_->apply(x, rhs, res);
    }

    template <class Operator>
    void
    FlexibleSolver<Operator>::
    apply(VectorType& x, VectorType& rhs, double reduction, Dune::InverseOperatorResult& res)
    {
        if (direct_solver_) {
            recreateDirectSolver();
        }
        linsolver_->apply(x, rhs, reduction, res);
    }

    /// Access the contained preconditioner.
    template <class Operator>
    auto
    FlexibleSolver<Operator>::
    preconditioner() -> AbstractPrecondType&
    {
        return *preconditioner_;
    }

    template <class Operator>
    Dune::SolverCategory::Category
    FlexibleSolver<Operator>::
    category() const
    {
        return linearoperator_for_solver_->category();
    }

    // Machinery for making sequential or parallel operators/preconditioners/scalar products.
    template <class Operator>
    template <class Comm>
    void
    FlexibleSolver<Operator>::
    initOpPrecSp(Operator& op,
                 const Opm::PropertyTree& prm,
                 const std::function<VectorType()> weightsCalculator,
                 const Comm& comm,
                 std::size_t pressureIndex)
    {
        // Parallel case.
        linearoperator_for_solver_ = &op;
        auto child = prm.get_child_optional("preconditioner");
        preconditioner_ = Opm::PreconditionerFactory<Operator, Comm>::create(op,
                                                                             child ? *child : Opm::PropertyTree(),
                                                                             weightsCalculator,
                                                                             comm,
                                                                             pressureIndex);
        scalarproduct_ = Dune::createScalarProduct<VectorType, Comm>(comm, op.category());
    }

    template <class Operator>
    void
    FlexibleSolver<Operator>::
    initOpPrecSp(Operator& op,
                 const Opm::PropertyTree& prm,
                 const std::function<VectorType()> weightsCalculator,
                 const Dune::Amg::SequentialInformation&,
                 std::size_t pressureIndex)
    {
        // Sequential case.
        linearoperator_for_solver_ = &op;
        auto child = prm.get_child_optional("preconditioner");
        preconditioner_ = Opm::PreconditionerFactory<Operator,Dune::Amg::SequentialInformation>::create(op,
                                                                       child ? *child : Opm::PropertyTree(),
                                                                       weightsCalculator,
                                                                       pressureIndex);
        scalarproduct_ = std::make_shared<Dune::SeqScalarProduct<VectorType>>();
    }


    template <class Operator>
    template <class Comm>
    void
    FlexibleSolver<Operator>::
    initSolver(const Opm::PropertyTree& prm, const Comm& comm)
    {
        const bool is_iorank = comm.communicator().rank() == 0;
        const double tol = prm.get<double>("tol", 1e-2);
        const int maxiter = prm.get<int>("maxiter", 200);
        const int verbosity = is_iorank ? prm.get<int>("verbosity", 0) : 0;
        const std::string solver_type = prm.get<std::string>("solver", "bicgstab");


        // make sure it is nullptr at the start (used for error checking in the end).
        // while the linSolver_ is initalized as a nullptr, we want to make sure it is reset here,
        // simply because we will check if it is at the end of this function and need to keep this invariant
        // (that it is nullptr at the start of this function).
        linsolver_.reset();
        if (solver_type == "bicgstab") {
            linsolver_ = std::make_shared<Dune::BiCGSTABSolver<VectorType>>(*linearoperator_for_solver_,
                                                                            *scalarproduct_,
                                                                            *preconditioner_,
                                                                            tol, // desired residual reduction factor
                                                                            maxiter, // maximum number of iterations
                                                                            verbosity);
        } else if (solver_type == "loopsolver") {
            linsolver_ = std::make_shared<Dune::LoopSolver<VectorType>>(*linearoperator_for_solver_,
                                                                        *scalarproduct_,
                                                                        *preconditioner_,
                                                                        tol, // desired residual reduction factor
                                                                        maxiter, // maximum number of iterations
                                                                        verbosity);
        } else if (solver_type == "gmres") {
            int restart = prm.get<int>("restart", 15);
            linsolver_ = std::make_shared<Dune::RestartedGMResSolver<VectorType>>(*linearoperator_for_solver_,
                                                                                  *scalarproduct_,
                                                                                  *preconditioner_,
                                                                                  tol,// desired residual reduction factor
                                                                                  restart,
                                                                                  maxiter, // maximum number of iterations
                                                                                  verbosity);

        } else {
            if constexpr (!Opm::is_gpu_operator_v<Operator>) {
                if (solver_type == "flexgmres") {
                    int restart = prm.get<int>("restart", 15);
                    linsolver_ = std::make_shared<Dune::RestartedFlexibleGMResSolver<VectorType>>(*linearoperator_for_solver_,
                                                                                            *scalarproduct_,
                                                                                            *preconditioner_,
                                                                                            tol,// desired residual reduction factor
                                                                                            restart,
                                                                                            maxiter, // maximum number of iterations
                                                                                            verbosity);
#if HAVE_SUITESPARSE_UMFPACK
                } else if (solver_type == "umfpack") {
                    if constexpr (std::is_same_v<typename VectorType::field_type,float>) {
                        OPM_THROW(std::invalid_argument, "UMFPack cannot be used with floats");
                    } else {
                        using MatrixType = std::remove_const_t<std::remove_reference_t<decltype(linearoperator_for_solver_->getmat())>>;
                        linsolver_ = std::make_shared<Dune::UMFPack<MatrixType>>(linearoperator_for_solver_->getmat(), verbosity, false);
                        direct_solver_ = true;
                    }
#endif
#if HAVE_CUDA
                } else if (solver_type == "gpubicgstab") {
                    linsolver_.reset(new Opm::gpuistl::SolverAdapter<Operator, Dune::BiCGSTABSolver, VectorType>(
                        *linearoperator_for_solver_,
                        *scalarproduct_,
                        preconditioner_,
                        tol, // desired residual reduction factor
                        maxiter, // maximum number of iterations
                        verbosity,
                        comm));
        #endif
               }
            }
        }
        if (!linsolver_) {
            OPM_THROW(std::invalid_argument,
                      "Properties: Solver " + solver_type + " not known.");
        }
    }


    // For now, the only direct solver we support is UMFPACK from SuiteSparse.
    // When the matrix is updated (keeping sparsity pattern) it is possible to
    // exploit separation of symbolic and numeric factorization, but we do not
    // do so at this point. For complete generality, the solver abstract class
    // Dune::InverseOperator<> should be extended with an update() function.
    template <class Operator>
    void
    FlexibleSolver<Operator>::
    recreateDirectSolver()
    {
#if HAVE_SUITESPARSE_UMFPACK
        if constexpr (!Opm::is_gpu_operator_v<Operator>) {
            if constexpr (std::is_same_v<typename VectorType::field_type, float>) {
                OPM_THROW(std::invalid_argument, "UMFPack cannot be used with floats");
            } else {
                using MatrixType = std::remove_const_t<std::remove_reference_t<decltype(linearoperator_for_solver_->getmat())>>;
                linsolver_ = std::make_shared<Dune::UMFPack<MatrixType>>(linearoperator_for_solver_->getmat(), 0, false);
            }
        }
#else
        OPM_THROW(std::logic_error, "Direct solver specified, but the FlexibleSolver class was not compiled with SuiteSparse support.");
#endif
    }


    // Main initialization routine.
    // Call with Comm == Dune::Amg::SequentialInformation to get a serial solver.
    template <class Operator>
    template <class Comm>
    void
    FlexibleSolver<Operator>::
    init(Operator& op,
         const Comm& comm,
         const Opm::PropertyTree& prm,
         const std::function<VectorType()> weightsCalculator,
         std::size_t pressureIndex)
    {
        initOpPrecSp(op, prm, weightsCalculator, comm, pressureIndex);
        initSolver(prm, comm);
    }

} // namespace Dune


// Macros to simplify explicit instantiation of FlexibleSolver for various block sizes.

// Vectors and matrices.
template<class Scalar, int N>
using BV = Dune::BlockVector<Dune::FieldVector<Scalar, N>>;
template<class Scalar, int N>
using OBM = Dune::BCRSMatrix<Opm::MatrixBlock<Scalar, N, N>>;

// Sequential operators.
template<class Scalar, int N>
using SeqOpM = Dune::MatrixAdapter<OBM<Scalar,N>, BV<Scalar,N>, BV<Scalar,N>>;
template<class Scalar, int N>
using SeqOpW = Opm::WellModelMatrixAdapter<OBM<Scalar,N>, BV<Scalar,N>, BV<Scalar,N>>;

#if HAVE_MPI

// Parallel communicator and operators.
using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;
template<class Scalar, int N>
using ParOpM = Opm::GhostLastMatrixAdapter<OBM<Scalar,N>, BV<Scalar,N>, BV<Scalar,N>, Comm>;
template<class Scalar, int N>
using ParOpW = Opm::WellModelGhostLastMatrixAdapter<OBM<Scalar,N>, BV<Scalar,N>, BV<Scalar,N>, true>;
template<class Scalar, int N>
using ParOpD = Dune::OverlappingSchwarzOperator<OBM<Scalar,N>, BV<Scalar,N>, BV<Scalar,N>, Comm>;

// Note: we must instantiate the constructor that is a template.
// This is only needed in the parallel case, since otherwise the Comm type is
// not a template argument but always SequentialInformation.

#define INSTANTIATE_FLEXIBLESOLVER_OP(...)                                                          \
    template class Dune::FlexibleSolver<__VA_ARGS__>;                                               \
    template Dune::FlexibleSolver<__VA_ARGS__>::                                                    \
        FlexibleSolver(__VA_ARGS__& op,                                                             \
                       const Comm& comm,                                                            \
                       const Opm::PropertyTree& prm,                                                \
                       const std::function<typename __VA_ARGS__::domain_type()>& weightsCalculator, \
                       std::size_t pressureIndex);

#define INSTANTIATE_FLEXIBLESOLVER(T,N)         \
    INSTANTIATE_FLEXIBLESOLVER_OP(SeqOpM<T,N>); \
    INSTANTIATE_FLEXIBLESOLVER_OP(SeqOpW<T,N>); \
    INSTANTIATE_FLEXIBLESOLVER_OP(ParOpM<T,N>); \
    INSTANTIATE_FLEXIBLESOLVER_OP(ParOpW<T,N>); \
    INSTANTIATE_FLEXIBLESOLVER_OP(ParOpD<T,N>);

#else // HAVE_MPI

#define INSTANTIATE_FLEXIBLESOLVER_OP(...) \
    template class Dune::FlexibleSolver<__VA_ARGS__>;

#define INSTANTIATE_FLEXIBLESOLVER(T,N)         \
    INSTANTIATE_FLEXIBLESOLVER_OP(SeqOpM<T,N>); \
    INSTANTIATE_FLEXIBLESOLVER_OP(SeqOpW<T,N>);

#endif // HAVE_MPI

#endif // OPM_FLEXIBLE_SOLVER_IMPL_HEADER_INCLUDED
