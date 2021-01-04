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

#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/ilufirstelement.hh>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/PreconditionerFactory.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>

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
        linsolver_->apply(x, rhs, res);
    }

    template <class Operator>
    void
    FlexibleSolver<Operator>::
    apply(VectorType& x, VectorType& rhs, double reduction, Dune::InverseOperatorResult& res)
    {
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
        preconditioner_ = Opm::PreconditionerFactory<Operator>::create(op,
                                                                       child ? *child : Opm::PropertyTree(),
                                                                       weightsCalculator,
                                                                       pressureIndex);
        scalarproduct_ = std::make_shared<Dune::SeqScalarProduct<VectorType>>();
    }

    template <class Operator>
    void
    FlexibleSolver<Operator>::
    initSolver(const Opm::PropertyTree& prm, const bool is_iorank)
    {
        const double tol = prm.get<double>("tol", 1e-2);
        const int maxiter = prm.get<int>("maxiter", 200);
        const int verbosity = is_iorank ? prm.get<int>("verbosity", 0) : 0;
        const std::string solver_type = prm.get<std::string>("solver", "bicgstab");
        if (solver_type == "bicgstab") {
            linsolver_.reset(new Dune::BiCGSTABSolver<VectorType>(*linearoperator_for_solver_,
                                                                  *scalarproduct_,
                                                                  *preconditioner_,
                                                                  tol, // desired residual reduction factor
                                                                  maxiter, // maximum number of iterations
                                                                  verbosity));
        } else if (solver_type == "loopsolver") {
            linsolver_.reset(new Dune::LoopSolver<VectorType>(*linearoperator_for_solver_,
                                                              *scalarproduct_,
                                                              *preconditioner_,
                                                              tol, // desired residual reduction factor
                                                              maxiter, // maximum number of iterations
                                                              verbosity));
        } else if (solver_type == "gmres") {
            int restart = prm.get<int>("restart", 15);
            linsolver_.reset(new Dune::RestartedGMResSolver<VectorType>(*linearoperator_for_solver_,
                                                                        *scalarproduct_,
                                                                        *preconditioner_,
                                                                        tol,
                                                                        restart, // desired residual reduction factor
                                                                        maxiter, // maximum number of iterations
                                                                        verbosity));
#if HAVE_SUITESPARSE_UMFPACK
        } else if (solver_type == "umfpack") {
            bool dummy = false;
            using MatrixType = std::remove_const_t<std::remove_reference_t<decltype(linearoperator_for_solver_->getmat())>>;
            linsolver_.reset(new Dune::UMFPack<MatrixType>(linearoperator_for_solver_->getmat(), verbosity, dummy));
#endif
        } else {
            OPM_THROW(std::invalid_argument, "Properties: Solver " << solver_type << " not known.");
        }
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
        initSolver(prm, comm.communicator().rank() == 0);
    }

} // namespace Dune


// Macros to simplify explicit instantiation of FlexibleSolver for various block sizes.
/*
template <int N>
using BV = Dune::BlockVector<Dune::FieldVector<double, N>>;
template <int N>
using BM = Dune::BCRSMatrix<Dune::FieldMatrix<double, N, N>>;
template <int N>
using OBM = Dune::BCRSMatrix<Opm::MatrixBlock<double, N, N>>;

#if HAVE_MPI

using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;

// Note: we must instantiate the constructor that is a template.
// This is only needed in the parallel case, since otherwise the Comm type is
// not a template argument but always SequentialInformation.

#define INSTANTIATE_FLEXIBLESOLVER(N)                  \
template class Dune::FlexibleSolver<BM<N>, BV<N>>;     \
template class Dune::FlexibleSolver<OBM<N>, BV<N>>;    \
template Dune::FlexibleSolver<BM<N>, BV<N>>::FlexibleSolver(AbstractOperatorType& op,                         \
                                                            const Comm& comm,                                 \
                                                            const Opm::PropertyTree& prm,                     \
                                                            const std::function<BV<N>()>& weightsCalculator, \
std::size_t);                                                                 \
template Dune::FlexibleSolver<OBM<N>, BV<N>>::FlexibleSolver(AbstractOperatorType& op,                        \
                                                             const Comm& comm,                                \
                                                             const Opm::PropertyTree& prm,                    \
                                                             const std::function<BV<N>()>& weightsCalculator, \
std::size_t);

#else // HAVE_MPI

#define INSTANTIATE_FLEXIBLESOLVER(N)                  \
template class Dune::FlexibleSolver<BM<N>, BV<N>>;     \
template class Dune::FlexibleSolver<OBM<N>, BV<N>>;

#endif // HAVE_MPI
*/

#define INSTANTIATE_FLEXIBLESOLVER(N) do(){}

#endif // OPM_FLEXIBLE_SOLVER_IMPL_HEADER_INCLUDED
