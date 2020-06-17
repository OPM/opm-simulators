/*
  Copyright 2019, 2020 SINTEF Digital, Mathematics and Cybernetics.

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

#include "config.h"

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <boost/property_tree/ptree.hpp>

namespace Dune
{
    /// Create a sequential solver.
    template <class MatrixType, class VectorType>
    FlexibleSolver<MatrixType, VectorType>::
    FlexibleSolver(const MatrixType& matrix,
                   const boost::property_tree::ptree& prm,
                   const std::function<VectorType()>& weightsCalculator)
    {
        init(matrix, Dune::Amg::SequentialInformation(), prm, weightsCalculator);
    }

    /// Create a parallel solver (if Comm is e.g. OwnerOverlapCommunication).
    template <class MatrixType, class VectorType>
    template <class Comm>
    FlexibleSolver<MatrixType, VectorType>::
    FlexibleSolver(const MatrixType& matrix,
                   const Comm& comm,
                   const boost::property_tree::ptree& prm,
                   const std::function<VectorType()>& weightsCalculator)
    {
        init(matrix, comm, prm, weightsCalculator);
    }

    template <class MatrixType, class VectorType>
    void
    FlexibleSolver<MatrixType, VectorType>::
    apply(VectorType& x, VectorType& rhs, Dune::InverseOperatorResult& res)
    {
        linsolver_->apply(x, rhs, res);
    }

    template <class MatrixType, class VectorType>
    void
    FlexibleSolver<MatrixType, VectorType>::
    apply(VectorType& x, VectorType& rhs, double reduction, Dune::InverseOperatorResult& res)
    {
        linsolver_->apply(x, rhs, reduction, res);
    }

    /// Access the contained preconditioner.
    template <class MatrixType, class VectorType>
    auto
    FlexibleSolver<MatrixType, VectorType>::
    preconditioner() -> AbstractPrecondType&
    {
        return *preconditioner_;
    }

    template <class MatrixType, class VectorType>
    Dune::SolverCategory::Category
    FlexibleSolver<MatrixType, VectorType>::
    category() const
    {
        return linearoperator_->category();
    }

    // Machinery for making sequential or parallel operators/preconditioners/scalar products.
    template <class MatrixType, class VectorType>
    template <class Comm>
    void
    FlexibleSolver<MatrixType, VectorType>::
    initOpPrecSp(const MatrixType& matrix, const boost::property_tree::ptree& prm,
                 const std::function<VectorType()> weightsCalculator, const Comm& comm)
    {
        // Parallel case.
        using ParOperatorType = Dune::OverlappingSchwarzOperator<MatrixType, VectorType, VectorType, Comm>;
        using pt = const boost::property_tree::ptree;
        auto linop = std::make_shared<ParOperatorType>(matrix, comm);
        linearoperator_ = linop;
        auto child = prm.get_child_optional("preconditioner");
        preconditioner_
            = Opm::PreconditionerFactory<ParOperatorType, Comm>::create(*linop, child? *child : pt(),
                                                                        weightsCalculator, comm);
        scalarproduct_ = Dune::createScalarProduct<VectorType, Comm>(comm, linearoperator_->category());
    }

    template <class MatrixType, class VectorType>
    void
    FlexibleSolver<MatrixType, VectorType>::
    initOpPrecSp(const MatrixType& matrix, const boost::property_tree::ptree& prm,
                 const std::function<VectorType()> weightsCalculator, const Dune::Amg::SequentialInformation&)
    {
        // Sequential case.
        using SeqOperatorType = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;
        using pt = const boost::property_tree::ptree;
        auto linop = std::make_shared<SeqOperatorType>(matrix);
        linearoperator_ = linop;
        auto child = prm.get_child_optional("preconditioner");
        preconditioner_ = Opm::PreconditionerFactory<SeqOperatorType>::create(*linop, child? *child : pt(),
                                                                              weightsCalculator);
        scalarproduct_ = std::make_shared<Dune::SeqScalarProduct<VectorType>>();
    }

    template <class MatrixType, class VectorType>
    void
    FlexibleSolver<MatrixType, VectorType>::
    initSolver(const boost::property_tree::ptree& prm, bool isMaster)
    {
        const double tol = prm.get<double>("tol", 1e-2);
        const int maxiter = prm.get<int>("maxiter", 200);
        const int verbosity = isMaster? prm.get<int>("verbosity", 0) : 0;
        const std::string solver_type = prm.get<std::string>("solver", "bicgstab");
        if (solver_type == "bicgstab") {
            linsolver_.reset(new Dune::BiCGSTABSolver<VectorType>(*linearoperator_,
                                                                  *scalarproduct_,
                                                                  *preconditioner_,
                                                                  tol, // desired residual reduction factor
                                                                  maxiter, // maximum number of iterations
                                                                  verbosity));
        } else if (solver_type == "loopsolver") {
            linsolver_.reset(new Dune::LoopSolver<VectorType>(*linearoperator_,
                                                              *scalarproduct_,
                                                              *preconditioner_,
                                                              tol, // desired residual reduction factor
                                                              maxiter, // maximum number of iterations
                                                              verbosity));
        } else if (solver_type == "gmres") {
            int restart = prm.get<int>("restart", 15);
            linsolver_.reset(new Dune::RestartedGMResSolver<VectorType>(*linearoperator_,
                                                                        *scalarproduct_,
                                                                        *preconditioner_,
                                                                        tol,
                                                                        restart, // desired residual reduction factor
                                                                        maxiter, // maximum number of iterations
                                                                        verbosity));
#if HAVE_SUITESPARSE_UMFPACK
        } else if (solver_type == "umfpack") {
            bool dummy = false;
            linsolver_.reset(new Dune::UMFPack<MatrixType>(linearoperator_->getmat(), verbosity, dummy));
#endif
        } else {
            OPM_THROW(std::invalid_argument, "Properties: Solver " << solver_type << " not known.");
        }
    }


    // Main initialization routine.
    // Call with Comm == Dune::Amg::SequentialInformation to get a serial solver.
    template <class MatrixType, class VectorType>
    template <class Comm>
    void
    FlexibleSolver<MatrixType, VectorType>::
    init(const MatrixType& matrix,
         const Comm& comm,
         const boost::property_tree::ptree& prm,
         const std::function<VectorType()> weightsCalculator)
    {
        initOpPrecSp(matrix, prm, weightsCalculator, comm);
        initSolver(prm, comm.communicator().rank()==0);
    }

} // namespace Dune


// Explicit instantiations of FlexibleSolver

template <int N>
using BV = Dune::BlockVector<Dune::FieldVector<double, N>>;
template <int N>
using BM = Dune::BCRSMatrix<Dune::FieldMatrix<double, N, N>>;
template <int N>
using OBM = Dune::BCRSMatrix<Opm::MatrixBlock<double, N, N>>;

// Variants using Dune::FieldMatrix blocks.
template class Dune::FlexibleSolver<BM<1>, BV<1>>;
template class Dune::FlexibleSolver<BM<2>, BV<2>>;
template class Dune::FlexibleSolver<BM<3>, BV<3>>;
template class Dune::FlexibleSolver<BM<4>, BV<4>>;

// Variants using Opm::MatrixBlock blocks.
template class Dune::FlexibleSolver<OBM<1>, BV<1>>;
template class Dune::FlexibleSolver<OBM<2>, BV<2>>;
template class Dune::FlexibleSolver<OBM<3>, BV<3>>;
template class Dune::FlexibleSolver<OBM<4>, BV<4>>;


using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;

template Dune::FlexibleSolver<OBM<1>, BV<1>>::FlexibleSolver(const MatrixType& matrix,
                                                             const Comm& comm,
                                                             const boost::property_tree::ptree& prm,
                                                             const std::function<BV<1>()>& weightsCalculator);

template Dune::FlexibleSolver<OBM<2>, BV<2>>::FlexibleSolver(const MatrixType& matrix,
                                                             const Comm& comm,
                                                             const boost::property_tree::ptree& prm,
                                                             const std::function<BV<2>()>& weightsCalculator);

template Dune::FlexibleSolver<OBM<3>, BV<3>>::FlexibleSolver(const MatrixType& matrix,
                                                             const Comm& comm,
                                                             const boost::property_tree::ptree& prm,
                                                             const std::function<BV<3>()>& weightsCalculator);

template Dune::FlexibleSolver<OBM<4>, BV<4>>::FlexibleSolver(const MatrixType& matrix,
                                                             const Comm& comm,
                                                             const boost::property_tree::ptree& prm,
                                                             const std::function<BV<4>()>& weightsCalculator);
