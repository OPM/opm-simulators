/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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


#ifndef OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED
#define OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED

#include <opm/simulators/linalg/PreconditionerFactory.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <boost/property_tree/ptree.hpp>

namespace Dune
{

template<class C>
struct IsComm : std::false_type
{};

template<>
struct IsComm<Dune::Amg::SequentialInformation> : std::true_type
{};

template<class Index>
struct IsComm<Dune::OwnerOverlapCopyCommunication<Index>> : std::true_type
{};

/// A solver class that encapsulates all needed objects for a linear solver
/// (operator, scalar product, iterative solver and preconditioner) and sets
/// them up based on runtime parameters, using the PreconditionerFactory for
/// setting up preconditioners.
template <class MatrixTypeT, class VectorTypeT>
class FlexibleSolver : public Dune::InverseOperator<VectorTypeT, VectorTypeT>
{
public:
    using MatrixType = MatrixTypeT;
    using VectorType = VectorTypeT;

    /// Create a sequential solver.
    FlexibleSolver(const boost::property_tree::ptree& prm, const MatrixType& matrix,
                   const std::function<VectorTypeT()>& weightsCalculator = std::function<VectorTypeT()>())
    {
        init(prm, matrix, weightsCalculator, Dune::Amg::SequentialInformation());
    }

    /// Create a parallel solver (if Comm is e.g. OwnerOverlapCommunication).
    template <class Comm>
    FlexibleSolver(const boost::property_tree::ptree& prm,
                   const MatrixType& matrix,
                   const typename std::enable_if<IsComm<Comm>::value, Comm>::type& comm)
    {
        init(prm, matrix, std::function<VectorTypeT()>(), comm);
    }

    /// Create a parallel solver (if Comm is e.g. OwnerOverlapCommunication).
    template <class Comm>
    FlexibleSolver(const boost::property_tree::ptree& prm, const MatrixType& matrix,
                   const std::function<VectorTypeT()>& weightsCalculator, const Comm& comm)
    {
        init(prm, matrix, weightsCalculator, comm);
    }

    virtual void apply(VectorType& x, VectorType& rhs, Dune::InverseOperatorResult& res) override
    {
        linsolver_->apply(x, rhs, res);
    }

    virtual void apply(VectorType& x, VectorType& rhs, double reduction, Dune::InverseOperatorResult& res) override
    {
        linsolver_->apply(x, rhs, reduction, res);
    }

    /// Type of the contained preconditioner.
    using AbstractPrecondType = Dune::PreconditionerWithUpdate<VectorType, VectorType>;

    /// Access the contained preconditioner.
    AbstractPrecondType& preconditioner()
    {
        return *preconditioner_;
    }

    virtual Dune::SolverCategory::Category category() const override
    {
        return linearoperator_->category();
    }

private:
    using AbstractOperatorType = Dune::AssembledLinearOperator<MatrixType, VectorType, VectorType>;
    using AbstractScalarProductType = Dune::ScalarProduct<VectorType>;
    using AbstractSolverType = Dune::InverseOperator<VectorType, VectorType>;

    // Machinery for making sequential or parallel operators/preconditioners/scalar products.
    template <class Comm>
    void initOpPrecSp(const MatrixType& matrix, const boost::property_tree::ptree& prm,
                      const std::function<VectorTypeT()> weightsCalculator, const Comm& comm)
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

    void initOpPrecSp(const MatrixType& matrix, const boost::property_tree::ptree& prm,
                      const std::function<VectorTypeT()> weightsCalculator, const Dune::Amg::SequentialInformation&)
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

    void initSolver(const boost::property_tree::ptree& prm)
    {
        const double tol = prm.get<double>("tol", 1e-2);
        const int maxiter = prm.get<int>("maxiter", 200);
        const int verbosity = prm.get<int>("verbosity", 0);
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
    template <class Comm>
    void init(const boost::property_tree::ptree& prm, const MatrixType& matrix,
              const std::function<VectorTypeT()> weightsCalculator, const Comm& comm)
    {
        initOpPrecSp(matrix, prm, weightsCalculator, comm);
        initSolver(prm);
    }

    std::shared_ptr<AbstractOperatorType> linearoperator_;
    std::shared_ptr<AbstractPrecondType> preconditioner_;
    std::shared_ptr<AbstractScalarProductType> scalarproduct_;
    std::shared_ptr<AbstractSolverType> linsolver_;
};

} // namespace Dune



#endif // OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED
