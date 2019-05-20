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

#include <opm/simulators/linalg/makePreconditioner.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>

#include <boost/property_tree/ptree.hpp>

namespace Dune
{



template <class MatrixTypeT, class VectorTypeT>
class FlexibleSolver : Dune::InverseOperator<VectorTypeT, VectorTypeT>
{
public:
    using MatrixType = MatrixTypeT;
    using VectorType = VectorTypeT;

    FlexibleSolver(const boost::property_tree::ptree& prm, const MatrixType& matrix)
    {
        makeSolver(prm, matrix);
    }

    virtual void apply(VectorType& x, VectorType& rhs, Dune::InverseOperatorResult& res) override
    {
        linsolver_->apply(x, rhs, res);
    }

    virtual void apply(VectorType& x, VectorType& rhs, double reduction, Dune::InverseOperatorResult& res) override
    {
        linsolver_->apply(x, rhs, reduction, res);
    }

    virtual Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    void makeSolver(const boost::property_tree::ptree& prm, const MatrixType& matrix)
    {
        const double tol = prm.get<double>("tol");
        const int maxiter = prm.get<int>("maxiter");
        linearoperator_.reset(new Dune::MatrixAdapter<MatrixType, VectorType, VectorType>(matrix));
        preconditioner_ = Dune::makePreconditioner<MatrixType, VectorType>(*linearoperator_, prm);
        int verbosity = prm.get<int>("verbosity");
        std::string solver_type = prm.get<std::string>("solver");
        if (solver_type == "bicgstab") {
            linsolver_.reset(new Dune::BiCGSTABSolver<VectorType>(*linearoperator_,
                                                                  *preconditioner_,
                                                                  tol, // desired residual reduction factor
                                                                  maxiter, // maximum number of iterations
                                                                  verbosity));
        } else if (solver_type == "loopsolver") {
            linsolver_.reset(new Dune::LoopSolver<VectorType>(*linearoperator_,
                                                              *preconditioner_,
                                                              tol, // desired residual reduction factor
                                                              maxiter, // maximum number of iterations
                                                              verbosity));
        } else if (solver_type == "gmres") {
            int restart = prm.get<int>("restart");
            linsolver_.reset(new Dune::RestartedGMResSolver<VectorType>(*linearoperator_,
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
            std::string msg("Solver not known ");
            msg += solver_type;
            throw std::runtime_error(msg);
        }
    }

    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> preconditioner_;
    std::shared_ptr<Dune::MatrixAdapter<MatrixType, VectorType, VectorType>> linearoperator_;
    std::shared_ptr<Dune::InverseOperator<VectorType, VectorType>> linsolver_;
};

} // namespace Dune



#endif // OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED
