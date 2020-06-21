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


#ifndef OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED
#define OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED

#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>

#include <dune/istl/solver.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <boost/property_tree/ptree.hpp>

namespace Dune
{

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
    FlexibleSolver(const MatrixType& matrix,
                   const boost::property_tree::ptree& prm,
                   const std::function<VectorTypeT()>& weightsCalculator = std::function<VectorTypeT()>());

    /// Create a parallel solver (if Comm is e.g. OwnerOverlapCommunication).
    template <class Comm>
    FlexibleSolver(const MatrixType& matrix,
                   const Comm& comm,
                   const boost::property_tree::ptree& prm,
                   const std::function<VectorTypeT()>& weightsCalculator = std::function<VectorTypeT()>());

    virtual void apply(VectorType& x, VectorType& rhs, Dune::InverseOperatorResult& res) override;

    virtual void apply(VectorType& x, VectorType& rhs, double reduction, Dune::InverseOperatorResult& res) override;

    /// Type of the contained preconditioner.
    using AbstractPrecondType = Dune::PreconditionerWithUpdate<VectorType, VectorType>;

    /// Access the contained preconditioner.
    AbstractPrecondType& preconditioner();

    virtual Dune::SolverCategory::Category category() const override;

private:
    using AbstractOperatorType = Dune::AssembledLinearOperator<MatrixType, VectorType, VectorType>;
    using AbstractScalarProductType = Dune::ScalarProduct<VectorType>;
    using AbstractSolverType = Dune::InverseOperator<VectorType, VectorType>;

    // Machinery for making sequential or parallel operators/preconditioners/scalar products.
    template <class Comm>
    void initOpPrecSp(const MatrixType& matrix, const boost::property_tree::ptree& prm,
                      const std::function<VectorTypeT()> weightsCalculator, const Comm& comm);

    void initOpPrecSp(const MatrixType& matrix, const boost::property_tree::ptree& prm,
                      const std::function<VectorTypeT()> weightsCalculator, const Dune::Amg::SequentialInformation&);

    void initSolver(const boost::property_tree::ptree& prm, bool isMaster);

    // Main initialization routine.
    // Call with Comm == Dune::Amg::SequentialInformation to get a serial solver.
    template <class Comm>
    void init(const MatrixType& matrix,
              const Comm& comm,
              const boost::property_tree::ptree& prm,
              const std::function<VectorTypeT()> weightsCalculator);

    std::shared_ptr<AbstractOperatorType> linearoperator_;
    std::shared_ptr<AbstractPrecondType> preconditioner_;
    std::shared_ptr<AbstractScalarProductType> scalarproduct_;
    std::shared_ptr<AbstractSolverType> linsolver_;
};

} // namespace Dune



#endif // OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED
