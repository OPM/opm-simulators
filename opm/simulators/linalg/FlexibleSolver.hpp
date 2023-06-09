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

namespace Opm
{
class PropertyTree;
}

namespace Dune
{

/// A solver class that encapsulates all needed objects for a linear solver
/// (operator, scalar product, iterative solver and preconditioner) and sets
/// them up based on runtime parameters, using the PreconditionerFactory for
/// setting up preconditioners.
template <class Operator>
class FlexibleSolver : public Dune::InverseOperator<typename Operator::domain_type,
                                                    typename Operator::range_type>
{
public:
    using VectorType = typename Operator::domain_type; // Assuming symmetry: domain == range

    /// Base class type of the contained preconditioner.
    using AbstractPrecondType = Dune::PreconditionerWithUpdate<VectorType, VectorType>;

    /// Create a sequential solver.
    FlexibleSolver(Operator& op,
                   const Opm::PropertyTree& prm,
                   const std::function<VectorType()>& weightsCalculator,
                   std::size_t pressureIndex);

    /// Create a parallel solver (if Comm is e.g. OwnerOverlapCommunication).
    template <class Comm>
    FlexibleSolver(Operator& op,
                   const Comm& comm,
                   const Opm::PropertyTree& prm,
                   const std::function<VectorType()>& weightsCalculator,
                   std::size_t pressureIndex);

    virtual void apply(VectorType& x, VectorType& rhs, Dune::InverseOperatorResult& res) override;

    virtual void apply(VectorType& x, VectorType& rhs, double reduction, Dune::InverseOperatorResult& res) override;

    /// Access the contained preconditioner.
    AbstractPrecondType& preconditioner();

    virtual Dune::SolverCategory::Category category() const override;

private:
    using AbstractScalarProductType = Dune::ScalarProduct<VectorType>;
    using AbstractSolverType = Dune::InverseOperator<VectorType, VectorType>;

    // Machinery for making sequential or parallel operators/preconditioners/scalar products.
    template <class Comm>
    void initOpPrecSp(Operator& op, const Opm::PropertyTree& prm,
                      const std::function<VectorType()> weightsCalculator, const Comm& comm,
                      std::size_t pressureIndex);

    void initOpPrecSp(Operator& op, const Opm::PropertyTree& prm,
                      const std::function<VectorType()> weightsCalculator, const Dune::Amg::SequentialInformation&,
                      std::size_t pressureIndex);

    void initSolver(const Opm::PropertyTree& prm, const bool is_iorank);

    void recreateDirectSolver();

    // Main initialization routine.
    // Call with Comm == Dune::Amg::SequentialInformation to get a serial solver.
    template <class Comm>
    void init(Operator& op,
              const Comm& comm,
              const Opm::PropertyTree& prm,
              const std::function<VectorType()> weightsCalculator,
              std::size_t pressureIndex);

    Operator* linearoperator_for_solver_;
    std::shared_ptr<AbstractPrecondType> preconditioner_;
    std::shared_ptr<AbstractScalarProductType> scalarproduct_;
    std::shared_ptr<AbstractSolverType> linsolver_;
    bool direct_solver_ = false;
};

} // namespace Dune

//#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>


#endif // OPM_FLEXIBLE_SOLVER_HEADER_INCLUDED
