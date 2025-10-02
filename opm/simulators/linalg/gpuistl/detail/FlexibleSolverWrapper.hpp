/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_FLEXIBLESOLVERWRAPPER_HEADER_INCLUDED
#define OPM_FLEXIBLESOLVERWRAPPER_HEADER_INCLUDED

#include <functional>
#include <memory>
#include <type_traits>

#include <dune/common/parallel/communication.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solver.hh>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/gpuistl/GpuOwnerOverlapCopy.hpp>

namespace Opm::gpuistl::detail
{

/**
 * @brief FlexibleSolverWrapper is compilational trick to reduce compile time overhead
 *
 * Essentially this class acts as a direct interface towards the FlexibleSolver class, but
 * we want to have a clear wall between the ISTLSolverGPUISTL and the FlexibleSolver to avoid
 * compilational overhead (the former is based on TypeTag the latter on operators).
 *
 * This wrapper handles the mapping from matrices to concrete operator types.
 */
template <class Matrix, class Vector, class Comm>
class FlexibleSolverWrapper
{
public:
    using AbstractSolverType = Dune::InverseOperator<Vector, Vector>;
    using AbstractPreconditionerType = Dune::PreconditionerWithUpdate<Vector, Vector>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;

    using AbstractOperatorPtrType = std::unique_ptr<AbstractOperatorType>;
    using AbstractSolverPtrType = std::unique_ptr<AbstractSolverType>;
    using GpuCommunicationType = std::conditional_t<
        std::is_same_v<Comm, Dune::Communication<int>>,
        int, // Dummy type for serial case
        GpuOwnerOverlapCopy<typename Matrix::field_type, Comm>>;

    FlexibleSolverWrapper(const Matrix& matrix,
                          bool parallel,
                          const PropertyTree& prm,
                          std::size_t pressureIndex,
                          const std::function<Vector()>& weightCalculator,
                          bool forceSerial,
                          const Comm* comm);

    void update();

    void apply(Vector& x, Vector& y, Dune::InverseOperatorResult& result);

private:
    AbstractOperatorPtrType m_operator;
    AbstractSolverPtrType m_solver;

    AbstractPreconditionerType& m_preconditioner;
    std::shared_ptr<GpuCommunicationType> m_gpuCommunication;

    FlexibleSolverWrapper(
        std::tuple<AbstractOperatorPtrType, AbstractSolverPtrType, std::reference_wrapper<AbstractPreconditionerType>, std::shared_ptr<GpuCommunicationType>>&&
            solverTuple);
};
} // namespace Opm::gpuistl::detail

#endif // OPM_FLEXIBLESOLVERWRAPPER_HEADER_INCLUDED
