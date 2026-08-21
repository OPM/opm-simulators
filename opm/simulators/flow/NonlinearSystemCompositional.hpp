/*
  Copyright 2026, SINTEF Digital

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

#ifndef OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_HEADER_INCLUDED
#define OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_HEADER_INCLUDED

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/flow/NonlinearSystem.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>

#include <flowexperimental/comp/wells/CompWellModel.hpp>

#include <filesystem>
#include <vector>

namespace Opm {

template <class TypeTag>
class NonlinearSystemCompositional : public NonlinearSystem<TypeTag>
{
public:
    using ParentType = NonlinearSystem<TypeTag>;
    using Simulator = typename ParentType::Simulator;
    using Grid = typename ParentType::Grid;
    using FluidSystem = typename ParentType::FluidSystem;
    using Indices = typename ParentType::Indices;
    using Scalar = typename ParentType::Scalar;
    using ComponentName = typename ParentType::ComponentName;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using ModelParameters = BlackoilModelParameters<Scalar>;

    static constexpr int numEq = Indices::numEq;

    using VectorBlockType = Dune::FieldVector<Scalar, numEq>;
    using BVector = Dune::BlockVector<VectorBlockType>;

    NonlinearSystemCompositional(Simulator& simulator,
                                 const ModelParameters& param,
                                 CompWellModel<TypeTag>& wellModel,
                                 bool terminalOutput);

    SimulatorReportSingle prepareStep(const SimulatorTimerInterface& timer);

    void initialLinearization(SimulatorReportSingle& report,
                              int minIter,
                              int maxIter,
                              const SimulatorTimerInterface& timer) override;

    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIteration(const SimulatorTimerInterface& timer,
                                             NonlinearSolverType& nonlinearSolver);

    template <class NonlinearSolverType>
    SimulatorReportSingle nonlinearIterationNewton(const SimulatorTimerInterface& timer,
                                                   NonlinearSolverType& nonlinearSolver);

    Scalar relativeChange() const;

    int linearIterationsLastSolve() const
    { return this->simulator_.model().newtonMethod().linearSolver().iterations(); }

    void solveJacobianSystem(BVector& x);

    void updateSolution(const BVector& dx);

    bool hasNlddSolver() const
    { return false; }

    const SimulatorReport& localAccumulatedReports() const
    {
      static const SimulatorReport emptyReport{};
      return emptyReport;
    }

    const std::vector<SimulatorReport>& domainAccumulatedReports() const
    {
      static const std::vector<SimulatorReport> emptyReports{};
      return emptyReports;
    }

    void writeNonlinearIterationsPerCell(const std::filesystem::path&) const {}

    template<class T>
    std::vector<std::vector<Scalar>> computeFluidInPlace(const T&, const std::vector<int>& fipnum) const
    { return computeFluidInPlace(fipnum); }

    std::vector<std::vector<Scalar>> computeFluidInPlace(const std::vector<int>&) const
    { return {}; }

    void writePartitions(const std::filesystem::path&) const {}

private:
    std::vector<Scalar> reservoirResidualMetrics() const;

    double linear_solve_setup_time_ = 0.0;
};

} // namespace Opm

#include <opm/simulators/flow/NonlinearSystemCompositional_impl.hpp>

#endif // OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_HEADER_INCLUDED
