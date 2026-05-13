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

#ifndef OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_IMPL_HEADER_INCLUDED
#define OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_IMPL_HEADER_INCLUDED

#ifndef OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/flow/NonlinearSystemCompositional.hpp>
#endif

#include <dune/common/timer.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <algorithm>
#include <array>
#include <cmath>

namespace Opm {

template <class TypeTag>
NonlinearSystemCompositional<TypeTag>::
NonlinearSystemCompositional(Simulator& simulator,
                             const ModelParameters& param,
                             CompWellModel<TypeTag>& wellModel,
                             const bool terminalOutput)
    : ParentType(simulator, terminalOutput)
    , param_(param)
    , well_model_(wellModel)
    , current_relaxation_(1.0)
    , dx_old_(this->simulator_.model().numGridDof())
{
    this->convergence_reports_.reserve(64);
}

template <class TypeTag>
SimulatorReportSingle
NonlinearSystemCompositional<TypeTag>::
prepareStep(const SimulatorTimerInterface& timer)
{
    return ParentType::prepareStep(timer);
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
initialLinearization(SimulatorReportSingle& report,
                     const int minIter,
                     const int,
                     const SimulatorTimerInterface& timer)
{
    ParentType::initialLinearization(report,
                                     timer,
                                     [this](const SimulatorTimerInterface& innerTimer)
                                     {
                                         return this->assembleReservoir(innerTimer);
                                     });

    Dune::Timer perfTimer;
    perfTimer.start();

    auto convrep = ConvergenceReport{timer.simulationTimeElapsed()};
    const auto residualMetrics = this->reservoirResidualMetrics();
    const auto tolerance = this->simulator_.model().newtonMethod().tolerance();

    for (int compIdx = 0; compIdx < numEq; ++compIdx) {
        const std::array<Scalar, 1> residual{residualMetrics[compIdx]};
        const std::array<ConvergenceReport::ReservoirFailure::Type, 1> types{
            ConvergenceReport::ReservoirFailure::Type::MassBalance
        };
        const std::array<Scalar, 1> tolerances{tolerance};

        this->addReservoirConvergenceMetrics(
            convrep,
            compIdx,
            this->compNames_.name(compIdx),
            residual,
            types,
            tolerances,
            param_.max_residual_allowed_,
            [this](const std::string& message)
            {
                if (this->terminal_output_) {
                    OpmLog::debug(message);
                }
            });
    }

    const auto severity = convrep.severityOfWorstFailure();
    const bool wellConverged = this->wellModel().getWellConvergence();
    report.converged = convrep.converged() && wellConverged &&
                       this->simulator_.problem().iterationContext().iteration() >= minIter;

    this->convergence_reports_.back().report.push_back(std::move(convrep));
    report.update_time += perfTimer.stop();
    this->residual_norms_history_.push_back(residualMetrics);

    if (severity == ConvergenceReport::Severity::NotANumber) {
        this->failureReport_ += report;
        OPM_THROW_PROBLEM(NumericalProblem, "NaN residual found!");
    }

    if (severity == ConvergenceReport::Severity::TooLarge) {
        this->failureReport_ += report;
        OPM_THROW_NOLOG(NumericalProblem, "Too large residual found!");
    }
}

template <class TypeTag>
template <class NonlinearSolverType>
SimulatorReportSingle
NonlinearSystemCompositional<TypeTag>::
nonlinearIteration(const SimulatorTimerInterface& timer,
                   NonlinearSolverType& nonlinearSolver)
{
    return ParentType::nonlinearIteration(
        timer,
        numEq,
        [this]()
        {
            this->residual_norms_history_.clear();
            this->current_relaxation_ = 1.0;
            this->dx_old_ = 0.0;
        },
        [this, &timer, &nonlinearSolver]() -> SimulatorReportSingle
        {
            return this->nonlinearIterationNewton(timer, nonlinearSolver);
        },
        []() {});
}

template <class TypeTag>
template <class NonlinearSolverType>
SimulatorReportSingle
NonlinearSystemCompositional<TypeTag>::
nonlinearIterationNewton(const SimulatorTimerInterface& timer,
                         NonlinearSolverType& nonlinearSolver)
{
    return ParentType::nonlinearIterationNewton(timer,
                                                nonlinearSolver,
                                                this->param_,
                                                this->wellModel(),
                                                this->residual_norms_history_,
                                                this->dx_old_,
                                                this->current_relaxation_,
                                                this->linear_solve_setup_time_,
                                                [this](SimulatorReportSingle& report,
                                                       const int minIter,
                                                       const int maxIter,
                                                       const SimulatorTimerInterface& innerTimer)
                                                {
                                                    this->initialLinearization(report, minIter, maxIter, innerTimer);
                                                },
                                                [this](BVector& x)
                                                {
                                                    this->solveJacobianSystem(x);
                                                },
                                                [this]()
                                                {
                                                    return this->linearIterationsLastSolve();
                                                },
                                                [this](const BVector& x)
                                                {
                                                    this->updateSolution(x);
                                                });
}

template <class TypeTag>
SimulatorReportSingle
NonlinearSystemCompositional<TypeTag>::
assembleReservoir(const SimulatorTimerInterface&)
{
    this->simulator_.problem().beginIteration();
    this->simulator_.model().linearizer().linearizeDomain();
    this->simulator_.problem().endIteration();
    return this->wellModel().lastReport();
}

template <class TypeTag>
typename NonlinearSystemCompositional<TypeTag>::Scalar
NonlinearSystemCompositional<TypeTag>::
relativeChange() const
{
    Scalar resultDelta = 0.0;
    Scalar resultDenom = 0.0;

    const auto& elemMapper = this->simulator_.model().elementMapper();
    const auto& gridView = this->simulator_.gridView();

    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        const unsigned globalElemIdx = elemMapper.index(elem);
        const auto& priVarsNew = this->simulator_.model().solution(/*timeIdx=*/0)[globalElemIdx];
        const auto& priVarsOld = this->simulator_.model().solution(/*timeIdx=*/1)[globalElemIdx];

        for (int pvIdx = 0; pvIdx < static_cast<int>(priVarsNew.size()); ++pvIdx) {
            const auto delta = priVarsNew[pvIdx] - priVarsOld[pvIdx];
            resultDelta += delta * delta;
            resultDenom += priVarsNew[pvIdx] * priVarsNew[pvIdx];
        }
    }

    resultDelta = gridView.comm().sum(resultDelta);
    resultDenom = gridView.comm().sum(resultDenom);

    return resultDenom > 0.0 ? resultDelta / resultDenom : 0.0;
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
solveJacobianSystem(BVector& x)
{
    auto& jacobian = this->simulator_.model().linearizer().jacobian();
    auto& residual = this->simulator_.model().linearizer().residual();
    auto& linSolver = this->simulator_.model().newtonMethod().linearSolver();

    x = 0.0;

    Dune::Timer perfTimer;
    perfTimer.start();
    linSolver.prepare(jacobian, residual);
    this->linear_solve_setup_time_ = perfTimer.stop();
    linSolver.setResidual(residual);
    linSolver.solve(x);
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
updateSolution(const BVector& dx)
{
    ParentType::updateSolution(dx,
                               false,
                               []() {},
                               [](const BVector&) {});
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
updateTUNING(const Tuning& tuning)
{
    this->param_.tolerance_cnv_ = tuning.TRGCNV;
    this->param_.tolerance_cnv_relaxed_ = tuning.XXXCNV;
    this->param_.tolerance_mb_ = tuning.TRGMBE;
    this->param_.tolerance_mb_relaxed_ = tuning.XXXMBE;
    this->param_.newton_max_iter_ = tuning.NEWTMX;
    this->param_.newton_min_iter_ = tuning.NEWTMN;
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
updateTUNINGDP(const TuningDp& tuning_dp)
{
    this->param_.tolerance_max_dp_ = tuning_dp.TRGDDP;
    this->param_.tolerance_max_ds_ = tuning_dp.TRGDDS;
    this->param_.tolerance_max_drs_ = tuning_dp.TRGDDRS;
    this->param_.tolerance_max_drv_ = tuning_dp.TRGDDRV;
}

template <class TypeTag>
std::vector<typename NonlinearSystemCompositional<TypeTag>::Scalar>
NonlinearSystemCompositional<TypeTag>::
reservoirResidualMetrics() const
{
    const auto& model = this->simulator_.model();
    const auto& residual = model.linearizer().residual();
    const auto& constraintsMap = model.linearizer().constraintsMap();

    std::vector<Scalar> residualMetrics(numEq, 0.0);

    for (unsigned dofIdx = 0; dofIdx < residual.size(); ++dofIdx) {
        if (dofIdx >= model.numGridDof() || model.dofTotalVolume(dofIdx) <= 0.0) {
            continue;
        }

        if (constraintsMap.count(dofIdx) > 0) {
            continue;
        }

        const auto& localResidual = residual[dofIdx];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            residualMetrics[eqIdx] = std::max(
                residualMetrics[eqIdx],
                std::abs(localResidual[eqIdx] * model.eqWeight(dofIdx, eqIdx)));
        }
    }

    if (this->grid_.comm().size() > 1 && !residualMetrics.empty()) {
        this->grid_.comm().max(residualMetrics.data(), residualMetrics.size());
    }

    return residualMetrics;
}

} // namespace Opm

#endif // OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_IMPL_HEADER_INCLUDED
