/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_GRAPH_MULTISEGMENTWELL_IMPL_HEADER_INCLUDED
#define OPM_GRAPH_MULTISEGMENTWELL_IMPL_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>

#include <opm/material/densead/Math.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/WellAssemble.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellConvergence.hpp>

#include <opm/models/utils/basicproperties.hh>

#include <cmath>
#include <functional>
#include <vector>

namespace Opm {

template<typename TypeTag>
GraphMultisegmentWell<TypeTag>::
GraphMultisegmentWell(const Well& well,
                      const ParallelWellInfo<Scalar>& pw_info,
                      const int time_step,
                      const ModelParameters& param,
                      const RateConverterType& rate_converter,
                      const int pvtRegionIdx,
                      const int num_conservation_quantities,
                      const int num_phases,
                      const int index_of_well,
                      const std::vector<PerforationData<Scalar>>& perf_data)
    : Base(well, pw_info, time_step, param, rate_converter, pvtRegionIdx,
           num_conservation_quantities, num_phases, index_of_well, perf_data)
{
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
init(const std::vector<Scalar>& depth_arg,
     const Scalar gravity_arg,
     const std::vector<Scalar>& B_avg,
     const bool changed_to_open_this_step)
{
    Base::init(depth_arg, gravity_arg, B_avg, changed_to_open_this_step);

    if (this->parallel_well_info_.communication().size() > 1) {
        OPM_THROW(std::runtime_error,
                  "GraphMultisegmentWell does not support wells distributed over "
                  "multiple MPI ranks (well " + this->name() + ").");
    }

    topo_ = GraphTopo::fromWellSegments(this->wellEcl().getSegments(),
                                        this->wellEcl().getConnections());
    graph_eqns_.init(topo_);

    // segment -> its outlet connection (the connection whose 'down' end is the segment)
    outlet_conn_.assign(topo_.numSegments(), -1);
    for (int c = 0; c < topo_.numConnections(); ++c)
        outlet_conn_[topo_.connection(c).down] = c;

    // MultisegmentWell primary-variable index -> GraphWell segment DOF
    using PVt = typename MSWEval::PrimaryVariables;
    var_to_gdof_.fill(-1);
    msw_wqtotal_ = PVt::WQTotal;
    var_to_gdof_[PVt::SPres] = 0;
    int gdof = 1;
    if constexpr (PVt::has_wfrac_variable)
        var_to_gdof_[PVt::WFrac] = gdof++;
    if constexpr (PVt::has_gfrac_variable)
        var_to_gdof_[PVt::GFrac] = gdof++;

    // GraphWell primary-variable phase map (component <-> fraction DOF)
    const int wc = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
        ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx) : -1;
    const int oc = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
        ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::oilCompIdx) : -1;
    const int gc = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
        ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx) : -1;
    std::array<int, NP> frac_comp{};
    if constexpr (PVt::has_wfrac_variable)
        frac_comp[var_to_gdof_[PVt::WFrac]] = wc;
    if constexpr (PVt::has_gfrac_variable)
        frac_comp[var_to_gdof_[PVt::GFrac]] = gc;
    const int ref_comp = (oc >= 0) ? oc : (wc >= 0 ? wc : gc);
    gpv_.setPhaseMap(ref_comp, frac_comp);
    ref_comp_ = ref_comp;
    gdof_comp_ = frac_comp;

    // Per-component rate scaling factors, identical to the production MultisegmentWell,
    // so the GraphWell unknowns map 1:1 onto the MultisegmentWell unknowns.
    std::array<Scalar, NP> scale{};
    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        if (!FluidSystem::phaseIsActive(phaseIdx))
            continue;
        const int comp = FluidSystem::canonicalToActiveCompIdx(FluidSystem::solventComponentIndex(phaseIdx));
        const int active_phase = FluidSystem::canonicalToActivePhaseIdx(phaseIdx);
        scale[comp] = this->scalingFactor(active_phase);
    }
    gpv_.setScale(scale);
    scale_comp_ = scale;
    gpv_.resize(topo_.numSegments(), topo_.numConnections());
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
setGraphStateFromBase()
{
    using PVt = typename MSWEval::PrimaryVariables;
    for (int seg = 0; seg < topo_.numSegments(); ++seg) {
        const auto& v = this->primary_variables_.value(seg);
        gpv_.segValue(seg, 0) = v[PVt::SPres];
        if constexpr (PVt::has_wfrac_variable)
            gpv_.segValue(seg, var_to_gdof_[PVt::WFrac]) = v[PVt::WFrac];
        if constexpr (PVt::has_gfrac_variable)
            gpv_.segValue(seg, var_to_gdof_[PVt::GFrac]) = v[PVt::GFrac];
    }
    for (int c = 0; c < topo_.numConnections(); ++c) {
        const int down = topo_.connection(c).down;
        gpv_.connValue(c) = -this->primary_variables_.value(down)[PVt::WQTotal];
    }
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
gatherPerforations(const Simulator& simulator,
                   const GroupStateHelperType& groupStateHelper,
                   std::vector<typename GAsm::PerforationInput>& perfs) const
{
    auto& deferred_logger = groupStateHelper.deferredLogger();
    perfs.assign(topo_.numPerforations(), typename GAsm::PerforationInput{});

    for (int seg = 0; seg < topo_.numSegments(); ++seg) {
        for (int perf : topo_.perforationsOfSegment(seg)) {
            const int cell_idx = this->cells()[perf];
            const auto& int_quants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/0);
            const auto& fs = int_quants.fluidState();

            auto& in = perfs[perf];
            in.pressure = fromRes(this->getPerfCellPressure(fs));
            in.rs = fromRes(Eval(fs.Rs()));
            in.rv = fromRes(Eval(fs.Rv()));

            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;
                const unsigned comp =
                    FluidSystem::canonicalToActiveCompIdx(FluidSystem::solventComponentIndex(phaseIdx));
                in.b[comp] = fromRes(Eval(fs.invB(phaseIdx)));
                in.mob[comp] = fromRes(Eval(int_quants.mobility(phaseIdx)));
            }
            static_cast<void>(deferred_logger);

            Scalar trans_mult = 1.0;
            this->getTransMult(trans_mult, simulator, cell_idx);
            in.Tw = this->wellIndex()[perf] * trans_mult;
            in.cell_perf_press_diff = this->cell_perforation_pressure_diffs_[perf];
        }
    }
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
assembleWellEqWithoutIteration(const Simulator& simulator,
                               const GroupStateHelperType& groupStateHelper,
                               const double dt,
                               const Well::InjectionControls& inj_controls,
                               const Well::ProductionControls& prod_controls,
                               WellStateType& well_state,
                               const bool solving_with_zero_rate)
{
    if (!this->isOperableAndSolvable() && !this->wellIsStopped())
        return;

    setGraphStateFromBase();

    // fluid properties from the GraphWell PVT port
    GFP fp(this->pvtRegionIdx());
    const auto fsinfo = this->getFirstPerforationFluidStateInfo(simulator);
    fp.update(gpv_, std::get<0>(fsinfo));

    std::vector<typename GAsm::PerforationInput> perfs;
    gatherPerforations(simulator, groupStateHelper, perfs);

    // initial (start-of-step) surface volumes for the accumulation term
    std::vector<std::array<Scalar, NP>> old_sv(topo_.numSegments(), std::array<Scalar, NP>{});
    for (int seg = 0; seg < topo_.numSegments(); ++seg)
        for (int comp = 0; comp < NP; ++comp)
            old_sv[seg][comp] = this->segment_fluid_initial_[seg][comp];

    typename GAsm::Controls ctrl;   // control eq assembled separately below
    ctrl.producer = this->isProducer();

    const bool allow_cf = this->getAllowCrossFlow()
        || this->openCrossFlowAvoidSingularity(simulator);

    const Scalar reg = this->regularize_ ? this->param_.regularization_factor_wells_ : Scalar{1};

    // Native AD assembly of mass conservation + momentum + perforation rates.
    assembler_.assemble(topo_, gpv_, fp, perfs, ctrl, static_cast<Scalar>(dt),
                        this->gravity(), allow_cf, old_sv, graph_eqns_,
                        /*make_solver=*/false, /*assemble_control=*/false, reg);

    assembleControl(simulator, groupStateHelper, well_state, inj_controls, prod_controls,
                    solving_with_zero_rate);

    updateWellStateRates(fp, perfs, allow_cf, well_state);
    graph_eqns_.createSolver();
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
updateWellStateRates(const GFP& fp,
                     const std::vector<typename GAsm::PerforationInput>& perfs,
                     bool allow_cf,
                     WellStateType& well_state)
{
    using Role = typename GPV::Role;
    auto& ws = well_state.well(this->index_of_well_);
    auto& perf_data = ws.perf_data;
    auto& perf_rates = perf_data.phase_rates;
    auto& perf_press = perf_data.pressure;
    ws.phase_mixing_rates.fill(0.0);

    const Scalar rhow = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
        ? FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, this->pvtRegionIdx()) : 0.0;
    const int watComp = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
        ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx) : -1;
    const bool producer = this->isProducer();

    for (int seg = 0; seg < topo_.numSegments(); ++seg) {
        const GEval p_seg = gpv_.segPressure(seg, Role::Self);
        const GEval rho_seg = fp.density(seg, Role::Self);
        for (int perf : topo_.perforationsOfSegment(seg)) {
            std::array<GEval, NP> cmix_s;
            for (int comp = 0; comp < NP; ++comp)
                cmix_s[comp] = gpv_.surfaceVolumeFraction(seg, comp, Role::Self);
            std::array<GEval, NP> cq_s;
            assembler_.computePerfRate(perfs[perf], p_seg, rho_seg, cmix_s,
                                       topo_.perforationDepthDiff(perf), this->gravity(),
                                       producer, allow_cf, cq_s);
            for (int comp = 0; comp < NP; ++comp) {
                const GEval cq_eff = cq_s[comp] * this->well_efficiency_factor_;
                this->connectionRates_[perf][comp] = this->restrictEval(cq_eff);
                perf_rates[perf * this->number_of_phases_
                           + FluidSystem::activeCompToActivePhaseIdx(comp)] = cq_s[comp].value();
            }
            const Scalar perf_seg_press_diff = this->gravity() * rho_seg.value()
                * topo_.perforationDepthDiff(perf);
            perf_press[perf] = p_seg.value() + perf_seg_press_diff;
            if (watComp >= 0)
                perf_data.wat_mass_rates[perf] = cq_s[watComp].value() * rhow;
        }
    }
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
assembleControl(const Simulator& simulator,
                const GroupStateHelperType& groupStateHelper,
                WellStateType& well_state,
                const Well::InjectionControls& inj_controls,
                const Well::ProductionControls& prod_controls,
                bool solving_with_zero_rate)
{
    // Compact AD type for the control equation: derivative slots are the top segment
    // DOFs [0,NP) plus the surface-connection flux at slot NP.
    static constexpr int QSlot = NP;
    using CEval = DenseAd::Evaluation<Scalar, NP + 1>;
    using CMath = MathToolbox<CEval>;
    static_cast<void>(CMath{});

    const int top = topo_.topSegment();
    const int sc = topo_.surfaceConnection();

    // bhp = top segment pressure (DOF 0); Q = surface connection flux (slot NP).
    CEval bhp(gpv_.segValue(top, 0));
    bhp.setDerivative(0, Scalar{1});
    CEval Q(gpv_.connValue(sc));
    Q.setDerivative(QSlot, Scalar{1});

    // scaled volume fraction F_p/g_p of the top segment (matches MSW volumeFractionScaled)
    auto volFracScaled = [&](int comp) -> CEval {
        const Scalar inv_scale = Scalar{1} / scale_comp_[comp];
        if (comp == ref_comp_) {
            Scalar v = Scalar{1};
            CEval f(Scalar{0});
            for (int k = 1; k < NP; ++k) {
                v -= gpv_.segValue(top, k);
                f.setDerivative(k, -inv_scale);
            }
            f.setValue(v * inv_scale);
            return f;
        }
        for (int k = 1; k < NP; ++k) {
            if (gdof_comp_[k] == comp) {
                CEval f(gpv_.segValue(top, k) * inv_scale);
                f.setDerivative(k, inv_scale);
                return f;
            }
        }
        return CEval(Scalar{0});
    };
    // surface phase rate getQs(comp) = WQTotal * volumeFractionScaled, WQTotal = -Q
    auto getQs = [&](int comp) -> CEval { return (-Q) * volFracScaled(comp); };

    // group-control bookkeeping, mirroring MultisegmentWell::assembleWellEqWithoutIteration
    const auto& summary_state = simulator.vanguard().summaryState();
    auto& deferred_logger = groupStateHelper.deferredLogger();
    const bool stopped_or_zero_target = this->stoppedOrZeroRateTarget(groupStateHelper);

    GroupState<Scalar> empty_group_state;
    auto& group_state = solving_with_zero_rate ? empty_group_state : groupStateHelper.groupState();
    GroupStateHelperType gsh = groupStateHelper;
    auto group_guard = gsh.pushGroupState(group_state);
    auto& ws = well_state.well(this->index_of_well_);
    if (this->wellUnderGroupControl(ws) && this->isProducer() && !stopped_or_zero_target)
        this->updateGroupTargetFallbackFlag(well_state, deferred_logger);

    static constexpr int Gas = FluidSystem::IndexTraitsType::gasPhaseIdx;
    static constexpr int Oil = FluidSystem::IndexTraitsType::oilPhaseIdx;
    static constexpr int Water = FluidSystem::IndexTraitsType::waterPhaseIdx;
    const auto& well = this->wellEcl();
    const Scalar rho = this->getRefDensity();

    auto getRates = [&]() {
        std::vector<CEval> rates(3, CEval(Scalar{0}));
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))
            rates[Water] = getQs(FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx));
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
            rates[Oil] = getQs(FluidSystem::canonicalToActiveCompIdx(FluidSystem::oilCompIdx));
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
            rates[Gas] = getQs(FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx));
        return rates;
    };

    CEval control_eq(Scalar{0});
    if (stopped_or_zero_target) {
        control_eq = -Q;   // WQTotal
    } else if (this->isInjector()) {
        const InjectorType injectorType = inj_controls.injector_type;
        int phase_pos = 0;
        switch (injectorType) {
        case InjectorType::WATER: phase_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx); break;
        case InjectorType::OIL:   phase_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx); break;
        case InjectorType::GAS:   phase_pos = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx); break;
        default: OPM_DEFLOG_THROW(std::runtime_error,
                     "Expected WATER, OIL or GAS as type for injector " + well.name(), deferred_logger);
        }
        const Scalar scaling = this->scalingFactor(phase_pos);
        const CEval injection_rate = (-Q) / scaling;   // WQTotal/scaling
        std::function<CEval()> bhp_from_thp = [&]() {
            const auto rates = getRates();
            return WellBhpThpCalculator(*this).template calculateBhpFromThp<CEval>(
                well_state, rates, well, summary_state, rho, deferred_logger);
        };
        WellAssemble(*this).template assembleControlEqInj<CEval>(
            gsh, inj_controls, bhp, injection_rate, bhp_from_thp, control_eq);
    } else {
        const auto rates = getRates();
        std::function<CEval()> bhp_from_thp = [&]() {
            return WellBhpThpCalculator(*this).template calculateBhpFromThp<CEval>(
                well_state, rates, well, summary_state, rho, deferred_logger);
        };
        WellAssemble(*this).template assembleControlEqProd<CEval>(
            gsh, prod_controls, bhp, rates, bhp_from_thp, control_eq);
    }

    // scatter into the surface-connection (control) row
    graph_eqns_.resConn()[sc][0] = control_eq.value();
    for (int k = 0; k < NP; ++k)
        graph_eqns_.Dcs()[sc][top][0][k] = control_eq.derivative(k);
    graph_eqns_.Dcc()[sc][sc][0][0] = control_eq.derivative(QSlot);
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
solveEqAndUpdateWellState(const Simulator& simulator,
                          const GroupStateHelperType& groupStateHelper,
                          WellStateType& well_state)
{
    if (!this->isOperableAndSolvable() && !this->wellIsStopped())
        return;
    try {
        const typename GraphEqns::BVectorWell gxw = graph_eqns_.solve();
        applyGraphSolution(simulator, gxw, groupStateHelper, well_state);
    }
    catch (const NumericalProblem& exp) {
        auto& deferred_logger = groupStateHelper.deferredLogger();
        deferred_logger.problem("In GraphMultisegmentWell::solveEqAndUpdateWellState for well "
                                + this->name() + ": " + exp.what());
        throw;
    }
}

template<typename TypeTag>
bool GraphMultisegmentWell<TypeTag>::
iterateWellEqWithControl(const Simulator& simulator,
                         const double dt,
                         const Well::InjectionControls& inj_controls,
                         const Well::ProductionControls& prod_controls,
                         const GroupStateHelperType& groupStateHelper,
                         WellStateType& well_state)
{
    if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return true;
    auto& deferred_logger = groupStateHelper.deferredLogger();
    const int max_iter_number = this->param_.max_inner_iter_ms_wells_;
    {
        const auto& [isFinite, residuals] = getFiniteWellResidualsGraph(Base::B_avg_, deferred_logger);
        static_cast<void>(residuals);
        if (!isFinite)
            return false;
    }
    this->updatePrimaryVariables(groupStateHelper);

    std::vector<std::vector<Scalar>> residual_history;
    std::vector<Scalar> measure_history;
    int it = 0;
    Scalar relaxation_factor = 1.;
    bool converged = false;
    bool relax_convergence = false;
    this->regularize_ = false;
    for (; it < max_iter_number; ++it) {
        if (it > this->param_.strict_inner_iter_wells_) {
            relax_convergence = true;
            this->regularize_ = true;
        }
        assembleWellEqWithoutIteration(simulator, groupStateHelper, dt, inj_controls, prod_controls,
                                       well_state, /*solving_with_zero_rate=*/false);
        const auto report = getWellConvergence(groupStateHelper, Base::B_avg_, relax_convergence);
        if (report.converged()) {
            converged = true;
            break;
        }
        {
            const auto& [isFinite, residuals] = getFiniteWellResidualsGraph(Base::B_avg_, deferred_logger);
            if (!isFinite)
                return false;
            residual_history.push_back(residuals);
            measure_history.push_back(this->getResidualMeasureValue(well_state, residual_history[it],
                this->param_.tolerance_wells_, this->param_.tolerance_pressure_ms_wells_, deferred_logger));
        }
        bool min_relaxation_reached =
            this->update_relaxation_factor(measure_history, relaxation_factor, this->regularize_, deferred_logger);
        if (min_relaxation_reached || this->repeatedStagnation(measure_history, this->regularize_, deferred_logger)) {
            const auto reportStag = getWellConvergence(groupStateHelper, Base::B_avg_, true);
            converged = reportStag.converged();
            break;
        }
        try {
            const auto gxw = graph_eqns_.solve();
            applyGraphSolution(simulator, gxw, groupStateHelper, well_state, relaxation_factor);
        } catch (const NumericalProblem& exp) {
            deferred_logger.problem("In GraphMultisegmentWell::iterateWellEqWithControl for well "
                                    + this->name() + ": " + exp.what());
            throw;
        }
    }
    return converged;
}

template<typename TypeTag>
bool GraphMultisegmentWell<TypeTag>::
iterateWellEqWithSwitching(const Simulator& simulator,
                           const double dt,
                           const Well::InjectionControls& inj_controls,
                           const Well::ProductionControls& prod_controls,
                           const GroupStateHelperType& groupStateHelper,
                           WellStateType& well_state,
                           const bool fixed_control,
                           const bool fixed_status,
                           const bool solving_with_zero_rate)
{
    auto& deferred_logger = groupStateHelper.deferredLogger();
    const int max_iter_number = this->param_.max_inner_iter_ms_wells_;
    {
        const auto& [isFinite, residuals] = getFiniteWellResidualsGraph(Base::B_avg_, deferred_logger);
        static_cast<void>(residuals);
        if (!isFinite)
            return false;
    }
    this->updatePrimaryVariables(groupStateHelper);

    std::vector<std::vector<Scalar>> residual_history;
    std::vector<Scalar> measure_history;
    int it = 0;
    Scalar relaxation_factor = 1.;
    bool converged = false;
    bool relax_convergence = false;
    this->regularize_ = false;
    const auto& summary_state = groupStateHelper.summaryState();

    const int min_its_after_switch = 3;
    const int max_status_switch = this->param_.max_well_status_switch_inner_iter_;
    int its_since_last_switch = min_its_after_switch;
    int switch_count = 0;
    int status_switch_count = 0;
    const auto well_status_orig = this->wellStatus_;
    const auto operability_orig = this->operability_status_;
    auto well_status_cur = well_status_orig;
    const bool allow_open = well_state.well(this->index_of_well_).status == WellStatus::OPEN;
    const bool allow_switching = !this->wellUnderZeroRateTarget(groupStateHelper) &&
                                 (!fixed_control || !fixed_status) && allow_open;
    bool final_check = false;
    this->operability_status_.resetOperability();
    this->operability_status_.solvable = true;

    for (; it < max_iter_number; ++it) {
        ++its_since_last_switch;
        if (allow_switching && its_since_last_switch >= min_its_after_switch && status_switch_count < max_status_switch) {
            const Scalar wqTotal = this->primary_variables_.getWQTotal().value();
            bool changed = this->updateWellControlAndStatusLocalIteration(
                simulator, groupStateHelper, inj_controls, prod_controls, wqTotal,
                well_state, fixed_control, fixed_status, solving_with_zero_rate);
            if (changed) {
                its_since_last_switch = 0;
                ++switch_count;
                if (well_status_cur != this->wellStatus_) {
                    well_status_cur = this->wellStatus_;
                    status_switch_count++;
                }
            }
            if (!changed && final_check)
                break;
            else
                final_check = false;
            if (status_switch_count == max_status_switch)
                this->wellStatus_ = well_status_orig;
        }

        if (it > this->param_.strict_inner_iter_wells_) {
            relax_convergence = true;
            this->regularize_ = true;
        }
        assembleWellEqWithoutIteration(simulator, groupStateHelper, dt, inj_controls, prod_controls,
                                       well_state, solving_with_zero_rate);
        const auto report = getWellConvergence(groupStateHelper, Base::B_avg_, relax_convergence);
        converged = report.converged();
        if (converged) {
            if (switch_count > 0 && its_since_last_switch < min_its_after_switch) {
                final_check = true;
                its_since_last_switch = min_its_after_switch;
            } else {
                break;
            }
        }
        {
            const auto& [isFinite, residuals] = getFiniteWellResidualsGraph(Base::B_avg_, deferred_logger);
            if (!isFinite) {
                converged = false;
                break;
            }
            residual_history.push_back(residuals);
        }
        if (!converged) {
            measure_history.push_back(this->getResidualMeasureValue(well_state, residual_history[it],
                this->param_.tolerance_wells_, this->param_.tolerance_pressure_ms_wells_, deferred_logger));
            bool min_relaxation_reached =
                this->update_relaxation_factor(measure_history, relaxation_factor, this->regularize_, deferred_logger);
            if (min_relaxation_reached || this->repeatedStagnation(measure_history, this->regularize_, deferred_logger)) {
                converged = false;
                break;
            }
        }
        try {
            const auto gxw = graph_eqns_.solve();
            applyGraphSolution(simulator, gxw, groupStateHelper, well_state, relaxation_factor);
        } catch (const NumericalProblem& exp) {
            deferred_logger.problem("In GraphMultisegmentWell::iterateWellEqWithSwitching for well "
                                    + this->name() + ": " + exp.what());
            throw;
        }
    }

    if (converged) {
        if (allow_switching) {
            const bool is_stopped = this->wellIsStopped();
            if (this->wellHasTHPConstraints(summary_state)) {
                this->operability_status_.can_obtain_bhp_with_thp_limit = !is_stopped;
                this->operability_status_.obey_thp_limit_under_bhp_limit = !is_stopped;
            } else {
                this->operability_status_.operable_under_only_bhp_limit = !is_stopped;
            }
        }
    } else {
        this->wellStatus_ = well_status_orig;
        this->operability_status_ = operability_orig;
        this->primary_variables_.outputLowLimitPressureSegments(deferred_logger);
    }
    return converged;
}

template<typename TypeTag>
ConvergenceReport GraphMultisegmentWell<TypeTag>::
getWellConvergence(const GroupStateHelperType& groupStateHelper,
                   const std::vector<Scalar>& B_avg,
                   const bool relax_tolerance) const
{
    using CR = ConvergenceReport;
    const auto& well_state = groupStateHelper.wellState();
    auto& deferred_logger = groupStateHelper.deferredLogger();
    const int numCons = this->numConservationQuantities();
    const bool well_is_stopped = this->wellIsStopped();

    const Scalar max_residual_allowed = this->param_.max_residual_allowed_;
    const Scalar tolerance_wells = this->param_.tolerance_wells_;
    const Scalar relaxed_flow = this->param_.relaxed_tolerance_flow_well_;
    const Scalar tolerance_pressure = this->param_.tolerance_pressure_ms_wells_;
    const Scalar relaxed_pressure = this->param_.relaxed_tolerance_pressure_ms_well_;

    const auto& resSeg = graph_eqns_.residual()[Dune::Indices::_0];
    const auto& resConn = graph_eqns_.residual()[Dune::Indices::_1];

    // mass-balance residuals (scaled by B_avg), maximised over segments
    std::vector<Scalar> max_mass(numCons, 0.0);
    for (int seg = 0; seg < topo_.numSegments(); ++seg)
        for (int comp = 0; comp < numCons; ++comp)
            max_mass[comp] = std::max(max_mass[comp], B_avg[comp] * std::abs(resSeg[seg][comp]));

    // pressure-equation residuals on the internal connections
    Scalar max_pressure = 0.0;
    for (int c = 0; c < topo_.numConnections(); ++c)
        if (c != topo_.surfaceConnection())
            max_pressure = std::max(max_pressure, std::abs(resConn[c][0]));

    ConvergenceReport report;
    for (int comp = 0; comp < numCons; ++comp) {
        const Scalar r = max_mass[comp];
        if (std::isnan(r))
            report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::NotANumber, comp, this->name()});
        else if (r > max_residual_allowed)
            report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::TooLarge, comp, this->name()});
        else if (!relax_tolerance && r > tolerance_wells)
            report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::Normal, comp, this->name()});
        else if (r > relaxed_flow)
            report.setWellFailed({CR::WellFailure::Type::MassBalance, CR::Severity::Normal, comp, this->name()});
    }
    {
        const Scalar r = max_pressure;
        const int dummy = -1;
        if (std::isnan(r))
            report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::NotANumber, dummy, this->name()});
        else if (std::isinf(r))
            report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::TooLarge, dummy, this->name()});
        else if (!relax_tolerance && r > tolerance_pressure)
            report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::Normal, dummy, this->name()});
        else if (r > relaxed_pressure)
            report.setWellFailed({CR::WellFailure::Type::Pressure, CR::Severity::Normal, dummy, this->name()});
    }

    // control equation residual (surface connection)
    WellConvergence(*this).
        checkConvergenceControlEq(well_state,
                                  {tolerance_pressure, tolerance_pressure,
                                   tolerance_wells, tolerance_wells, max_residual_allowed},
                                  std::abs(resConn[topo_.surfaceConnection()][0]),
                                  well_is_stopped, report, deferred_logger);
    return report;
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
apply(const BVector& x, BVector& Ax) const
{
    if (!this->isOperableAndSolvable() && !this->wellIsStopped())
        return;
    if (this->param_.matrix_add_well_contributions_)
        return;
    graph_eqns_.apply(x, Ax);
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
addWellPressureEquations(PressureMatrix& /*mat*/,
                         const BVector& /*x*/,
                         const int /*pressureVarIndex*/,
                         const bool /*use_well_weights*/,
                         const WellStateType& /*well_state*/) const
{
    // No-op: the GraphWell does not contribute to the CPRW well-pressure system.
    // The well is fully accounted for through the Schur apply() path.
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
apply(BVector& r) const
{
    if (!this->isOperableAndSolvable() && !this->wellIsStopped())
        return;
    graph_eqns_.apply(r);
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
recoverWellSolutionAndUpdateWellState(const Simulator& simulator,
                                      const BVector& x,
                                      const GroupStateHelperType& groupStateHelper,
                                      WellStateType& well_state)
{
    if (!this->isOperableAndSolvable() && !this->wellIsStopped())
        return;

    typename GraphEqns::BVectorWell gxw;
    graph_eqns_.recoverSolutionWell(x, gxw);
    applyGraphSolution(simulator, gxw, groupStateHelper, well_state);
}

template<typename TypeTag>
void GraphMultisegmentWell<TypeTag>::
applyGraphSolution(const Simulator& simulator,
                   const typename GraphEqns::BVectorWell& gxw,
                   const GroupStateHelperType& groupStateHelper,
                   WellStateType& well_state,
                   Scalar relaxation_factor)
{
    using namespace Dune::Indices;
    using PVt = typename MSWEval::PrimaryVariables;
    constexpr int numWellEq = MSWEval::numWellEq;

    typename MSWEval::BVectorWell dwells(topo_.numSegments());
    for (int seg = 0; seg < topo_.numSegments(); ++seg) {
        for (int e = 0; e < numWellEq; ++e)
            dwells[seg][e] = 0.0;
        dwells[seg][PVt::SPres] = gxw[_0][seg][0];
        if constexpr (PVt::has_wfrac_variable)
            dwells[seg][PVt::WFrac] = gxw[_0][seg][var_to_gdof_[PVt::WFrac]];
        if constexpr (PVt::has_gfrac_variable)
            dwells[seg][PVt::GFrac] = gxw[_0][seg][var_to_gdof_[PVt::GFrac]];
        // WQTotal = -Q  =>  d(WQTotal) = -dQ
        dwells[seg][PVt::WQTotal] = -gxw[_1][outlet_conn_[seg]][0];
    }

    this->updateWellState(simulator, dwells, groupStateHelper, well_state, relaxation_factor);
}

template<typename TypeTag>
std::pair<bool, std::vector<typename GraphMultisegmentWell<TypeTag>::Scalar>>
GraphMultisegmentWell<TypeTag>::
getFiniteWellResidualsGraph(const std::vector<Scalar>& B_avg,
                            DeferredLogger& deferred_logger) const
{
    constexpr int numWellEq = MSWEval::numWellEq;
    const int numCons = this->numConservationQuantities();
    std::vector<Scalar> residuals(numWellEq + 1, 0.0);

    const auto& resSeg = graph_eqns_.residual()[Dune::Indices::_0];
    const auto& resConn = graph_eqns_.residual()[Dune::Indices::_1];

    // mass-balance residuals (indices [0, numCons))
    for (int seg = 0; seg < topo_.numSegments(); ++seg) {
        for (int comp = 0; comp < numCons; ++comp) {
            const Scalar r = std::abs(resSeg[seg][comp]) * B_avg[comp];
            if (std::isnan(r) || std::isinf(r)) {
                deferred_logger.debug("nan/inf residual in GraphMultisegmentWell " + this->name());
                return {false, residuals};
            }
            residuals[comp] = std::max(residuals[comp], r);
        }
    }
    // pressure-equation residuals (index numWellEq-1) over internal connections
    for (int c = 0; c < topo_.numConnections(); ++c) {
        if (c == topo_.surfaceConnection())
            continue;
        const Scalar r = std::abs(resConn[c][0]);
        if (std::isnan(r) || std::isinf(r))
            return {false, residuals};
        residuals[numWellEq - 1] = std::max(residuals[numWellEq - 1], r);
    }
    // control-equation residual (index numWellEq)
    {
        const Scalar r = std::abs(resConn[topo_.surfaceConnection()][0]);
        if (std::isnan(r) || std::isinf(r))
            return {false, residuals};
        residuals[numWellEq] = r;
    }
    return {true, residuals};
}

} // namespace Opm

#endif // OPM_GRAPH_MULTISEGMENTWELL_IMPL_HEADER_INCLUDED
