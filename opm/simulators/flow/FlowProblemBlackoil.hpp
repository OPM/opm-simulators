// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2023 INRIA
  Copyright 2024 SINTEF Digital

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::FlowProblemBlackoil
 */
#ifndef OPM_FLOW_PROBLEM_BLACK_HPP
#define OPM_FLOW_PROBLEM_BLACK_HPP

#include <opm/models/blackoil/blackoilmodules.hpp>

#include <opm/output/eclipse/EclipseIO.hpp>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQSet.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/io/eclipse/SummaryNode.hpp>

#include <opm/simulators/flow/ActionHandler.hpp>
#include <opm/grid/LookUpData.hh>

#include <opm/simulators/aquifers/NumericalAquiferAuxCells.hpp>
#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/flow/MixingRateControls.hpp>
#include <opm/simulators/flow/OutputBlackoilModule.hpp>
#include <opm/simulators/flow/VtkTracerModule.hpp>
#include <opm/simulators/flow/HybridNewton.hpp>
#include <opm/simulators/flow/HybridNewtonConfig.hpp>

#include <opm/simulators/utils/satfunc/SatfuncConsistencyCheckManager.hpp>

#if HAVE_DAMARIS
#include <opm/simulators/flow/DamarisWriter.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace Opm {

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief This problem simulates an input file given in the data format used by the
 *        commercial ECLiPSE simulator.
 */
template <class TypeTag>
class FlowProblemBlackoil : public FlowProblem<TypeTag>
{
    // TODO: the naming of the Types might be able to be adjusted
public:
    using FlowProblemType = FlowProblem<TypeTag>;

private:
    using typename FlowProblemType::Scalar;
    using typename FlowProblemType::Simulator;
    using typename FlowProblemType::GridView;
    using typename FlowProblemType::FluidSystem;
    using typename FlowProblemType::Vanguard;
    using typename FlowProblemType::GlobalEqVector;
    using typename FlowProblemType::EqVector;
    using FlowProblemType::dim;
    using FlowProblemType::dimWorld;
    using FlowProblemType::numEq;
    using FlowProblemType::numPhases;
    using FlowProblemType::numComponents;

    // TODO: potentially some cleaning up depending on the usage later here
    using FlowProblemType::enableBioeffects;
    using FlowProblemType::enableBrine;
    using FlowProblemType::enableConvectiveMixing;
    using FlowProblemType::enableDiffusion;
    using FlowProblemType::enableDispersion;
    using FlowProblemType::energyModuleType;
    using FlowProblemType::enableExperiments;
    using FlowProblemType::enableExtbo;
    using FlowProblemType::enableFoam;
    using FlowProblemType::enableMICP;
    using FlowProblemType::enablePolymer;
    using FlowProblemType::enablePolymerMolarWeight;
    using FlowProblemType::enableSaltPrecipitation;
    using FlowProblemType::enableSolvent;
    using FlowProblemType::enableThermalFluxBoundaries;

    using FlowProblemType::gasPhaseIdx;
    using FlowProblemType::oilPhaseIdx;
    using FlowProblemType::waterPhaseIdx;

    using FlowProblemType::waterCompIdx;
    using FlowProblemType::oilCompIdx;
    using FlowProblemType::gasCompIdx;

    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using typename FlowProblemType::RateVector;
    using typename FlowProblemType::PrimaryVariables;
    using typename FlowProblemType::Indices;
    using typename FlowProblemType::IntensiveQuantities;
    using typename FlowProblemType::ElementContext;

    using typename FlowProblemType::MaterialLaw;
    using typename FlowProblemType::DimMatrix;

    static constexpr bool enableDissolvedGas =
        Indices::compositionSwitchIdx != std::numeric_limits<unsigned>::max();
    enum { enableVapwat = getPropValue<TypeTag, Properties::EnableVapwat>() };
    enum { enableDisgasInWater = getPropValue<TypeTag, Properties::EnableDisgasInWater>() };
    enum { enableGeochemistry = getPropValue<TypeTag, Properties::EnableGeochemistry>() };
    enum { enableMech = getPropValue<TypeTag, Properties::EnableMech>() };

    using BioeffectsModule = BlackOilBioeffectsModule<TypeTag, enableBioeffects>;
    using BrineModule = BlackOilBrineModule<TypeTag, enableBrine>;
    using ConvectiveMixingModule = BlackOilConvectiveMixingModule<TypeTag, enableConvectiveMixing>;
    using DiffusionModule = BlackOilDiffusionModule<TypeTag, enableDiffusion>;
    using DispersionModule = BlackOilDispersionModule<TypeTag, enableDispersion>;
    using ExtboModule = BlackOilExtboModule<TypeTag, enableExtbo>;
    using FoamModule = BlackOilFoamModule<TypeTag, enableFoam>;
    using PolymerModule = BlackOilPolymerModule<TypeTag, enablePolymer>;
    using SolventModule = BlackOilSolventModule<TypeTag, enableSolvent>;

    using EclWriterType = EclWriter<TypeTag, OutputBlackOilModule<TypeTag> >;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using InitialFluidState = typename EquilInitializer<TypeTag>::ScalarFluidState;
    using HybridNewton = BlackOilHybridNewton<TypeTag>;

#if HAVE_DAMARIS
    using DamarisWriterType = DamarisWriter<TypeTag>;
#endif

public:
    using FlowProblemType::porosity;
    using FlowProblemType::pvtRegionIndex;

    /*!
     * \copydoc FvBaseProblem::registerParameters
     */
    static void registerParameters()
    {
        FlowProblemType::registerParameters();

        EclWriterType::registerParameters();
#if HAVE_DAMARIS
        DamarisWriterType::registerParameters();
#endif
        VtkTracerModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    explicit FlowProblemBlackoil(Simulator& simulator)
        : FlowProblemType(simulator)
        , mixControls_(simulator.vanguard().schedule())
        , prev_timestep_state_(simulator.vanguard().schedule())
        , actionHandler_(simulator.vanguard().eclState(),
                         simulator.vanguard().schedule(),
                         simulator.vanguard().actionState(),
                         simulator.vanguard().summaryState(),
                         this->wellModel_,
                         simulator.vanguard().grid().comm())
        , hybridNewton_(simulator)
    {
        this->model().addOutputModule(std::make_unique<VtkTracerModule<TypeTag>>(simulator));

        // Tell the black-oil extensions to initialize their internal data structures
        const auto& vanguard = simulator.vanguard();

        if constexpr (enableBrine) {
            BlackOilBrineParams<Scalar> brineParams;
            brineParams.template initFromState<enableBrine,
                                               enableSaltPrecipitation>(vanguard.eclState());
            BrineModule::setParams(std::move(brineParams));
        }

        if constexpr (enableDiffusion) {
            DiffusionModule::initFromState(vanguard.eclState());
        }

        if constexpr (enableDispersion) {
            DispersionModule::initFromState(vanguard.eclState());
        }

        if constexpr (enableExtbo) {
            BlackOilExtboParams<Scalar> extboParams;
            extboParams.template initFromState<enableExtbo>(vanguard.eclState());
            ExtboModule::setParams(std::move(extboParams));
        }

        if constexpr (enableFoam) {
            BlackOilFoamParams<Scalar> foamParams;
            foamParams.template initFromState<enableFoam>(vanguard.eclState());
            FoamModule::setParams(std::move(foamParams));
        }

        if constexpr (enableBioeffects) {
            BlackOilBioeffectsParams<Scalar> bioeffectsParams;
            bioeffectsParams.template initFromState<enableBioeffects, enableMICP>(vanguard.eclState());
            BioeffectsModule::setParams(std::move(bioeffectsParams));
        }

        if constexpr (enablePolymer) {
            BlackOilPolymerParams<Scalar> polymerParams;
            polymerParams.template initFromState<enablePolymer, enablePolymerMolarWeight>(vanguard.eclState());
            PolymerModule::setParams(std::move(polymerParams));
        }

        if constexpr (enableSolvent) {
            BlackOilSolventParams<Scalar> solventParams;
            solventParams.template initFromState<enableSolvent>(vanguard.eclState(), vanguard.schedule());
            SolventModule::setParams(std::move(solventParams));
        }

        // create the ECL writer
        eclWriter_ = std::make_unique<EclWriterType>(simulator);
        enableEclOutput_ = Parameters::Get<Parameters::EnableEclOutput>();

        // Safeguard against geochemistry since it exsist in a separate module with a separate problem class
        if constexpr (!enableGeochemistry) {
            if (vanguard.eclState().runspec().geochem().enabled()) {
                throw std::runtime_error("GEOCHEM keyword in the deck but geochemistry module "
                                         "was disabled at compile time!");
            }
        }

        // Safeguard against TPSA-geomechanics since it requires FlowProblemTPSA
        if constexpr (!enableMech) {
            const auto& rspec = vanguard.eclState().runspec();
            if (rspec.mech() && rspec.mechSolver().tpsa()) {
                throw std::runtime_error("TPSA solver enabled in the deck, but geomechanics "
                                         "module was disabled at compile time!");
            }
        }

#if HAVE_DAMARIS
        // create Damaris writer
        damarisWriter_ = std::make_unique<DamarisWriterType>(simulator);
        enableDamarisOutput_ = Parameters::Get<Parameters::EnableDamarisOutput>();
#endif
    }

    /*!
     * \brief Called by the simulator before an episode begins.
     */
    void beginEpisode() override
    {
        // Before the well model sets up the report step's group controls, which
        // need the wells outside the region to be counted already. Until the
        // first time step is chosen this is the latest interval the parent has
        // an average for, which is the best there is.
        {
            const auto episodeStart = static_cast<double>(this->simulator().episodeStartTime())
                - static_cast<double>(this->simulator().startTime());
            this->updateFluxFixedWellRates_(this->simulator().episodeIndex(),
                                            episodeStart, episodeStart);
        }

        FlowProblemType::beginEpisode();

        auto& simulator = this->simulator();

        const int episodeIdx = simulator.episodeIndex();
        auto& vanguard = simulator.vanguard();
        const auto& schedule = vanguard.schedule();

        // Time here is the start of the episode. simulator.time() has not been
        // advanced yet at this point and timeStepSize() still holds the WHOLE
        // episode length, so the usual time + timeStepSize would land at the
        // end of the episode and pick up the parent's state from days later.
        // That is what produced a one-step excursion at every report step: the
        // boundary was briefly driven by a far more depleted exterior state.
        //
        // The boundary record itself is deliberately not refreshed here.
        // beginTimeStep runs immediately afterwards for the first step of the
        // episode and selects it with the proper step size.
        const auto episodeStart = static_cast<double>(simulator.episodeStartTime())
                                - static_cast<double>(simulator.startTime());

        if (episodeIdx >= 0) {
            this->seedParentSummaryState_(episodeStart, episodeStart);
        }

        // Evaluate UDQ assign statements to make sure the settings are
        // available as UDA controls for the current report step.
        this->actionHandler_
            .evalUDQAssignments(episodeIdx, vanguard.udqState());

        this->restoreFluxGroupTargetUDQs_(episodeIdx);

        if (episodeIdx >= 0) {
            const auto& oilVap = schedule[episodeIdx].oilvap();
            if (oilVap.getType() == OilVaporizationProperties::OilVaporization::VAPPARS) {
                FluidSystem::setVapPars(oilVap.vap1(), oilVap.vap2());
            }
            else {
                FluidSystem::setVapPars(0.0, 0.0);
            }

            if constexpr (enableConvectiveMixing) {
                ConvectiveMixingModule::beginEpisode(simulator.vanguard().eclState(), schedule, episodeIdx,
                                                     this->moduleParams_.convectiveMixingModuleParam);
            }
        }
    }

    //! \brief Seed the SummaryState from the parent run's summary vectors.
    //!
    //! \details Called at the start of every time step, before UDQ assignments
    //!          are evaluated, so that quantities the reduced run cannot
    //!          compute for itself are nevertheless available to UDQ, UDA and
    //!          group control through the standard code paths.
    //!
    //!          Values are assigned verbatim -- no delta arithmetic -- so the
    //!          state starts out as a clean snapshot of the parent at this
    //!          point in time. Everything the reduced run does compute is
    //!          subsequently overwritten by Summary::eval().
    //!
    //!          Cumulative (Total) quantities are deliberately NOT seeded.
    //!          SummaryState::update() accumulates into those keys rather than
    //!          assigning, so seeding a cumulative and then letting
    //!          Summary::eval() add this step's contribution on top would count
    //!          the same production twice. Cumulatives are instead built up by
    //!          the normal accumulation path from the rates, including the
    //!          rates injected for wells outside the region.
    void seedParentSummaryState_(const double start, const double time)
    {
        const auto* parent = this->fluxParentSummary_();
        if (parent == nullptr) {
            return;
        }

        auto& summaryState = this->simulator().vanguard().summaryState();

        for (const auto& key : parent->keys()) {
            if (ParentSummary::keyType(key) == SummaryConfigNode::Type::Total) {
                continue;
            }

            const auto value = parent->valueOver(key, start, time);
            if (!std::isfinite(value)) {
                continue;
            }

            ParentSummary::assign(summaryState, key, value);
        }
    }

    //! \brief Names of the wells that have at least one connection inside the
    //!        USEFLUX region.
    //! \details The consumer deck is the parent deck, so the schedule still
    //!          describes every well. Wells completed outside the region have
    //!          no active connection and take no part in the sector run.
    //!
    //!          A well belongs to the sector if ANY rank holds one of its
    //!          connections, which is why this is a collective and has to be
    //!          reached on every rank. Answering it rank by rank gives each a
    //!          different set of absent wells, and the group target correction
    //!          then subtracts a different amount on each; the ranks proceed to
    //!          disagree about how the group is controlled, and the first
    //!          collective that decision reaches no longer matches across them.
    std::vector<std::string> fluxLocalWellNames_(const std::size_t stepIdx) const
    {
        const auto& vanguard = this->simulator().vanguard();
        const auto& schedule = vanguard.schedule();

        const auto& wellNames = schedule.wellNames(stepIdx);

        auto connected = std::vector<int>(wellNames.size(), 0);
        std::transform(wellNames.begin(), wellNames.end(), connected.begin(),
                       [&vanguard, &schedule, stepIdx](const std::string& wellName)
                       {
                           const auto& well = schedule.getWell(wellName, stepIdx);

                           return std::any_of(well.getConnections().begin(),
                                              well.getConnections().end(),
                                              [&vanguard](const auto& conn)
                                              {
                                                  return vanguard.compressedIndex(conn.global_index()) >= 0;
                                              })
                               ? 1 : 0;
                       });

        // Logical OR across the ranks. The well name order comes from the
        // schedule, so every rank indexes the same wells.
        if (!connected.empty()) {
            vanguard.grid().comm().max(connected.data(), connected.size());
        }

        std::vector<std::string> localWells;
        for (std::size_t i = 0; i < wellNames.size(); ++i) {
            if (connected[i] != 0) {
                localWells.push_back(wellNames[i]);
            }
        }

        return localWells;
    }

    //! \brief Names of the wells that take no part in the sector run.
    std::vector<std::string> fluxAbsentWellNames_(const std::size_t stepIdx) const
    {
        const auto localWells = this->fluxLocalWellNames_(stepIdx);
        const auto& schedule = this->simulator().vanguard().schedule();

        std::vector<std::string> absentWells;
        for (const auto& wellName : schedule.wellNames(stepIdx)) {
            const auto isLocal = std::find(localWells.begin(), localWells.end(), wellName)
                != localWells.end();

            if (!isLocal) {
                absentWells.push_back(wellName);
            }
        }

        return absentWells;
    }

    //! \brief Hand the well model the rates of the wells outside the region.
    //!
    //! \details A group in the parent deck shares its target among all of its
    //!          members, but only the members inside the region are simulated
    //!          here. The others are still in the schedule, and what they
    //!          produced in the parent run is known, so they take part in group
    //!          control as wells on individual control at those rates: they
    //!          count towards every group sum and come off every group target,
    //!          and are never handed a share of one. Leaving them out entirely
    //!          lets the sector's own wells take over the whole target, and any
    //!          group limit, voidage replacement or sales target then sees a
    //!          fraction of what the parent's did.
    //!
    //!          The parent reports only a total reservoir volume per well. Every
    //!          reader of a producer's reservoir rates sums them over phases, so
    //!          it is booked against oil, or the first active phase without oil;
    //!          an injector's goes to the phase it injects.
    //!
    //!          Collective, because deciding which wells are absent is.
    void updateFluxFixedWellRates_(const int episodeIdx, const double start, const double time)
    {
        const auto* parent = this->fluxParentSummary_();
        if ((parent == nullptr) || (episodeIdx < 0)) {
            return;
        }

        using Helper = std::remove_reference_t<decltype(this->wellModel_.groupStateHelper())>;
        using FixedRates = typename Helper::FixedWellRates;

        const auto& schedule = this->simulator().vanguard().schedule();
        const auto& units = this->simulator().vanguard().eclState().getUnits();
        const auto stepIdx = static_cast<std::size_t>(episodeIdx);

        // The parent's summary is in the deck's output units; the well model
        // works in SI.
        const auto rate = [parent, start, time, &units](const std::string& key,
                                                        const UnitSystem::measure measure)
        {
            const auto value = parent->valueOver(key, start, time);
            return std::isfinite(value) ? units.to_si(measure, value) : 0.0;
        };

        const auto firstActive = [](std::initializer_list<unsigned> phases)
        {
            for (const auto phaseIdx : phases) {
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    return phaseIdx;
                }
            }
            return *phases.begin();
        };

        std::unordered_map<std::string, FixedRates> fixed;
        for (const auto& wellName : this->fluxAbsentWellNames_(stepIdx)) {
            const auto& well = schedule.getWell(wellName, stepIdx);
            const auto suffix = std::string{well.isInjector() ? "IR:" : "PR:"} + wellName;

            using M = UnitSystem::measure;

            FixedRates rates{};
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                rates.surface[FluidSystem::oilPhaseIdx] = rate("WO" + suffix, M::liquid_surface_rate);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                rates.surface[FluidSystem::waterPhaseIdx] = rate("WW" + suffix, M::liquid_surface_rate);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                rates.surface[FluidSystem::gasPhaseIdx] = rate("WG" + suffix, M::gas_surface_rate);
            }

            const auto bookedTo = [&well, &firstActive]()
            {
                if (!well.isInjector()) {
                    return firstActive({FluidSystem::oilPhaseIdx,
                                        FluidSystem::waterPhaseIdx,
                                        FluidSystem::gasPhaseIdx});
                }

                switch (well.injectorType()) {
                case InjectorType::GAS:
                    return static_cast<unsigned>(FluidSystem::gasPhaseIdx);
                case InjectorType::OIL:
                    return static_cast<unsigned>(FluidSystem::oilPhaseIdx);
                default:
                    return static_cast<unsigned>(FluidSystem::waterPhaseIdx);
                }
            }();

            rates.reservoir[bookedTo] = rate("WV" + suffix, M::rate);

            fixed.emplace(wellName, rates);
        }

        this->wellModel_.groupStateHelper().setFixedWellRates(std::move(fixed));
    }

    //! \brief Put back this run's own value of every UDQ that sets a group
    //!        target, over the one the parent summary seeded.
    //!
    //! \details The seed is taken at the end of the time step, which is right
    //!          for a rate -- the parent's rate over the step -- but not for a
    //!          UDQ read as a target: the parent's group control ran on the
    //!          value from the start of the step, and a DEFINE evaluated once
    //!          per step has moved on by the end of it. The UDQ state holds the
    //!          start-of-step value, as evaluated here.
    //!
    //!          Only the I/O rank has it, though. evalSummary() evaluates the
    //!          UDQs there and shares only the summary state afterwards, which
    //!          leaves every other rank with whatever it was handed when the
    //!          deck was distributed, and the group would then be driven to a
    //!          different target on each. Hence the broadcast, which is why this
    //!          has to be reached on every rank.
    void restoreFluxGroupTargetUDQs_(const int episodeIdx)
    {
        if ((this->fluxParentSummary_() == nullptr) || (episodeIdx < 0)) {
            return;
        }

        auto& vanguard = this->simulator().vanguard();
        const auto& schedule = vanguard.schedule();
        const auto stepIdx = static_cast<std::size_t>(episodeIdx);

        // In the order the groups come out of the schedule, so that every rank
        // builds the same list.
        std::vector<std::string> udqNames;
        std::unordered_set<std::string> seen;
        const auto note = [&udqNames, &seen](const UDAValue& target)
        {
            if (target.is<std::string>() && seen.insert(target.get<std::string>()).second) {
                udqNames.push_back(target.get<std::string>());
            }
        };

        for (const auto& groupName : schedule.groupNames(stepIdx)) {
            const auto& group = schedule.getGroup(groupName, stepIdx);

            if (group.isProductionGroup()) {
                const auto& production = group.productionProperties();
                note(production.oil_target);
                note(production.water_target);
                note(production.gas_target);
                note(production.liquid_target);
                note(production.resv_target);
            }

            if (group.isInjectionGroup()) {
                for (const auto& [phase, injection] : group.injectionProperties()) {
                    note(injection.surface_max_rate);
                    note(injection.resv_max_rate);
                    note(injection.target_reinj_fraction);
                    note(injection.target_void_fraction);
                }
            }
        }

        if (udqNames.empty()) {
            return;
        }

        const auto& comm = vanguard.grid().comm();
        auto values = std::vector<double>(udqNames.size(),
                                          std::numeric_limits<double>::quiet_NaN());

        if (comm.rank() == 0) {
            const auto& udqState = vanguard.udqState();
            std::transform(udqNames.begin(), udqNames.end(), values.begin(),
                           [&udqState](const std::string& name)
                           {
                               return udqState.has(name)
                                   ? udqState.get(name)
                                   : std::numeric_limits<double>::quiet_NaN();
                           });
        }

        comm.broadcast(values.data(), static_cast<int>(values.size()), 0);

        auto& summaryState = vanguard.summaryState();
        for (std::size_t i = 0; i < udqNames.size(); ++i) {
            if (std::isfinite(values[i])) {
                summaryState.set(udqNames[i], values[i]);
            }
        }
    }

    /*!
     * \brief Called by the simulator before each time integration.
     */
    void beginTimeStep() override
    {
        // The boundary record for this step is the one covering it, which is
        // found from the step's start: a record ending exactly there describes
        // the step before. And it is chosen ahead of the base class, which
        // caches the boundary state and mass rates from whatever record is
        // current and starts the well model's step, whose rate converter reads
        // the record's sums for the cells outside the sector.
        this->refreshFluxBoundaryRecord_(static_cast<double>(this->simulator().time()),
                                         this->fluxStepEndTime_());

        // The parent's rates are averages over its own steps, which need not
        // be this run's, so everything taken from its summary for this step is
        // taken over the step: a rate as its average across (start, end], all
        // else at the end.
        const auto stepStart = static_cast<double>(this->simulator().time());
        const auto time = this->fluxStepEndTime_();

        // Ahead of the base class, which starts the well model's time step and
        // with it the group controls. Refreshed per time step, not per report
        // step.
        this->updateFluxFixedWellRates_(this->simulator().episodeIndex(), stepStart, time);

        FlowProblemType::beginTimeStep();

        // Seed the summary state from the parent run before anything reads it.
        // UDA group targets may be driven by UDQs that evolve on every time
        // step, so this has to follow the time step cadence rather than being
        // applied once per report step.
        this->seedParentSummaryState_(stepStart, time);
        this->restoreFluxGroupTargetUDQs_(this->simulator().episodeIndex());

        hybridNewton_.tryApplyHybridNewton();
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void registerAuxiliaryCellModules()
    {
        const auto& eclState = this->simulator().vanguard().eclState();

        if ((eclState.numericalAquiferMode() == NumericalAquiferMode::AuxiliaryCells) &&
            eclState.aquifer().hasNumericalAquifer())
        {
            const auto& aquifers = this->registerAuxCellModule_(
                std::make_unique<NumericalAquiferAuxCells<TypeTag>>(this->simulator()));

            if (this->simulator().gridView().comm().rank() == 0) {
                OpmLog::info(fmt::format("Numerical aquifers represented as {} auxiliary "
                                         "cells outside the grid", aquifers.numDofs()));
            }
        }
    }

    void finishInit()
    {
        FlowProblemType::finishInit();

        auto& simulator = this->simulator();
        const bool transmissibilitiesFinished =
            this->prepareTransmissibilityOutput_(*eclWriter_, enableEclOutput_);

        const auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();
        this->initializeSimulatorTime_();

        // conserve inner energy instead of enthalpy if TEMP is used
        // or THERMAL and parameter ConserveInnerEnergyThermal is true (default false)
        bool isThermal = eclState.getSimulationConfig().isThermal();
        bool isTemp = eclState.getSimulationConfig().isTemp();
        bool conserveInnerEnergy = isTemp || (isThermal && Parameters::Get<Parameters::ConserveInnerEnergyThermal>());
        FluidSystem::setEnergyEqualEnthalpy(conserveInnerEnergy);

        this->initializeModelProperties_();

        // write the static output files (EGRID, INIT)
        if (enableEclOutput_) {
            this->eclWriter_->writeInit();
        }

        if (!transmissibilitiesFinished) {
            this->finishTransmissibilities_();
        }

        this->applyAuxCellTransmissibilities_();

        const auto& initconfig = eclState.getInitConfig();

        // The restart file is per grid element; auxiliary state would come back undefined.
        if (initconfig.restartRequested() && !this->auxCellModules_.empty()) {
            OPM_THROW(std::runtime_error,
                      "Restart is not supported together with auxiliary cells "
                      "(numerical aquifers represented outside the grid): their state "
                      "is not written to the restart file.");
        }

        this->tracerModel_.init(initconfig.restartRequested());
        if (initconfig.restartRequested()) {
            this->readEclRestartSolution_();
        }
        else {
            this->readInitialCondition_();
        }
        this->temperatureModel_.init();
        this->tracerModel_.prepareTracerBatches();

        this->updatePffDofData_();

        if constexpr (getPropValue<TypeTag, Properties::EnablePolymer>()) {
            const auto& vanguard = this->simulator().vanguard();
            const auto& gridView = vanguard.gridView();
            const int numElements = gridView.size(/*codim=*/0);
            this->polymer_.maxAdsorption.resize(numElements, 0.0);
        }

        this->readBoundaryConditions_();

        // compute and set eq weights based on initial b values
        this->computeAndSetEqWeights_();

        if (this->enableDriftCompensation_ || this->enableDriftCompensationTemp_) {
            this->drift_.resize(this->model().numTotalDof());
            this->drift_ = 0.0;
        }

        // after finishing the initialization and writing the initial solution, we move
        // to the first "real" episode/report step
        // for restart the episode index and start is already set
        if (!initconfig.restartRequested() && !eclState.getIOConfig().initOnly()) {
            simulator.startNextEpisode(schedule.seconds(1));
            simulator.setEpisodeIndex(0);
            simulator.setTimeStepIndex(0);
        }

        if (Parameters::Get<Parameters::CheckSatfuncConsistency>() &&
            ! this->satfuncConsistencyRequirementsMet())
        {
            // User requested that saturation functions be checked for
            // consistency and essential/critical requirements are not met.
            // Abort simulation run.
            //
            // Note: We need synchronisation here lest ranks other than the
            // I/O rank throw exceptions too early thereby risking an
            // incomplete failure report being shown to the user.
            this->simulator().vanguard().grid().comm().barrier();

            throw std::domain_error {
                "Saturation function end-points do not "
                "meet requisite consistency conditions"
            };
        }

        // TODO: move to the end for later refactoring of the function finishInit()
        //
        // deal with DRSDT
        this->mixControls_.init(this->model().numTotalDof(),
                                this->episodeIndex(),
                                eclState.runspec().tabdims().getNumPVTTables());

        // Seed the DRSDT/DRVDT history from the initial composition.  Left at
        // zero, the limiter caps Rs at DRSDT * dt on the first step and boils
        // the dissolved gas out of an undersaturated reservoir.  A restarted
        // run seeds it from the restart solution instead.
        if (!initconfig.restartRequested()) {
            for (std::size_t elemIdx = 0; elemIdx < this->initialFluidStates_.size(); ++elemIdx) {
                const auto& fs = this->initialFluidStates_[elemIdx];
                this->mixControls_.updateLastValues(elemIdx, fs.Rs(), fs.Rv());
            }
        }

        if (this->enableVtkOutput_() && eclState.getIOConfig().initOnly()) {
            simulator.setTimeStepSize(0.0);
            simulator.model().applyInitialSolution();
            FlowProblemType::writeOutput(true);
        }

        if (!eclState.getIOConfig().initOnly()) {
            if (!this->enableTuning_ && eclState.getSimulationConfig().anyTUNING()) {
                OpmLog::info("\nThe deck has TUNING in the SCHEDULE section, but "
                             "it is ignored due\nto the flag --enable-tuning=false. "
                             "Set this flag to true to activate it.\n"
                             "Manually tuning the simulator with the TUNING keyword may "
                             "increase run time.\nIt is recommended using the simulator's "
                             "default tuning (--enable-tuning=false).");
            }
        }
    }

    /*!
     * \brief Called by the simulator after each time integration.
     */
    void endTimeStep() override
    {
        FlowProblemType::endTimeStep();
        this->endStepApplyAction();
    }

    void endStepApplyAction()
    {
        // Sample the DUMPFLUX boundary rates while the linearizer still holds
        // the converged fluxes of the time step that just finished.
        this->eclWriter_->sampleFluxDumperRates(this->simulator().timeStepSize());

        // After the solution is updated, the values in output module needs
        // also updated.
        this->eclWriter().mutableOutputModule().invalidateLocalData();

        // Evaluate the summary state for CpGrid with LGRs as well.  In serial
        // the CollectDataOnIORank maps are the identity and the leaf-grid data
        // is written directly; in parallel the name-keyed summary data (wells,
        // groups, region values) is gathered to the I/O rank while cell-based
        // restart/block output on refined grids is still being wired up.
        this->eclWriter_->evalSummaryState(!this->episodeWillBeOver());

        {
            OPM_TIMEBLOCK(applyActions);

            const int episodeIdx = this->episodeIndex();
            auto& simulator = this->simulator();

            // Clear out any existing events as these have already been
            // processed when we're running an action block
            this->simulator().vanguard().schedule().clearEvents(episodeIdx);

            // Re-ordering in case of Alugrid
            this->actionHandler_
                .applyActions(episodeIdx, simulator.time() + simulator.timeStepSize(),
                              [this](const bool global)
            {
                using TransUpdateQuantities = typename
                    Vanguard::TransmissibilityType::TransUpdateQuantities;

                this->transmissibilities_
                    .update(global, TransUpdateQuantities::All,
                            [&vg = this->simulator().vanguard()]
                            (const unsigned int i)
                    {
                        return vg.gridIdxToEquilGridIdx(i);
                    });
            });
        }
    }

    /*!
     * \brief Called by the simulator after the end of an episode.
     */
    void endEpisode() override
    {
        OPM_TIMEBLOCK(endEpisode);

        // Rerun UDQ assignents following action processing on the final
        // time step of this episode to make sure that any UDQ ASSIGN
        // operations triggered in action blocks take effect.  This is
        // mainly to work around a shortcoming of the ScheduleState copy
        // constructor which clears pending UDQ assignments under the
        // assumption that all such assignments have been processed.  If an
        // action block happens to trigger on the final time step of an
        // episode and that action block runs a UDQ assignment, then that
        // assignment would be dropped and the rest of the simulator will
        // never see its effect without this hack.
        this->actionHandler_
            .evalUDQAssignments(this->episodeIndex(), this->simulator().vanguard().udqState());

        FlowProblemType::endEpisode();
    }

    void writeReports(const SimulatorTimer& timer)
    {
        if (this->enableEclOutput_) {
            this->eclWriter_->writeReports(timer);
        }
    }


    /*!
     * \brief Write the requested quantities of the current solution into the output
     *        files.
     */
    void writeOutput(const bool verbose) override
    {
        FlowProblemType::writeOutput(verbose);

        const auto isSubStep = !this->episodeWillBeOver();

        auto localCellData = data::Solution {};

#if HAVE_DAMARIS
        // N.B. the Damaris output has to be done before the ECL output as the ECL one
        // does all kinds of std::move() relocation of data
        if (this->enableDamarisOutput_ && (this->damarisWriter_ != nullptr)) {
            this->damarisWriter_->writeOutput(localCellData, isSubStep);
        }
#endif

        if (this->enableEclOutput_ && (this->eclWriter_ != nullptr)) {
            this->eclWriter_->writeOutput(std::move(localCellData), isSubStep,
                                          this->simulator().vanguard().schedule()
                                         .exitStatus().has_value());
        }
    }

    void finalizeOutput()
    {
        OPM_TIMEBLOCK(finalizeOutput);
        // this will write all pending output to disk
        // to avoid corruption of output files
        eclWriter_.reset();
    }


    /*!
     * \copydoc FvBaseProblem::initialSolutionApplied()
     */
    void initialSolutionApplied() override
    {
        FlowProblemType::initialSolutionApplied();

        // let the object for threshold pressures initialize itself. this is done only at
        // this point, because determining the threshold pressures may require to access
        // the initial solution.
        this->thresholdPressures_.finishInit();
        this->adoptFluxThresholdPressures_();

        // Compute the initial FIP report (also for CpGrid with LGRs) so the
        // in-place reference values are available for the summary balance.
        if (this->simulator().episodeIndex() == 0) {
            eclWriter_->writeInitialFIPReport();
        }
    }

    void addToSourceDense(RateVector& rate,
                          unsigned globalDofIdx,
                          unsigned timeIdx) const override
    {
        this->aquiferModel_.addToSource(rate, globalDofIdx, timeIdx);

        this->addFluxNncSource_(rate, globalDofIdx, timeIdx);

        // Add source term from deck
        const auto& source = this->simulator().vanguard().schedule()[this->episodeIndex()].source();
        std::array<int,3> ijk{};
        // Auxiliary DOFs have no Cartesian position and no deck source.
        const bool isGridDof = globalDofIdx < this->model().numGridDof();
        if (isGridDof) {
            this->simulator().vanguard().cartesianCoordinate(globalDofIdx, ijk);
        }

        if (isGridDof && source.hasSource(ijk)) {
            const int pvtRegionIdx = this->pvtRegionIndex(globalDofIdx);
            static std::array<SourceComponent, 3> sc_map = {SourceComponent::WATER, SourceComponent::OIL, SourceComponent::GAS};
            static std::array<int, 3> phidx_map = {FluidSystem::waterPhaseIdx, FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx};
            static std::array<int, 3> cidx_map = {waterCompIdx, oilCompIdx, gasCompIdx};

            for (unsigned i = 0; i < phidx_map.size(); ++i) {
                const auto phaseIdx = phidx_map[i];
                const auto sourceComp = sc_map[i];
                const auto compIdx = cidx_map[i];
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }
                Scalar mass_rate = source.rate(ijk, sourceComp) / this->model().dofTotalVolume(globalDofIdx);
                if constexpr (getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>()) {
                    mass_rate /= FluidSystem::referenceDensity(phaseIdx, pvtRegionIdx);
                }
                rate[FluidSystem::canonicalToActiveCompIdx(compIdx)] += mass_rate;
            }

            if constexpr (enableSolvent) {
                Scalar mass_rate = source.rate(ijk, SourceComponent::SOLVENT) / this->model().dofTotalVolume(globalDofIdx);
                if constexpr (getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>()) {
                    const auto& solventPvt = SolventModule::solventPvt();
                    mass_rate /= solventPvt.referenceDensity(pvtRegionIdx);
                }
                rate[Indices::contiSolventEqIdx] += mass_rate;
            }
            if constexpr (enablePolymer) {
                rate[Indices::polymerConcentrationIdx] += source.rate(ijk, SourceComponent::POLYMER) / this->model().dofTotalVolume(globalDofIdx);
            }
            if constexpr (enableMICP) {
                rate[Indices::microbialConcentrationIdx] += source.rate(ijk, SourceComponent::MICR) / this->model().dofTotalVolume(globalDofIdx);
                rate[Indices::oxygenConcentrationIdx] += source.rate(ijk, SourceComponent::OXYG) / this->model().dofTotalVolume(globalDofIdx);
                rate[Indices::ureaConcentrationIdx] += source.rate(ijk, SourceComponent::UREA) / (this->model().dofTotalVolume(globalDofIdx));
            }
            if constexpr (energyModuleType == EnergyModules::FullyImplicitThermal) {
                for (unsigned i = 0; i < phidx_map.size(); ++i) {
                    const auto phaseIdx = phidx_map[i];
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }
                    const auto sourceComp = sc_map[i];
                    const auto source_hrate = source.hrate(ijk, sourceComp);
                    if (source_hrate) {
                        rate[Indices::contiEnergyEqIdx] += source_hrate.value() / this->model().dofTotalVolume(globalDofIdx);
                    } else {
                        const auto& intQuants = this->simulator().model().intensiveQuantities(globalDofIdx, /*timeIdx*/ 0);
                        auto fs = intQuants.fluidState();
                        // if temperature is not set, use cell temperature as default
                        const auto source_temp = source.temperature(ijk, sourceComp);
                        if (source_temp) {
                            Scalar temperature = source_temp.value();
                            fs.setTemperature(temperature);
                        }
                        const auto& h = FluidSystem::enthalpy(fs, phaseIdx, pvtRegionIdx);
                        Scalar mass_rate = source.rate(ijk, sourceComp)/ this->model().dofTotalVolume(globalDofIdx);
                        Scalar energy_rate = getValue(h)*mass_rate;
                        rate[Indices::contiEnergyEqIdx] += energy_rate;
                    }
                }
            }
        }

        // if requested, compensate systematic mass loss for cells which were "well
        // behaved" in the last time step
        if (this->enableDriftCompensation_) {
            const auto& simulator = this->simulator();
            const auto& model = this->model();

            // we use a lower tolerance for the compensation too
            // assure the added drift from the last step does not
            // cause convergence issues on the current step
            Scalar maxCompensation = model.newtonMethod().tolerance()/10;
            Scalar poro = this->porosity(globalDofIdx, timeIdx);
            Scalar dt = simulator.timeStepSize();
            EqVector dofDriftRate = this->drift_[globalDofIdx];
            dofDriftRate /= dt*model.dofTotalVolume(globalDofIdx);

            // restrict drift compensation to the CNV tolerance
            for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
                Scalar cnv = std::abs(dofDriftRate[eqIdx])*dt*model.eqWeight(globalDofIdx, eqIdx)/poro;
                if (cnv > maxCompensation) {
                    dofDriftRate[eqIdx] *= maxCompensation/cnv;
                }
            }

            for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx)
                rate[eqIdx] -= dofDriftRate[eqIdx];
        }
    }

    /*!
     * \brief Calculate the transmissibility multiplier due to porosity reduction.
     */
    template <class LhsEval, class Callback>
    LhsEval permFactTransMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx, Callback& obtain) const
    {
        OPM_TIMEBLOCK_LOCAL(permFactTransMultiplier, Subsystem::PvtProps);
        if constexpr (enableSaltPrecipitation) {
            const auto& fs = intQuants.fluidState();
            unsigned tableIdx = this->simulator().problem().satnumRegionIndex(elementIdx);
            LhsEval porosityFactor = obtain(1. - fs.saltSaturation());
            porosityFactor = min(porosityFactor, 1.0);
            const auto& permfactTable = BrineModule::permfactTable(tableIdx);
            return permfactTable.eval(porosityFactor, /*extrapolation=*/true);
        }
        else if constexpr (enableBioeffects) {
            return obtain(intQuants.permFactor());
        }
        else {
            return 1.0;
        }
    }

    // temporary solution to facilitate output of initial state from flow
    const InitialFluidState& initialFluidState(unsigned globalDofIdx) const
    { return initialFluidStates_[globalDofIdx]; }

    std::vector<InitialFluidState>& initialFluidStates()
    { return initialFluidStates_; }

    const std::vector<InitialFluidState>& initialFluidStates() const
    { return initialFluidStates_; }

    const EclipseIO& eclIO() const
    { return eclWriter_->eclIO(); }

    void setSubStepReport(const SimulatorReportSingle& report)
    { return eclWriter_->setSubStepReport(report); }

    void setSimulationReport(const SimulatorReport& report)
    { return eclWriter_->setSimulationReport(report); }

    //! \brief The state of the parent's cell on the far side of a sector
    //!        boundary face, as recorded in the .FLUX file.
    //!
    //! \details Pressures come carried from that cell's centre to \p targetDepth,
    //!   using its own densities; pass NaN to leave them at the centre.
    InitialFluidState fluxExteriorFluidState_(const EclIO::FluxFile::Data& fluxData,
                                              const EclIO::FluxFile::ReportStep& fluxStep,
                                              const FluxBoundary::Face& fluxFace,
                                              const std::size_t fluxFaceIndex,
                                              const unsigned globalDofIdx,
                                              const Scalar targetDepth) const
    {
        InitialFluidState fluidState;

        // The stream on the far side belongs to the parent's cell, not
        // to this one, so its PVT region is the one the producing run
        // recorded. Falling back to the interior region would evaluate
        // the density, formation volume factor and viscosity of the
        // inflow against the wrong tables.
        const int pvtRegionIdx = (fluxFace.exteriorPvtRegion >= 0)
            ? fluxFace.exteriorPvtRegion
            : this->pvtRegionIndex(globalDofIdx);
        fluidState.setPvtRegionIndex(pvtRegionIdx);

        const auto& initialState = initialFluidStates_[globalDofIdx];
        double sw = FluidSystem::phaseIsActive(waterPhaseIdx)
                  ? initialState.saturation(waterPhaseIdx)
                  : 0.0;
        double sg = FluidSystem::phaseIsActive(gasPhaseIdx)
                  ? initialState.saturation(gasPhaseIdx)
                  : 0.0;

        if (fluxFaceIndex < fluxStep.swat.size()) {
            sw = fluxStep.swat[fluxFaceIndex];
        }
        if (fluxFaceIndex < fluxStep.sgas.size()) {
            sg = fluxStep.sgas[fluxFaceIndex];
        }

        sw = std::clamp(sw, 0.0, 1.0);
        sg = std::clamp(sg, 0.0, 1.0);
        const double so = std::clamp(1.0 - sw - sg, 0.0, 1.0);

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            fluidState.setSaturation(FluidSystem::waterPhaseIdx, sw);
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            fluidState.setSaturation(FluidSystem::gasPhaseIdx, sg);
        }
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            fluidState.setSaturation(FluidSystem::oilPhaseIdx, so);
        }
        fluidState.setTotalSaturation(1.0);

        const double pressure = fluxStep.pressures[fluxFaceIndex];

        // Capillary pressure of the exterior cell, relative to the
        // reference phase. The producing run recorded it because
        // recomputing it here would use this cell's saturation-function
        // region, scaled end points and hysteresis state, none of which
        // describe the cell the pressure came from.
        std::array<Scalar, numPhases> pc = {0};
        const auto numPhaseSlots =
            static_cast<std::size_t>(fluxData.header.numPhases);
        const auto pcBase = fluxFaceIndex * numPhaseSlots;
        if (fluxStep.capPressure.size() >= pcBase + numPhaseSlots) {
            std::size_t slot = 0;
            const auto load = [&](const unsigned phaseIdx)
            {
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    pc[phaseIdx] = fluxStep.capPressure[pcBase + slot];
                    ++slot;
                }
            };
            load(oilPhaseIdx);
            load(waterPhaseIdx);
            load(gasPhaseIdx);
        }
        else {
            const auto& matParams = this->materialLawParams(globalDofIdx);
            MaterialLaw::capillaryPressures(pc, matParams, fluidState);
        }

        Valgrind::CheckDefined(pressure);
        Valgrind::CheckDefined(pc);
        for (unsigned activePhaseIdx = 0; activePhaseIdx < FluidSystem::numActivePhases(); ++activePhaseIdx) {
            const auto phaseIdx = FluidSystem::activeToCanonicalPhaseIdx(activePhaseIdx);
            if (Indices::oilEnabled)
                fluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
            else if (Indices::gasEnabled)
                fluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[gasPhaseIdx]));
            else if (Indices::waterEnabled)
                fluidState.setPressure(phaseIdx, pressure);
        }

        if constexpr (energyModuleType != EnergyModules::NoTemperature) {
            double temperature = initialState.temperature(0);
            if (fluxFaceIndex < fluxStep.temperature.size()) {
                temperature = fluxStep.temperature[fluxFaceIndex];
            }
            fluidState.setTemperature(temperature);
        }

        if constexpr (enableDissolvedGas) {
            if (FluidSystem::enableDissolvedGas()) {
                if (fluxFaceIndex < fluxStep.rs.size()) {
                    fluidState.setRs(fluxStep.rs[fluxFaceIndex]);
                }
                else {
                    fluidState.setRs(0.0);
                }
                if (fluxFaceIndex < fluxStep.rv.size()) {
                    fluidState.setRv(fluxStep.rv[fluxFaceIndex]);
                }
                else {
                    fluidState.setRv(0.0);
                }
            }
        }
        if constexpr (enableDisgasInWater) {
            if (FluidSystem::enableDissolvedGasInWater()) {
                fluidState.setRsw(0.0);
            }
        }
        if constexpr (enableVapwat) {
            if (FluidSystem::enableVaporizedWater()) {
                fluidState.setRvw(0.0);
            }
        }

        // The pressures above are the exterior cell's, taken at its
        // centre, but a boundary condition is imposed at the face and
        // the discretisation carries it from there to this cell's
        // centre. The leg from the exterior cell up to the face is
        // nobody else's to walk, and skipping it imposes the exterior
        // cell's pressure at a depth where the exterior cell does not
        // have it. In a dipping layer that is metres of head: on one
        // field it left an equilibrated sector with no wells open
        // gaining a quarter of a bar everywhere in five days, and more
        // than a bar against the steeper parts of the boundary.
        //
        // The correction wants densities, which are not known until the
        // pressure is settled, so take them at the uncorrected pressure
        // and let the loop below recompute them once it has moved. Over
        // a few metres of head the second pass changes nothing that
        // matters.
        if (std::isfinite(fluxFace.exteriorDepth) && std::isfinite(targetDepth)) {
            const auto dz = targetDepth - fluxFace.exteriorDepth;

            if (dz != 0.0) {
                const auto g = this->gravity()[dimWorld - 1];

                for (unsigned activePhaseIdx = 0;
                     activePhaseIdx < FluidSystem::numActivePhases(); ++activePhaseIdx)
                {
                    const auto phaseIdx =
                        FluidSystem::activeToCanonicalPhaseIdx(activePhaseIdx);

                    const auto rho =
                        FluidSystem::density(fluidState, phaseIdx, pvtRegionIdx);

                    fluidState.setPressure(phaseIdx,
                                           fluidState.pressure(phaseIdx) + rho * g * dz);
                }
            }
        }

        for (unsigned activePhaseIdx = 0; activePhaseIdx < FluidSystem::numActivePhases(); ++activePhaseIdx) {
            const auto phaseIdx = FluidSystem::activeToCanonicalPhaseIdx(activePhaseIdx);

            const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState, phaseIdx, pvtRegionIdx);
            fluidState.setInvB(phaseIdx, b);

            const auto& rho = FluidSystem::density(fluidState, phaseIdx, pvtRegionIdx);
            fluidState.setDensity(phaseIdx, rho);
            if constexpr (energyModuleType == EnergyModules::FullyImplicitThermal) {
                const auto& h = FluidSystem::enthalpy(fluidState, phaseIdx, pvtRegionIdx);
                fluidState.setEnthalpy(phaseIdx, h);
            }
        }

        fluidState.checkDefined();
        return fluidState;
    }

    InitialFluidState boundaryFluidState(unsigned globalDofIdx, const int directionId) const
    {
        OPM_TIMEBLOCK_LOCAL(boundaryFluidState, Subsystem::Assembly);
        const FaceDir::DirEnum dir = FaceDir::FromIntersectionIndex(directionId);

        const auto* fluxData = this->fluxBoundaryData_();
        const auto* fluxFace = this->fluxBoundaryFace_(globalDofIdx, dir);
        const auto* fluxStep = this->fluxBoundaryReportStep_();
        if (fluxData && fluxFace && fluxStep
            && (fluxData->header.mode == EclIO::FluxFile::Mode::Pressure
                || fluxData->header.mode == EclIO::FluxFile::Mode::Both)) {
            const auto fluxSlot = this->fluxBoundaryFaceSlot_(globalDofIdx, dir);
            const auto fluxFaceIndex = static_cast<std::size_t>(fluxSlot - 1);
            if (fluxSlot > 0 && fluxFaceIndex < fluxStep->pressures.size()) {
                return this->fluxExteriorFluidState_(*fluxData, *fluxStep, *fluxFace,
                                                     fluxFaceIndex, globalDofIdx,
                                                     this->fluxBoundaryFaceDepthAt_(globalDofIdx, dir));
            }
        }

        const auto& bcstate = this->simulator().vanguard().schedule()[this->episodeIndex()].bcstate;
        if (bcstate.size() > 0) {
            // index == 0: no boundary conditions for this
            // global cell and direction
            if (this->bcindex_(dir)[globalDofIdx] == 0)
                return initialFluidStates_[globalDofIdx];

            const auto& bc = bcstate[this->bcindex_(dir)[globalDofIdx]];
            if (bc.bctype == BCType::DIRICHLET )
            {
                InitialFluidState fluidState;
                const int pvtRegionIdx = this->pvtRegionIndex(globalDofIdx);
                fluidState.setPvtRegionIndex(pvtRegionIdx);

                switch (bc.component) {
                    case BCComponent::OIL:
                        if (!FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
                            throw std::logic_error("oil is not active and you're trying to add oil BC");

                        fluidState.setSaturation(FluidSystem::oilPhaseIdx, 1.0);
                        break;
                    case BCComponent::GAS:
                        if (!FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
                            throw std::logic_error("gas is not active and you're trying to add gas BC");

                        fluidState.setSaturation(FluidSystem::gasPhaseIdx, 1.0);
                        break;
                    case BCComponent::WATER:
                        if (!FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))
                            throw std::logic_error("water is not active and you're trying to add water BC");

                        fluidState.setSaturation(FluidSystem::waterPhaseIdx, 1.0);
                        break;
                    case BCComponent::SOLVENT:
                    case BCComponent::POLYMER:
                    case BCComponent::MICR:
                    case BCComponent::OXYG:
                    case BCComponent::UREA:
                    case BCComponent::NONE:
                        throw std::logic_error("you need to specify a valid component (OIL, WATER or GAS) when DIRICHLET type is set in BC");
                }
                fluidState.setTotalSaturation(1.0);
                double pressure = initialFluidStates_[globalDofIdx].pressure(this->refPressurePhaseIdx_());
                const auto pressure_input = bc.pressure;
                if (pressure_input) {
                    pressure = *pressure_input;
                }

                std::array<Scalar, numPhases> pc = {0};
                const auto& matParams = this->materialLawParams(globalDofIdx);
                MaterialLaw::capillaryPressures(pc, matParams, fluidState);
                Valgrind::CheckDefined(pressure);
                Valgrind::CheckDefined(pc);
                for (unsigned activePhaseIdx = 0; activePhaseIdx < FluidSystem::numActivePhases(); ++activePhaseIdx) {
                    const auto phaseIdx = FluidSystem::activeToCanonicalPhaseIdx(activePhaseIdx);
                    if (Indices::oilEnabled)
                        fluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
                    else if (Indices::gasEnabled)
                        fluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[gasPhaseIdx]));
                    else if (Indices::waterEnabled)
                        //single (water) phase
                        fluidState.setPressure(phaseIdx, pressure);
                }
                if constexpr (energyModuleType != EnergyModules::NoTemperature) {
                    double temperature = initialFluidStates_[globalDofIdx].temperature(0); // we only have one temperature
                    const auto temperature_input = bc.temperature;
                    if(temperature_input)
                        temperature = *temperature_input;
                    fluidState.setTemperature(temperature);
                }

                if constexpr (enableDissolvedGas) {
                    if (FluidSystem::enableDissolvedGas()) {
                        fluidState.setRs(0.0);
                        fluidState.setRv(0.0);
                    }
                }
                if constexpr (enableDisgasInWater) {
                    if (FluidSystem::enableDissolvedGasInWater()) {
                        fluidState.setRsw(0.0);
                    }
                }
                if constexpr (enableVapwat) {
                    if (FluidSystem::enableVaporizedWater()) {
                        fluidState.setRvw(0.0);
                    }
                }

                for (unsigned activePhaseIdx = 0; activePhaseIdx < FluidSystem::numActivePhases(); ++activePhaseIdx) {
                    const auto phaseIdx = FluidSystem::activeToCanonicalPhaseIdx(activePhaseIdx);

                    const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState, phaseIdx, pvtRegionIdx);
                    fluidState.setInvB(phaseIdx, b);

                    const auto& rho = FluidSystem::density(fluidState, phaseIdx, pvtRegionIdx);
                    fluidState.setDensity(phaseIdx, rho);
                    if constexpr (energyModuleType == EnergyModules::FullyImplicitThermal) {
                        const auto& h = FluidSystem::enthalpy(fluidState, phaseIdx, pvtRegionIdx);
                        fluidState.setEnthalpy(phaseIdx, h);
                    }
                }
                fluidState.checkDefined();
                return fluidState;
            }
        }
        return initialFluidStates_[globalDofIdx];
    }


    const EclWriterType& eclWriter() const
    { return *eclWriter_; }

    EclWriterType& eclWriter()
    { return *eclWriter_; }

    /*!
     * \brief Returns the maximum value of the gas dissolution factor at the current time
     *        for a given degree of freedom.
     */
    Scalar maxGasDissolutionFactor(unsigned timeIdx, unsigned globalDofIdx) const
    {
        return this->mixControls_.maxGasDissolutionFactor(timeIdx, globalDofIdx,
                                                          this->episodeIndex(),
                                                          this->pvtRegionIndex(globalDofIdx));
    }

    /*!
     * \brief Returns the maximum value of the oil vaporization factor at the current
     *        time for a given degree of freedom.
     */
    Scalar maxOilVaporizationFactor(unsigned timeIdx, unsigned globalDofIdx) const
    {
        return this->mixControls_.maxOilVaporizationFactor(timeIdx, globalDofIdx,
                                                           this->episodeIndex(),
                                                           this->pvtRegionIndex(globalDofIdx));
    }

    /*!
     * \brief Return if the storage term of the first iteration is identical to the storage
     *        term for the solution of the previous time step.
     *
     * Storage recycling is disabled for configurations where the first iteration storage
     * can differ from the previous timestep storage:
     * 1. DRSDT/DRVDT update dissolved/vaporized composition limits explicitly.
     * 2. Rock compaction multipliers make pore volume state-dependent across timesteps.
     * 3. TPSA geomechanics can update mechanics state between coupled Flow/TPSA solves,
     *    which changes porosity/pore-volume terms used by Flow.
     */
    bool recycleFirstIterationStorage() const
    {
        const auto& rspec = this->simulator().vanguard().eclState().runspec();
        const bool tpsaActive = rspec.mech() && rspec.mechSolver().tpsa();
        if (tpsaActive) {
            return false;
        }

        int episodeIdx = this->episodeIndex();
        return !this->mixControls_.drsdtActive(episodeIdx) &&
               !this->mixControls_.drvdtActive(episodeIdx) &&
               this->rockCompPoroMultWc_.empty() &&
               this->rockCompPoroMult_.empty();
    }

    /*!
     * \copydoc FvBaseProblem::initial
     *
     * The reservoir problem uses a constant boundary condition for
     * the whole domain.
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        values.setPvtRegionIndex(pvtRegionIndex(context, spaceIdx, timeIdx));
        values.assignNaive(initialFluidStates_[globalDofIdx]);

        if constexpr (enableSolvent) {
            SolventModule::assignPrimaryVars(values,
                                             this->solventSaturation_[globalDofIdx],
                                             this->solventRsw_[globalDofIdx]);
        }

        if constexpr (enablePolymer) {
            values[Indices::polymerConcentrationIdx] = this->polymer_.concentration[globalDofIdx];
        }

        if constexpr (enablePolymerMolarWeight) {
            values[Indices::polymerMoleWeightIdx]= this->polymer_.moleWeight[globalDofIdx];
        }

        if constexpr (enableBrine) {
            if (enableSaltPrecipitation && values.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp) {
                values[Indices::saltConcentrationIdx] = initialFluidStates_[globalDofIdx].saltSaturation();
            }
            else {
                values[Indices::saltConcentrationIdx] = initialFluidStates_[globalDofIdx].saltConcentration();
            }
        }

        if constexpr (enableBioeffects) {
            values[Indices::microbialConcentrationIdx] = this->bioeffects_.microbialConcentration[globalDofIdx];
            values[Indices::biofilmVolumeFractionIdx] = this->bioeffects_.biofilmVolumeFraction[globalDofIdx];
            if constexpr (enableMICP) {
                values[Indices::oxygenConcentrationIdx] = this->bioeffects_.oxygenConcentration[globalDofIdx];
                values[Indices::ureaConcentrationIdx] = this->bioeffects_.ureaConcentration[globalDofIdx];
                values[Indices::calciteVolumeFractionIdx] = this->bioeffects_.calciteVolumeFraction[globalDofIdx];
            }
        }

        values.checkDefined();
    }


    Scalar drsdtcon(unsigned elemIdx, int episodeIdx) const
    {
        return this->mixControls_.drsdtcon(elemIdx, episodeIdx,
                                           this->pvtRegionIndex(elemIdx));
    }

    bool drsdtconIsActive(unsigned elemIdx, int episodeIdx) const
    {
        return this->mixControls_.drsdtConvective(episodeIdx, this->pvtRegionIndex(elemIdx));
    }

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * Reservoir simulation uses no-flow conditions as default for all boundaries.
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(eclProblemBoundary, Subsystem::Assembly);
        if (!context.intersection(spaceIdx).boundary())
            return;

        if constexpr (energyModuleType != EnergyModules::FullyImplicitThermal || !enableThermalFluxBoundaries)
            values.setNoFlow();
        else {
            // in the energy case we need to specify a non-trivial boundary condition
            // because the geothermal gradient needs to be maintained. for this, we
            // simply assume the initial temperature at the boundary and specify the
            // thermal flow accordingly. in this context, "thermal flow" means energy
            // flow due to a temerature gradient while assuming no-flow for mass
            unsigned interiorDofIdx = context.interiorScvIndex(spaceIdx, timeIdx);
            unsigned globalDofIdx = context.globalSpaceIndex(interiorDofIdx, timeIdx);
            values.setThermalFlow(context, spaceIdx, timeIdx, this->initialFluidStates_[globalDofIdx] );
        }

        if (this->nonTrivialBoundaryConditions()) {
            unsigned indexInInside  = context.intersection(spaceIdx).indexInInside();
            unsigned interiorDofIdx = context.interiorScvIndex(spaceIdx, timeIdx);
            unsigned globalDofIdx = context.globalSpaceIndex(interiorDofIdx, timeIdx);
            unsigned pvtRegionIdx = pvtRegionIndex(context, spaceIdx, timeIdx);

            RateVector phaseVolRate = 0.0;
            const auto dir = FaceDir::FromIntersectionIndex(indexInInside);
            bool isMassRate = false;
            if (this->fluxBoundaryRate_(globalDofIdx, dir, phaseVolRate, isMassRate)) {
                if (isMassRate) {
                    // Already a mass rate, formed by the producing run with the
                    // density of the cell the flow comes from.
                    values.setMassRate(phaseVolRate, pvtRegionIdx);
                }
                else {
                    const auto& insideFs =
                        context.intensiveQuantities(interiorDofIdx, timeIdx).fluidState();
                    values.setMassRate(this->fluxPhaseVolToMassRate_(phaseVolRate,
                                                                     insideFs,
                                                                     globalDofIdx,
                                                                     indexInInside),
                                       pvtRegionIdx);
                }
                return;
            }

            const auto [type, massrate] = this->boundaryCondition(globalDofIdx, indexInInside);
            if (type == BCType::THERMAL)
                values.setThermalFlow(context, spaceIdx, timeIdx, this->boundaryFluidState(globalDofIdx, indexInInside));
            else if (type == BCType::FREE || type == BCType::DIRICHLET)
                values.setFreeFlow(context, spaceIdx, timeIdx, this->boundaryFluidState(globalDofIdx, indexInInside));
            else if (type == BCType::RATE)
                values.setMassRate(massrate, pvtRegionIdx);
        }
    }

private:
    //! \brief Add what crosses the sector boundary through a non-neighbour
    //!   connection as a source in the interior cell.
    //!
    //! \details The ordinary boundary faces are applied through the grid face
    //!   they sit on. An NNC has no such face here: the cell on the far side is
    //!   outside the sector and so absent from this grid, and the connection
    //!   went with it. Left alone, everything the parent run recorded across
    //!   one of these would simply vanish.
    //!
    //!   The file stores a total component mass rate, positive into the sector.
    //!   A source is expressed per unit volume and, when the model conserves
    //!   surface volume, in surface volume rather than mass, so both
    //!   conversions happen here. The reference density is the interior cell's,
    //!   which is the one the producing run used to form the mass, so the round
    //!   trip is exact.
    void addFluxNncSource_(RateVector& rate,
                           const unsigned globalDofIdx,
                           const unsigned timeIdx) const
    {
        const auto& faces = this->fluxNncFacesAt_(globalDofIdx);
        if (faces.empty()) {
            return;
        }

        const auto* fluxData = this->fluxBoundaryData_();
        const auto* fluxStep = this->fluxBoundaryReportStep_();
        if ((fluxData == nullptr) || (fluxStep == nullptr)) {
            return;
        }

        const auto fluxEnabled = (static_cast<int>(fluxData->header.mode)
                                  & static_cast<int>(EclIO::FluxFile::Mode::Flux)) != 0;
        if (!fluxEnabled) {
            this->addFluxNncPressureSource_(rate, globalDofIdx, timeIdx,
                                            faces, *fluxData, *fluxStep);
            return;
        }

        const auto faceCount = static_cast<std::size_t>(fluxData->header.numBoundaryFaces);
        const auto phaseCount = static_cast<std::size_t>(fluxData->header.numPhases);
        if ((phaseCount == 0)
            || (fluxStep->massRates.size() != faceCount * phaseCount))
        {
            return;
        }

        const auto volume = this->model().dofTotalVolume(globalDofIdx);
        if (!(volume > 0.0)) {
            return;
        }

        const auto pvtRegionIdx = this->pvtRegionIndex(globalDofIdx);
        const auto& intQuants = this->model().intensiveQuantities(globalDofIdx, timeIdx);
        const auto& fsIn = intQuants.fluidState();
        int faceIndex = -1;

        const auto add = [&](const std::size_t slot,
                             const unsigned phaseIdx,
                             const unsigned compIdx)
        {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                return;
            }

            if constexpr (energyModuleType == EnergyModules::FullyImplicitThermal) {
                // Inflow brings the exterior cell's temperature, outflow this cell's enthalpy.
                const Scalar massIn = fluxStep->massRates[slot];
                const auto f = static_cast<std::size_t>(faceIndex);
                if (massIn > 0.0 && f < fluxStep->temperature.size()) {
                    auto fs = fsIn;
                    fs.setTemperature(fluxStep->temperature[f]);
                    rate[Indices::contiEnergyEqIdx] += massIn / volume
                        * FluidSystem::enthalpy(fs, phaseIdx, pvtRegionIdx);
                }
                else {
                    rate[Indices::contiEnergyEqIdx] += massIn / volume * fsIn.enthalpy(phaseIdx);
                }
            }

            auto massRate = fluxStep->massRates[slot] / volume;
            if constexpr (getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>()) {
                massRate /= FluidSystem::referenceDensity(phaseIdx, pvtRegionIdx);
            }

            rate[FluidSystem::canonicalToActiveCompIdx(compIdx)] += massRate;
        };

        for (const auto face : faces) {
            faceIndex = face;
            if ((faceIndex < 0) || (static_cast<std::size_t>(faceIndex) >= faceCount)) {
                continue;
            }

            if constexpr (energyModuleType == EnergyModules::FullyImplicitThermal) {
                const auto f = static_cast<std::size_t>(faceIndex);
                if (this->fluxBoundary_ && f < this->fluxBoundary_->faces().size()
                    && f < fluxStep->temperature.size())
                {
                    const auto alpha = this->fluxBoundary_->faces()[f].thermalTransmissibility;
                    if (alpha > 0.0) {
                        rate[Indices::contiEnergyEqIdx] += intQuants.totalThermalConductivity() * alpha
                            * (fluxStep->temperature[f] - fsIn.temperature(0)) / volume;
                    }
                }
            }

            const auto base = static_cast<std::size_t>(faceIndex) * phaseCount;
            std::size_t phaseSlot = 0;

            if (fluxData->header.hasPhase(EclIO::FluxFile::Phase::Oil)) {
                add(base + phaseSlot, oilPhaseIdx, oilCompIdx);
                ++phaseSlot;
            }
            if (fluxData->header.hasPhase(EclIO::FluxFile::Phase::Water)) {
                add(base + phaseSlot, waterPhaseIdx, waterCompIdx);
                ++phaseSlot;
            }
            if (fluxData->header.hasPhase(EclIO::FluxFile::Phase::Gas)) {
                add(base + phaseSlot, gasPhaseIdx, gasCompIdx);
            }
        }
    }

    //! \brief Flow through the NNC boundary faces of a pressure-mode sector.
    //!
    //! \details A directional boundary face is imposed through the grid face
    //!   it sits on, and the discretisation works out the flux across it. An
    //!   NNC has no such face here, the cell on the far side having gone with
    //!   the rest of the outside, so the flux is worked out here instead and
    //!   enters as a source. It is the two-point flux the parent would have
    //!   computed across the same connection: its transmissibility, the
    //!   exterior cell's recorded state as the far end, the average density
    //!   for the gravity head between the two cell centres, the threshold
    //!   pressure if there is one, and the upstream cell's mobility and
    //!   formation volume factors. The interior end is this cell's current
    //!   state, derivatives and all, so the connection is as implicit as the
    //!   rest of the boundary.
    //!
    //!   Black-oil components, the enthalpy they carry and conduction through
    //!   the recorded FLXTHTR; water-borne gas, solvent and polymer are not
    //!   carried across.
    void addFluxNncPressureSource_(RateVector& rate,
                                   const unsigned globalDofIdx,
                                   const unsigned timeIdx,
                                   const std::vector<int>& faces,
                                   const EclIO::FluxFile::Data& fluxData,
                                   const EclIO::FluxFile::ReportStep& fluxStep) const
    {
        using Evaluation = typename FlowProblemType::Evaluation;

        const auto pressureEnabled = (static_cast<int>(fluxData.header.mode)
                                      & static_cast<int>(EclIO::FluxFile::Mode::Pressure)) != 0;
        if (!pressureEnabled || !this->fluxBoundary_) {
            return;
        }

        const auto& allFaces = this->fluxBoundary_->faces();
        const auto volume = this->model().dofTotalVolume(globalDofIdx);
        if (!(volume > 0.0)) {
            return;
        }

        const auto& intQuants = this->model().intensiveQuantities(globalDofIdx, timeIdx);
        const auto& fsIn = intQuants.fluidState();
        const auto pvtRegionIn = this->pvtRegionIndex(globalDofIdx);
        const auto zIn = this->dofCenterDepth(globalDofIdx);
        const auto g = this->gravity()[dimWorld - 1];
        const auto numPhaseSlots = static_cast<std::size_t>(fluxData.header.numPhases);

        for (const auto faceIndexInt : faces) {
            const auto faceIndex = static_cast<std::size_t>(faceIndexInt);
            if ((faceIndexInt < 0) || (faceIndex >= allFaces.size())
                || (faceIndex >= fluxStep.pressures.size()))
            {
                continue;
            }

            const auto& face = allFaces[faceIndex];
            if (!std::isfinite(face.transmissibility) || !(face.transmissibility > 0.0)) {
                continue;
            }

            const auto exFs = this->fluxExteriorFluidState_(fluxData, fluxStep, face, faceIndex,
                                                            globalDofIdx,
                                                            std::numeric_limits<Scalar>::quiet_NaN());
            const auto zEx = std::isfinite(face.exteriorDepth) ? face.exteriorDepth : zIn;

            // The exterior cell's relative permeabilities, as recorded, or this
            // cell's curves at the exterior saturations when they were not.
            std::array<Scalar, numPhases> krEx{};
            const auto krBase = faceIndex * numPhaseSlots;
            if (fluxStep.relPerm.size() >= krBase + numPhaseSlots) {
                std::size_t slot = 0;
                for (const unsigned phaseIdx : {unsigned{oilPhaseIdx},
                                                unsigned{waterPhaseIdx},
                                                unsigned{gasPhaseIdx}}) {
                    if (FluidSystem::phaseIsActive(phaseIdx)) {
                        krEx[phaseIdx] = fluxStep.relPerm[krBase + slot];
                        ++slot;
                    }
                }
            }
            else {
                MaterialLaw::relativePermeabilities(krEx, this->materialLawParams(globalDofIdx), exFs);
            }

            std::array<Evaluation, numPhases> surfaceRate{};
            Evaluation energyRate = 0.0;
            const auto upstreamRates = [&](const unsigned phaseIdx,
                                           const Evaluation& reservoirRate,
                                           const auto& fs)
            {
                const auto b = fs.invB(phaseIdx);
                surfaceRate[phaseIdx] += reservoirRate * b;

                if (!FluidSystem::phaseIsActive(oilPhaseIdx)
                    || !FluidSystem::phaseIsActive(gasPhaseIdx))
                {
                    return;
                }
                if (phaseIdx == oilPhaseIdx && FluidSystem::enableDissolvedGas()) {
                    surfaceRate[gasPhaseIdx] += reservoirRate * b * fs.Rs();
                }
                if (phaseIdx == gasPhaseIdx && FluidSystem::enableVaporizedOil()) {
                    surfaceRate[oilPhaseIdx] += reservoirRate * b * fs.Rv();
                }
            };

            for (unsigned activePhaseIdx = 0;
                 activePhaseIdx < FluidSystem::numActivePhases(); ++activePhaseIdx)
            {
                const auto phaseIdx = FluidSystem::activeToCanonicalPhaseIdx(activePhaseIdx);

                const Evaluation rhoAvg = (fsIn.density(phaseIdx) + exFs.density(phaseIdx)) / 2;
                Evaluation dp = exFs.pressure(phaseIdx) + rhoAvg * (g * (zIn - zEx))
                    - fsIn.pressure(phaseIdx);

                const bool inflow = dp > 0.0;
                const Scalar thpres =
                    this->fluxFaceThresholdPressure_(globalDofIdx, face, /*interiorToExterior=*/!inflow);
                if (thpres > 0.0) {
                    if (std::abs(getValue(dp)) > thpres) {
                        dp += inflow ? -thpres : thpres;
                    }
                    else {
                        dp = 0.0;
                    }
                }

                if (inflow) {
                    const Scalar mobility = krEx[phaseIdx] / exFs.viscosity(phaseIdx);
                    const Evaluation q = face.transmissibility * mobility * dp;
                    upstreamRates(phaseIdx, q, exFs);
                    if constexpr (energyModuleType == EnergyModules::FullyImplicitThermal) {
                        energyRate += q * exFs.density(phaseIdx)
                            * FluidSystem::enthalpy(exFs, phaseIdx, exFs.pvtRegionIndex());
                    }
                }
                else {
                    const auto trans = face.transmissibility
                        * getValue(intQuants.rockCompTransMultiplier());
                    const Evaluation q = trans * intQuants.mobility(phaseIdx) * dp;
                    upstreamRates(phaseIdx, q, fsIn);
                    if constexpr (energyModuleType == EnergyModules::FullyImplicitThermal) {
                        energyRate += q * fsIn.density(phaseIdx) * fsIn.enthalpy(phaseIdx);
                    }
                }
            }

            if constexpr (energyModuleType == EnergyModules::FullyImplicitThermal) {
                if (face.thermalTransmissibility > 0.0) {
                    energyRate += intQuants.totalThermalConductivity() * face.thermalTransmissibility
                        * (exFs.temperature(0) - fsIn.temperature(0));
                }
                rate[Indices::contiEnergyEqIdx] += energyRate / volume;
            }

            const auto add = [&](const unsigned phaseIdx, const unsigned compIdx)
            {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    return;
                }

                auto componentRate = surfaceRate[phaseIdx] / volume;
                if constexpr (!getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>()) {
                    componentRate *= FluidSystem::referenceDensity(phaseIdx, pvtRegionIn);
                }

                rate[FluidSystem::canonicalToActiveCompIdx(compIdx)] += componentRate;
            };

            add(oilPhaseIdx, oilCompIdx);
            add(waterPhaseIdx, waterCompIdx);
            add(gasPhaseIdx, gasCompIdx);
        }
    }

    template <class ScalarFluidState>
    RateVector fluxPhaseVolToMassRate_(const RateVector& phaseVolRate,
                                       const ScalarFluidState& insideFs,
                                       const unsigned globalDofIdx,
                                       const unsigned boundaryFaceIdx) const
    {
        auto checkedDensity = [&](const unsigned phaseIdx,
                                  const char* phaseName)
        {
            const Scalar rho = getValue(insideFs.density(phaseIdx));
            if (std::isfinite(rho) && rho > 0.0) {
                return rho;
            }

            std::ostringstream os;
            os << "Invalid live phase density for USEFLUX FLUX-mode conversion at dof "
               << globalDofIdx << ", boundary face " << boundaryFaceIdx
               << ", phase '" << phaseName << "': " << rho;
            throw std::runtime_error(os.str());
        };

        RateVector massRate = 0.0;
        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            const auto comp = FluidSystem::canonicalToActiveCompIdx(oilCompIdx);
            const auto volRate = phaseVolRate[comp];
            const Scalar volRateValue = getValue(volRate);
            if (!std::isfinite(volRateValue)) {
                std::ostringstream os;
                os << "Invalid USEFLUX FLUX-mode volumetric rate at dof "
                   << globalDofIdx << ", boundary face " << boundaryFaceIdx
                   << ", phase 'OIL': " << volRateValue;
                throw std::runtime_error(os.str());
            }
            massRate[comp] = volRate * checkedDensity(oilPhaseIdx, "OIL");
        }
        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            const auto comp = FluidSystem::canonicalToActiveCompIdx(waterCompIdx);
            const auto volRate = phaseVolRate[comp];
            const Scalar volRateValue = getValue(volRate);
            if (!std::isfinite(volRateValue)) {
                std::ostringstream os;
                os << "Invalid USEFLUX FLUX-mode volumetric rate at dof "
                   << globalDofIdx << ", boundary face " << boundaryFaceIdx
                   << ", phase 'WATER': " << volRateValue;
                throw std::runtime_error(os.str());
            }
            massRate[comp] = volRate * checkedDensity(waterPhaseIdx, "WATER");
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            const auto comp = FluidSystem::canonicalToActiveCompIdx(gasCompIdx);
            const auto volRate = phaseVolRate[comp];
            const Scalar volRateValue = getValue(volRate);
            if (!std::isfinite(volRateValue)) {
                std::ostringstream os;
                os << "Invalid USEFLUX FLUX-mode volumetric rate at dof "
                   << globalDofIdx << ", boundary face " << boundaryFaceIdx
                   << ", phase 'GAS': " << volRateValue;
                throw std::runtime_error(os.str());
            }
            massRate[comp] = volRate * checkedDensity(gasPhaseIdx, "GAS");
        }

        return massRate;
    }

public:

    //!\brief Read simulator solution state from the outputmodule (used with restart)
    //! \param restart_step Step to read at
    //! \param fip_init True to do limited simulator initialization
    //! \details \a fip_init is used when calculating original FIP from restart state
    void readSolutionFromOutputModule(const int restart_step, bool fip_init)
    {
        auto& simulator = this->simulator();
        const auto& eclState = simulator.vanguard().eclState();

        std::size_t numElems = this->model().numGridDof();
        this->initialFluidStates_.resize(numElems);
        if constexpr (enableSolvent) {
            this->solventSaturation_.resize(numElems, 0.0);
            this->solventRsw_.resize(numElems, 0.0);
        }

        if constexpr (enablePolymer)
            this->polymer_.concentration.resize(numElems, 0.0);

        if constexpr (enablePolymerMolarWeight) {
            const std::string msg {"Support of the RESTART for polymer molecular weight "
                                   "is not implemented yet. The polymer weight value will be "
                                   "zero when RESTART begins"};
            OpmLog::warning("NO_POLYMW_RESTART", msg);
            this->polymer_.moleWeight.resize(numElems, 0.0);
        }

        if constexpr (enableBioeffects) {
            this->bioeffects_.resize(numElems);
        }

        // Initialize mixing controls before trying to set any lastRx valuesx
        this->mixControls_.init(numElems, restart_step, eclState.runspec().tabdims().getNumPVTTables());

        if constexpr (enableBioeffects) {
            this->bioeffects_ = this->eclWriter_->outputModule().getBioeffects().getSolution();
        }

        // TEMPI is given on the input grid while the loop below indexes by leaf
        // cell; map it across once (identity without LGRs) rather than reading
        // past the end of the input-grid array for every refined cell.
        std::vector<double> tempiOnLeaf{};
        if constexpr (energyModuleType != EnergyModules::NoTemperature) {
            const auto& rspec = eclState.runspec();
            if (rspec.co2Storage() || rspec.h2Storage() || rspec.co2Sol() || rspec.h2Sol()) {
                using LeafGrid = GetPropType<TypeTag, Properties::Grid>;
                const LookUpData<LeafGrid, GridView> lookUpData(simulator.gridView());
                tempiOnLeaf = lookUpData
                    .assignFieldPropsDoubleOnLeaf(eclState.fieldProps(), "TEMPI");
            }
        }

        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = this->initialFluidStates_[elemIdx];
            elemFluidState.setPvtRegionIndex(pvtRegionIndex(elemIdx));
            this->eclWriter_->outputModule().initHysteresisParams(simulator, elemIdx);
            this->eclWriter_->outputModule().assignToFluidState(elemFluidState, elemIdx);

            // Note: Function processRestartSaturations_() mutates the
            // 'ssol' argument--the value from the restart file--if solvent
            // is enabled.  Then, store the updated solvent saturation into
            // 'solventSaturation_'.  Otherwise, just pass a dummy value to
            // the function and discard the unchanged result.  Do not index
            // into 'solventSaturation_' unless solvent is enabled.
            {
                auto ssol = enableSolvent
                            ? this->eclWriter_->outputModule().getSolventSaturation(elemIdx)
                            : Scalar(0);

                this->processRestartSaturations_(elemFluidState, ssol);

                if constexpr (enableSolvent) {
                    this->solventSaturation_[elemIdx] = ssol;
                    this->solventRsw_[elemIdx] = this->eclWriter_->outputModule().getSolventRsw(elemIdx);
                }
            }

            // For CO2STORE, H2STORE, CO2SOL and H2SOL we need to set the initial temperature
            // for isothermal simulations
            if constexpr (energyModuleType != EnergyModules::NoTemperature) {
                const auto& rspec = eclState.runspec();
                bool needTemperature = rspec.co2Storage() || rspec.h2Storage()
                    || rspec.co2Sol() || rspec.h2Sol();
                if (needTemperature) {
                    elemFluidState.setTemperature(tempiOnLeaf[elemIdx]);
                }
            }

            this->mixControls_.updateLastValues(elemIdx, elemFluidState.Rs(), elemFluidState.Rv());

            if constexpr (enablePolymer)
                this->polymer_.concentration[elemIdx] = this->eclWriter_->outputModule().getPolymerConcentration(elemIdx);
            // if we need to restart for polymer molecular weight simulation, we need to add related here
        }

        const int episodeIdx = this->episodeIndex();
        this->mixControls_.updateMaxValues(episodeIdx, simulator.timeStepSize());

        // assign the restart solution to the current solution. note that we still need
        // to compute real initial solution after this because the initial fluid states
        // need to be correct for stuff like boundary conditions.
        auto& sol = this->model().solution(/*timeIdx=*/0);
        const auto& gridView = this->gridView();
        ElementContext elemCtx(simulator);
        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            int elemIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            this->initial(sol[elemIdx], elemCtx, /*spaceIdx=*/0, /*timeIdx=*/0);
        }

        // make sure that the ghost and overlap entities exhibit the correct
        // solution. alternatively, this could be done in the loop above by also
        // considering non-interior elements. Since the initial() method might not work
        // 100% correctly for such elements, let's play safe and explicitly synchronize
        // using message passing.
        this->model().syncOverlap();

        if (fip_init) {
            this->updateReferencePorosity_();
            this->mixControls_.init(this->model().numTotalDof(),
                                    this->episodeIndex(),
                                    eclState.runspec().tabdims().getNumPVTTables());
        }
    }

    /*!
     * \copydoc FlowProblem::thresholdPressureBoundary
     */
    Scalar thresholdPressureBoundary(unsigned globalSpaceIdx,
                                     unsigned boundaryFaceIdx,
                                     bool interiorToExterior) const
    {
        const auto* face = this->fluxBoundaryFaceAtOrdinal_(globalSpaceIdx, boundaryFaceIdx);
        return (face == nullptr)
            ? Scalar{0}
            : this->fluxFaceThresholdPressure_(globalSpaceIdx, *face, interiorToExterior);
    }

    //! \brief Threshold pressure across a sector boundary face, directional or
    //!        NNC alike.
    Scalar fluxFaceThresholdPressure_(const unsigned globalSpaceIdx,
                                      const FluxBoundary::Face& face,
                                      const bool interiorToExterior) const
    {
        if (!this->thresholdPressures_.enableThresholdPressure()) {
            return 0.0;
        }

        if (face.exteriorEquilRegion < 0) {
            return 0.0;
        }

        const auto interiorRegion = this->thresholdPressures_.equilRegionIndex(globalSpaceIdx);
        const auto exteriorRegion = static_cast<unsigned>(face.exteriorEquilRegion);
        const auto numRegions = this->thresholdPressures_.numEquilRegions();

        if ((interiorRegion >= numRegions) || (exteriorRegion >= numRegions)) {
            return 0.0;
        }

        const auto& table = this->thresholdPressures_.data();
        const auto offset = interiorToExterior
            ? interiorRegion * numRegions + exteriorRegion
            : exteriorRegion * numRegions + interiorRegion;

        return (offset < table.size()) ? table[offset] : Scalar{0};
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<FlowProblemType&>(*this));
        serializer(mixControls_);
        serializer(*eclWriter_);
    }

protected:
    //! \brief Snapshot the mixing rate controls before the timestep runs.
    void captureBeginTimeStepState_() override
    {
        FlowProblemType::captureBeginTimeStepState_();
        prev_timestep_state_.mixControls = mixControls_;
    }

    //! \brief Restore mixing rate controls after a failed timestep.
    void restoreBeginTimeStepState_() override
    {
        FlowProblemType::restoreBeginTimeStepState_();
        mixControls_ = prev_timestep_state_.mixControls;
    }

    //! \brief Take the producing run's threshold pressures in place of the ones
    //!        this run worked out for itself.
    //!
    //! \details A THPRES entry given as defaulted is not a number in the deck.
    //!   It is the largest initial potential difference found anywhere along
    //!   that region boundary, so arriving at it means equilibrating the whole
    //!   field. A sector holds only part of each region boundary and so arrives
    //!   at a smaller number -- on one field, 0.01 bar where the parent had
    //!   0.81 -- which lets flow through faces the parent held shut. The
    //!   producing run recorded its table for exactly this reason, and it
    //!   covers the sector's interior region boundaries as well as the ones at
    //!   its edge.
    //!
    //!   Deliberately silent about a table of the wrong size: that means the
    //!   two runs disagree about how many equilibration regions there are, and
    //!   entries picked out of it would be meaningless rather than merely
    //!   imprecise, so the locally computed table is kept and the mismatch
    //!   reported.
    void adoptFluxThresholdPressures_()
    {
        const auto* fluxData = this->fluxBoundaryData_();
        if (fluxData == nullptr) {
            return;
        }

        const auto& table = fluxData->thresholdPressure;
        if (table.empty()) {
            return;
        }

        const auto numRegions = this->thresholdPressures_.numEquilRegions();
        const auto expected = static_cast<std::size_t>(numRegions) * numRegions;

        if (!this->thresholdPressures_.enableThresholdPressure()) {
            OpmLog::warning("USEFLUX: the producing run had threshold pressures but this "
                            "one does not; add EQLOPTS/THPRES to honour them");
            return;
        }

        if (table.size() != expected) {
            OpmLog::warning(fmt::format("USEFLUX: ignoring the producing run's threshold "
                                        "pressures: {} entries for {} equilibration "
                                        "regions, expected {}",
                                        table.size(), numRegions, expected));
            return;
        }

        this->thresholdPressures_.setFromRestart(table);

        OpmLog::info(fmt::format("USEFLUX: adopted the producing run's threshold pressures "
                                 "over {} equilibration region(s)", numRegions));

        // The table printed during finishInit() was the one this run worked out
        // for itself, which is no longer the one it will use.
        this->thresholdPressures_.logPressures();
    }

    void updateExplicitQuantities_(int episodeIdx, int timeStepSize, const bool first_step_after_restart) override
    {
        this->updateExplicitQuantities_(first_step_after_restart);

        if constexpr (getPropValue<TypeTag, Properties::EnablePolymer>())
            updateMaxPolymerAdsorption_();

        mixControls_.updateExplicitQuantities(episodeIdx, timeStepSize);
    }

    void updateMaxPolymerAdsorption_()
    {
        // we need to update the max polymer adsoption data for all elements
        this->updateProperty_("FlowProblemBlackoil::updateMaxPolymerAdsorption_() failed:",
                              [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                              {
                                  this->updateMaxPolymerAdsorption_(compressedDofIdx,iq);
                              });
    }

    bool updateMaxPolymerAdsorption_(unsigned compressedDofIdx, const IntensiveQuantities& iq)
    {
        const Scalar pa = scalarValue(iq.polymerAdsorption());
        auto& mpa = this->polymer_.maxAdsorption;
        if (mpa[compressedDofIdx] < pa) {
            mpa[compressedDofIdx] = pa;
            return true;
        } else {
            return false;
        }
    }

    void computeAndSetEqWeights_()
    {
        std::vector<Scalar> sumInvB(numPhases, 0.0);
        const auto& gridView = this->gridView();
        ElementContext elemCtx(this->simulator());
        for(const auto& elem: elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            int elemIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& dofFluidState = this->initialFluidStates_[elemIdx];
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                sumInvB[phaseIdx] += dofFluidState.invB(phaseIdx);
            }
        }

        std::size_t numDof = this->model().numGridDof();
        const auto& comm = this->simulator().vanguard().grid().comm();
        comm.sum(sumInvB.data(),sumInvB.size());
        Scalar numTotalDof = comm.sum(numDof);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            Scalar avgB = numTotalDof / sumInvB[phaseIdx];
            const unsigned solventCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
            const unsigned activeSolventCompIdx = FluidSystem::canonicalToActiveCompIdx(solventCompIdx);
            this->model().setEqWeight(activeSolventCompIdx, avgB);
        }
    }

    // update the parameters needed for DRSDT and DRVDT
    bool updateCompositionChangeLimits_()
    {
        OPM_TIMEBLOCK(updateCompositionChangeLimits);
        // update the "last Rs" values for all elements, including the ones in the ghost
        // and overlap regions
        int episodeIdx = this->episodeIndex();
        std::array<bool,3> active{this->mixControls_.drsdtConvective(episodeIdx),
                                  this->mixControls_.drsdtActive(episodeIdx),
                                  this->mixControls_.drvdtActive(episodeIdx)};
        if (!active[0] && !active[1] && !active[2]) {
            return false;
        }

        this->updateProperty_("FlowProblemBlackoil::updateCompositionChangeLimits_()) failed:",
                              [this,episodeIdx,active](unsigned compressedDofIdx,
                                                       const IntensiveQuantities& iq)
                              {
                                  const DimMatrix& perm = this->intrinsicPermeability(compressedDofIdx);
                                  const Scalar distZ = active[0] ? this->simulator().vanguard().cellThickness(compressedDofIdx) : 0.0;
                                  const int pvtRegionIdx = this->pvtRegionIndex(compressedDofIdx);
                                  this->mixControls_.update(compressedDofIdx,
                                                            iq,
                                                            episodeIdx,
                                                            this->gravity_[dim - 1],
                                                            perm[dim - 1][dim - 1],
                                                            distZ,
                                                            pvtRegionIdx);
                              }
        );

        return true;
    }

    void readEclRestartSolution_()
    {
        // Restarting a refined run reads the solution section each level was
        // written to and puts the leaf back together. In parallel the reference
        // grid holding that leaf ordering exists only on the I/O rank, and
        // handing the assembled solution to the others is not solved yet, so
        // say so rather than fail somewhere further in.
        if ((this->simulator().vanguard().grid().maxLevel() > 0) &&
            (this->simulator().vanguard().grid().comm().size() > 1))
        {
            throw std::invalid_argument {
                "Restarting a run with LGRs is supported on a single MPI process "
                "only. Run the restart in serial, or drop the refinement."
            };
        }

        // Set the start time of the simulation
        auto& simulator = this->simulator();
        const auto& schedule = simulator.vanguard().schedule();
        const auto& eclState = simulator.vanguard().eclState();
        const auto& initconfig = eclState.getInitConfig();
        const int restart_step = initconfig.getRestartStep();
        {
            simulator.setTime(schedule.seconds(restart_step));

            simulator.startNextEpisode(simulator.startTime() + simulator.time(),
                                       schedule.stepLength(restart_step));
            simulator.setEpisodeIndex(restart_step);
        }
        this->eclWriter_->beginRestart();

        Scalar dt = std::min(this->eclWriter_->restartTimeStepSize(), simulator.episodeLength());
        simulator.setTimeStepSize(dt);

        this->readSolutionFromOutputModule(restart_step, false);

        this->eclWriter_->endRestart();
    }

    void readEquilInitialCondition_() override
    {
        const auto& simulator = this->simulator();

        // initial condition corresponds to hydrostatic conditions.
        EquilInitializer<TypeTag> equilInitializer(simulator, *(this->materialLawManager_));

        std::size_t numElems = this->model().numGridDof();
        this->initialFluidStates_.resize(numElems);
        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = this->initialFluidStates_[elemIdx];
            elemFluidState.assign(equilInitializer.initialFluidState(elemIdx));
        }
    }

    void readExplicitInitialCondition_() override
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& fp = eclState.fieldProps();
        bool has_swat     = fp.has_double("SWAT");
        bool has_sgas     = fp.has_double("SGAS");
        bool has_rs       = fp.has_double("RS");
        bool has_rsw      = fp.has_double("RSW");
        bool has_rv       = fp.has_double("RV");
        bool has_rvw      = fp.has_double("RVW");
        bool has_pressure = fp.has_double("PRESSURE");
        bool has_salt     = fp.has_double("SALT");
        bool has_saltp    = fp.has_double("SALTP");

        // make sure all required quantities are enables
        if (Indices::numPhases > 1) {
            if (FluidSystem::phaseIsActive(waterPhaseIdx) && !has_swat)
                throw std::runtime_error("The ECL input file requires the presence of the SWAT keyword if "
                                         "the water phase is active");
            if (FluidSystem::phaseIsActive(gasPhaseIdx) && !has_sgas && FluidSystem::phaseIsActive(oilPhaseIdx))
                throw std::runtime_error("The ECL input file requires the presence of the SGAS keyword if "
                                         "the gas phase is active");
        }
        if (!has_pressure)
            throw std::runtime_error("The ECL input file requires the presence of the PRESSURE "
                                     "keyword if the model is initialized explicitly");
        if (FluidSystem::enableDissolvedGas() && !has_rs)
            throw std::runtime_error("The ECL input file requires the RS keyword to be present if"
                                     " dissolved gas is enabled and the model is initialized explicitly");
        if (FluidSystem::enableDissolvedGasInWater() && !has_rsw)
            OpmLog::warning("The model is initialized explicitly and the RSW keyword is not present in the"
                            " ECL input file. The RSW values are set equal to 0");
        if (FluidSystem::enableVaporizedOil() && !has_rv)
            throw std::runtime_error("The ECL input file requires the RV keyword to be present if"
                                     " vaporized oil is enabled and the model is initialized explicitly");
        if (FluidSystem::enableVaporizedWater() && !has_rvw)
            throw std::runtime_error("The ECL input file requires the RVW keyword to be present if"
                                     " vaporized water is enabled and the model is initialized explicitly");
        if (enableBrine && !has_salt)
            throw std::runtime_error("The ECL input file requires the SALT keyword to be present if"
                                     " brine is enabled and the model is initialized explicitly");
        if (enableSaltPrecipitation && !has_saltp)
            throw std::runtime_error("The ECL input file requires the SALTP keyword to be present if"
                                     " salt precipitation is enabled and the model is initialized explicitly");

        std::size_t numDof = this->model().numGridDof();

        initialFluidStates_.resize(numDof);

        // The arrays below are given on the (unrefined) input grid but are read
        // per leaf cell.  Refinement makes the leaf longer, so map them across --
        // a refined cell inherits its parent cell's value.  Reading the input-grid
        // array by leaf index instead runs off its end, and the refined cells come
        // up initialised from whatever follows in memory; the run then fails to
        // converge on step one, which reads as a solver problem rather than an
        // initialisation one.  LookUpData is the identity without LGRs.
        using LeafGrid = GetPropType<TypeTag, Properties::Grid>;
        const LookUpData<LeafGrid, GridView> lookUpData(this->simulator().gridView());
        const auto onLeaf = [&lookUpData, &fp](const std::string& kw)
        { return lookUpData.assignFieldPropsDoubleOnLeaf(fp, kw); };

        std::vector<double> waterSaturationData;
        std::vector<double> gasSaturationData;
        std::vector<double> pressureData;
        std::vector<double> rsData;
        std::vector<double> rswData;
        std::vector<double> rvData;
        std::vector<double> rvwData;
        std::vector<double> tempiData;
        std::vector<double> saltData;
        std::vector<double> saltpData;

        if (FluidSystem::phaseIsActive(waterPhaseIdx) && Indices::numPhases > 1)
            waterSaturationData = onLeaf("SWAT");
        else
            waterSaturationData.resize(numDof);

        if (FluidSystem::phaseIsActive(gasPhaseIdx) && FluidSystem::phaseIsActive(oilPhaseIdx))
            gasSaturationData = onLeaf("SGAS");
        else
            gasSaturationData.resize(numDof);

        pressureData = onLeaf("PRESSURE");
        if (FluidSystem::enableDissolvedGas())
            rsData = onLeaf("RS");

        if (FluidSystem::enableDissolvedGasInWater() && has_rsw)
            rswData = onLeaf("RSW");

        if (FluidSystem::enableVaporizedOil())
            rvData = onLeaf("RV");

        if (FluidSystem::enableVaporizedWater())
            rvwData = onLeaf("RVW");

        // initial reservoir temperature
        tempiData = onLeaf("TEMPI");

        // initial salt concentration data
        if constexpr (enableBrine)
            saltData = onLeaf("SALT");

        // initial precipitated salt saturation data
        if constexpr (enableSaltPrecipitation)
            saltpData = onLeaf("SALTP");

        // calculate the initial fluid states
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofFluidState = initialFluidStates_[dofIdx];

            dofFluidState.setPvtRegionIndex(pvtRegionIndex(dofIdx));

            //////
            // set temperature
            //////
            if constexpr (energyModuleType != EnergyModules::NoTemperature) {
                Scalar temperatureLoc = tempiData[dofIdx];
                if (!std::isfinite(temperatureLoc) || temperatureLoc <= 0)
                    temperatureLoc = FluidSystem::surfaceTemperature;
                dofFluidState.setTemperature(temperatureLoc);
            }

            //////
            // set salt concentration
            //////
            if constexpr (enableBrine)
                dofFluidState.setSaltConcentration(saltData[dofIdx]);

            //////
            // set precipitated salt saturation
            //////
            if constexpr (enableSaltPrecipitation)
                dofFluidState.setSaltSaturation(saltpData[dofIdx]);

            //////
            // set saturations
            //////
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))
                dofFluidState.setSaturation(FluidSystem::waterPhaseIdx,
                                            waterSaturationData[dofIdx]);

            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)){
                if (!FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)){
                    dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                                1.0
                                                - waterSaturationData[dofIdx]);
                }
                else
                    dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                                gasSaturationData[dofIdx]);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const Scalar soil = 1.0 - waterSaturationData[dofIdx] - gasSaturationData[dofIdx];
                if (soil < smallSaturationTolerance_) {
                    dofFluidState.setSaturation(FluidSystem::oilPhaseIdx, 0.0);
                }
                else {
                    dofFluidState.setSaturation(FluidSystem::oilPhaseIdx, soil);
                }
            }

            //////
            // set phase pressures
            //////
            Scalar pressure = pressureData[dofIdx]; // oil pressure (or gas pressure for water-gas system or water pressure for single phase)

            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            std::array<Scalar, numPhases> pc = {0};
            const auto& matParams = this->materialLawParams(dofIdx);
            MaterialLaw::capillaryPressures(pc, matParams, dofFluidState);
            Valgrind::CheckDefined(pressure);
            Valgrind::CheckDefined(pc);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                if (Indices::oilEnabled)
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
                else if (Indices::gasEnabled)
                    dofFluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[gasPhaseIdx]));
                else if (Indices::waterEnabled)
                    //single (water) phase
                    dofFluidState.setPressure(phaseIdx, pressure);
            }

            if constexpr (enableDissolvedGas) {
                if (FluidSystem::enableDissolvedGas())
                    dofFluidState.setRs(rsData[dofIdx]);
                else if (Indices::gasEnabled && Indices::oilEnabled)
                    dofFluidState.setRs(0.0);
                if (FluidSystem::enableVaporizedOil())
                    dofFluidState.setRv(rvData[dofIdx]);
                else if (Indices::gasEnabled && Indices::oilEnabled)
                    dofFluidState.setRv(0.0);
            }

            if constexpr (enableDisgasInWater) {
                if (FluidSystem::enableDissolvedGasInWater() && has_rsw)
                    dofFluidState.setRsw(rswData[dofIdx]);
            }

            if constexpr (enableVapwat) {
                if (FluidSystem::enableVaporizedWater())
                    dofFluidState.setRvw(rvwData[dofIdx]);
            }

            //////
            // set invB_
            //////
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                const auto& b = FluidSystem::inverseFormationVolumeFactor(dofFluidState, phaseIdx, pvtRegionIndex(dofIdx));
                dofFluidState.setInvB(phaseIdx, b);

                const auto& rho = FluidSystem::density(dofFluidState, phaseIdx, pvtRegionIndex(dofIdx));
                dofFluidState.setDensity(phaseIdx, rho);

            }
        }
    }


    void processRestartSaturations_(InitialFluidState& elemFluidState, Scalar& solventSaturation)
    {
        // each phase needs to be above certain value to be claimed to be existing
        // this is used to recover some RESTART running with the defaulted single-precision format
        Scalar sumSaturation = 0.0;
        for (std::size_t phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                if (elemFluidState.saturation(phaseIdx) < smallSaturationTolerance_)
                    elemFluidState.setSaturation(phaseIdx, 0.0);

                sumSaturation += elemFluidState.saturation(phaseIdx);
            }

        }
        if constexpr (enableSolvent) {
            if (solventSaturation < smallSaturationTolerance_)
                solventSaturation = 0.0;

            sumSaturation += solventSaturation;
        }

        assert(sumSaturation > 0.0);

        for (std::size_t phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                const Scalar saturation = elemFluidState.saturation(phaseIdx) / sumSaturation;
                elemFluidState.setSaturation(phaseIdx, saturation);
            }
        }
        if constexpr (enableSolvent) {
            solventSaturation = solventSaturation / sumSaturation;
        }
    }

    void readInitialCondition_() override
    {
        FlowProblemType::readInitialCondition_();

        if constexpr (enableSolvent || enablePolymer || enablePolymerMolarWeight || enableBioeffects)
            this->readBlackoilExtentionsInitialConditions_(this->model().numGridDof(),
                                                           enableSolvent,
                                                           enablePolymer,
                                                           enablePolymerMolarWeight,
                                                           enableBioeffects,
                                                           enableMICP);

    }

    void handleSolventBC(const BCState::BCFace& bc, RateVector& rate) const override
    {
        if constexpr (!enableSolvent)
            throw std::logic_error("solvent is disabled and you're trying to add solvent to BC");

        rate[Indices::solventSaturationIdx] = bc.rate;
    }

    void handlePolymerBC(const BCState::BCFace& bc, RateVector& rate) const override
    {
        if constexpr (!enablePolymer)
            throw std::logic_error("polymer is disabled and you're trying to add polymer to BC");

        rate[Indices::polymerConcentrationIdx] = bc.rate;
    }

    void handleMicrBC(const BCState::BCFace& bc, RateVector& rate) const override
    {
        if constexpr (!enableMICP)
            throw std::logic_error("MICP is disabled and you're trying to add microbes to BC");

        rate[Indices::microbialConcentrationIdx] = bc.rate;
    }

    void handleOxygBC(const BCState::BCFace& bc, RateVector& rate) const override
    {
        if constexpr (!enableMICP)
            throw std::logic_error("MICP is disabled and you're trying to add oxygen to BC");

        rate[Indices::oxygenConcentrationIdx] = bc.rate;
    }

    void handleUreaBC(const BCState::BCFace& bc, RateVector& rate) const override
    {
        if constexpr (!enableMICP)
            throw std::logic_error("MICP is disabled and you're trying to add urea to BC");

        rate[Indices::ureaConcentrationIdx] = bc.rate;
        // since the urea concentration can be much larger than 1, then we apply a scaling factor
        rate[Indices::ureaConcentrationIdx] *= getPropValue<TypeTag, Properties::BlackOilUreaScalingFactor>();
    }

    void updateExplicitQuantities_(const bool first_step_after_restart)
    {
        OPM_TIMEBLOCK(updateExplicitQuantities);
        const bool invalidateFromMaxWaterSat = this->updateMaxWaterSaturation_();
        const bool invalidateFromMinPressure = this->updateMinPressure_();

        // update hysteresis and max oil saturation used in vappars
        const bool invalidateFromHyst = this->updateHysteresis_();
        const bool invalidateFromMaxOilSat = this->updateMaxOilSaturation_();

        // deal with DRSDT and DRVDT
        const bool invalidateDRDT = !first_step_after_restart && this->updateCompositionChangeLimits_();

        // the derivatives may have changed
        const bool invalidateIntensiveQuantities
            = invalidateFromMaxWaterSat || invalidateFromMinPressure || invalidateFromHyst || invalidateFromMaxOilSat || invalidateDRDT;
        if (invalidateIntensiveQuantities) {
            OPM_TIMEBLOCK(beginTimeStepInvalidateIntensiveQuantities);
            this->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        }

        this->updateRockCompTransMultVal_();
    }

    bool satfuncConsistencyRequirementsMet() const
    {
        if (const auto nph = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
            + FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
            + FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
            nph < 2)
        {
            // Single phase runs don't need saturation functions and there's
            // nothing to do here.  Return 'true' to tell caller that the
            // consistency requirements are Met.
            return true;
        }

        const auto numSamplePoints = static_cast<std::size_t>
            (Parameters::Get<Parameters::NumSatfuncConsistencySamplePoints>());

        auto sfuncConsistencyChecks =
            Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar> {
            numSamplePoints, this->simulator().vanguard().eclState(),
            [&cmap = this->simulator().vanguard().cartesianIndexMapper()](const int elemIdx)
            { return cmap.cartesianIndex(elemIdx); }
        };

        const auto ioRank = 0;
        const auto isIoRank = this->simulator().vanguard()
            .grid().comm().rank() == ioRank;

        // Note: Run saturation function consistency checks on main grid
        // only (i.e., levelGridView(0)).  These checks are not supported
        // for LGRs at this time.
        sfuncConsistencyChecks.collectFailuresTo(ioRank)
            .run(this->simulator().vanguard().grid().levelGridView(0),
                 [&vg   = this->simulator().vanguard(),
                  &emap = this->simulator().model().elementMapper()]
                 (const auto& elem)
                 { return vg.gridIdxToEquilGridIdx(emap.index(elem)); });

        using ViolationLevel = typename Satfunc::PhaseChecks::
            SatfuncConsistencyCheckManager<Scalar>::ViolationLevel;

        auto reportFailures = [&sfuncConsistencyChecks]
            (const ViolationLevel level)
        {
            sfuncConsistencyChecks.reportFailures
                (level, [](std::string_view record)
                { OpmLog::info(std::string { record }); });
        };

        if (sfuncConsistencyChecks.anyFailedStandardChecks()) {
            if (isIoRank) {
                OpmLog::warning("Saturation Function "
                                "End-point Consistency Problems");

                reportFailures(ViolationLevel::Standard);
            }
        }

        if (sfuncConsistencyChecks.anyFailedCriticalChecks()) {
            if (isIoRank) {
                OpmLog::error("Saturation Function "
                              "End-point Consistency Failures");

                reportFailures(ViolationLevel::Critical);
            }

            // There are "critical" check failures.  Report that consistency
            // requirements are not Met.
            return false;
        }

        // If we get here then there are no critical failures.  Report
        // Met = true, i.e., that the consistency requirements ARE met.
        return true;
    }

    std::vector<InitialFluidState> initialFluidStates_;

    bool enableEclOutput_;
    std::unique_ptr<EclWriterType> eclWriter_;

    const Scalar smallSaturationTolerance_ = 1.e-6;
#if HAVE_DAMARIS
    bool enableDamarisOutput_ = false ;
    std::unique_ptr<DamarisWriterType> damarisWriter_;
#endif
    MixingRateControls<FluidSystem> mixControls_;

    //! \brief Blackoil part of the begin-of-timestep snapshot; see
    //!        FlowProblem::PrevTimestepState.
    struct PrevTimestepState
    {
        explicit PrevTimestepState(const Schedule& schedule)
            : mixControls(schedule)
        {}

        MixingRateControls<FluidSystem> mixControls; //!< DRSDT / DRVDT
    };

    PrevTimestepState prev_timestep_state_;

    ActionHandler<Scalar, IndexTraits> actionHandler_;

    HybridNewton hybridNewton_;

private:
    /// Whether ECL summary/FIP evaluation is wired up for this grid
    /// configuration.
    ///
    /// The only unsupported case is a *distributed* CpGrid with LGRs:
    /// CollectDataOnIORank does not build its index maps for a distributed
    /// refined grid yet.  In serial the index maps are the identity and the
    /// leaf-grid well/cell data is written directly, so refined serial runs
    /// are fine.  Kept in one place so no call site can miss a condition.
    bool eclOutputEvalSupported_() const
    {
        const auto& grid = this->simulator().vanguard().gridView().grid();

        using GridType = std::remove_cv_t<std::remove_reference_t<decltype(grid)>>;
        constexpr bool isCpGrid = std::is_same_v<GridType, Dune::CpGrid>;

        return !isCpGrid
            || (grid.maxLevel() == 0)
            || (grid.comm().size() == 1);
    }

    /// Whether or not the current epsiode will end at the end of the
    /// current time step.
    ///
    /// Custom implementation for black-oil cases.
    ///
    /// [2025-08-14] This is arguably a hack.  We intentionally do not use
    /// the 'simulator()'s episodeWillBeOver() function, as that has been
    /// shown to have non-trivial false negatives in the context of certain
    /// regression tests.  It is likely that we will need to revisit this
    /// predicate in the future.
    bool episodeWillBeOver() const override
    {
        const auto currTime = this->simulator().time()
            + this->simulator().timeStepSize();

        const auto nextReportStep =
            this->simulator().vanguard().schedule()
            .seconds(this->simulator().episodeIndex() + 1);

        const auto isSubStep = (nextReportStep - currTime)
            > (2 * std::numeric_limits<float>::epsilon()) * nextReportStep;

        return !isSubStep;
    }
};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_BLACK_HPP
