// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
/**
 * \file
 *
 * \copydoc Opm::EclWellManager
 */
#ifndef EWOMS_ECL_WELL_MANAGER_HH
#define EWOMS_ECL_WELL_MANAGER_HH

#include "eclpeacemanwell.hh"

#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Events.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/WGState.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellConnections.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/output/data/Wells.hpp>
#include <opm/output/data/Groups.hpp>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/parallel/threadedentityiterator.hh>

#include <dune/grid/common/gridenums.hh>

#include <map>
#include <unordered_set>
#include <string>
#include <vector>

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief A class which handles well controls as specified by an
 *        Eclipse deck
 */
template <class TypeTag>
class EclWellManager
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = FluidSystem::numPhases };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    typedef typename GridView::template Codim<0>::Entity Element;

    using Well = EclPeacemanWell<TypeTag>;

    typedef std::map<int, std::pair<const Connection*, std::shared_ptr<Well> > > WellConnectionsMap;

    typedef Dune::FieldVector<Evaluation, numEq> EvalEqVector;

public:
    EclWellManager(Simulator& simulator)
        : simulator_(simulator)
    { }

    /*!
     * \brief This sets up the basic properties of all wells.
     *
     * I.e., well positions, names etc...
     */
    void init()
    {
        const Schedule& deckSchedule = simulator_.vanguard().schedule();
        const auto& summaryState = simulator_.vanguard().summaryState();
        // create the wells which intersect with the current process' grid
        const auto wellsatEnd = deckSchedule.getWellsatEnd();
        for (size_t deckWellIdx = 0; deckWellIdx < deckSchedule.numWells(); ++deckWellIdx)
        {
            const auto& deckWell = wellsatEnd[deckWellIdx];
            const std::string& wellName = deckWell.name();
            Scalar wellTemperature = 273.15 + 15.56; // [K]
            if (deckWell.isInjector())
                wellTemperature = deckWell.injectionControls(summaryState).temperature;

            // set the name of the well but not much else. (i.e., if it is not completed,
            // the well primarily serves as a placeholder.) The big rest of the well is
            // specified by the updateWellCompletions_() method
            auto well = std::make_shared<Well>(simulator_);
            well->setName(wellName);
            well->setWellStatus(Well::Shut);
            well->setTemperature(wellTemperature);

            wells_.push_back(well);
            wellNameToIndex_[well->name()] = wells_.size() - 1;
        }
    }

    /*!
     * \brief This should be called the problem before each simulation
     *        episode to adapt the well controls.
     */
    void beginEpisode(bool wasRestarted=false)
    {
        const EclipseState& eclState = simulator_.vanguard().eclState();
        const Schedule& deckSchedule = simulator_.vanguard().schedule();
        const auto& summaryState = simulator_.vanguard().summaryState();
        unsigned episodeIdx = simulator_.episodeIndex();
        WellConnectionsMap wellCompMap;
        computeWellConnectionsMap_(episodeIdx, wellCompMap);

        if (wasRestarted || wellTopologyChanged_(eclState, deckSchedule, episodeIdx))
            updateWellTopology_(episodeIdx, wellCompMap, gridDofIsPenetrated_);

        // set those parameters of the wells which do not change the topology of the
        // linearized system of equations
        updateWellParameters_(episodeIdx, wellCompMap);

        const auto& deckWells = deckSchedule.getWells(episodeIdx);
        // set the injection data for the respective wells.
        for (const auto& deckWell : deckWells) {
            if (!hasWell(deckWell.name()))
                continue;

            auto well = this->well(deckWell.name());

            if (deckWell.isInjector( ))
                well->setTemperature(deckWell.injectionControls(summaryState).temperature);

            auto deckWellStatus = deckWell.getStatus( );
            switch (deckWellStatus) {
            case ::Opm::Well::Status::AUTO:
                // TODO: for now, auto means open...
            case ::Opm::Well::Status::OPEN:
                well->setWellStatus(Well::Open);
                break;
            case ::Opm::Well::Status::STOP:
                well->setWellStatus(Well::Closed);
                break;
            case ::Opm::Well::Status::SHUT:
                well->setWellStatus(Well::Shut);
                break;
            }

            // make sure that the well is either an injector or a
            // producer for the current episode. (it is not allowed to
            // be neither or to be both...)
            assert((deckWell.isInjector( )?1:0) +
                   (deckWell.isProducer( )?1:0) == 1);

            if (deckWell.isInjector( )) {
                well->setWellType(Well::Injector);
                const auto controls = deckWell.injectionControls(summaryState);
                switch (controls.injector_type) {
                case InjectorType::WATER:
                    well->setInjectedPhaseIndex(FluidSystem::waterPhaseIdx);
                    break;
                case InjectorType::GAS:
                    well->setInjectedPhaseIndex(FluidSystem::gasPhaseIdx);
                    break;
                case InjectorType::OIL:
                    well->setInjectedPhaseIndex(FluidSystem::oilPhaseIdx);
                    break;
                case InjectorType::MULTI:
                    throw std::runtime_error("Not implemented: Multi-phase injector wells");
                }

                using InjectorCMode = ::Opm::Well::InjectorCMode;
                switch (controls.cmode) {
                case InjectorCMode::RATE:
                    well->setControlMode(Well::ControlMode::VolumetricSurfaceRate);
                    break;

                case InjectorCMode::RESV:
                    well->setControlMode(Well::ControlMode::VolumetricReservoirRate);
                    break;

                case InjectorCMode::BHP:
                    well->setControlMode(Well::ControlMode::BottomHolePressure);
                    break;

                case InjectorCMode::THP:
                    well->setControlMode(Well::ControlMode::TubingHeadPressure);
                    break;

                case InjectorCMode::GRUP:
                    throw std::runtime_error("Not implemented: Well groups");

                case InjectorCMode::CMODE_UNDEFINED:
                    std::cout << "Warning: Control mode of injection well " << well->name()
                              << " is undefined. Assuming well to be shut.\n";
                    well->setWellStatus(Well::WellStatus::Shut);
                    continue;
                }

                switch (controls.injector_type) {
                case InjectorType::WATER:
                    well->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/0.0, /*water=*/1.0);
                    break;

                case InjectorType::OIL:
                    well->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/0.0, /*water=*/0.0);
                    break;

                case InjectorType::GAS:
                    well->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/1.0, /*water=*/0.0);
                    break;

                case InjectorType::MULTI:
                    throw std::runtime_error("Not implemented: Multi-phase injection wells");
                }

                well->setMaximumSurfaceRate(controls.surface_rate);
                well->setMaximumReservoirRate(controls.reservoir_rate);
                well->setTargetBottomHolePressure(controls.bhp_limit);

                // TODO
                well->setTargetTubingHeadPressure(1e30);
                //well->setTargetTubingHeadPressure(controls.thp_limit);
            }

            if (deckWell.isProducer( )) {
                well->setWellType(Well::Producer);
                const auto controls = deckWell.productionControls(summaryState);

                using ProducerCMode = ::Opm::Well::ProducerCMode;
                switch (controls.cmode) {
                case ProducerCMode::ORAT:
                    well->setControlMode(Well::ControlMode::VolumetricSurfaceRate);
                    well->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/0.0, /*water=*/0.0);
                    well->setMaximumSurfaceRate(controls.oil_rate);
                    break;

                case ProducerCMode::GRAT:
                    well->setControlMode(Well::ControlMode::VolumetricSurfaceRate);
                    well->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/1.0, /*water=*/0.0);
                    well->setMaximumSurfaceRate(controls.gas_rate);
                    break;

                case ProducerCMode::WRAT:
                    well->setControlMode(Well::ControlMode::VolumetricSurfaceRate);
                    well->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/0.0, /*water=*/1.0);
                    well->setMaximumSurfaceRate(controls.water_rate);
                    break;

                case ProducerCMode::LRAT:
                    well->setControlMode(Well::ControlMode::VolumetricSurfaceRate);
                    well->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/0.0, /*water=*/1.0);
                    well->setMaximumSurfaceRate(controls.liquid_rate);
                    break;

                case ProducerCMode::CRAT:
                    throw std::runtime_error("Not implemented: Linearly combined rates");

                case ProducerCMode::RESV:
                    well->setControlMode(Well::ControlMode::VolumetricReservoirRate);
                    well->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/1.0, /*water=*/1.0);
                    well->setMaximumSurfaceRate(controls.resv_rate);
                    break;

                case ProducerCMode::BHP:
                    well->setControlMode(Well::ControlMode::BottomHolePressure);
                    break;

                case ProducerCMode::THP:
                    well->setControlMode(Well::ControlMode::TubingHeadPressure);
                    break;

                case ProducerCMode::GRUP:
                    throw std::runtime_error("Not implemented: Well groups");

                case ProducerCMode::NONE:
                    // fall-through
                case ProducerCMode::CMODE_UNDEFINED:
                    std::cout << "Warning: Control mode of production well " << well->name()
                              << " is undefined. Assuming well to be shut.";
                    well->setWellStatus(Well::WellStatus::Shut);
                    continue;
                }

                well->setTargetBottomHolePressure(controls.bhp_limit);

                // TODO
                well->setTargetTubingHeadPressure(-1e30);
                //well->setTargetTubingHeadPressure(controls.thp_limit);
            }
        }
    }

    /*!
     * \brief Return the number of wells considered by the EclWellManager.
     */
    unsigned numWells() const
    { return wells_.size(); }

    /*!
     * \brief Return if a given well name is known to the wells manager
     */
    bool hasWell(const std::string& wellName) const
    {
        return wellNameToIndex_.find( wellName ) != wellNameToIndex_.end();
    }

    /*!
     * \brief Returns true iff a given degree of freedom is currently penetrated by any well.
     */
    bool gridDofIsPenetrated(unsigned globalDofIdx) const
    { return gridDofIsPenetrated_[globalDofIdx]; }

    /*!
     * \brief Given a well name, return the corresponding index.
     *
     * A std::runtime_error will be thrown if the well name is unknown.
     */
    unsigned wellIndex(const std::string& wellName) const
    {
        assert( hasWell( wellName ) );
        const auto& it = wellNameToIndex_.find(wellName);
        if (it == wellNameToIndex_.end())
            throw std::runtime_error("No well called '"+wellName+"'found");
        return it->second;
    }

    /*!
     * \brief Given a well name, return the corresponding well.
     *
     * A std::runtime_error will be thrown if the well name is unknown.
     */
    std::shared_ptr<const Well> well(const std::string& wellName) const
    { return wells_[wellIndex(wellName)]; }

    /*!
     * \brief Given a well name, return the corresponding well.
     *
     * A std::runtime_error will be thrown if the well name is unknown.
     */
    std::shared_ptr<Well> well(const std::string& wellName)
    { return wells_[wellIndex(wellName)]; }

    /*!
     * \brief Given a well index, return the corresponding well.
     */
    std::shared_ptr<const Well> well(size_t wellIdx) const
    { return wells_[wellIdx]; }

    /*!
     * \brief Given a well index, return the corresponding well.
     */
    std::shared_ptr<Well> well(size_t wellIdx)
    { return wells_[wellIdx]; }

    /*!
     * \brief Informs the well manager that a time step has just begun.
     */
    void beginTimeStep()
    {
        // iterate over all wells and notify them individually
        for (size_t wellIdx = 0; wellIdx < wells_.size(); ++wellIdx)
            wells_[wellIdx]->beginTimeStep();
    }

    /*!
     * \brief Informs the well that an iteration has just begun.
     *
     * In this method, the well calculates the bottom hole and tubing head pressures, the
     * actual unconstraint production and injection rates, etc.
     */
    void beginIteration()
    {
        // call the preprocessing routines
        const size_t wellSize = wells_.size();
        for (size_t wellIdx = 0; wellIdx < wellSize; ++wellIdx)
            wells_[wellIdx]->beginIterationPreProcess();

        // call the accumulation routines
        ElementContext elemCtx(simulator_);
        const auto gridView = simulator_.vanguard().gridView();
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            for (size_t wellIdx = 0; wellIdx < wellSize; ++wellIdx)
                wells_[wellIdx]->beginIterationAccumulate(elemCtx, /*timeIdx=*/0);
        }
        OPM_END_PARALLEL_TRY_CATCH("EclWellManager::beginIteration() failed: ");

        // call the postprocessing routines
        for (size_t wellIdx = 0; wellIdx < wellSize; ++wellIdx)
            wells_[wellIdx]->beginIterationPostProcess();
    }

    /*!
     * \brief Informs the well manager that an iteration has just been finished.
     */
    void endIteration()
    {
        // iterate over all wells and notify them individually
        const size_t wellSize = wells_.size();
        for (size_t wellIdx = 0; wellIdx < wellSize; ++wellIdx)
            wells_[wellIdx]->endIteration();
    }

    /*!
     * \brief Informs the well manager that a time step has just been finished.
     */
    void endTimeStep()
    {
        Scalar dt = simulator_.timeStepSize();

        // iterate over all wells and notify them individually. also, update the
        // production/injection totals for the active wells.
        const size_t wellSize = wells_.size();
        for (size_t wellIdx = 0; wellIdx < wellSize; ++wellIdx) {
            auto well = wells_[wellIdx];
            well->endTimeStep();

            // update the surface volumes of the produced/injected fluids
            std::array<Scalar, numPhases>* injectedVolume;
            if (wellTotalInjectedVolume_.count(well->name()) == 0) {
                injectedVolume = &wellTotalInjectedVolume_[well->name()];
                std::fill(injectedVolume->begin(), injectedVolume->end(), 0.0);
            }
            else
                injectedVolume = &wellTotalInjectedVolume_[well->name()];

            std::array<Scalar, numPhases>* producedVolume;
            if (wellTotalProducedVolume_.count(well->name()) == 0) {
                producedVolume = &wellTotalProducedVolume_[well->name()];
                std::fill(producedVolume->begin(), producedVolume->end(), 0.0);
            }
            else
                producedVolume = &wellTotalProducedVolume_[well->name()];

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // this assumes that the implicit Euler method is used for time
                // integration. TODO: Once the time discretization becomes pluggable,
                // this integration needs to be done by the time discretization code!
                Scalar vol = dt * well->surfaceRate(phaseIdx);

                if (vol < 0)
                    (*producedVolume)[phaseIdx] += -vol;
                else
                    (*injectedVolume)[phaseIdx] += vol;
            }
        }
    }

    /*!
     * \brief Informs the well manager that a simulation episode has just been finished.
     */
    void endEpisode()
    { }

    /*!
     * \brief Returns the surface volume of a fluid phase produced by a well.
     */
    Scalar totalProducedVolume(const std::string& wellName, unsigned phaseIdx) const
    {
        if (wellTotalProducedVolume_.count(wellName) == 0)
            return 0.0; // well not yet seen
        return wellTotalProducedVolume_.at(wellName)[phaseIdx];
    }

    /*!
     * \brief Returns the surface volume of a fluid phase injected by a well.
     */
    Scalar totalInjectedVolume(const std::string& wellName, unsigned phaseIdx) const
    {
        if (wellTotalInjectedVolume_.count(wellName) == 0)
            return 0.0; // well not yet seen
        return wellTotalInjectedVolume_.at(wellName)[phaseIdx];
    }

    /*!
     * \brief Computes the source term due to wells for a degree of
     *        freedom.
     */
    template <class Context>
    void computeTotalRatesForDof(EvalEqVector& q,
                                 const Context& context,
                                 unsigned dofIdx,
                                 unsigned timeIdx) const
    {
        q = 0.0;

        if (!gridDofIsPenetrated(context.globalSpaceIndex(dofIdx, timeIdx)))
            return;

        RateVector wellRate;

        // iterate over all wells and add up their individual rates
        const size_t wellSize = wells_.size();
        for (size_t wellIdx = 0; wellIdx < wellSize; ++wellIdx) {
            wellRate = 0.0;
            wells_[wellIdx]->computeTotalRatesForDof(wellRate, context, dofIdx, timeIdx);
            for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx)
                q[eqIdx] += wellRate[eqIdx];
        }
    }

    data::Wells wellData() const
    {
        data::Wells wellDat;

        using rt = data::Rates::opt;
        for (unsigned wellIdx = 0; wellIdx < numWells(); ++wellIdx) {
            const auto& ebosWell = well(wellIdx);
            auto& wellOut = wellDat[ebosWell->name()];

            wellOut.bhp = ebosWell->bottomHolePressure();
            wellOut.thp = ebosWell->tubingHeadPressure();
            wellOut.temperature = 0;
            wellOut.rates.set( rt::wat, ebosWell->surfaceRate(waterPhaseIdx) );
            wellOut.rates.set( rt::oil, ebosWell->surfaceRate(oilPhaseIdx) );
            wellOut.rates.set( rt::gas, ebosWell->surfaceRate(gasPhaseIdx) );

            const int numConnections = ebosWell->numConnections();
            wellOut.connections.resize(numConnections);

            for( int i = 0; i < numConnections; ++i ) {
                auto& connection = wellOut.connections[ i ];
                connection.index = 0;
                connection.pressure = 0.0;
                connection.reservoir_rate = 0.0;
                connection.rates.set( rt::wat, 0.0 );
                connection.rates.set( rt::oil, 0.0 );
                connection.rates.set( rt::gas, 0.0 );
            }
        }

        return wellDat;
    }

    data::GroupAndNetworkValues
    groupAndNetworkData(const int /* reportStepIdx */)
    {
        return {};
    }

    /*!
     * \brief This method writes the complete state of all wells
     *        to the hard disk.
     */
    template <class Restarter>
    void serialize(Restarter&)
    {
        /* do nothing: Everything which we need here is provided by the deck->.. */
    }

    /*!
     * \brief This method restores the complete state of the all wells
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     */
    template <class Restarter>
    void deserialize(Restarter&)
    {
        // initialize the wells for the current episode
        beginEpisode(/*wasRestarted=*/true);
    }

    /*!
     * \brief Returns true if something in a well changed compared to the previous report
     *        step.
     *
     * "Something" can either be the well topology (i.e., which grid blocks are contained
     * in which well) or it can be a well parameter like the bottom hole pressure...
     */
    bool wellsChanged(const EclipseState& eclState,
                      const Schedule& schedule,
                      unsigned reportStepIdx) const
    {
        if (wellTopologyChanged_(eclState, reportStepIdx))
            return true;

        if ((schedule.size() - 1) <= (unsigned) reportStepIdx)
            // for the "until the universe dies" episode, the wells don't change
            return false;

        const Events& events = schedule[reportStepIdx].events();
        return events.hasEvent(ScheduleEvents::PRODUCTION_UPDATE |
                               ScheduleEvents::INJECTION_UPDATE |
                               ScheduleEvents::WELL_STATUS_CHANGE );
    }

    void initFromRestartFile(const RestartValue&) {
        throw std::logic_error("initFromRestartFile() method not implemented for class eclwellmanager");
    }

    const WellState& wellState() const
    {
        throw std::logic_error("wellState() method not implemented for class eclwellmanager");
    }

    WellState& wellState()
    {
        throw std::logic_error("wellState() method not implemented for class eclwellmanager");
    }

    void commitWGState()
    {
        throw std::logic_error("commitWellState() method not implemented for class eclwellmanager");
    }

    void commitWGState(WGState)
    {
        throw std::logic_error("commitWGState() method not implemented for class eclwellmanager");
    }

    WellTestState& wellTestState() {
        throw std::logic_error("wellTestState() method not implemented for class eclwellmanager");
    }



    void resetWGState()
    {
        throw std::logic_error("resetWGState() method not implemented for class eclwellmanager");
    }

    void updateNupcolWGState()
    {
        throw std::logic_error("updateNupcolWGState() method not implemented for class eclwellmanager");
    }

    void
    updateEclWell(int, int)
    {
        throw std::logic_error("wellPI() method not implemented for class eclwellmanager");
    }


    void
    updateEclWells(int, const std::unordered_set<std::string>&) {
        throw std::logic_error("wellPI() method not implemented for class eclwellmanager");
    }


    double
    wellPI(int ) const
    {
        throw std::logic_error("wellPI() method not implemented for class eclwellmanager");
    }

    double
    wellPI(const std::string& ) const
    {
        throw std::logic_error("wellPI() method not implemented for class eclwellmanager");
    }





protected:
    bool wellTopologyChanged_(const EclipseState&,
                              const Schedule& schedule,
                              unsigned reportStepIdx) const
    {
        if (reportStepIdx == 0) {
            // the well topology has always been changed relative to before the
            // simulation is started...
            return true;
        }

        if ((schedule.size() - 1) <= (unsigned) reportStepIdx)
            // for the "until the universe dies" episode, the wells don't change
            return false;

        const Events& events = schedule[reportStepIdx].events();
        return events.hasEvent(ScheduleEvents::NEW_WELL |
                               ScheduleEvents::COMPLETION_CHANGE);
    }

    void updateWellTopology_(unsigned,
                             const WellConnectionsMap& wellConnections,
                             std::vector<bool>& gridDofIsPenetrated) const
    {
        auto& model = simulator_.model();
        const auto& vanguard = simulator_.vanguard();

        // first, remove all wells from the reservoir
        model.clearAuxiliaryModules();
        auto wellIt = wells_.begin();
        const auto& wellEndIt = wells_.end();
        for (; wellIt != wellEndIt; ++wellIt) {
            (*wellIt)->clear();
            (*wellIt)->beginSpec();
        }

        //////
        // tell the active wells which DOFs they contain
        const auto gridView = simulator_.vanguard().gridView();

        gridDofIsPenetrated.resize(model.numGridDof());
        std::fill(gridDofIsPenetrated.begin(), gridDofIsPenetrated.end(), false);

        ElementContext elemCtx(simulator_);
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto elemEndIt = gridView.template end</*codim=*/0>();
        std::set<std::shared_ptr<Well> > wells;
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue; // non-local entities need to be skipped

            elemCtx.updateStencil(elem);
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++ dofIdx) {
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                unsigned cartesianDofIdx = vanguard.cartesianIndex(globalDofIdx);

                if (wellConnections.count(cartesianDofIdx) == 0)
                    // the current DOF is not contained in any well, so we must skip
                    // it...
                    continue;

                gridDofIsPenetrated[globalDofIdx] = true;

                auto eclWell = wellConnections.at(cartesianDofIdx).second;
                eclWell->addDof(elemCtx, dofIdx);

                wells.insert(eclWell);
            }
            //////
        }

        // register all wells at the model as auxiliary equations
        wellIt = wells_.begin();
        for (; wellIt != wellEndIt; ++wellIt) {
            (*wellIt)->endSpec();
            model.addAuxiliaryModule(wellIt->get());
        }
    }

    void computeWellConnectionsMap_(unsigned reportStepIdx,
                                    WellConnectionsMap& cartesianIdxToConnectionMap)
    {
        const auto& deckSchedule = simulator_.vanguard().schedule();

#ifndef NDEBUG
        const auto& eclState = simulator_.vanguard().eclState();
        const auto& eclGrid = eclState.getInputGrid();
        assert( int(eclGrid.getNX()) == simulator_.vanguard().cartesianDimensions()[ 0 ] );
        assert( int(eclGrid.getNY()) == simulator_.vanguard().cartesianDimensions()[ 1 ] );
        assert( int(eclGrid.getNZ()) == simulator_.vanguard().cartesianDimensions()[ 2 ] );
#endif

        // compute the mapping from logically Cartesian indices to the well the
        // respective connection.
        const auto deckWells = deckSchedule.getWells(reportStepIdx);
        for (const auto& deckWell : deckWells) {
            const std::string& wellName = deckWell.name();

            if (!hasWell(wellName))
            {
#ifndef NDEBUG
                if( simulator_.vanguard().grid().comm().size() == 1 )
                {
                    std::cout << "Well '" << wellName << "' suddenly appears in the connection "
                              << "for the report step, but has not been previously specified. "
                              << "Ignoring.\n";
                }
#endif
                continue;
            }

            std::array<int, 3> cartesianCoordinate;
            // set the well parameters defined by the current set of connections
            const auto& connectionSet = deckWell.getConnections();
            for (size_t connIdx = 0; connIdx < connectionSet.size(); connIdx ++) {
                const auto& connection = connectionSet.get(connIdx);
                cartesianCoordinate[ 0 ] = connection.getI();
                cartesianCoordinate[ 1 ] = connection.getJ();
                cartesianCoordinate[ 2 ] = connection.getK();
                unsigned cartIdx = simulator_.vanguard().cartesianIndex( cartesianCoordinate );

                // in this code we only support each cell to be part of at most a single
                // well. TODO (?) change this?
                assert(cartesianIdxToConnectionMap.count(cartIdx) == 0);

                auto eclWell = wells_[wellIndex(wellName)];
                cartesianIdxToConnectionMap[cartIdx] = std::make_pair(&connection, eclWell);
            }
        }
    }


    void updateWellParameters_(unsigned reportStepIdx, const WellConnectionsMap& wellConnections)
    {
        const auto& deckSchedule = simulator_.vanguard().schedule();
        const auto deckWells = deckSchedule.getWells(reportStepIdx);

        // set the reference depth for all wells
        for (const auto& deckWell : deckWells) {
            const std::string& wellName = deckWell.name();

            if( hasWell( wellName ) )
            {
                wells_[wellIndex(wellName)]->clear();
                wells_[wellIndex(wellName)]->setReferenceDepth(deckWell.getRefDepth());
            }
        }

        // associate the well connections with grid cells and register them in the
        // Peaceman well object
        const auto& vanguard = simulator_.vanguard();
        const GridView gridView = vanguard.gridView();

        ElementContext elemCtx(simulator_);
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto elemEndIt = gridView.template end</*codim=*/0>();

        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue; // non-local entities need to be skipped

            elemCtx.updateStencil(elem);
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++ dofIdx)
            {
                assert( elemCtx.numPrimaryDof(/*timeIdx=*/0) == 1 );
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                unsigned cartesianDofIdx = vanguard.cartesianIndex(globalDofIdx);

                if (wellConnections.count(cartesianDofIdx) == 0)
                    // the current DOF is not contained in any well, so we must skip
                    // it...
                    continue;

                const auto& connInfo = wellConnections.at(cartesianDofIdx);
                const Connection* connection = connInfo.first;
                std::shared_ptr<Well> eclWell = connInfo.second;
                eclWell->addDof(elemCtx, dofIdx);
                eclWell->setConnectionTransmissibilityFactor(elemCtx, dofIdx, connection->CF());
                eclWell->setRadius(elemCtx, dofIdx, connection->rw());
                //eclWell->setEffectivePermeability(elemCtx, dofIdx, connection->Kh());
            }
        }
    }

    Simulator& simulator_;

    std::vector<std::shared_ptr<Well> > wells_;
    std::vector<bool> gridDofIsPenetrated_;
    std::map<std::string, int> wellNameToIndex_;
    std::map<std::string, std::array<Scalar, numPhases> > wellTotalInjectedVolume_;
    std::map<std::string, std::array<Scalar, numPhases> > wellTotalProducedVolume_;
};
} // namespace Opm

#endif
