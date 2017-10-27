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
 * \copydoc Ewoms::EclWellManager
 */
#ifndef EWOMS_ECL_WELL_MANAGER_HH
#define EWOMS_ECL_WELL_MANAGER_HH

#include "eclpeacemanwell.hh"

#include <ewoms/disc/common/fvbaseproperties.hh>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Events.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/CompletionSet.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/TimeMap.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <ewoms/common/propertysystem.hh>
#include <ewoms/parallel/threadedentityiterator.hh>

#include <dune/grid/common/gridenums.hh>

#include <map>
#include <string>
#include <vector>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Grid);
}

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief A class which handles well controls as specified by an
 *        Eclipse deck
 */
template <class TypeTag>
class EclWellManager
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = FluidSystem::numPhases };

    typedef typename GridView::template Codim<0>::Entity Element;

    typedef Ewoms::EclPeacemanWell<TypeTag> Well;

    typedef std::map<int, std::pair<const Opm::Completion*, std::shared_ptr<Well> > > WellCompletionsMap;

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
    void init(const Opm::EclipseState& eclState,
              const Opm::Schedule& deckSchedule)
    {

        // create the wells which intersect with the current process' grid
        for (size_t deckWellIdx = 0; deckWellIdx < deckSchedule.numWells(); ++deckWellIdx)
        {
            const Opm::Well* deckWell = deckSchedule.getWells()[deckWellIdx];
            const std::string& wellName = deckWell->name();

            // set the name of the well but not much else. (i.e., if it is not completed,
            // the well primarily serves as a placeholder.) The big rest of the well is
            // specified by the updateWellCompletions_() method
            auto well = std::make_shared<Well>(simulator_);
            well->setName(wellName);
            well->setWellStatus(Well::Shut);

            wells_.push_back(well);
            wellNameToIndex_[well->name()] = wells_.size() - 1;
        }
    }

    /*!
     * \brief This should be called the problem before each simulation
     *        episode to adapt the well controls.
     */
    void beginEpisode(const Opm::EclipseState& eclState, const Opm::Schedule& deckSchedule, bool wasRestarted=false)
    {
        unsigned episodeIdx = simulator_.episodeIndex();

        WellCompletionsMap wellCompMap;
        computeWellCompletionsMap_(episodeIdx, wellCompMap);

        if (wasRestarted || wellTopologyChanged_(eclState, deckSchedule, episodeIdx))
            updateWellTopology_(episodeIdx, wellCompMap, gridDofIsPenetrated_);

        // set those parameters of the wells which do not change the topology of the
        // linearized system of equations
        updateWellParameters_(episodeIdx, wellCompMap);

        const std::vector<const Opm::Well*>& deckWells = deckSchedule.getWells(episodeIdx);
        // set the injection data for the respective wells.
        for (size_t deckWellIdx = 0; deckWellIdx < deckWells.size(); ++deckWellIdx) {
            const Opm::Well* deckWell = deckWells[deckWellIdx];

            if (!hasWell(deckWell->name()))
                continue;

            auto well = this->well(deckWell->name());

            Opm::WellCommon::StatusEnum deckWellStatus = deckWell->getStatus(episodeIdx);
            switch (deckWellStatus) {
            case Opm::WellCommon::AUTO:
                // TODO: for now, auto means open...
            case Opm::WellCommon::OPEN:
                well->setWellStatus(Well::Open);
                break;
            case Opm::WellCommon::STOP:
                well->setWellStatus(Well::Closed);
                break;
            case Opm::WellCommon::SHUT:
                well->setWellStatus(Well::Shut);
                break;
            }

            // make sure that the well is either an injector or a
            // producer for the current episode. (it is not allowed to
            // be neither or to be both...)
            assert((deckWell->isInjector(episodeIdx)?1:0) +
                   (deckWell->isProducer(episodeIdx)?1:0) == 1);

            if (deckWell->isInjector(episodeIdx)) {
                well->setWellType(Well::Injector);

                const Opm::WellInjectionProperties& injectProperties =
                    deckWell->getInjectionProperties(episodeIdx);

                switch (injectProperties.injectorType) {
                case Opm::WellInjector::WATER:
                    well->setInjectedPhaseIndex(FluidSystem::waterPhaseIdx);
                    break;
                case Opm::WellInjector::GAS:
                    well->setInjectedPhaseIndex(FluidSystem::gasPhaseIdx);
                    break;
                case Opm::WellInjector::OIL:
                    well->setInjectedPhaseIndex(FluidSystem::oilPhaseIdx);
                    break;
                case Opm::WellInjector::MULTI:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Multi-phase injector wells");
                }

                switch (injectProperties.controlMode) {
                case Opm::WellInjector::RATE:
                    well->setControlMode(Well::ControlMode::VolumetricSurfaceRate);
                    break;

                case Opm::WellInjector::RESV:
                    well->setControlMode(Well::ControlMode::VolumetricReservoirRate);
                    break;

                case Opm::WellInjector::BHP:
                    well->setControlMode(Well::ControlMode::BottomHolePressure);
                    break;

                case Opm::WellInjector::THP:
                    well->setControlMode(Well::ControlMode::TubingHeadPressure);
                    break;

                case Opm::WellInjector::GRUP:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Well groups");

                case Opm::WellInjector::CMODE_UNDEFINED:
                    std::cout << "Warning: Control mode of injection well " << well->name()
                              << " is undefined. Assuming well to be shut.\n";
                    well->setWellStatus(Well::WellStatus::Shut);
                    continue;
                }

                switch (injectProperties.injectorType) {
                case Opm::WellInjector::WATER:
                    well->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/0.0, /*water=*/1.0);
                    break;

                case Opm::WellInjector::OIL:
                    well->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/0.0, /*water=*/0.0);
                    break;

                case Opm::WellInjector::GAS:
                    well->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/1.0, /*water=*/0.0);
                    break;

                case Opm::WellInjector::MULTI:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Multi-phase injection wells");
                }

                well->setMaximumSurfaceRate(injectProperties.surfaceInjectionRate);
                well->setMaximumReservoirRate(injectProperties.reservoirInjectionRate);
                well->setTargetBottomHolePressure(injectProperties.BHPLimit);

                // TODO
                well->setTargetTubingHeadPressure(1e30);
                //well->setTargetTubingHeadPressure(injectProperties.THPLimit);
            }

            if (deckWell->isProducer(episodeIdx)) {
                well->setWellType(Well::Producer);

                const Opm::WellProductionProperties& producerProperties =
                    deckWell->getProductionProperties(episodeIdx);

                switch (producerProperties.controlMode) {
                case Opm::WellProducer::ORAT:
                    well->setControlMode(Well::ControlMode::VolumetricSurfaceRate);
                    well->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/0.0, /*water=*/0.0);
                    well->setMaximumSurfaceRate(producerProperties.OilRate);
                    break;

                case Opm::WellProducer::GRAT:
                    well->setControlMode(Well::ControlMode::VolumetricSurfaceRate);
                    well->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/1.0, /*water=*/0.0);
                    well->setMaximumSurfaceRate(producerProperties.GasRate);
                    break;

                case Opm::WellProducer::WRAT:
                    well->setControlMode(Well::ControlMode::VolumetricSurfaceRate);
                    well->setVolumetricPhaseWeights(/*oil=*/0.0, /*gas=*/0.0, /*water=*/1.0);
                    well->setMaximumSurfaceRate(producerProperties.WaterRate);
                    break;

                case Opm::WellProducer::LRAT:
                    well->setControlMode(Well::ControlMode::VolumetricSurfaceRate);
                    well->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/0.0, /*water=*/1.0);
                    well->setMaximumSurfaceRate(producerProperties.LiquidRate);
                    break;

                case Opm::WellProducer::CRAT:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Linearly combined rates");

                case Opm::WellProducer::RESV:
                    well->setControlMode(Well::ControlMode::VolumetricReservoirRate);
                    well->setVolumetricPhaseWeights(/*oil=*/1.0, /*gas=*/1.0, /*water=*/1.0);
                    well->setMaximumSurfaceRate(producerProperties.ResVRate);
                    break;

                case Opm::WellProducer::BHP:
                    well->setControlMode(Well::ControlMode::BottomHolePressure);
                    break;

                case Opm::WellProducer::THP:
                    well->setControlMode(Well::ControlMode::TubingHeadPressure);
                    break;

                case Opm::WellProducer::GRUP:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Well groups");

                case Opm::WellProducer::NONE:
                    // fall-through
                case Opm::WellProducer::CMODE_UNDEFINED:
                    std::cout << "Warning: Control mode of production well " << well->name()
                              << " is undefined. Assuming well to be shut.";
                    well->setWellStatus(Well::WellStatus::Shut);
                    continue;
                }

                well->setTargetBottomHolePressure(producerProperties.BHPLimit);

                // TODO
                well->setTargetTubingHeadPressure(-1e30);
                //well->setTargetTubingHeadPressure(producerProperties.THPLimit);
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
        {
            OPM_THROW(std::runtime_error,
                      "No well called '" << wellName << "'found");
        }
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
        ThreadedEntityIterator<GridView, /*codim=*/0> threadedElemIt(simulator_.gridManager().gridView());
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            ElementContext elemCtx(simulator_);
            auto elemIt = threadedElemIt.beginParallel();
            for (; !threadedElemIt.isFinished(elemIt); elemIt = threadedElemIt.increment()) {
                const Element& elem = *elemIt;
                if (elem.partitionType() != Dune::InteriorEntity)
                    continue;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                for (size_t wellIdx = 0; wellIdx < wellSize; ++wellIdx)
                    wells_[wellIdx]->beginIterationAccumulate(elemCtx, /*timeIdx=*/0);
            }
        }

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

    /*!
     * \brief This method writes the complete state of all wells
     *        to the hard disk.
     */
    template <class Restarter>
    void serialize(Restarter& res OPM_UNUSED)
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
    void deserialize(Restarter& res OPM_UNUSED)
    {
        // initialize the wells for the current episode
        beginEpisode(simulator_.gridManager().eclState(), simulator_.gridManager().schedule(), /*wasRestarted=*/true);
    }

    /*!
     * \brief Returns true if something in a well changed compared to the previous report
     *        step.
     *
     * "Something" can either be the well topology (i.e., which grid blocks are contained
     * in which well) or it can be a well parameter like the bottom hole pressure...
     */
    bool wellsChanged(const Opm::EclipseState& eclState,
                      const Opm::Schedule& schedule,
                      unsigned reportStepIdx) const
    {
        if (wellTopologyChanged_(eclState, reportStepIdx))
            return true;

        if (schedule.getTimeMap().numTimesteps() <= (unsigned) reportStepIdx)
            // for the "until the universe dies" episode, the wells don't change
            return false;

        const Opm::Events& events = schedule.getEvents();
        return events.hasEvent(Opm::ScheduleEvents::PRODUCTION_UPDATE |
                               Opm::ScheduleEvents::INJECTION_UPDATE |
                               Opm::ScheduleEvents::WELL_STATUS_CHANGE,
                               reportStepIdx);
    }

protected:
    bool wellTopologyChanged_(const Opm::EclipseState& eclState,
                              const Opm::Schedule& schedule,
                              unsigned reportStepIdx) const
    {
        if (reportStepIdx == 0) {
            // the well topology has always been changed relative to before the
            // simulation is started...
            return true;
        }

        if (schedule.getTimeMap().numTimesteps() <= (unsigned) reportStepIdx)
            // for the "until the universe dies" episode, the wells don't change
            return false;

        const Opm::Events& events = schedule.getEvents();
        return events.hasEvent(Opm::ScheduleEvents::NEW_WELL |
                               Opm::ScheduleEvents::COMPLETION_CHANGE,
                               reportStepIdx);
    }

    void updateWellTopology_(unsigned reportStepIdx OPM_UNUSED,
                             const WellCompletionsMap& wellCompletions,
                             std::vector<bool>& gridDofIsPenetrated) const
    {
        auto& model = simulator_.model();
        const auto& gridManager = simulator_.gridManager();

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
        const auto gridView = simulator_.gridManager().gridView();

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
                unsigned cartesianDofIdx = gridManager.cartesianIndex(globalDofIdx);

                if (wellCompletions.count(cartesianDofIdx) == 0)
                    // the current DOF is not contained in any well, so we must skip
                    // it...
                    continue;

                gridDofIsPenetrated[globalDofIdx] = true;

                auto eclWell = wellCompletions.at(cartesianDofIdx).second;
                eclWell->addDof(elemCtx, dofIdx);

                wells.insert(eclWell);
            }
            //////
        }

        // register all wells at the model as auxiliary equations
        wellIt = wells_.begin();
        for (; wellIt != wellEndIt; ++wellIt) {
            (*wellIt)->endSpec();
            model.addAuxiliaryModule(*wellIt);
        }
    }

    void computeWellCompletionsMap_(unsigned reportStepIdx OPM_UNUSED, WellCompletionsMap& cartesianIdxToCompletionMap)
    {
        const auto& deckSchedule = simulator_.gridManager().schedule();

#ifndef NDEBUG
        const auto& eclState = simulator_.gridManager().eclState();
        const auto& eclGrid = eclState.getInputGrid();
        assert( int(eclGrid.getNX()) == simulator_.gridManager().cartesianDimensions()[ 0 ] );
        assert( int(eclGrid.getNY()) == simulator_.gridManager().cartesianDimensions()[ 1 ] );
        assert( int(eclGrid.getNZ()) == simulator_.gridManager().cartesianDimensions()[ 2 ] );
#endif

        // compute the mapping from logically Cartesian indices to the well the
        // respective completion.
        const std::vector<const Opm::Well*>& deckWells = deckSchedule.getWells(reportStepIdx);
        for (size_t deckWellIdx = 0; deckWellIdx < deckWells.size(); ++deckWellIdx) {
            const Opm::Well* deckWell = deckWells[deckWellIdx];
            const std::string& wellName = deckWell->name();

            if (!hasWell(wellName))
            {
#ifndef NDEBUG
                if( simulator_.gridManager().grid().comm().size() == 1 )
                {
                    std::cout << "Well '" << wellName << "' suddenly appears in the completions "
                              << "for the report step, but has not been previously specified. "
                              << "Ignoring.\n";
                }
#endif
                continue;
            }

            std::array<int, 3> cartesianCoordinate;
            // set the well parameters defined by the current set of completions
            const auto& completionSet = deckWell->getCompletions(reportStepIdx);
            for (size_t complIdx = 0; complIdx < completionSet.size(); complIdx ++) {
                const auto& completion = completionSet.get(complIdx);
                cartesianCoordinate[ 0 ] = completion.getI();
                cartesianCoordinate[ 1 ] = completion.getJ();
                cartesianCoordinate[ 2 ] = completion.getK();
                unsigned cartIdx = simulator_.gridManager().cartesianIndex( cartesianCoordinate );

                // in this code we only support each cell to be part of at most a single
                // well. TODO (?) change this?
                assert(cartesianIdxToCompletionMap.count(cartIdx) == 0);

                auto eclWell = wells_[wellIndex(wellName)];
                cartesianIdxToCompletionMap[cartIdx] = std::make_pair(&completion, eclWell);
            }
        }
    }

    void updateWellParameters_(unsigned reportStepIdx, const WellCompletionsMap& wellCompletions)
    {
        const auto& deckSchedule = simulator_.gridManager().schedule();
        const std::vector<const Opm::Well*>& deckWells = deckSchedule.getWells(reportStepIdx);

        // set the reference depth for all wells
        for (size_t deckWellIdx = 0; deckWellIdx < deckWells.size(); ++deckWellIdx) {
            const Opm::Well* deckWell = deckWells[deckWellIdx];
            const std::string& wellName = deckWell->name();

            if( hasWell( wellName ) )
            {
                wells_[wellIndex(wellName)]->clear();
                wells_[wellIndex(wellName)]->setReferenceDepth(deckWell->getRefDepth());
            }
        }

        // associate the well completions with grid cells and register them in the
        // Peaceman well object
        const auto& gridManager = simulator_.gridManager();
        const GridView gridView = gridManager.gridView();

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
                unsigned cartesianDofIdx = gridManager.cartesianIndex(globalDofIdx);

                if (wellCompletions.count(cartesianDofIdx) == 0)
                    // the current DOF is not contained in any well, so we must skip
                    // it...
                    continue;

                const auto& compInfo = wellCompletions.at(cartesianDofIdx);
                const Opm::Completion* completion = compInfo.first;
                std::shared_ptr<Well> eclWell = compInfo.second;
                eclWell->addDof(elemCtx, dofIdx);

                // the catch is a hack for a ideosyncrasy of opm-parser with regard to
                // defaults handling: if the deck did not specify a radius for the
                // completion, there seems to be no other way to detect this except for
                // catching the exception
                try {
                    eclWell->setRadius(elemCtx, dofIdx, 0.5*completion->getDiameter());
                }
                catch (const std::logic_error& e)
                {}

                // overwrite the automatically computed effective
                // permeability by the one specified in the deck-> Note: this
                // is not implemented by opm-parser yet...
                /*
                  Scalar Kh = completion->getEffectivePermeability();
                  if (std::isfinite(Kh) && Kh > 0.0)
                      eclWell->setEffectivePermeability(elemCtx, dofIdx, Kh);
                */

                // overwrite the automatically computed connection
                // transmissibilty factor by the one specified in the deck->
                const auto& ctf = completion->getConnectionTransmissibilityFactorAsValueObject();
                if (ctf.hasValue() && ctf.getValue() > 0.0)
                    eclWell->setConnectionTransmissibilityFactor(elemCtx, dofIdx, ctf.getValue());
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
} // namespace Ewoms

#endif
