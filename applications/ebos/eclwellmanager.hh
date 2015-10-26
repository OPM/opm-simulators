// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2014-2015 by Andreas Lauser

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
#include <opm/parser/eclipse/EclipseState/Schedule/CompletionSet.hpp>

#include <ewoms/common/propertysystem.hh>

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
    EclWellManager(Simulator &simulator)
        : simulator_(simulator)
    { }

    /*!
     * \brief This sets up the basic properties of all wells.
     *
     * I.e., well positions, names etc...
     */
    void init(Opm::EclipseStateConstPtr eclState)
    {
        const auto &deckSchedule = eclState->getSchedule();

        typedef std::array<int, Grid::dimension> CartesianCoordinate;

        // vector containing wells used for numerics
        // this list can differ from the overall number of wells
        // due to partitioning issues
        std::map< const std::string, std::pair< std::shared_ptr<Well>, CartesianCoordinate> > wells;

        // create the wells
        for (size_t deckWellIdx = 0; deckWellIdx < deckSchedule->numWells(); ++deckWellIdx)
        {
            Opm::WellConstPtr deckWell = deckSchedule->getWells()[deckWellIdx];
            const std::string &wellName = deckWell->name();
            CartesianCoordinate cartesianCoord;
            cartesianCoord.fill( 0 );
            cartesianCoord[ 0 ] = deckWell->getHeadI();
            cartesianCoord[ 1 ] = deckWell->getHeadJ();

            std::shared_ptr<Well> well(new Well(simulator_));
            // insert well into well map
            wells[ wellName ] = std::make_pair( well, cartesianCoord );

            // set the name of the well but not much else. (i.e., if it is not completed,
            // the well primarily serves as a placeholder.) The big rest of the well is
            // specified by the updateWellCompletions_() method
            well->beginSpec();
            well->setName(wellName);
            well->setWellStatus(Well::Shut);
            well->endSpec();
        }

        const auto gridView = simulator_.gridManager().gridView();
        ElementContext elemCtx(simulator_);
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto elemEndIt = gridView.template end</*codim=*/0>();

        // loop over elements and insert those wells that appear on the current process
        for (; elemIt != elemEndIt; ++elemIt)
        {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue; // non-local entities need to be skipped

            elemCtx.updateStencil(elem);
            for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++ dofIdx)
            {
                const int globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                CartesianCoordinate cartCoord;
                simulator_.gridManager().cartesianCoordinate( globalDofIdx, cartCoord );

                for( auto wellIt = wells.begin(); wellIt != wells.end(); ++wellIt )
                {
                    std::shared_ptr< Well > well = wellIt->second.first;
                    const CartesianCoordinate& wellCoord = wellIt->second.second;
                    if( wellCoord[ 0 ] == cartCoord[ 0 ] && wellCoord[ 1 ] == cartCoord[ 1 ] )
                    {
                        wellNameToIndex_[ well->name() ] = wells_.size();
                        wells_.push_back( well );

                        // only insert a well once and therefore erase it
                        wells.erase( wellIt++ );
                        // if all wells inserted, return
                        if( wells.empty() )
                            return ;
                    }
                }
            }

        }
    }

    /*!
     * \brief This should be called the problem before each simulation
     *        episode to adapt the well controls.
     */
    void beginEpisode(Opm::EclipseStateConstPtr eclState, bool wasRestarted=false)
    {
        int episodeIdx = simulator_.episodeIndex();

        const auto &deckSchedule = eclState->getSchedule();
        WellCompletionsMap wellCompMap;
        computeWellCompletionsMap_(episodeIdx, wellCompMap);

        if (wasRestarted || wellTopologyChanged_(eclState, episodeIdx)) {
            updateWellTopology_(episodeIdx, wellCompMap);
        }

        // set those parameters of the wells which do not change the topology of the
        // linearized system of equations
        updateWellParameters_(episodeIdx, wellCompMap);

        const std::vector<Opm::WellConstPtr>& deckWells = deckSchedule->getWells(episodeIdx);
        // set the injection data for the respective wells.
        for (size_t deckWellIdx = 0; deckWellIdx < deckWells.size(); ++deckWellIdx) {
            Opm::WellConstPtr deckWell = deckWells[deckWellIdx];

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

                const Opm::WellInjectionProperties &injectProperties =
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
                    OPM_THROW(std::runtime_error,
                              "Control mode of well " << well->name() << " is undefined.");
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
                well->setTargetTubingHeadPressure(1e100);
                //well->setTargetTubingHeadPressure(injectProperties.THPLimit);
            }

            if (deckWell->isProducer(episodeIdx)) {
                well->setWellType(Well::Producer);

                const Opm::WellProductionProperties &producerProperties =
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

                case Opm::WellProducer::CMODE_UNDEFINED:
                    OPM_THROW(std::runtime_error,
                              "Control mode of well " << well->name() << " is undefined.");
                }

                well->setTargetBottomHolePressure(producerProperties.BHPLimit);

                // TODO
                well->setTargetTubingHeadPressure(-1e100);
                //well->setTargetTubingHeadPressure(producerProperties.THPLimit);
            }
        }
    }

    /*!
     * \brief Return the number of wells considered by the EclWellManager.
     */
    int numWells() const
    { return wells_.size(); }

    /*!
     * \brief Return if a given well name is known to the wells manager
     */
    bool hasWell(const std::string &wellName) const
    {
        return wellNameToIndex_.find( wellName ) != wellNameToIndex_.end();
    }

    /*!
     * \brief Given a well name, return the corresponding index.
     *
     * A std::runtime_error will be thrown if the well name is unknown.
     */
    int wellIndex(const std::string &wellName) const
    {
        assert( hasWell( wellName ) );
        const auto &it = wellNameToIndex_.find(wellName);
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
    std::shared_ptr<const Well> well(const std::string &wellName) const
    { return wells_[wellIndex(wellName)]; }

    /*!
     * \brief Given a well name, return the corresponding well.
     *
     * A std::runtime_error will be thrown if the well name is unknown.
     */
    std::shared_ptr<Well> well(const std::string &wellName)
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
        auto elemIt = simulator_.gridManager().gridView().template begin</*codim=*/0>();
        const auto &elemEndIt = simulator_.gridManager().gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updateStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            for (size_t wellIdx = 0; wellIdx < wellSize; ++wellIdx)
                wells_[wellIdx]->beginIterationAccumulate(elemCtx, /*timeIdx=*/0);
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

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
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
    Scalar totalProducedVolume(const std::string& wellName, int phaseIdx) const
    {
        if (wellTotalProducedVolume_.count(wellName) == 0)
            return 0.0; // well not yet seen
        return wellTotalProducedVolume_.at(wellName)[phaseIdx];
    }

    /*!
     * \brief Returns the surface volume of a fluid phase injected by a well.
     */
    Scalar totalInjectedVolume(const std::string& wellName, int phaseIdx) const
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
    void computeTotalRatesForDof(EvalEqVector &q,
                                 const Context &context,
                                 int dofIdx,
                                 int timeIdx) const
    {
        for (int eqIdx = 0; eqIdx < numEq; ++ eqIdx)
            q[eqIdx] = 0.0;

        RateVector wellRate;

        // iterate over all wells and add up their individual rates
        const size_t wellSize = wells_.size();
        for (size_t wellIdx = 0; wellIdx < wellSize; ++wellIdx) {
            wellRate = 0.0;
            wells_[wellIdx]->computeTotalRatesForDof(wellRate, context, dofIdx, timeIdx);
            q += wellRate;
        }
    }

    /*!
     * \brief This method writes the complete state of all wells
     *        to the hard disk.
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        /* do nothing: Everything which we need here is provided by the deck... */
    }

    /*!
     * \brief This method restores the complete state of the all wells
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        // initialize the wells for the current episode
        beginEpisode(simulator_.gridManager().eclState(), /*wasRestarted=*/true);
    }

    /*!
     * \brief Returns true if something in a well changed compared to the previous report
     *        step.
     *
     * "Something" can either be the well topology (i.e., which grid blocks are contained
     * in which well) or it can be a well parameter like the bottom hole pressure...
     */
    bool wellsChanged(Opm::EclipseStateConstPtr eclState, int reportStepIdx) const
    {
        if (wellTopologyChanged_(eclState, reportStepIdx))
            return true;

        // this is slightly hacky because it assumes that the object which stores the set
        // of wells which are relevant for a report step does not change if there are no
        // changed well parameters. opm-parser does not guarantee this, but so far it
        // seems to adhere to it...
        auto deckSchedule = eclState->getSchedule();

        if (deckSchedule->getTimeMap()->numTimesteps() <= (unsigned) reportStepIdx)
            // for the "until the universe dies" episode, the wells don't change
            return false;

        const auto& curDeckWells = deckSchedule->getWells(reportStepIdx);
        const auto& prevDeckWells = deckSchedule->getWells(reportStepIdx - 1);

        for (unsigned i = 0; i < curDeckWells.size(); ++i) {
            if (curDeckWells[i] != prevDeckWells[i])
                return true;
        }
        return false;
    }

protected:
    bool wellTopologyChanged_(Opm::EclipseStateConstPtr eclState, int reportStepIdx) const
    {
        if (reportStepIdx == 0) {
            // the well topology has always been changed relative to before the
            // simulation is started
            return true;
        }

        auto deckSchedule = eclState->getSchedule();
        if (deckSchedule->getTimeMap()->numTimesteps() <= (unsigned) reportStepIdx)
            // for the "until the universe dies" episode, the wells don't change
            return false;

        const auto& curDeckWells = deckSchedule->getWells(reportStepIdx);
        const auto& prevDeckWells = deckSchedule->getWells(reportStepIdx - 1);

        if (curDeckWells.size() != prevDeckWells.size())
            // the number of wells changed
            return true;

        auto curWellIt = curDeckWells.begin();
        const auto& curWellEndIt = curDeckWells.end();
        for (; curWellIt != curWellEndIt; ++curWellIt) {
            // find the well in the previous time step
            auto prevWellIt = prevDeckWells.begin();
            const auto& prevWellEndIt = prevDeckWells.end();
            for (; ; ++prevWellIt) {
                if (prevWellIt == prevWellEndIt)
                    // current well has not been featured in previous report step, i.e.,
                    // the well topology has changed...
                    return true;

                if ((*prevWellIt)->name() == (*curWellIt)->name())
                    // the previous report step had a well with the same name as the
                    // current one!
                    break;
            }

            // make sure that the wells exhibit the same completions!
            const auto curCompletionSet = (*curWellIt)->getCompletions(reportStepIdx);
            const auto prevCompletionSet = (*prevWellIt)->getCompletions(reportStepIdx);

            if (curCompletionSet->size() != prevCompletionSet->size())
                // number of completions of the well has changed!
                return true;

            for (size_t curWellComplIdx = 0;
                 curWellComplIdx < curCompletionSet->size();
                 ++ curWellComplIdx)
            {
                Opm::CompletionConstPtr curCompletion = curCompletionSet->get(curWellComplIdx);

                for (size_t prevWellComplIdx = 0;; ++ prevWellComplIdx)
                {
                    if (prevWellComplIdx == prevCompletionSet->size())
                        // a new completion has appeared in the current report step
                        return true;

                    Opm::CompletionConstPtr prevCompletion = prevCompletionSet->get(curWellComplIdx);

                    if (curCompletion->getI() == prevCompletion->getI()
                        && curCompletion->getJ() == prevCompletion->getJ()
                        && curCompletion->getK() == prevCompletion->getK())
                        // completion is present in both wells, look at next completion!
                        break;
                }
            }
        }

        return false;
    }

    void updateWellTopology_(int reportStepIdx, const WellCompletionsMap& wellCompletions)
    {
        auto& model = simulator_.model();
        const auto& gridManager = simulator_.gridManager();

        // first, remove all wells from the reservoir
        model.clearAuxiliaryModules();
        auto wellIt = wells_.begin();
        const auto wellEndIt = wells_.end();
        for (; wellIt != wellEndIt; ++wellIt)
            (*wellIt)->clear();

        // tell the active wells which DOFs they contain
        const auto gridView = simulator_.gridManager().gridView();
        ElementContext elemCtx(simulator_);
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto elemEndIt = gridView.template end</*codim=*/0>();
        std::set<std::shared_ptr<Well> > wells;
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue; // non-local entities need to be skipped

            elemCtx.updateStencil(elem);
            for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++ dofIdx) {
                int globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                int cartesianDofIdx = gridManager.cartesianIndex(globalDofIdx);

                if (wellCompletions.count(cartesianDofIdx) == 0)
                    // the current DOF is not contained in any well, so we must skip
                    // it...
                    continue;

                auto eclWell = wellCompletions.at(cartesianDofIdx).second;
                eclWell->addDof(elemCtx, dofIdx);

                wells.insert(eclWell);
            }
        }

        // register all wells at the model as auxiliary equations
        auto wellIt2 = wells.begin();
        const auto& wellEndIt2 = wells.end();
        for (; wellIt2 != wellEndIt2; ++wellIt2)
        {
            model.addAuxiliaryModule(*wellIt2);
        }
    }

    void computeWellCompletionsMap_(int reportStepIdx, WellCompletionsMap& cartesianIdxToCompletionMap)
    {
        auto eclStatePtr = simulator_.gridManager().eclState();
        auto deckSchedule = eclStatePtr->getSchedule();
        auto eclGrid = eclStatePtr->getEclipseGrid();

        assert( int(eclGrid->getNX()) == simulator_.gridManager().cartesianDimensions()[ 0 ] );
        assert( int(eclGrid->getNY()) == simulator_.gridManager().cartesianDimensions()[ 1 ] );

        // compute the mapping from logically Cartesian indices to the well the
        // respective completion.
        const std::vector<Opm::WellConstPtr>& deckWells = deckSchedule->getWells(reportStepIdx);
        for (size_t deckWellIdx = 0; deckWellIdx < deckWells.size(); ++deckWellIdx) {
            Opm::WellConstPtr deckWell = deckWells[deckWellIdx];
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
            Opm::CompletionSetConstPtr completionSet = deckWell->getCompletions(reportStepIdx);
            for (size_t complIdx = 0; complIdx < completionSet->size(); complIdx ++) {
                Opm::CompletionConstPtr completion = completionSet->get(complIdx);
                cartesianCoordinate[ 0 ] = completion->getI();
                cartesianCoordinate[ 1 ] = completion->getJ();
                cartesianCoordinate[ 2 ] = completion->getK();
                const int cartIdx = simulator_.gridManager().cartesianIndex( cartesianCoordinate );
                assert( cartIdx == (completion->getI() + completion->getJ()*eclGrid->getNX()
                                      + completion->getK()*eclGrid->getNX()*eclGrid->getNY() ) );

                // in this code we only support each cell to be part of at most a single
                // well. TODO (?) change this?
                assert(cartesianIdxToCompletionMap.count(cartIdx) == 0);

                auto eclWell = wells_[wellIndex(wellName)];
                cartesianIdxToCompletionMap[cartIdx] = std::make_pair(&(*completion), eclWell);
            }
        }
    }

    void updateWellParameters_(int reportStepIdx, const WellCompletionsMap& wellCompletions)
    {
        auto eclStatePtr = simulator_.gridManager().eclState();
        auto deckSchedule = eclStatePtr->getSchedule();
        const std::vector<Opm::WellConstPtr>& deckWells = deckSchedule->getWells(reportStepIdx);

        // set the reference depth for all wells
        for (size_t deckWellIdx = 0; deckWellIdx < deckWells.size(); ++deckWellIdx) {
            Opm::WellConstPtr deckWell = deckWells[deckWellIdx];
            const std::string& wellName = deckWell->name();

            if( hasWell( wellName ) )
            {
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
            for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++ dofIdx)
            {
                assert( elemCtx.numPrimaryDof(/*timeIdx=*/0) == 1 );
                int globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                int cartesianDofIdx = gridManager.cartesianIndex(globalDofIdx);

                if (wellCompletions.count(cartesianDofIdx) == 0)
                    // the current DOF is not contained in any well, so we must skip
                    // it...
                    continue;

                const auto& compInfo = wellCompletions.at(cartesianDofIdx);
                const Opm::Completion* completion = compInfo.first;
                std::shared_ptr<Well> eclWell = compInfo.second;

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
                // permeability by the one specified in the deck. Note: this
                // is not implemented by opm-parser yet...
                /*
                  Scalar Kh = completion->getEffectivePermeability();
                  if (std::isfinite(Kh) && Kh > 0.0)
                      eclWell->setEffectivePermeability(elemCtx, dofIdx, Kh);
                */

                // overwrite the automatically computed connection
                // transmissibilty factor by the one specified in the deck.
                const auto& ctf = completion->getConnectionTransmissibilityFactorAsValueObject();
                if (ctf.hasValue() && ctf.getValue() > 0.0)
                    eclWell->setConnectionTransmissibilityFactor(elemCtx, dofIdx, ctf.getValue());
            }
        }
    }

    Simulator &simulator_;

    std::vector<std::shared_ptr<Well> > wells_;
    std::map<std::string, int> wellNameToIndex_;
    std::map<std::string, std::array<Scalar, numPhases> > wellTotalInjectedVolume_;
    std::map<std::string, std::array<Scalar, numPhases> > wellTotalProducedVolume_;
};
} // namespace Ewoms

#endif
