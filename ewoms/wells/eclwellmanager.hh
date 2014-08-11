/*
  Copyright (C) 2014 by Andreas Lauser

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

#include <ewoms/wells/eclpeacemanwell.hh>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/CompletionSet.hpp>

#include <opm/core/utility/PropertySystem.hpp>

#include <dune/grid/common/gridenums.hh>

#include <map>
#include <string>
#include <vector>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(Grid);
}}

namespace Ewoms {
/*!
 * \brief A class which handles well controls as specified by an
 *        Eclipse deck
 */
template <class TypeTag>
class EclWellManager
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

    typedef Ewoms::EclPeacemanWell<TypeTag> Well;

public:
    EclWellManager(const Simulator &simulator)
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
        const Grid &grid = simulator_.gridManager().grid();
        const GridView gridView = simulator_.gridManager().gridView();

        // create the wells
        for (size_t deckWellIdx = 0; deckWellIdx < deckSchedule->numWells(); ++deckWellIdx) {
            Opm::WellConstPtr deckWell = deckSchedule->getWells()[deckWellIdx];
            const std::string &wellName = deckWell->name();

            std::shared_ptr<Well> well(new Well(simulator_));
            wellNameToIndex_[wellName] = wells_.size();
            wells_.push_back(well);

            // specify the DOFs directly affected by the
            // well. Probably this could be done quite a bit more
            // efficiently, but for now it should be Fast Enough (TM).
            well->beginSpec();

            well->setName(wellName);

            ElementContext elemCtx(simulator_);
            auto elemIt = gridView.template begin</*codim=*/0>();
            const auto elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                if (elemIt->partitionType() != Dune::InteriorEntity)
                    continue;

                elemCtx.updateStencil(elemIt);
                for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++ dofIdx) {
                    int globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                    std::array<int,3> ijk;
                    // if the compiler complains here, you're not
                    // using Dune::CpGrid. Other grids are not
                    // supported by the EclWellsManager, sorry.
                    grid.getIJK(globalDofIdx, ijk);

                    // TODO: time dependent wells (i.e. move this code into the
                    // beginEpisode() method!?)
                    Opm::CompletionSetConstPtr completionSet =
                        deckWell->getCompletions(/*timeStepIdx=*/0);
                    for (size_t complIdx = 0; complIdx < completionSet->size(); complIdx ++) {
                        Opm::CompletionConstPtr completion =
                            completionSet->get(complIdx);
                        if (ijk[0] == completion->getI()
                            && ijk[1] == completion->getJ()
                            && ijk[2] == completion->getK())
                        {
                            well->addDof(elemCtx, dofIdx);
                            well->setRadius(elemCtx, dofIdx, 0.5*completion->getDiameter());
                        }
                    }
                }
            }

            well->endSpec();
        }
    }

    /*!
     * \brief This should be called the problem before each simulation
     *        episode to adapt the well controls.
     */
    void beginEpisode(Opm::EclipseStateConstPtr eclState)
    {
        int episodeIdx = simulator_.episodeIndex();

        const auto &deckSchedule = eclState->getSchedule();

        // set the injection data for the respective wells.
        for (size_t deckWellIdx = 0; deckWellIdx < deckSchedule->numWells(); ++deckWellIdx) {
            Opm::WellConstPtr deckWell = deckSchedule->getWells()[deckWellIdx];

            if (!hasWell(deckWell->name()))
                continue;

            auto well = this->well(deckWell->name());

            Opm::WellCommon::StatusEnum deckWellStatus = deckWell->getStatus(episodeIdx);
            switch (deckWellStatus) {
            case Opm::WellCommon::AUTO:
                // TODO: for now, auto means open...
            case Opm::WellCommon::OPEN:
                well->setOpen(true);
                break;
            case Opm::WellCommon::STOP:
                // TODO: cross flow
            case Opm::WellCommon::SHUT:
                well->setOpen(false);
                break;
            }

            // make sure that the well is either an injector or a
            // producer for the current episode. (it is not allowed to
            // be neither or to be both...)
            assert( (deckWell->isInjector(episodeIdx)?1:0) +
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
                    well->setControlMode(Well::ControlMode::TopHolePressure);
                    break;

                case Opm::WellInjector::GRUP:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Well groups");
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
                well->setTargetTopHolePressure(injectProperties.THPLimit);
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
                    well->setControlMode(Well::ControlMode::TopHolePressure);
                    break;

                case Opm::WellProducer::GRUP:
                    OPM_THROW(std::runtime_error,
                              "Not implemented: Well groups");
                }

                well->setTargetBottomHolePressure(producerProperties.BHPLimit);
                well->setTargetTopHolePressure(producerProperties.THPLimit);
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
    { return wellNameToIndex_.count(wellName) > 0; }

    /*!
     * \brief Given a well name, return the corresponding index.
     *
     * A std::runtime_error will be thrown if the well name is unknown.
     */
    int wellIndex(const std::string &wellName) const
    {
        const auto &it = wellNameToIndex_.find(wellName);
        if (it == wellNameToIndex_.end())
            OPM_THROW(std::runtime_error,
                      "No well called '" << wellName << "'found");
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

        // TODO: adapt well controls
    }

    /*!
     * \brief Informs the well that an iteration has just begun.
     *
     * In this method, the well calculates the bottom and top hole
     * pressures, the actual unconstraint production and injection
     * rates, etc.
     */
    void beginIteration()
    {
        // call the preprocessing routines
        for (size_t wellIdx = 0; wellIdx < wells_.size(); ++wellIdx)
            wells_[wellIdx]->beginIterationPreProcess();

        // call the accumulation routines
        ElementContext elemCtx(simulator_);
        auto elemIt = simulator_.gridManager().gridView().template begin</*codim=*/0>();
        const auto &elemEndIt = simulator_.gridManager().gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updateStencil(*elemIt);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            for (size_t wellIdx = 0; wellIdx < wells_.size(); ++wellIdx)
                wells_[wellIdx]->beginIterationAccumulate(elemCtx, /*timeIdx=*/0);
        }

        // call the postprocessing routines
        for (size_t wellIdx = 0; wellIdx < wells_.size(); ++wellIdx)
            wells_[wellIdx]->beginIterationPostProcess();
    }

    /*!
     * \brief Informs the well manager that an iteration has just been finished.
     */
    void endIteration()
    {
        // iterate over all wells and notify them individually
        for (size_t wellIdx = 0; wellIdx < wells_.size(); ++wellIdx)
            wells_[wellIdx]->endIteration();
    }

    /*!
     * \brief Informs the well manager that a time step has just been finished.
     */
    void endTimeStep()
    {
        // iterate over all wells and notify them individually
        for (size_t wellIdx = 0; wellIdx < wells_.size(); ++wellIdx)
            wells_[wellIdx]->endTimeStep();
    }

    /*!
     * \brief Informs the well manager that a simulation episode has just been finished.
     */
    void endEpisode()
    { }

    /*!
     * \brief Computes the source term due to wells for a degree of
     *        freedom.
     */
    template <class Context>
    void computeTotalRatesForDof(RateVector &q,
                                 const Context &context,
                                 int dofIdx,
                                 int timeIdx) const
    {
        q = 0.0;

        RateVector wellRate;

        // iterate over all wells and add up their individual rates
        for (size_t wellIdx = 0; wellIdx < wells_.size(); ++wellIdx) {
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
        // iterate over all wells and serialize them individually
        for (int wellIdx = 0; wellIdx < wells_.size(); ++wellIdx)
            wells_[wellIdx]->serialize(res);
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
        // iterate over all wells and deserialize them individually
        for (int wellIdx = 0; wellIdx < wells_.size(); ++wellIdx) {
            std::shared_ptr<Well> well(new Well(simulator_));

            well->deserialize(res);

            wells_.push_back(well);
            wellNameToIndex_[well->name()] = wells_.size() - 1;
        }
    }

protected:
    const Simulator &simulator_;

    std::vector<std::shared_ptr<Well> > wells_;
    std::map<std::string, int> wellNameToIndex_;
};
} // namespace Ewoms

#endif
