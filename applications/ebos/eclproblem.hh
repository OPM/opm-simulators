// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
/*!
 * \file
 *
 * \copydoc Ewoms::EclProblem
 */
#ifndef EWOMS_ECL_PROBLEM_HH
#define EWOMS_ECL_PROBLEM_HH

#include <opm/material/localad/Evaluation.hpp>

#include "eclgridmanager.hh"
#include "eclwellmanager.hh"
#include "eclequilinitializer.hh"
#include "eclwriter.hh"
#include "eclsummarywriter.hh"
#include "ecloutputblackoilmodule.hh"
#include "ecltransmissibility.hh"
#include "ecldummygradientcalculator.hh"
#include "eclfluxmodule.hh"
#include "ecldeckunits.hh"

#include <ewoms/models/blackoil/blackoilmodel.hh>
#include <ewoms/disc/ecfv/ecfvdiscretization.hh>
#include <ewoms/io/polyhedralgridconverter.hh>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WetGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DeadOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityWaterPvt.hpp>

// for this simulator to make sense, dune-cornerpoint and opm-parser
// must be available
#include <dune/grid/CpGrid.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <boost/date_time.hpp>

#include <vector>
#include <string>

namespace Ewoms {
template <class TypeTag>
class EclProblem;

namespace Properties {
NEW_TYPE_TAG(EclBaseProblem, INHERITS_FROM(EclGridManager, EclOutputBlackOil));

// Write all solutions for visualization, not just the ones for the
// report steps...
NEW_PROP_TAG(EnableWriteAllSolutions);

// The number of time steps skipped between writing two consequtive restart files
NEW_PROP_TAG(RestartWritingInterval);

// Set the problem property
SET_TYPE_PROP(EclBaseProblem, Problem, Ewoms::EclProblem<TypeTag>);

// Select the element centered finite volume method as spatial discretization
SET_TAG_PROP(EclBaseProblem, SpatialDiscretizationSplice, EcfvDiscretization);

//! for ebos, use automatic differentiation to linearize the system of PDEs
SET_TAG_PROP(EclBaseProblem, LocalLinearizerSplice, AutoDiffLocalLinearizer);

// Set the material Law
SET_PROP(EclBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Opm::ThreePhaseMaterialTraits<Scalar,
                                          /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                          /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                          /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

public:
    typedef Opm::EclMaterialLawManager<Traits> EclMaterialLawManager;

    typedef typename EclMaterialLawManager::MaterialLaw type;
};

// Enable gravity
SET_BOOL_PROP(EclBaseProblem, EnableGravity, true);

// Reuse the last linearization if possible?
SET_BOOL_PROP(EclBaseProblem, EnableLinearizationRecycling, false);

// Only relinearize the parts where the current solution is sufficiently "bad"
SET_BOOL_PROP(EclBaseProblem, EnablePartialRelinearization, false);

// only write the solutions for the report steps to disk
SET_BOOL_PROP(EclBaseProblem, EnableWriteAllSolutions, false);

// The default for the end time of the simulation [s]
//
// By default, stop it after the universe will probably have stopped
// to exist. (the ECL problem will finish the simulation explicitly
// after it simulated the last episode specified in the deck.)
SET_SCALAR_PROP(EclBaseProblem, EndTime, 1e100);

// The default for the initial time step size of the simulation [s].
//
// The chosen value means that the size of the first time step is the
// one of the initial episode (if the length of the initial episode is
// not millions of trillions of years, that is...)
SET_SCALAR_PROP(EclBaseProblem, InitialTimeStepSize, 1e100);

// increase the default raw tolerance for the newton solver to 10^-4 because this is what
// everone else seems to be doing...
SET_SCALAR_PROP(EclBaseProblem, NewtonRawTolerance, 1e-4);

// Disable the VTK output by default for this problem ...
SET_BOOL_PROP(EclBaseProblem, EnableVtkOutput, false);

// ... but enable the ECL output by default
SET_BOOL_PROP(EclBaseProblem, EnableEclOutput, true);

// also enable the summary output.
SET_BOOL_PROP(EclBaseProblem, EnableEclSummaryOutput, true);

// the cache for intensive quantities can be used for ECL problems and also yields a
// decent speedup...
SET_BOOL_PROP(EclBaseProblem, EnableIntensiveQuantityCache, true);

// the cache for the storage term can also be used and also yields a decent speedup
SET_BOOL_PROP(EclBaseProblem, EnableStorageCache, true);

// Use the "velocity module" which uses the Eclipse "NEWTRAN" transmissibilities
SET_TYPE_PROP(EclBaseProblem, FluxModule, Ewoms::EclTransFluxModule<TypeTag>);

// Use the dummy gradient calculator in order not to do unnecessary work.
SET_TYPE_PROP(EclBaseProblem, GradientCalculator, Ewoms::EclDummyGradientCalculator<TypeTag>);

// The default name of the data file to load
SET_STRING_PROP(EclBaseProblem, GridFile, "data/ecl.DATA");

// The frequency of writing restart (*.ers) files. This is the number of time steps
// between writing restart files
SET_INT_PROP(EclBaseProblem, RestartWritingInterval, 0xffffff); // disable

} // namespace Properties

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This problem simulates an input file given in the data format used by the
 *        commercial ECLiPSE simulator.
 */
template <class TypeTag>
class EclProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // Grid and world dimension
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP(TypeTag, MaterialLaw)::EclMaterialLawManager EclMaterialLawManager;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

    typedef Opm::CompositionalFluidState<Scalar, FluidSystem> ScalarFluidState;
    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Ewoms::EclSummaryWriter<TypeTag> EclSummaryWriter;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    struct RockParams {
        Scalar referencePressure;
        Scalar compressibility;
    };

public:
    /*!
     * \copydoc FvBaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        Ewoms::EclOutputBlackOilModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableWriteAllSolutions,
                             "Write all solutions to disk instead of only the ones for the "
                             "report steps");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableEclOutput,
                             "Write binary output which is compatible with the commercial "
                             "Eclipse simulator");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, RestartWritingInterval,
                             "The frequencies of which time steps are serialized to disk");
    }

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    EclProblem(Simulator &simulator)
        : ParentType(simulator)
        , transmissibilities_(simulator)
        , wellManager_(simulator)
        , deckUnits_(simulator)
        , eclWriter_(simulator)
        , summaryWriter_(simulator)
    {
        // add the output module for the Ecl binary output
        simulator.model().addOutputModule(new Ewoms::EclOutputBlackOilModule<TypeTag>(simulator));
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        auto& simulator = this->simulator();

        // set the value of the gravity constant to the one used by the FLOW simulator
        this->gravity_ = 0.0;

        // the "NOGRAV" keyword from Frontsim disables gravity...
        const auto& deck = simulator.gridManager().deck();
        if (!deck->hasKeyword("NOGRAV") && EWOMS_GET_PARAM(TypeTag, bool, EnableGravity))
            this->gravity_[dim - 1] = 9.80665;

        initFluidSystem_();
        readRockParameters_();
        readMaterialParameters_();
        transmissibilities_.finishInit();
        readInitialCondition_();

        // initialize the wells. Note that this needs to be done after initializing the
        // intrinsic permeabilities because the well model uses them...
        wellManager_.init(simulator.gridManager().eclState());

        // Set the start time of the simulation
        Opm::TimeMapConstPtr timeMap = simulator.gridManager().schedule()->getTimeMap();
        tm curTime = boost::posix_time::to_tm(timeMap->getStartTime(/*timeStepIdx=*/0));

        Scalar startTime = std::mktime(&curTime);
        simulator.setStartTime(startTime);

        // We want the episode index to be the same as the report step index to make
        // things simpler, so we have to set the episode index to -1 because it is
        // incremented inside beginEpisode()...
        simulator.setEpisodeIndex(-1);
    }

    /*!
     * \brief This method restores the complete state of the well
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     *
     * \tparam Restarter The deserializer type
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        // reload the current episode/report step from the deck
        beginEpisode(/*isOnRestart=*/true);

        // deserialize the wells
        wellManager_.deserialize(res);
    }

    /*!
     * \brief This method writes the complete state of the well
     *        to the harddisk.
     */
    template <class Restarter>
    void serialize(Restarter &res)
    { wellManager_.serialize(res); }

    /*!
     * \brief Called by the simulator before an episode begins.
     */
    void beginEpisode(bool isOnRestart = false)
    {
        // Proceed to the next report step
        Simulator &simulator = this->simulator();
        Opm::EclipseStateConstPtr eclState = this->simulator().gridManager().eclState();
        Opm::TimeMapConstPtr timeMap = eclState->getSchedule()->getTimeMap();

        // Opm::TimeMap deals with points in time, so the number of time intervals (i.e.,
        // report steps) is one less!
        int numReportSteps = timeMap->size() - 1;

        // start the next episode if there are additional report steps, else finish the
        // simulation
        int nextEpisodeIdx = simulator.episodeIndex();
        while (nextEpisodeIdx < numReportSteps &&
               simulator.time() >= timeMap->getTimePassedUntil(nextEpisodeIdx + 1)*(1 - 1e-10))
        {
            ++ nextEpisodeIdx;
        }

        Scalar episodeLength = timeMap->getTimeStepLength(nextEpisodeIdx);
        Scalar dt = episodeLength;
        if (nextEpisodeIdx == 0) {
            // allow the size of the initial time step to be set via an external parameter
            Scalar initialDt = EWOMS_GET_PARAM(TypeTag, Scalar, InitialTimeStepSize);
            dt = std::min(dt, initialDt);
        }

        if (nextEpisodeIdx < numReportSteps) {
            simulator.startNextEpisode(episodeLength);
            simulator.setTimeStepSize(dt);
        }

        // set up the wells
        wellManager_.beginEpisode(this->simulator().gridManager().eclState(), isOnRestart);
    }

    /*!
     * \brief Called by the simulator before each time integration.
     */
    void beginTimeStep()
    {
        wellManager_.beginTimeStep();

        // this is a little hack to write the initial condition, which we need to do
        // before the first time step has finished.
        static bool initialWritten = false;
        if (this->simulator().episodeIndex() == 0 && !initialWritten) {
            summaryWriter_.write(wellManager_, /*isInitial=*/true);
            initialWritten = true;
        }
    }

    /*!
     * \brief Called by the simulator before each Newton-Raphson iteration.
     */
    void beginIteration()
    { wellManager_.beginIteration(); }

    /*!
     * \brief Called by the simulator after each Newton-Raphson iteration.
     */
    void endIteration()
    { wellManager_.endIteration(); }

    /*!
     * \brief Called by the simulator after each time integration.
     */
    void endTimeStep()
    {
        wellManager_.endTimeStep();

        // write the summary information after each time step
        summaryWriter_.write(wellManager_);

        updateHysteresis_();

#ifndef NDEBUG
        // in debug mode, we don't care about performance, so we check if the model does
        // the right thing (i.e., the mass change inside the whole reservoir must be
        // equivalent to the fluxes over the grid's boundaries plus the source rates
        // specified by the problem)
        this->model().checkConservativeness(/*tolerance=*/-1, /*verbose=*/true);
#endif // NDEBUG
    }

    /*!
     * \brief Called by the simulator after the end of an episode.
     */
    void endEpisode()
    {
        auto& simulator = this->simulator();
        const auto& eclState = simulator.gridManager().eclState();
        auto& linearizer = this->model().linearizer();
        int episodeIdx = simulator.episodeIndex();

        bool wellsWillChange = wellManager_.wellsChanged(eclState, episodeIdx + 1);
        linearizer.setLinearizationReusable(!wellsWillChange);

        Opm::TimeMapConstPtr timeMap = eclState->getSchedule()->getTimeMap();
        int numReportSteps = timeMap->size() - 1;
        if (episodeIdx + 1 >= numReportSteps) {
            simulator.setFinished(true);
            return;
        }
    }

    /*!
     * \brief Returns true if the current solution should be written
     *        to disk for visualization.
     *
     * For the ECL simulator we only write at the end of
     * episodes/report steps...
     */
    bool shouldWriteOutput() const
    {
        if (this->simulator().timeStepIndex() < 0)
            // always write the initial solution
            return true;

        if (EWOMS_GET_PARAM(TypeTag, bool, EnableWriteAllSolutions))
            return true;

        return this->simulator().episodeWillBeOver();
    }

    /*!
     * \brief Returns true if an eWoms restart file should be written to disk.
     */
    bool shouldWriteRestartFile() const
    {
        unsigned n = EWOMS_GET_PARAM(TypeTag, unsigned, RestartWritingInterval);
        unsigned i = this->simulator().timeStepIndex();
        if (i > 0 && (i%n) == 0)
            return true; // we don't write a restart file for the initial condition
        return false;
    }

    /*!
     * \brief Write the requested quantities of the current solution into the output
     *        files.
     */
    void writeOutput(bool verbose = true)
    {
        // calculate the time _after_ the time was updated
        Scalar t = this->simulator().time() + this->simulator().timeStepSize();

        // prepare the ECL and the VTK writers
        if (enableEclOutput_())
            eclWriter_.beginWrite(t);

        // use the generic code to prepare the output fields and to
        // write the desired VTK files.
        ParentType::writeOutput(verbose);

        if (enableEclOutput_()) {
            this->model().appendOutputFields(eclWriter_);
            eclWriter_.endWrite();
        }
    }

    /*!
     * \brief Returns the object which converts between SI and deck units.
     */
    const EclDeckUnits<TypeTag>& deckUnits() const
    { return deckUnits_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context,
                                           unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return intrinsicPermeability_[globalSpaceIdx];
    }

    /*!
     * \brief This method returns the intrinsic permeability tensor
     *        given a global element index.
     *
     * Its main (only?) usage is the ECL transmissibility calculation code...
     */
    const DimMatrix &intrinsicPermeability(unsigned globalElemIdx) const
    { return intrinsicPermeability_[globalElemIdx]; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::transmissibility
     */
    Scalar transmissibility(unsigned elem1Idx, unsigned elem2Idx) const
    { return transmissibilities_.transmissibility(elem1Idx, elem2Idx); }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return porosity_[globalSpaceIdx];
    }

    /*!
     * \copydoc BlackoilProblem::rockCompressibility
     */
    template <class Context>
    Scalar rockCompressibility(const Context &context, unsigned spaceIdx, unsigned timeIdx) const
    {
        if (rockParams_.empty())
            return 0.0;

        unsigned tableIdx = 0;
        if (!rockTableIdx_.empty()) {
            unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
            tableIdx = rockTableIdx_[globalSpaceIdx];
        }

        return rockParams_[tableIdx].compressibility;
    }

    /*!
     * \copydoc BlackoilProblem::rockReferencePressure
     */
    template <class Context>
    Scalar rockReferencePressure(const Context &context, unsigned spaceIdx, unsigned timeIdx) const
    {
        if (rockParams_.empty())
            return 1e5;

        unsigned tableIdx = 0;
        if (!rockTableIdx_.empty()) {
            unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
            tableIdx = rockTableIdx_[globalSpaceIdx];
        }

        return rockParams_[tableIdx].referencePressure;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return materialLawParams(globalSpaceIdx);
    }

    const MaterialLawParams& materialLawParams(unsigned globalDofIdx) const
    { return materialLawManager_->materialLawParams(globalDofIdx); }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned pvtRegionIndex(const Context &context, unsigned spaceIdx, unsigned timeIdx) const
    { return pvtRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the index the relevant PVT region given a cell index
     */
    unsigned pvtRegionIndex(unsigned elemIdx) const
    {
        Opm::DeckConstPtr deck = this->simulator().gridManager().deck();

        if (!deck->hasKeyword("PVTNUM"))
            return 0;

        const auto& gridManager = this->simulator().gridManager();

        unsigned cartesianDofIdx = gridManager.cartesianIndex(elemIdx);
        return deck->getKeyword("PVTNUM")->getIntData()[cartesianDofIdx] - 1;
    }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return this->simulator().gridManager().caseName(); }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context &context, unsigned spaceIdx, unsigned timeIdx) const
    {
        // use the temporally constant temperature, i.e. use the initial temperature of
        // the DOF
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return initialFluidStates_[globalDofIdx].temperature(/*phaseIdx=*/0);
    }

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * ECLiPSE uses no-flow conditions for all boundaries. \todo really?
     */
    template <class Context>
    void boundary(BoundaryRateVector &values,
                  const Context &context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    { values.setNoFlow(); }

    /*!
     * \copydoc FvBaseProblem::initial
     *
     * The reservoir problem uses a constant boundary condition for
     * the whole domain.
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        values.setPvtRegionIndex(pvtRegionIndex(context, spaceIdx, timeIdx));

        if (useMassConservativeInitialCondition_) {
            const auto& matParams = materialLawParams(context, spaceIdx, timeIdx);
            values.assignMassConservative(initialFluidStates_[globalDofIdx], matParams);
        }
        else
            values.assignNaive(initialFluidStates_[globalDofIdx]);
    }

    void initialSolutionApplied()
    {
        updateHysteresis_();
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0 everywhere.
     */
    template <class Context>
    void source(RateVector &rate,
                const Context &context,
                unsigned spaceIdx,
                unsigned timeIdx) const
    {
        rate = Toolbox::createConstant(0);

        for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx)
            rate[eqIdx] = Toolbox::createConstant(0.0);

        wellManager_.computeTotalRatesForDof(rate, context, spaceIdx, timeIdx);

        // convert the source term from the total mass rate of the
        // cell to the one per unit of volume as used by the model.
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx)
            rate[eqIdx] /= this->model().dofTotalVolume(globalDofIdx);
    }

    /*!
     * \brief Returns a reference to the ECL well manager used by the problem.
     *
     * This can be used for inspecting wells outside of the problem.
     */
    const EclWellManager<TypeTag>& wellManager() const
    { return wellManager_; }

private:
    static bool enableEclOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput); }

    void readRockParameters_()
    {
        auto deck = this->simulator().gridManager().deck();
        auto eclState = this->simulator().gridManager().eclState();
        const auto& gridManager = this->simulator().gridManager();

        // the ROCK keyword has not been specified, so we don't need
        // to read rock parameters
        if (!deck->hasKeyword("ROCK"))
            return;

        const auto rockKeyword = deck->getKeyword("ROCK");
        rockParams_.resize(rockKeyword->size());
        for (size_t rockRecordIdx = 0; rockRecordIdx < rockKeyword->size(); ++ rockRecordIdx) {
            const auto rockRecord = rockKeyword->getRecord(rockRecordIdx);
            rockParams_[rockRecordIdx].referencePressure =
                rockRecord->getItem("PREF")->getSIDouble(0);
            rockParams_[rockRecordIdx].compressibility =
                rockRecord->getItem("COMPRESSIBILITY")->getSIDouble(0);
        }

        // PVTNUM has not been specified, so everything is in the first region and we
        // don't need to care...
        if (!eclState->hasIntGridProperty("PVTNUM"))
            return;

        const std::vector<int>& pvtnumData =
            eclState->getIntGridProperty("PVTNUM")->getData();
        rockTableIdx_.resize(gridManager.gridView().size(/*codim=*/0));
        for (size_t elemIdx = 0; elemIdx < rockTableIdx_.size(); ++ elemIdx) {
            unsigned cartElemIdx = gridManager.cartesianIndex(elemIdx);

            // reminder: Eclipse uses FORTRAN-style indices
            rockTableIdx_[elemIdx] = pvtnumData[cartElemIdx] - 1;
        }
    }

    void readMaterialParameters_()
    {
        const auto &gridManager = this->simulator().gridManager();
        auto deck = gridManager.deck();
        auto eclState = gridManager.eclState();

        size_t numDof = this->model().numGridDof();

        intrinsicPermeability_.resize(numDof);
        porosity_.resize(numDof);

        ////////////////////////////////
        // permeability

        // read the intrinsic permeabilities from the eclState. Note that all arrays
        // provided by eclState are one-per-cell of "uncompressed" grid, whereas the
        // dune-cornerpoint grid object might remove a few elements...
        if (eclState->hasDoubleGridProperty("PERMX")) {
            const std::vector<double> &permxData =
                eclState->getDoubleGridProperty("PERMX")->getData();
            std::vector<double> permyData(permxData);
            if (eclState->hasDoubleGridProperty("PERMY"))
                permyData = eclState->getDoubleGridProperty("PERMY")->getData();
            std::vector<double> permzData(permxData);
            if (eclState->hasDoubleGridProperty("PERMZ"))
                permzData = eclState->getDoubleGridProperty("PERMZ")->getData();

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
                unsigned cartesianElemIdx = gridManager.cartesianIndex(dofIdx);
                intrinsicPermeability_[dofIdx] = 0.0;
                intrinsicPermeability_[dofIdx][0][0] = permxData[cartesianElemIdx];
                intrinsicPermeability_[dofIdx][1][1] = permyData[cartesianElemIdx];
                intrinsicPermeability_[dofIdx][2][2] = permzData[cartesianElemIdx];
            }

            // for now we don't care about non-diagonal entries
        }
        else
            OPM_THROW(std::logic_error,
                      "Can't read the intrinsic permeability from the ecl state. "
                      "(The PERM{X,Y,Z} keywords are missing)");
        ////////////////////////////////


        ////////////////////////////////
        // compute the porosity
        if (!eclState->hasDoubleGridProperty("PORO") && !eclState->hasDoubleGridProperty("PORV"))
            OPM_THROW(std::runtime_error,
                      "Can't read the porosity from the ECL state object. "
                      "(The PORO and PORV keywords are missing)");

        if (eclState->hasDoubleGridProperty("PORO")) {
            const std::vector<double> &poroData =
                eclState->getDoubleGridProperty("PORO")->getData();

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
                unsigned cartesianElemIdx = gridManager.cartesianIndex(dofIdx);
                porosity_[dofIdx] = poroData[cartesianElemIdx];
            }
        }

        // overwrite the porosity using the PORV keyword for the elements for which PORV
        // is defined...
        if (eclState->hasDoubleGridProperty("PORV")) {
            const std::vector<double> &porvData =
                eclState->getDoubleGridProperty("PORV")->getData();

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
                unsigned cartesianElemIdx = gridManager.cartesianIndex(dofIdx);
                if (std::isfinite(porvData[cartesianElemIdx])) {
                    Scalar dofVolume = this->simulator().model().dofTotalVolume(dofIdx);
                    porosity_[dofIdx] = porvData[cartesianElemIdx]/dofVolume;
                }
            }
        }

        // apply the NTG keyword to the porosity
        if (eclState->hasDoubleGridProperty("NTG")) {
            const std::vector<double> &ntgData =
                eclState->getDoubleGridProperty("NTG")->getData();

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx)
                porosity_[dofIdx] *= ntgData[gridManager.cartesianIndex(dofIdx)];
        }

        // apply the MULTPV keyword to the porosity
        if (eclState->hasDoubleGridProperty("MULTPV")) {
            const std::vector<double> &multpvData =
                eclState->getDoubleGridProperty("MULTPV")->getData();

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx)
                porosity_[dofIdx] *= multpvData[gridManager.cartesianIndex(dofIdx)];
        }

        // the fluid-matrix interactions for ECL problems are dealt with by a separate class
        std::vector<int> compressedToCartesianElemIdx(numDof);
        for (unsigned elemIdx = 0; elemIdx < numDof; ++elemIdx)
            compressedToCartesianElemIdx[elemIdx] = gridManager.cartesianIndex(elemIdx);

        materialLawManager_ = std::make_shared<EclMaterialLawManager>();
        materialLawManager_->initFromDeck(deck, eclState, compressedToCartesianElemIdx);
    }

    void initFluidSystem_()
    {
        const auto deck = this->simulator().gridManager().deck();
        const auto eclState = this->simulator().gridManager().eclState();

        FluidSystem::initFromDeck(deck, eclState);
   }

    void readInitialCondition_()
    {
        const auto &gridManager = this->simulator().gridManager();
        const auto deck = gridManager.deck();

        if (!deck->hasKeyword("EQUIL"))
            readExplicitInitialCondition_();
        else
            readEquilInitialCondition_();
    }

    void readEquilInitialCondition_()
    {
        // The EQUIL initializer also modifies the material law manager according to
        // SWATINIT (although it does not belong there strictly speaking)
        typedef Ewoms::EclEquilInitializer<TypeTag> EquilInitializer;
        EquilInitializer equilInitializer(this->simulator(), materialLawManager_);

        // since the EquilInitializer provides fluid states that are consistent with the
        // black-oil model, we can use naive instead of mass conservative determination
        // of the primary variables.
        useMassConservativeInitialCondition_ = false;

        size_t numElems = this->model().numGridDof();
        initialFluidStates_.resize(numElems);
        for (size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto &elemFluidState = initialFluidStates_[elemIdx];
            elemFluidState.assign(equilInitializer.initialFluidState(elemIdx));
        }

        // release the equil grid pointer since it's no longer needed.
        this->simulator().gridManager().releaseEquilGrid();
    }

    void readExplicitInitialCondition_()
    {
        const auto &gridManager = this->simulator().gridManager();
        const auto deck = gridManager.deck();
        const auto eclState = gridManager.eclState();

        // since the values specified in the deck do not need to be consistent, we use an
        // initial condition that conserves the total mass specified by these values.
        useMassConservativeInitialCondition_ = true;

        bool enableDisgas = deck->hasKeyword("DISGAS");
        bool enableVapoil = deck->hasKeyword("VAPOIL");

        // make sure all required quantities are enables
         if (!deck->hasKeyword("SWAT") ||
             !deck->hasKeyword("SGAS"))
             OPM_THROW(std::runtime_error,
                      "The ECL input file requires the presence of the SWAT "
                      "and SGAS keywords if the model is initialized explicitly");
         if (!deck->hasKeyword("PRESSURE"))
             OPM_THROW(std::runtime_error,
                      "The ECL input file requires the presence of the PRESSURE "
                      "keyword if the model is initialized explicitly");
         if (enableDisgas && !deck->hasKeyword("RS"))
             OPM_THROW(std::runtime_error,
                      "The ECL input file requires the RS keyword to be present if"
                      " dissolved gas is enabled");
         if (enableVapoil && !deck->hasKeyword("RV"))
             OPM_THROW(std::runtime_error,
                      "The ECL input file requires the RV keyword to be present if"
                      " vaporized oil is enabled");

        size_t numDof = this->model().numGridDof();

        initialFluidStates_.resize(numDof);

        const std::vector<double> &waterSaturationData =
            deck->getKeyword("SWAT")->getSIDoubleData();
        const std::vector<double> &gasSaturationData =
            deck->getKeyword("SGAS")->getSIDoubleData();
        const std::vector<double> &pressureData =
            deck->getKeyword("PRESSURE")->getSIDoubleData();
        const std::vector<double> *rsData = 0;
        if (enableDisgas)
            rsData = &deck->getKeyword("RS")->getSIDoubleData();
        const std::vector<double> *rvData = 0;
        if (enableVapoil)
            rvData = &deck->getKeyword("RV")->getSIDoubleData();
        // initial reservoir temperature
        const std::vector<double> &tempiData =
            eclState->getDoubleGridProperty("TEMPI")->getData();

        // make sure that the size of the data arrays is correct
#ifndef NDEBUG
        const auto &cartSize = this->simulator().gridManager().cartesianDimensions();
        size_t numCartesianCells = cartSize[0] * cartSize[1] * cartSize[2];
        assert(waterSaturationData.size() == numCartesianCells);
        assert(gasSaturationData.size() == numCartesianCells);
        assert(pressureData.size() == numCartesianCells);
        if (enableDisgas)
            assert(rsData->size() == numCartesianCells);
        if (enableVapoil)
            assert(rvData->size() == numCartesianCells);
#endif

        // calculate the initial fluid states
        for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto &dofFluidState = initialFluidStates_[dofIdx];

            size_t cartesianDofIdx = gridManager.cartesianIndex(dofIdx);
            assert(0 <= cartesianDofIdx);
            assert(cartesianDofIdx <= numCartesianCells);

            //////
            // set temperature
            //////
            Scalar temperature = tempiData[cartesianDofIdx];
            if (!std::isfinite(temperature) || temperature <= 0)
                temperature = FluidSystem::surfaceTemperature;
            dofFluidState.setTemperature(temperature);

            //////
            // set saturations
            //////
            dofFluidState.setSaturation(FluidSystem::waterPhaseIdx,
                                        waterSaturationData[cartesianDofIdx]);
            dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                        gasSaturationData[cartesianDofIdx]);
            dofFluidState.setSaturation(FluidSystem::oilPhaseIdx,
                                        1
                                        - waterSaturationData[cartesianDofIdx]
                                        - gasSaturationData[cartesianDofIdx]);

            //////
            // set phase pressures
            //////
            Scalar oilPressure = pressureData[cartesianDofIdx];

            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            Dune::FieldVector< Scalar, numPhases > pc( 0 );
            const auto& matParams = materialLawParams(dofIdx);
            MaterialLaw::capillaryPressures(pc, matParams, dofFluidState);
            Valgrind::CheckDefined(oilPressure);
            Valgrind::CheckDefined(pc);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                dofFluidState.setPressure(phaseIdx, oilPressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
            Scalar gasPressure = dofFluidState.pressure(gasPhaseIdx);

            //////
            // set compositions
            //////

            // reset all mole fractions to 0
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                    dofFluidState.setMoleFraction(phaseIdx, compIdx, 0.0);

            // by default, assume immiscibility for all phases
            dofFluidState.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);
            dofFluidState.setMoleFraction(gasPhaseIdx, gasCompIdx, 1.0);
            dofFluidState.setMoleFraction(oilPhaseIdx, oilCompIdx, 1.0);

            if (enableDisgas) {
                // set the composition of the oil phase:
                //
                // first, retrieve the relevant black-oil parameters from
                // the fluid system.
                //
                // note that we use the gas pressure here. this is because the primary
                // varibles and the intensive quantities of the black oil model also do
                // this...
                Scalar RsSat = FluidSystem::gasDissolutionFactor(temperature,
                                                                 gasPressure,
                                                                 /*regionIdx=*/0);
                Scalar RsReal = (*rsData)[cartesianDofIdx];

                if (RsReal > RsSat) {
                    std::array<int, 3> ijk;
                    gridManager.cartesianCoordinate(dofIdx, ijk);
                    std::cerr << "Warning: The specified amount gas (R_s = " << RsReal << ") is more"
                              << " than the maximium\n"
                              << "         amount which can be dissolved in oil"
                              << " (R_s,max=" << RsSat << ")"
                              << " for cell (" << ijk[0] << ", " << ijk[1] << ", " << ijk[2] << ")."
                              << " Ignoring.\n";
                    RsReal = RsSat;
                }

                // calculate composition of the real and the saturated oil phase in terms of
                // mass fractions.
                Scalar rhooRef = FluidSystem::referenceDensity(oilPhaseIdx, /*regionIdx=*/0);
                Scalar rhogRef = FluidSystem::referenceDensity(gasPhaseIdx, /*regionIdx=*/0);
                Scalar XoGReal = RsReal/(RsReal + rhooRef/rhogRef);

                // convert mass to mole fractions
                Scalar MG = FluidSystem::molarMass(gasCompIdx);
                Scalar MO = FluidSystem::molarMass(oilCompIdx);

                Scalar xoGReal = XoGReal * MO / ((MO - MG) * XoGReal + MG);
                Scalar xoOReal = 1 - xoGReal;

                // finally, set the oil-phase composition
                dofFluidState.setMoleFraction(oilPhaseIdx, gasCompIdx, xoGReal);
                dofFluidState.setMoleFraction(oilPhaseIdx, oilCompIdx, xoOReal);
            }

            if (enableVapoil) {
                // set the composition of the gas phase:
                //
                // first, retrieve the relevant black-gas parameters from
                // the fluid system.
                Scalar RvSat = FluidSystem::oilVaporizationFactor(temperature,
                                                                  gasPressure,
                                                                  /*regionIdx=*/0);
                Scalar RvReal = (*rvData)[cartesianDofIdx];

                if (RvReal > RvSat) {
                    std::array<int, 3> ijk;
                    gridManager.cartesianCoordinate(dofIdx, ijk);
                    std::cerr << "Warning: The specified amount oil (R_v = " << RvReal << ") is more"
                              << " than the maximium\n"
                              << "         amount which can be dissolved in gas"
                              << " (R_v,max=" << RvSat << ")"
                              << " for cell (" << ijk[0] << ", " << ijk[1] << ", " << ijk[2] << ")."
                              << " Ignoring.\n";
                    RvReal = RvSat;
                }

                // calculate composition of the real and the saturated gas phase in terms of
                // mass fractions.
                Scalar rhooRef = FluidSystem::referenceDensity(oilPhaseIdx, /*regionIdx=*/0);
                Scalar rhogRef = FluidSystem::referenceDensity(gasPhaseIdx, /*regionIdx=*/0);
                Scalar XgOReal = RvReal/(RvReal + rhogRef/rhooRef);

                // convert mass to mole fractions
                Scalar MG = FluidSystem::molarMass(gasCompIdx);
                Scalar MO = FluidSystem::molarMass(oilCompIdx);

                Scalar xgOReal = XgOReal * MG / ((MG - MO) * XgOReal + MO);
                Scalar xgGReal = 1 - xgOReal;

                // finally, set the gas-phase composition
                dofFluidState.setMoleFraction(gasPhaseIdx, oilCompIdx, xgOReal);
                dofFluidState.setMoleFraction(gasPhaseIdx, gasCompIdx, xgGReal);
            }
        }
    }

    // update the hysteresis parameters of the material laws for the whole grid
    void updateHysteresis_()
    {
        if (!materialLawManager_->enableHysteresis())
            return;

        ElementContext elemCtx(this->simulator());
        const auto& gridManager = this->simulator().gridManager();
        auto elemIt = gridManager.gridView().template begin</*codim=*/0>();
        const auto &elemEndIt = gridManager.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updateStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            materialLawManager_->updateHysteresis(intQuants.fluidState(), compressedDofIdx);
        }
    }

    std::vector<Scalar> porosity_;
    std::vector<DimMatrix> intrinsicPermeability_;
    EclTransmissibility<TypeTag> transmissibilities_;

    std::shared_ptr<EclMaterialLawManager> materialLawManager_;

    std::vector<unsigned short> rockTableIdx_;
    std::vector<RockParams> rockParams_;

    bool useMassConservativeInitialCondition_;
    std::vector<ScalarFluidState> initialFluidStates_;

    EclWellManager<TypeTag> wellManager_;

    EclDeckUnits<TypeTag> deckUnits_;

    EclWriter<TypeTag> eclWriter_;
    EclSummaryWriter summaryWriter_;
};
} // namespace Ewoms

#endif
