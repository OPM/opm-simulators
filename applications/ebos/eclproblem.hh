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

#include "eclgridmanager.hh"
#include "eclwellmanager.hh"
#include "eclwriter.hh"
#include "eclsummarywriter.hh"
#include "ecloutputblackoilmodule.hh"
#include "ecltransmissibility.hh"
#include "ecldummygradientcalculator.hh"
#include "eclfluxmodule.hh"

#include <ewoms/models/blackoil/blackoilmodel.hh>
#include <ewoms/disc/ecfv/ecfvdiscretization.hh>

#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/SplineTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EclDefaultMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/core/utility/Average.hpp>

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
}

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(EclBaseProblem, INHERITS_FROM(EclGridManager, EclOutputBlackOil));

// The temperature inside the reservoir
NEW_PROP_TAG(Temperature);

// Write all solutions for visualization, not just the ones for the
// report steps...
NEW_PROP_TAG(EnableWriteAllSolutions);

// Set the problem property
SET_TYPE_PROP(EclBaseProblem, Problem, Ewoms::EclProblem<TypeTag>);

// Select the element centered finite volume method as spatial discretization
SET_TAG_PROP(EclBaseProblem, SpatialDiscretizationSplice, EcfvDiscretization);

// Set the material Law
SET_PROP(EclBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx> OilWaterTraits;

    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx> GasOilTraits;

    typedef Opm::ThreePhaseMaterialTraits<Scalar,
                                          /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                          /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                          /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

    typedef typename Opm::PiecewiseLinearTwoPhaseMaterial<OilWaterTraits> OilWaterLaw;
    typedef typename Opm::PiecewiseLinearTwoPhaseMaterial<GasOilTraits> GasOilLaw;

//    typedef typename Opm::SplineTwoPhaseMaterial<OilWaterTraits> OilWaterLaw;
//    typedef typename Opm::SplineTwoPhaseMaterial<GasOilTraits> GasOilLaw;

public:
    typedef Opm::EclDefaultMaterial<Traits, GasOilLaw, OilWaterLaw> type;
};

// Enable gravity
SET_BOOL_PROP(EclBaseProblem, EnableGravity, true);

// Reuse the last linearization if possible?
SET_BOOL_PROP(EclBaseProblem, EnableLinearizationRecycling, true);

// Only relinearize the parts where the current solution is sufficiently "bad"
SET_BOOL_PROP(EclBaseProblem, EnablePartialRelinearization, true);

// only write the solutions for the report steps to disk
SET_BOOL_PROP(EclBaseProblem, EnableWriteAllSolutions, false);

// set the defaults for some problem specific properties
SET_SCALAR_PROP(EclBaseProblem, Temperature, 293.15);

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

// Use the "velocity module" which uses the Eclipse "NEWTRAN" transmissibilities
SET_TYPE_PROP(EclBaseProblem, FluxModule, Ewoms::EclTransFluxModule<TypeTag>);

// Use the dummy gradient calculator in order not to do unnecessary work.
SET_TYPE_PROP(EclBaseProblem, GradientCalculator, Ewoms::EclDummyGradientCalculator<TypeTag>);

// The default name of the data file to load
SET_STRING_PROP(EclBaseProblem, GridFile, "data/ecl.DATA");
}} // namespace Properties, Opm

namespace Ewoms {
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

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // Grid and world dimension
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
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
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, BlackOilFluidState) BlackOilFluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

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

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableWriteAllSolutions,
                             "Write all solutions to disk instead of only the ones for the "
                             "report steps");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableEclOutput,
                             "Write binary output which is compatible with the commercial "
                             "Eclipse simulator");
    }

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    EclProblem(Simulator &simulator)
        : ParentType(simulator)
        , transmissibilities_(simulator)
        , wellManager_(simulator)
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

        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);

        // invert the direction of the gravity vector for ECL problems
        // (z coodinates represent depth, not height.)
        this->gravity_[dim - 1] *= -1;

        // the "NOGRAV" keyword from Frontsim disables gravity...
        const auto& deck = simulator.gridManager().deck();
        if (deck->hasKeyword("NOGRAV"))
            this->gravity_ = 0.0;

        initFluidSystem_();
        readRockParameters_();
        readMaterialParameters_();
        transmissibilities_.finishInit();
        readInitialCondition_();

        // initialize the wells. Note that this needs to be done after initializing the
        // intrinsic permeabilities because the well model uses them...
        wellManager_.init(simulator.gridManager().eclState());

        // Start the first episode. For this, ask the ECL schedule.
        Opm::TimeMapConstPtr timeMap = simulator.gridManager().schedule()->getTimeMap();
        tm curTime = boost::posix_time::to_tm(timeMap->getStartTime(/*timeStepIdx=*/0));

        Scalar startTime = std::mktime(&curTime);
        simulator.setStartTime(startTime);
        simulator.startNextEpisode(/*startTime=*/startTime,
                                   /*length=*/timeMap->getTimeStepLength(/*timeStepIdx=*/0));

        // we want the episode index to be the same as the report step
        // index to make things simpler...
        simulator.setEpisodeIndex(0);
    }

    /*!
     * \brief Called by the simulator before an episode begins.
     */
    void beginEpisode()
    { wellManager_.beginEpisode(this->simulator().gridManager().eclState()); }

    /*!
     * \brief Called by the simulator before each time integration.
     */
    void beginTimeStep()
    { wellManager_.beginTimeStep(); }
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

#ifndef NDEBUG
        this->model().checkConservativeness(/*tolerance=*/-1, /*verbose=*/true);
#endif // NDEBUG
    }

    /*!
     * \brief Called by the simulator after the end of an episode.
     */
    void endEpisode()
    {
        // first, write the summary information ...
        summaryWriter_.write(wellManager_);

        // ... then proceed to the next report step
        Simulator &simulator = this->simulator();
        Opm::EclipseStateConstPtr eclState = this->simulator().gridManager().eclState();
        Opm::TimeMapConstPtr timeMap = eclState->getSchedule()->getTimeMap();

        // TimeMap deals with points in time, so the number of time
        // intervals (i.e., report steps) is one less!
        int numReportSteps = timeMap->size() - 1;

        // start the next episode if there are additional report
        // steps, else finish the simulation
        int nextEpisodeIdx = simulator.episodeIndex() + 1;
        if (nextEpisodeIdx < numReportSteps) {
            simulator.startNextEpisode(timeMap->getTimeStepLength(nextEpisodeIdx));
            simulator.setTimeStepSize(timeMap->getTimeStepLength(nextEpisodeIdx));
        }
        else
            simulator.setFinished(true);
    }

    /*!
     * \brief Returns true if the current solution should be written
     *        to disk for visualization.
     *
     * For the ECL simulator we only write at the end of
     * episodes/report steps...
     */
    bool shouldWriteOutput()
    {
        if (this->simulator().timeStepIndex() < 0)
            // always write the initial solution
            return true;

        if (EWOMS_GET_PARAM(TypeTag, bool, EnableWriteAllSolutions))
            return true;

        return this->simulator().episodeWillBeOver();
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
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context,
                                           int spaceIdx,
                                           int timeIdx) const
    {
        int globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return intrinsicPermeability_[globalSpaceIdx];
    }

    /*!
     * \brief This method returns the intrinsic permeability tensor
     *        given a global element index.
     *
     * Its main (only?) usage is the ECL transmissibility calculation code...
     */
    const DimMatrix &intrinsicPermeability(int globalElemIdx) const
    { return intrinsicPermeability_[globalElemIdx]; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::transmissibility
     */
    Scalar transmissibility(int elem1Idx, int elem2Idx) const
    { return transmissibilities_.transmissibility(elem1Idx, elem2Idx); }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        int globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return porosity_[globalSpaceIdx];
    }

    /*!
     * \copydoc BlackoilProblem::rockCompressibility
     */
    template <class Context>
    Scalar rockCompressibility(const Context &context, int spaceIdx, int timeIdx) const
    {
        if (rockParams_.empty())
            return 0.0;

        int tableIdx = 0;
        if (!rockTableIdx_.empty()) {
            int globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
            tableIdx = rockTableIdx_[globalSpaceIdx];
        }

        return rockParams_[tableIdx].compressibility;
    }

    /*!
     * \copydoc BlackoilProblem::rockReferencePressure
     */
    template <class Context>
    Scalar rockReferencePressure(const Context &context, int spaceIdx, int timeIdx) const
    {
        if (rockParams_.empty())
            return 1e5;

        int tableIdx = 0;
        if (!rockTableIdx_.empty()) {
            int globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
            tableIdx = rockTableIdx_[globalSpaceIdx];
        }

        return rockParams_[tableIdx].referencePressure;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context,
                                               int spaceIdx, int timeIdx) const
    {
        int tableIdx = 0;
        if (materialParamTableIdx_.size() > 0) {
            int globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
            tableIdx = materialParamTableIdx_[globalSpaceIdx];
        }
        return materialParams_[tableIdx];
    }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    int pvtRegionIndex(const Context &context, int spaceIdx, int timeIdx) const
    {
        Opm::DeckConstPtr deck = this->simulator().gridManager().deck();

        if (!deck->hasKeyword("PVTNUM"))
            return 0;

        const auto &cartesianCellId = this->simulator().gridManager().cartesianCellId();

        // this is quite specific to the ECFV discretization. But so is everything in an
        // ECL deck, i.e., we don't need to care here...
        int compressedDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        int cartesianDofIdx = cartesianCellId[compressedDofIdx];

        return deck->getKeyword("PVTNUM")->getIntData()[cartesianDofIdx] - 1;
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return this->simulator().gridManager().caseName(); }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     *
     * The black-oil model assumes constant temperature to define its
     * parameters. Although temperature is thus not really used by the
     * model, it gets written to the VTK output. Who nows, maybe we
     * will need it one day?
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; }

    // \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * ECLiPSE uses no-flow conditions for all boundaries. \todo really?
     */
    template <class Context>
    void boundary(BoundaryRateVector &values,
                  const Context &context,
                  int spaceIdx,
                  int timeIdx) const
    { values.setNoFlow(); }

    //! \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     *
     * The reservoir problem uses a constant boundary condition for
     * the whole domain.
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        int globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        values.assignNaive(initialFluidStates_[globalDofIdx]);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0 everywhere.
     */
    template <class Context>
    void source(RateVector &rate,
                const Context &context,
                int spaceIdx,
                int timeIdx) const
    {
        rate = 0;
        wellManager_.computeTotalRatesForDof(rate, context, spaceIdx, timeIdx);

        // convert the source term from the total mass rate of the
        // cell to the one per unit of volume as used by the model.
        int globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        rate /= this->model().dofTotalVolume(globalDofIdx);
    }

    //! \}

private:
    static bool enableEclOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput); }

    void readRockParameters_()
    {
        auto deck = this->simulator().gridManager().deck();
        auto eclState = this->simulator().gridManager().eclState();

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

        // ROCKTAB has not been specified, so everything is in the
        // first region and we don't need to care...
        if (!eclState->hasIntGridProperty("ROCKTAB"))
            return;

        const std::vector<int>& rocktabData =
            eclState->getIntGridProperty("ROCKTAB")->getData();
        for (size_t elemIdx = 0; elemIdx < rocktabData.size(); ++ elemIdx)
            // reminder: Eclipse uses FORTRAN indices
            rockTableIdx_[elemIdx] = rocktabData[elemIdx] - 1;
    }

    void readMaterialParameters_()
    {
        auto deck = this->simulator().gridManager().deck();
        auto eclState = this->simulator().gridManager().eclState();
        const auto &cartesianCellId = this->simulator().gridManager().cartesianCellId();

        size_t numDof = this->model().numGridDof();

        intrinsicPermeability_.resize(numDof);
        porosity_.resize(numDof);
        materialParams_.resize(numDof);

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
                int cartesianElemIdx = cartesianCellId[dofIdx];
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
        if (eclState->hasDoubleGridProperty("PORO")) {
            const std::vector<double> &poroData =
                eclState->getDoubleGridProperty("PORO")->getData();

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
                int cartesianElemIdx = cartesianCellId[dofIdx];
                porosity_[dofIdx] = poroData[cartesianElemIdx];
            }
        }
        else
            OPM_THROW(std::runtime_error,
                      "Can't read the porosity from the ECL state object. "
                      "(The PORO keyword is missing)");

        // apply the NTG keyword to the porosity
        if (eclState->hasDoubleGridProperty("NTG")) {
            const std::vector<double> &ntgData =
                eclState->getDoubleGridProperty("NTG")->getData();

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
                int cartesianElemIdx = cartesianCellId[dofIdx];
                porosity_[dofIdx] *= ntgData[cartesianElemIdx];
            }
        }

        // apply the MULTPV keyword to the porosity
        if (eclState->hasDoubleGridProperty("MULTPV")) {
            const std::vector<double> &multpvData =
                eclState->getDoubleGridProperty("MULTPV")->getData();

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
                int cartesianElemIdx = cartesianCellId[dofIdx];
                porosity_[dofIdx] *= multpvData[cartesianElemIdx];
            }
        }

        ////////////////////////////////
        // fluid parameters
        const auto& swofTables = eclState->getSwofTables();
        const auto& sgofTables = eclState->getSgofTables();

        // the number of tables for the SWOF and the SGOF keywords
        // must be identical
        assert(swofTables.size() == sgofTables.size());

        size_t numSatfuncTables = swofTables.size();
        materialParams_.resize(numSatfuncTables);

        typedef typename MaterialLawParams::GasOilParams GasOilParams;
        typedef typename MaterialLawParams::OilWaterParams OilWaterParams;

        for (size_t tableIdx = 0; tableIdx < numSatfuncTables; ++ tableIdx) {
            // set the parameters of the material law for a given table
            OilWaterParams owParams;
            GasOilParams goParams;

            const auto& swofTable = swofTables[tableIdx];
            const auto& sgofTable = sgofTables[tableIdx];

            const auto &SwColumn = swofTable.getSwColumn();

            owParams.setKrwSamples(SwColumn, swofTable.getKrwColumn());
            owParams.setKrnSamples(SwColumn, swofTable.getKrowColumn());
            owParams.setPcnwSamples(SwColumn, swofTable.getPcowColumn());

            // convert the saturations of the SGOF keyword from gas to oil saturations
            std::vector<double> SoSamples(sgofTable.numRows());
            for (size_t sampleIdx = 0; sampleIdx < sgofTable.numRows(); ++ sampleIdx)
                SoSamples[sampleIdx] = 1 - sgofTable.getSgColumn()[sampleIdx];

            goParams.setKrwSamples(SoSamples, sgofTable.getKrogColumn());
            goParams.setKrnSamples(SoSamples, sgofTable.getKrgColumn());
            goParams.setPcnwSamples(SoSamples, sgofTable.getPcogColumn());

            owParams.finalize();
            goParams.finalize();

            // compute the connate water saturation. In ECL decks that is defined as
            // the first saturation value of the SWOF keyword.
            Scalar Swco = SwColumn.front();
            materialParams_[tableIdx].setConnateWaterSaturation(Swco);

            materialParams_[tableIdx].setOilWaterParams(owParams);
            materialParams_[tableIdx].setGasOilParams(goParams);

            materialParams_[tableIdx].finalize();
        }

        // set the index of the table to be used
        if (eclState->hasIntGridProperty("SATNUM")) {
            const std::vector<int> &satnumData =
                eclState->getIntGridProperty("SATNUM")->getData();

            materialParamTableIdx_.resize(numDof);
            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
                int cartesianElemIdx = cartesianCellId[dofIdx];

                // make sure that all values are in the correct range
                assert(1 <= satnumData[dofIdx]);
                assert(satnumData[dofIdx] <= static_cast<int>(numSatfuncTables));

                // ECL uses Fortran-style indices which start at
                // 1, but this here is C++...
                materialParamTableIdx_[dofIdx] = satnumData[cartesianElemIdx] - 1;
            }
        }
        else
            materialParamTableIdx_.clear();
        ////////////////////////////////
    }

    void initFluidSystem_()
    {
        const auto deck = this->simulator().gridManager().deck();
        const auto eclState = this->simulator().gridManager().eclState();

        FluidSystem::initBegin();

        int numRegions = deck->getKeyword("DENSITY")->size();
        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            // set the reference densities
            Opm::DeckRecordConstPtr densityRecord =
                deck->getKeyword("DENSITY")->getRecord(regionIdx);
            FluidSystem::setReferenceDensities(densityRecord->getItem("OIL")->getSIDouble(0),
                                             densityRecord->getItem("WATER")->getSIDouble(0),
                                             densityRecord->getItem("GAS")->getSIDouble(0),
                                             regionIdx);

            // so far, we require the presence of the PVTO, PVTW and PVDG
            // keywords...
            FluidSystem::setPvtoTable(eclState->getPvtoTables()[regionIdx], regionIdx);
            FluidSystem::setPvtw(deck->getKeyword("PVTW"), regionIdx);
            FluidSystem::setPvdgTable(eclState->getPvdgTables()[regionIdx], regionIdx);
        }

        FluidSystem::initEnd();
   }

    void readInitialCondition_()
    {
        const auto deck = this->simulator().gridManager().deck();
        const auto &gridManager = this->simulator().gridManager();
        const auto &cartesianCellId = gridManager.cartesianCellId();

        if (!deck->hasKeyword("SWAT") ||
            !deck->hasKeyword("SGAS"))
            OPM_THROW(std::runtime_error,
                      "So far, the ECL input file requires the presence of the SWAT "
                      "and SGAS keywords");
        if (!deck->hasKeyword("PRESSURE"))
            OPM_THROW(std::runtime_error,
                      "So far, the ECL input file requires the presence of the PRESSURE "
                      "keyword");
        if (!deck->hasKeyword("DISGAS"))
            OPM_THROW(std::runtime_error,
                      "The deck must exhibit gas dissolved in the oil phase"
                      " (DISGAS keyword is missing)");
        if (!deck->hasKeyword("RS"))
            OPM_THROW(std::runtime_error,
                      "The ECL input file requires the presence of the RS keyword");

        if (deck->hasKeyword("VAPOIL"))
            OPM_THROW(std::runtime_error,
                      "The deck must _not_ exhibit vaporized oil"
                      " (The VAPOIL keyword is unsupported)");
        if (deck->hasKeyword("RV"))
            OPM_THROW(std::runtime_error,
                      "The ECL input file requires the RV keyword to be non-present");

        size_t numDof = this->model().numGridDof();

        initialFluidStates_.resize(numDof);

        const std::vector<double> &waterSaturationData =
            deck->getKeyword("SWAT")->getSIDoubleData();
        const std::vector<double> &gasSaturationData =
            deck->getKeyword("SGAS")->getSIDoubleData();
        const std::vector<double> &pressureData =
            deck->getKeyword("PRESSURE")->getSIDoubleData();
        const std::vector<double> &rsData =
            deck->getKeyword("RS")->getSIDoubleData();

        // make sure that the size of the data arrays is correct
#ifndef NDEBUG
        const auto &cartSize = this->simulator().gridManager().logicalCartesianSize();
        size_t numCartesianCells = cartSize[0] * cartSize[1] * cartSize[2];
        assert(waterSaturationData.size() == numCartesianCells);
        assert(gasSaturationData.size() == numCartesianCells);
        assert(pressureData.size() == numCartesianCells);
        assert(rsData.size() == numCartesianCells);
#endif

        // calculate the initial fluid states
        for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto &dofFluidState = initialFluidStates_[dofIdx];

            size_t cartesianDofIdx = cartesianCellId[dofIdx];
            assert(0 <= cartesianDofIdx);
            assert(cartesianDofIdx <= numCartesianCells);

            //////
            // set temperatures
            //////
            dofFluidState.setTemperature(temperature_);

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
            // set pressures
            //////
            Scalar oilPressure = pressureData[cartesianDofIdx];
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                dofFluidState.setPressure(phaseIdx, oilPressure);
            }

            //////
            // set compositions
            //////

            // reset all mole fractions to 0
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    dofFluidState.setMoleFraction(phaseIdx, compIdx, 0.0);

            // set compositions of the gas and water phases
            dofFluidState.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);
            dofFluidState.setMoleFraction(gasPhaseIdx, gasCompIdx, 1.0);


            // set the composition of the oil phase:
            //
            // first, retrieve the relevant black-oil parameters from
            // the fluid system.
            Scalar RsSat = FluidSystem::gasDissolutionFactor(oilPressure, /*regionIdx=*/0);
            Scalar RsReal = rsData[cartesianDofIdx];

            if (RsReal > RsSat) {
                std::array<int, 3> ijk;
                gridManager.getIJK(dofIdx, ijk);
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
    }

    std::vector<Scalar> porosity_;
    std::vector<DimMatrix> intrinsicPermeability_;
    EclTransmissibility<TypeTag> transmissibilities_;

    std::vector<unsigned short> materialParamTableIdx_;
    std::vector<MaterialLawParams> materialParams_;

    std::vector<unsigned short> rockTableIdx_;
    std::vector<RockParams> rockParams_;

    std::vector<BlackOilFluidState> initialFluidStates_;

    Scalar temperature_;

    EclWellManager<TypeTag> wellManager_;

    EclWriter<TypeTag> eclWriter_;
    EclSummaryWriter summaryWriter_;
};
} // namespace Ewoms

#endif
