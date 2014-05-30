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

#include <ewoms/models/blackoil/blackoilmodel.hh>
#include <ewoms/disc/ecfv/ecfvdiscretization.hh>

#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EclDefaultMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

// for this simulator to make sense, dune-cornerpoint and opm-parser
// must be available
#include <dune/grid/CpGrid.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Utility/SgofTable.hpp>
#include <opm/parser/eclipse/Utility/SwofTable.hpp>
#include <opm/parser/eclipse/Utility/PvtoTable.hpp>
#include <opm/parser/eclipse/Utility/PvtwTable.hpp>
#include <opm/parser/eclipse/Utility/PvdgTable.hpp>

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

NEW_TYPE_TAG(EclBaseProblem, INHERITS_FROM(EclGridManager));

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

public:
    typedef Opm::EclDefaultMaterial<Traits, GasOilLaw, OilWaterLaw> type;
};

// Enable gravity
SET_BOOL_PROP(EclBaseProblem, EnableGravity, true);

// Reuse the last linearization if possible?
SET_BOOL_PROP(EclBaseProblem, EnableLinearizationRecycling, true);

// Re-assemble the linearization only for the cells which have changed?
SET_BOOL_PROP(EclBaseProblem, EnablePartialRelinearization, true);

// only write the solutions for the report steps to disk
SET_BOOL_PROP(EclBaseProblem, EnableWriteAllSolutions, false);

// set the defaults for some problem specific properties
SET_SCALAR_PROP(EclBaseProblem, Temperature, 293.15);

// The default for the end time of the simulation [s]
//
// By default, stop after the first year...
SET_SCALAR_PROP(EclBaseProblem, EndTime, 1*365*24*60*60);

// The default for the initial time step size of the simulation [s].
//
// The chosen value means that the size of the first time step is the
// one of the initial episode (if the length of the initial episode is
// not millions of trillions of years, that is...)
SET_SCALAR_PROP(EclBaseProblem, InitialTimeStepSize, 1e100);

// Disable the VTK output by default for this problem ...
SET_BOOL_PROP(EclBaseProblem, EnableVtkOutput, false);

// ... but enable the Eclipse output by default
SET_BOOL_PROP(EclBaseProblem, EnableEclipseOutput, true);

// The default DGF file to load
SET_STRING_PROP(EclBaseProblem, GridFile, "data/ecl.DATA");
}} // namespace Properties, Opm

namespace Ewoms {
/*!
 * \ingroup TestProblems
 *
 * \brief This problem uses a deck in the format of the Eclipse
 *        simulator.
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
    typedef typename GET_PROP_TYPE(TypeTag, GridManager) GridManager;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:

    /*!
     * \copydoc FvBaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableWriteAllSolutions,
                             "Write all solutions to disk instead of only the ones for the report steps");
    }

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    EclProblem(Simulator &simulator)
        : ParentType(simulator)
    {
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);

        // invert the direction of the gravity vector for ECL problems
        // (z coodinates represent depth, not height.)
        this->gravity_[dim - 1] *= -1;

        const auto deck = this->simulator().gridManager().deck();

        initFluidSystem_(deck);
        readMaterialParameters_(deck);
        readInitialCondition_(deck);

        // Start the first episode. For this, ask the Eclipse schedule.
        Opm::TimeMapConstPtr timeMap = simulator.gridManager().schedule()->getTimeMap();
        tm curTime = boost::posix_time::to_tm(timeMap->getStartTime(/*timeStepIdx=*/0));

        Scalar startTime = std::mktime(&curTime);
        simulator.setStartTime(startTime);
        simulator.startNextEpisode(/*startTime=*/startTime,
                                   /*length=*/timeMap->getTimeStepLength(/*timeStepIdx=*/0));

        // we want the episode index to be the same as the report step
        // index to make things simpler...
        simulator.setEpisodeIndex(0);

        // the user-specified initial time step can be shorter than
        // the initial report step size given in the deck, but it
        // can't be longer.
        Scalar dt = simulator.timeStepSize();
        if (dt > simulator.episodeLength())
            simulator.setTimeStepSize(simulator.episodeLength());
    }

    /*!
     * \brief Called by the time manager after the end of an episode.
     */
    void episodeEnd()
    {
        Simulator &simulator = this->simulator();
        Opm::TimeMapConstPtr timeMap = simulator.gridManager().schedule()->getTimeMap();

        // TimeMap deals with points in time, so the number of time
        // intervalls (i.e., report steps) is one less!
        int numReportSteps = timeMap->size() - 1;

        // start the next episode if there are additional report
        // steps, else finish the simulation
        int episodeIdx = simulator.episodeIndex();
        if (episodeIdx < numReportSteps)
            simulator.startNextEpisode(timeMap->getTimeStepLength(episodeIdx + 1));
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
        if (this->simulator().timeStepIndex() == 0)
            // always write the initial solution
            return true;

        if (EWOMS_GET_PARAM(TypeTag, bool, EnableWriteAllSolutions))
            return true;

        return this->simulator().episodeWillBeOver();
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
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        int globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return porosity_[globalSpaceIdx];
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
     * Eclipse uses no-flow conditions for all boundaries. \todo really?
     */
    template <class Context>
    void boundary(BoundaryRateVector &values,
                  const Context &context,
                  int spaceIdx,
                  int timeIdx) const
    { values.setNoFlow(); }

    //! \}

    /*!
     * \name Volume terms
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
    void source(RateVector &rate, const Context &context, int spaceIdx,
                int timeIdx) const
    {
#warning "TODO: wells"
        rate = Scalar(0.0);
    }

    //! \}

private:
    void readMaterialParameters_(Opm::DeckConstPtr deck)
    {
        size_t numDof = this->model().numDof();

        intrinsicPermeability_.resize(numDof);
        porosity_.resize(numDof);
        materialParams_.resize(numDof);

        // read the intrinsic permeabilities from the deck
        if (deck->hasKeyword("PERM")) {
            // the PERM and PERM{X,Y,Z,{X,Y,Z}{X,Y,Z}} keywords are
            // mutually exclusive, but if the deck does shit, it is
            // not our fault!
            const std::vector<double> &permData =
                deck->getKeyword("PERM")->getSIDoubleData();

            assert(permData.size() == numDof);

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx)
                intrinsicPermeability_[dofIdx] = this->toDimMatrix_(permData[dofIdx]);
        }
        else if (deck->hasKeyword("PERMX")) {
            const std::vector<double> &permxData =
                deck->getKeyword("PERMX")->getSIDoubleData();
            std::vector<double> permyData(permxData);
            if (deck->hasKeyword("PERMY"))
                permyData = deck->getKeyword("PERMY")->getSIDoubleData();
            std::vector<double> permzData(permxData);
            if (deck->hasKeyword("PERMZ"))
                permzData = deck->getKeyword("PERMZ")->getSIDoubleData();

            assert(permxData.size() == numDof);
            assert(permyData.size() == numDof);
            assert(permzData.size() == numDof);

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
                intrinsicPermeability_[dofIdx] = 0.0;
                intrinsicPermeability_[dofIdx][0][0] = permxData[dofIdx];
                intrinsicPermeability_[dofIdx][1][1] = permyData[dofIdx];
                intrinsicPermeability_[dofIdx][2][2] = permzData[dofIdx];
            }

            // we don't care about non-diagonal entries
        }
        else
            OPM_THROW(std::logic_error,
                      "Can't read the intrinsic permeability from the deck. "
                      "(The PERM* keywords are missing)");

        if (deck->hasKeyword("PORO")) {
            const std::vector<double> &poroData =
                deck->getKeyword("PORO")->getSIDoubleData();

            assert(poroData.size() == numDof);

            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx)
                porosity_[dofIdx] = poroData[dofIdx];
        }
        else
            OPM_THROW(std::logic_error,
                      "Can't read the porosity from the deck. "
                      "(The PORO keyword is missing)");

        Opm::DeckKeywordConstPtr swofKeyword = deck->getKeyword("SWOF");
        Opm::DeckKeywordConstPtr sgofKeyword = deck->getKeyword("SGOF");

        // the number of tables for the SWOF and the SGOF keywords
        // must be identical
        assert(Opm::SwofTable::numTables(swofKeyword) == Opm::SgofTable::numTables(sgofKeyword));

        typedef typename MaterialLawParams::GasOilParams GasOilParams;
        typedef typename MaterialLawParams::OilWaterParams OilWaterParams;

        size_t numSatfuncTables = Opm::SwofTable::numTables(swofKeyword);
        materialParams_.resize(numSatfuncTables);
        for (size_t tableIdx = 0; tableIdx < numSatfuncTables; ++ tableIdx) {
            // set the parameters of the material law for a given table
            OilWaterParams owParams;
            GasOilParams goParams;

            Opm::SwofTable swofTable(swofKeyword, tableIdx);
            Opm::SgofTable sgofTable(sgofKeyword, tableIdx);

            owParams.setSwSamples(swofTable.getSwColumn());
            owParams.setKrwSamples(swofTable.getKrwColumn());
            owParams.setKrnSamples(swofTable.getKrowColumn());
            owParams.setPcnwSamples(swofTable.getPcowColumn());

            // convert the saturations from gas to oil saturations
            auto SoSamples = sgofTable.getSgColumn();
            for (size_t sampleIdx = 0; sampleIdx < SoSamples.size(); ++ sampleIdx) {
                SoSamples[sampleIdx] = 1 - SoSamples[sampleIdx];
            }
            goParams.setSwSamples(SoSamples);
            goParams.setKrwSamples(sgofTable.getKrogColumn());
            goParams.setKrnSamples(sgofTable.getKrgColumn());
            goParams.setPcnwSamples(sgofTable.getPcogColumn());

            owParams.finalize();
            goParams.finalize();

            materialParams_[tableIdx].setOilWaterParams(owParams);
            materialParams_[tableIdx].setGasOilParams(goParams);

            materialParams_[tableIdx].finalize();
        }

        // set the index of the table to be used
        if (deck->hasKeyword("SATNUM")) {
            const std::vector<int> &satnumData =
                deck->getKeyword("SATNUM")->getRecord(0)->getItem(0)->getIntData();
            assert(satnumData.size() == numDof);

            materialParamTableIdx_.resize(numDof);
            for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
                // make sure that all values are in the correct range
                assert(1 <= satnumData[dofIdx]);
                assert(satnumData[dofIdx] <= static_cast<int>(numSatfuncTables));

                // Eclipse uses Fortran-style indices which start at
                // 1, but this here is C++...
                materialParamTableIdx_[dofIdx] = satnumData[dofIdx] - 1;
            }
        }
        else
            materialParamTableIdx_.clear();
    }

    void initFluidSystem_(Opm::DeckConstPtr deck)
    {
        FluidSystem::initBegin();

        // so far, we require the presence of the PVTO, PVTW and PVDG
        // keywords...
        Opm::PvtoTable pvtoTable(deck->getKeyword("PVTO"), /*tableIdx=*/0);
        Opm::PvtwTable pvtwTable(deck->getKeyword("PVTW"));
        Opm::PvdgTable pvdgTable(deck->getKeyword("PVDG"));

        FluidSystem::setPvtoTable(pvtoTable);
        FluidSystem::setPvtwTable(pvtwTable);
        FluidSystem::setPvdgTable(pvdgTable);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            FluidSystem::setReferenceVolumeFactor(phaseIdx, 1.0);

        // set the reference densities
        Opm::DeckRecordConstPtr densityRecord = deck->getKeyword("DENSITY")->getRecord(0);
        FluidSystem::setSurfaceDensities(densityRecord->getItem("OIL")->getSIDouble(0),
                                         densityRecord->getItem("WATER")->getSIDouble(0),
                                         densityRecord->getItem("GAS")->getSIDouble(0));

        FluidSystem::initEnd();
   }

    void readInitialCondition_(Opm::DeckConstPtr deck)
    {
        size_t numDof = this->model().numDof();

        initialFluidStates_.resize(numDof);

        if (!deck->hasKeyword("SWAT") ||
            !deck->hasKeyword("SGAS"))
            OPM_THROW(std::runtime_error,
                      "So far, the Eclipse input file requires the presence of the SWAT "
                      "and SGAS keywords");
        if (!deck->hasKeyword("PRESSURE"))
            OPM_THROW(std::runtime_error,
                      "So far, the Eclipse input file requires the presence of the PRESSURE "
                      "keyword");

        const std::vector<double> &waterSaturationData =
            deck->getKeyword("SWAT")->getSIDoubleData();
        const std::vector<double> &gasSaturationData =
            deck->getKeyword("SGAS")->getSIDoubleData();
        const std::vector<double> &pressureData =
            deck->getKeyword("PRESSURE")->getSIDoubleData();

        // make sure that the size of the data arrays is correct
        assert(waterSaturationData.size() == numDof);
        assert(gasSaturationData.size() == numDof);
        assert(pressureData.size() == numDof);

        // calculate the initial fluid states
        for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto &dofFluidState = initialFluidStates_[dofIdx];

            //////
            // set temperatures
            //////
            dofFluidState.setTemperature(temperature_);

            //////
            // set saturations
            //////
            dofFluidState.setSaturation(FluidSystem::waterPhaseIdx,
                                        waterSaturationData[dofIdx]);
            dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                        gasSaturationData[dofIdx]);
            dofFluidState.setSaturation(FluidSystem::oilPhaseIdx,
                                        1
                                        - waterSaturationData[dofIdx]
                                        - gasSaturationData[dofIdx]);

            //////
            // set pressures
            //////
            Scalar oilPressure = pressureData[dofIdx];
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
            Scalar Bo = FluidSystem::oilFormationVolumeFactor(oilPressure);
            Scalar Rs = FluidSystem::gasDissolutionFactor(oilPressure);
            Scalar rhoo = FluidSystem::surfaceDensity(oilPhaseIdx) / Bo;
            Scalar rhogref = FluidSystem::surfaceDensity(gasPhaseIdx);

            // calculate composition of oil phase in terms of mass
            // fractions.
            Scalar XoG = Rs * rhogref / rhoo;

            // convert mass to mole fractions
            Scalar MG = FluidSystem::molarMass(gasCompIdx);
            Scalar MO = FluidSystem::molarMass(oilCompIdx);

            Scalar xoG = XoG * MO / ((MO - MG) * XoG + MG);
            Scalar xoO = 1 - xoG;

            // finally set the oil-phase composition
            dofFluidState.setMoleFraction(oilPhaseIdx, gasCompIdx, xoG);
            dofFluidState.setMoleFraction(oilPhaseIdx, oilCompIdx, xoO);
        }
    }

    std::vector<Scalar> porosity_;
    std::vector<DimMatrix> intrinsicPermeability_;

    std::vector<unsigned short> materialParamTableIdx_;
    std::vector<MaterialLawParams> materialParams_;

    std::vector<BlackOilFluidState> initialFluidStates_;

    Scalar temperature_;
};
} // namespace Ewoms

#endif
