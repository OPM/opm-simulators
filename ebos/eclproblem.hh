// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2023 INRIA
  
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
 * \copydoc Opm::EclProblem
 */
#ifndef EWOMS_ECL_PROBLEM_HH
#define EWOMS_ECL_PROBLEM_HH

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <ebos/eclbaseaquifermodel.hh>
#include <ebos/eclcpgridvanguard.hh>
#include <ebos/ecldummygradientcalculator.hh>
#include <ebos/eclequilinitializer.hh>
#include <ebos/eclfluxmodule.hh>
#include <ebos/eclgenericproblem.hh>
#include <ebos/eclnewtonmethod.hh>
#include <ebos/ecloutputblackoilmodule.hh>
#include <ebos/eclproblem_properties.hh>
#include <ebos/eclthresholdpressure.hh>
#include <ebos/ecltransmissibility.hh>
#include <ebos/eclwriter.hh>
#if HAVE_DAMARIS
#include <ebos/damariswriter.hh>
#endif
#include <ebos/ecltracermodel.hh>
#include <ebos/FIBlackOilModel.hpp>
#include <ebos/vtkecltracermodule.hh>

#include <opm/common/utility/TimeService.hpp>

#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/E.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/material/common/ConditionalStorage.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WetGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DeadOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityWaterPvt.hpp>
#include <opm/material/thermal/EclThermalLawManager.hpp>

#include <opm/models/common/directionalmobility.hh>
#include <opm/models/utils/pffgridvector.hh>
#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>

#include <opm/output/eclipse/EclipseIO.hpp>

#include <opm/simulators/flow/EclActionHandler.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <opm/utility/CopyablePtr.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <algorithm>
#include <functional>
#include <set>
#include <string>
#include <vector>

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This problem simulates an input file given in the data format used by the
 *        commercial ECLiPSE simulator.
 */
template <class TypeTag>
class EclProblem : public GetPropType<TypeTag, Properties::BaseProblem>
                 , public EclGenericProblem<GetPropType<TypeTag, Properties::GridView>,
                                            GetPropType<TypeTag, Properties::FluidSystem>,
                                            GetPropType<TypeTag, Properties::Scalar>>
{
    using BaseType = EclGenericProblem<GetPropType<TypeTag, Properties::GridView>,
                                       GetPropType<TypeTag, Properties::FluidSystem>,
                                       GetPropType<TypeTag, Properties::Scalar>>;
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;

    // Grid and world dimension
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { enableExperiments = getPropValue<TypeTag, Properties::EnableExperiments>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };
    enum { enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>() };
    enum { enableBrine = getPropValue<TypeTag, Properties::EnableBrine>() };
    enum { enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>() };
    enum { enablePolymerMolarWeight = getPropValue<TypeTag, Properties::EnablePolymerMW>() };
    enum { enableFoam = getPropValue<TypeTag, Properties::EnableFoam>() };
    enum { enableExtbo = getPropValue<TypeTag, Properties::EnableExtbo>() };
    enum { enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableDispersion = getPropValue<TypeTag, Properties::EnableDispersion>() };
    enum { enableThermalFluxBoundaries = getPropValue<TypeTag, Properties::EnableThermalFluxBoundaries>() };
    enum { enableApiTracking = getPropValue<TypeTag, Properties::EnableApiTracking>() };
    enum { enableMICP = getPropValue<TypeTag, Properties::EnableMICP>() };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using EclMaterialLawManager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
    using EclThermalLawManager = typename GetProp<TypeTag, Properties::SolidEnergyLaw>::EclThermalLawManager;
    using MaterialLawParams = typename EclMaterialLawManager::MaterialLawParams;
    using SolidEnergyLawParams = typename EclThermalLawManager::SolidEnergyLawParams;
    using ThermalConductionLawParams = typename EclThermalLawManager::ThermalConductionLawParams;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using DofMapper = GetPropType<TypeTag, Properties::DofMapper>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using EclWellModel = GetPropType<TypeTag, Properties::EclWellModel>;
    using EclAquiferModel = GetPropType<TypeTag, Properties::EclAquiferModel>;

    using SolventModule = BlackOilSolventModule<TypeTag>;
    using PolymerModule = BlackOilPolymerModule<TypeTag>;
    using FoamModule = BlackOilFoamModule<TypeTag>;
    using BrineModule = BlackOilBrineModule<TypeTag>;
    using ExtboModule = BlackOilExtboModule<TypeTag>;
    using MICPModule = BlackOilMICPModule<TypeTag>;

    using InitialFluidState = typename EclEquilInitializer<TypeTag>::ScalarFluidState;

    using Toolbox = MathToolbox<Evaluation>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using EclWriterType = EclWriter<TypeTag>;
#if HAVE_DAMARIS
    using DamarisWriterType = DamarisWriter<TypeTag>;
#endif

    using TracerModel = EclTracerModel<TypeTag>;
    using DirectionalMobilityPtr = Opm::Utility::CopyablePtr<DirectionalMobility<TypeTag, Evaluation>>;

public:
    using EclGenericProblem<GridView,FluidSystem,Scalar>::briefDescription;
    using EclGenericProblem<GridView,FluidSystem,Scalar>::helpPreamble;
    using EclGenericProblem<GridView,FluidSystem,Scalar>::shouldWriteOutput;
    using EclGenericProblem<GridView,FluidSystem,Scalar>::shouldWriteRestartFile;
    using EclGenericProblem<GridView,FluidSystem,Scalar>::rockCompressibility;
    using EclGenericProblem<GridView,FluidSystem,Scalar>::rockReferencePressure;
    using EclGenericProblem<GridView,FluidSystem,Scalar>::porosity;

    /*!
     * \copydoc FvBaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
        EclWriterType::registerParameters();
#if HAVE_DAMARIS
        DamarisWriterType::registerParameters();
#endif

        VtkEclTracerModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableWriteAllSolutions,
                             "Write all solutions to disk instead of only the ones for the "
                             "report steps");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableEclOutput,
                             "Write binary output which is compatible with the commercial "
                             "Eclipse simulator");
#if HAVE_DAMARIS
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableDamarisOutput,
                             "Write a specific variable using Damaris in a separate core");
#endif
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclOutputDoublePrecision,
                             "Tell the output writer to use double precision. Useful for 'perfect' restarts");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, RestartWritingInterval,
                             "The frequencies of which time steps are serialized to disk");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclEnableDriftCompensation,
                             "Enable partial compensation of systematic mass losses via the source term of the next time step");
        if constexpr (enableExperiments)
            EWOMS_REGISTER_PARAM(TypeTag, bool, EclEnableAquifers,
                                 "Enable analytic and numeric aquifer models");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclEnableTuning,
                             "Honor some aspects of the TUNING keyword from the ECL deck.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, OutputMode,
                             "Specify which messages are going to be printed. Valid values are: none, log, all (default)");
        EWOMS_REGISTER_PARAM(TypeTag, int, NumPressurePointsEquil,
                             "Number of pressure points (in each direction) in tables used for equilibration");
        EWOMS_HIDE_PARAM(TypeTag, NumPressurePointsEquil); // Users will typically not need to modify this parameter..
    }


    /*!
     * \copydoc FvBaseProblem::handlePositionalParameter
     */
    static int handlePositionalParameter(std::set<std::string>& seenParams,
                                         std::string& errorMsg,
                                         int,
                                         const char** argv,
                                         int paramIdx,
                                         int)
    {
        using ParamsMeta = GetProp<TypeTag, Properties::ParameterMetaData>;
        Dune::ParameterTree& tree = ParamsMeta::tree();
        return eclPositionalParameter(tree,
                                      seenParams,
                                      errorMsg,
                                      argv,
                                      paramIdx);
    }

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    EclProblem(Simulator& simulator)
        : ParentType(simulator)
        , EclGenericProblem<GridView,FluidSystem,Scalar>(simulator.vanguard().eclState(),
                                                         simulator.vanguard().schedule(),
                                                         simulator.vanguard().gridView())
        , transmissibilities_(simulator.vanguard().eclState(),
                              simulator.vanguard().gridView(),
                              simulator.vanguard().cartesianIndexMapper(),
                              simulator.vanguard().grid(),
                              simulator.vanguard().cellCentroids(),
                              enableEnergy,
                              enableDiffusion,
                              enableDispersion)
        , thresholdPressures_(simulator)
        , wellModel_(simulator)
        , aquiferModel_(simulator)
        , pffDofData_(simulator.gridView(), this->elementMapper())
        , tracerModel_(simulator)
        , actionHandler_(simulator.vanguard().eclState(),
                         simulator.vanguard().schedule(),
                         simulator.vanguard().actionState(),
                         simulator.vanguard().summaryState(),
                         wellModel_,
                         simulator.vanguard().grid().comm())
    {
        this->model().addOutputModule(new VtkEclTracerModule<TypeTag>(simulator));
        // Tell the black-oil extensions to initialize their internal data structures
        const auto& vanguard = simulator.vanguard();
        SolventModule::initFromState(vanguard.eclState(), vanguard.schedule());
        PolymerModule::initFromState(vanguard.eclState());
        FoamModule::initFromState(vanguard.eclState());
        BrineModule::initFromState(vanguard.eclState());
        ExtboModule::initFromState(vanguard.eclState());
        MICPModule::initFromState(vanguard.eclState());

        // create the ECL writer
        eclWriter_ = std::make_unique<EclWriterType>(simulator);
#if HAVE_DAMARIS
        // create Damaris writer
        damarisWriter_ = std::make_unique<DamarisWriterType>(simulator);
        enableDamarisOutput_ = EWOMS_GET_PARAM(TypeTag, bool, EnableDamarisOutput) ;
#endif
        enableDriftCompensation_ = EWOMS_GET_PARAM(TypeTag, bool, EclEnableDriftCompensation);

        enableEclOutput_ = EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput);

        if constexpr (enableExperiments)
            enableAquifers_ = EWOMS_GET_PARAM(TypeTag, bool, EclEnableAquifers);
        else
            enableAquifers_ = true;

        this->enableTuning_ = EWOMS_GET_PARAM(TypeTag, bool, EclEnableTuning);
        this->initialTimeStepSize_ = EWOMS_GET_PARAM(TypeTag, Scalar, InitialTimeStepSize);
        this->maxTimeStepAfterWellEvent_ = EWOMS_GET_PARAM(TypeTag, double, TimeStepAfterEventInDays)*24*60*60;

        // The value N for this parameter is defined in the following order of presedence:
        // 1. Command line value (--num-pressure-points-equil=N)
        // 2. EQLDIMS item 2
        // Default value is defined in opm-common/src/opm/input/eclipse/share/keywords/000_Eclipse100/E/EQLDIMS
        if (EWOMS_PARAM_IS_SET(TypeTag, int, NumPressurePointsEquil))
        {
            this->numPressurePointsEquil_ = EWOMS_GET_PARAM(TypeTag, int, NumPressurePointsEquil);
        } else {
            this->numPressurePointsEquil_ = simulator.vanguard().eclState().getTableManager().getEqldims().getNumDepthNodesP();
        }

        RelpermDiagnostics relpermDiagnostics;
        relpermDiagnostics.diagnosis(vanguard.eclState(), vanguard.cartesianIndexMapper());
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        auto& simulator = this->simulator();
        const auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();

        // Set the start time of the simulation
        simulator.setStartTime(schedule.getStartTime());
        simulator.setEndTime(schedule.simTime(schedule.size() - 1));

        // We want the episode index to be the same as the report step index to make
        // things simpler, so we have to set the episode index to -1 because it is
        // incremented by endEpisode(). The size of the initial time step and
        // length of the initial episode is set to zero for the same reason.
        simulator.setEpisodeIndex(-1);
        simulator.setEpisodeLength(0.0);

        // the "NOGRAV" keyword from Frontsim or setting the EnableGravity to false
        // disables gravity, else the standard value of the gravity constant at sea level
        // on earth is used
        this->gravity_ = 0.0;
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity))
            this->gravity_[dim - 1] = 9.80665;
        if (!eclState.getInitConfig().hasGravity())
            this->gravity_[dim - 1] = 0.0;

        if (this->enableTuning_) {
            // if support for the TUNING keyword is enabled, we get the initial time
            // steping parameters from it instead of from command line parameters
            const auto& tuning = schedule[0].tuning();
            this->initialTimeStepSize_ = tuning.TSINIT.has_value() ? tuning.TSINIT.value() : -1.0;
            this->maxTimeStepAfterWellEvent_ = tuning.TMAXWC;
        }

        this->initFluidSystem_();

        // deal with DRSDT
        this->mixControls_.init(this->model().numGridDof(),
                                this->episodeIndex(),
                                eclState.runspec().tabdims().getNumPVTTables());

        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
            FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            this->maxOilSaturation_.resize(this->model().numGridDof(), 0.0);
        }

        this->readRockParameters_(simulator.vanguard().cellCenterDepths(),
                                  [&simulator](const unsigned idx)
                                  {
                                      std::array<int,dim> coords;
                                      simulator.vanguard().cartesianCoordinate(idx, coords);
                                      for (auto& c : coords) {
                                          ++c;
                                      }
                                      return coords;
                                  });
        readMaterialParameters_();
        readThermalParameters_();

        // Re-ordering in case of ALUGrid
        std::function<unsigned int(unsigned int)> gridToEquilGrid = [&simulator](unsigned int i) {
            return simulator.vanguard().gridIdxToEquilGridIdx(i);
        };
        transmissibilities_.finishInit(gridToEquilGrid);

        const auto& initconfig = eclState.getInitConfig();
        tracerModel_.init(initconfig.restartRequested());
        if (initconfig.restartRequested())
            readEclRestartSolution_();
        else
            readInitialCondition_();

        tracerModel_.prepareTracerBatches();

        updatePffDofData_();

        if constexpr (getPropValue<TypeTag, Properties::EnablePolymer>()) {
            const auto& vanguard = this->simulator().vanguard();
            const auto& gridView = vanguard.gridView();
            int numElements = gridView.size(/*codim=*/0);
            this->polymer_.maxAdsorption.resize(numElements, 0.0);
        }

        readBoundaryConditions_();

        // compute and set eq weights based on initial b values
        computeAndSetEqWeights_();

        if (enableDriftCompensation_) {
            drift_.resize(this->model().numGridDof());
            drift_ = 0.0;
        }

        // write the static output files (EGRID, INIT, SMSPEC, etc.)
        if (enableEclOutput_) {
            if (simulator.vanguard().grid().comm().size() > 1) {
                if (simulator.vanguard().grid().comm().rank() == 0)
                    eclWriter_->setTransmissibilities(&simulator.vanguard().globalTransmissibility());
            } else
                eclWriter_->setTransmissibilities(&simulator.problem().eclTransmissibilities());

            // Re-ordering in case of ALUGrid
            std::function<unsigned int(unsigned int)> equilGridToGrid = [&simulator](unsigned int i) {
                return simulator.vanguard().gridEquilIdxToGridIdx(i);
            };
            eclWriter_->writeInit(equilGridToGrid);
        }

        simulator.vanguard().releaseGlobalTransmissibilities();

        // after finishing the initialization and writing the initial solution, we move
        // to the first "real" episode/report step
        // for restart the episode index and start is already set
        if (!initconfig.restartRequested()) {
            simulator.startNextEpisode(schedule.seconds(1));
            simulator.setEpisodeIndex(0);
        }
    }

    void prefetch(const Element& elem) const
    { pffDofData_.prefetch(elem); }

    /*!
     * \brief This method restores the complete state of the problem and its sub-objects
     *        from disk.
     *
     * The serialization format used by this method is ad-hoc. It is the inverse of the
     * serialize() method.
     *
     * \tparam Restarter The deserializer type
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        // reload the current episode/report step from the deck
        beginEpisode();

        // deserialize the wells
        wellModel_.deserialize(res);

        if (enableAquifers_)
            // deserialize the aquifer
            aquiferModel_.deserialize(res);
    }

    /*!
     * \brief This method writes the complete state of the problem and its subobjects to
     *        disk.
     *
     * The file format used here is ad-hoc.
     */
    template <class Restarter>
    void serialize(Restarter& res)
    {
        wellModel_.serialize(res);

        if (enableAquifers_)
            aquiferModel_.serialize(res);
    }

    int episodeIndex() const
    {
        return std::max(this->simulator().episodeIndex(), 0);
    }

    /*!
     * \brief Called by the simulator before an episode begins.
     */
    void beginEpisode()
    {
        OPM_TIMEBLOCK(beginEpisode);
        // Proceed to the next report step
        auto& simulator = this->simulator();
        int episodeIdx = simulator.episodeIndex();
        auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();
        const auto& events = schedule[episodeIdx].events();

        if (episodeIdx >= 0 && events.hasEvent(ScheduleEvents::GEO_MODIFIER)) {
            // bring the contents of the keywords to the current state of the SCHEDULE
            // section.
            //
            // TODO (?): make grid topology changes possible (depending on what exactly
            // has changed, the grid may need be re-created which has some serious
            // implications on e.g., the solution of the simulation.)
            const auto& miniDeck = schedule[episodeIdx].geo_keywords();
            const auto& cc = simulator.vanguard().grid().comm();
            eclState.apply_schedule_keywords( miniDeck );
            eclBroadcast(cc, eclState.getTransMult() );

            // Re-ordering in case of ALUGrid
            std::function<unsigned int(unsigned int)> equilGridToGrid = [&simulator](unsigned int i) {
                  return simulator.vanguard().gridEquilIdxToGridIdx(i);
            };

            // re-compute all quantities which may possibly be affected.
            transmissibilities_.update(true, equilGridToGrid);
            this->referencePorosity_[1] = this->referencePorosity_[0];
            updateReferencePorosity_();
            updatePffDofData_();
            this->model().linearizer().updateDiscretizationParameters();
        }

        bool tuningEvent = this->beginEpisode_(enableExperiments, this->episodeIndex());

        // set up the wells for the next episode.
        wellModel_.beginEpisode();

        // set up the aquifers for the next episode.
        if (enableAquifers_)
            // set up the aquifers for the next episode.
            aquiferModel_.beginEpisode();

        // set the size of the initial time step of the episode
        Scalar dt = limitNextTimeStepSize_(simulator.episodeLength());
        // negative value of initialTimeStepSize_ indicates no active limit from TSINIT or NEXTSTEP
        if ( (episodeIdx == 0 || tuningEvent) && this->initialTimeStepSize_ > 0)
            // allow the size of the initial time step to be set via an external parameter
            // if TUNING is enabled, also limit the time step size after a tuning event to TSINIT
            dt = std::min(dt, this->initialTimeStepSize_);
        simulator.setTimeStepSize(dt);

        // Evaluate UDQ assign statements to make sure the settings are
        // available as UDA controls for the current report step.
        actionHandler_.evalUDQAssignments(episodeIdx, simulator.vanguard().udqState());

        if (episodeIdx >= 0) {
            const auto& oilVap = schedule[episodeIdx].oilvap();
            if (oilVap.getType() == OilVaporizationProperties::OilVaporization::VAPPARS) {
                FluidSystem::setVapPars(oilVap.vap1(), oilVap.vap2());
            } else {
                FluidSystem::setVapPars(0.0, 0.0);
            }
        }
    }

    /*!
     * \brief Called by the simulator before each time integration.
     */
    void beginTimeStep()
    {
        OPM_TIMEBLOCK(beginTimeStep);
        int episodeIdx = this->episodeIndex();

        this->beginTimeStep_(enableExperiments,
                             episodeIdx,
                             this->simulator().timeStepIndex(),
                             this->simulator().startTime(),
                             this->simulator().time(),
                             this->simulator().timeStepSize(),
                             this->simulator().endTime());

        // update maximum water saturation and minimum pressure
        // used when ROCKCOMP is activated
        asImp_().updateExplicitQuantities_();

        if (nonTrivialBoundaryConditions()) {
            this->model().linearizer().updateBoundaryConditionData();
        }

        wellModel_.beginTimeStep();
        if (enableAquifers_)
            aquiferModel_.beginTimeStep();
        tracerModel_.beginTimeStep();

    }

    /*!
     * \brief Called by the simulator before each Newton-Raphson iteration.
     */
    void beginIteration()
    {
        OPM_TIMEBLOCK(beginIteration);
        wellModel_.beginIteration();
        if (enableAquifers_)
            aquiferModel_.beginIteration();
    }

    /*!
     * \brief Called by the simulator after each Newton-Raphson iteration.
     */
    void endIteration()
    {
        OPM_TIMEBLOCK(endIteration);
        wellModel_.endIteration();
        if (enableAquifers_)
            aquiferModel_.endIteration();
    }

    /*!
     * \brief Called by the simulator after each time integration.
     */
    void endTimeStep()
    {
        OPM_TIMEBLOCK(endTimeStep);
#ifndef NDEBUG
        if constexpr (getPropValue<TypeTag, Properties::EnableDebuggingChecks>()) {
            // in debug mode, we don't care about performance, so we check if the model does
            // the right thing (i.e., the mass change inside the whole reservoir must be
            // equivalent to the fluxes over the grid's boundaries plus the source rates
            // specified by the problem)
            int rank = this->simulator().gridView().comm().rank();
            if (rank == 0)
                std::cout << "checking conservativeness of solution\n";
            this->model().checkConservativeness(/*tolerance=*/-1, /*verbose=*/true);
            if (rank == 0)
                std::cout << "solution is sufficiently conservative\n";
        }
#endif // NDEBUG

        auto& simulator = this->simulator();
        wellModel_.endTimeStep();
        if (enableAquifers_)
            aquiferModel_.endTimeStep();
        tracerModel_.endTimeStep();


        // Compute flux for output
        this->model().linearizer().updateFlowsInfo();

        // deal with DRSDT and DRVDT
        asImp_().updateCompositionChangeLimits_();
        {
        OPM_TIMEBLOCK(driftCompansation);
        if (enableDriftCompensation_) {
            const auto& residual = this->model().linearizer().residual();
            for (unsigned globalDofIdx = 0; globalDofIdx < residual.size(); globalDofIdx ++) {
                drift_[globalDofIdx] = residual[globalDofIdx];
                drift_[globalDofIdx] *= simulator.timeStepSize();
                if constexpr (getPropValue<TypeTag, Properties::UseVolumetricResidual>())
                    drift_[globalDofIdx] *= this->model().dofTotalVolume(globalDofIdx);
            }
        }
        }
        bool isSubStep = !EWOMS_GET_PARAM(TypeTag, bool, EnableWriteAllSolutions) && !this->simulator().episodeWillBeOver();
        eclWriter_->evalSummaryState(isSubStep);

        int episodeIdx = this->episodeIndex();

        // Re-ordering in case of Alugrid
        std::function<unsigned int(unsigned int)> gridToEquilGrid = [&simulator](unsigned int i) {
            return simulator.vanguard().gridIdxToEquilGridIdx(i);
        };

        std::function<void(bool)> transUp =
            [this,gridToEquilGrid](bool global) {
                this->transmissibilities_.update(global,gridToEquilGrid);
            };
        {
        OPM_TIMEBLOCK(applyActions);
        actionHandler_.applyActions(episodeIdx,
                                    simulator.time() + simulator.timeStepSize(),
                                    transUp);
        }
        // deal with "clogging" for the MICP model
        if constexpr (enableMICP){
          auto& model = this->model();
          const auto& residual = this->model().linearizer().residual();
          for (unsigned globalDofIdx = 0; globalDofIdx < residual.size(); globalDofIdx ++) {
            auto& phi = this->referencePorosity_[/*timeIdx=*/1][globalDofIdx];
            MICPModule::checkCloggingMICP(model, phi, globalDofIdx);
        }
      }
    }

    /*!
     * \brief Called by the simulator after the end of an episode.
     */
    void endEpisode()
    {
        OPM_TIMEBLOCK(endEpisode);
        auto& simulator = this->simulator();
        auto& schedule = simulator.vanguard().schedule();

        wellModel_.endEpisode();
        if (enableAquifers_)
            aquiferModel_.endEpisode();

        int episodeIdx = this->episodeIndex();
        // check if we're finished ...
        if (episodeIdx + 1 >= static_cast<int>(schedule.size() - 1)) {
            simulator.setFinished(true);
            return;
        }

        // .. if we're not yet done, start the next episode (report step)
        simulator.startNextEpisode(schedule.stepLength(episodeIdx + 1));
    }

    /*!
     * \brief Write the requested quantities of the current solution into the output
     *        files.
     */
    void writeOutput(const SimulatorTimer& timer, bool verbose = true)
    {
        OPM_TIMEBLOCK(problemWriteOutput);
        // use the generic code to prepare the output fields and to
        // write the desired VTK files.
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableWriteAllSolutions) || this->simulator().episodeWillBeOver()){
            ParentType::writeOutput(verbose);
        }

        bool isSubStep = !EWOMS_GET_PARAM(TypeTag, bool, EnableWriteAllSolutions) && !this->simulator().episodeWillBeOver();
        
        data::Solution localCellData = {};
#if HAVE_DAMARIS
        // N.B. the Damaris output has to be done before the ECL output as the ECL one 
        // does all kinds of std::move() relocation of data
        if (enableDamarisOutput_) {
            damarisWriter_->writeOutput(localCellData, isSubStep) ;
        }
#endif 
        if (enableEclOutput_){
            eclWriter_->writeOutput(std::move(localCellData), timer, isSubStep);
        }
    }

    void finalizeOutput() {
        OPM_TIMEBLOCK(finalizeOutput);
        // this will write all pending output to disk
        // to avoid corruption of output files
        eclWriter_.reset();
    }


    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context,
                                           unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return transmissibilities_.permeability(globalSpaceIdx);
    }

    /*!
     * \brief This method returns the intrinsic permeability tensor
     *        given a global element index.
     *
     * Its main (only?) usage is the ECL transmissibility calculation code...
     */
    const DimMatrix& intrinsicPermeability(unsigned globalElemIdx) const
    { return transmissibilities_.permeability(globalElemIdx); }

    /*!
     * \copydoc EclTransmissiblity::transmissibility
     */
    template <class Context>
    Scalar transmissibility(const Context& context,
                            [[maybe_unused]] unsigned fromDofLocalIdx,
                            unsigned toDofLocalIdx) const
    {
        assert(fromDofLocalIdx == 0);
        return pffDofData_.get(context.element(), toDofLocalIdx).transmissibility;
    }

    /*!
     * \brief Direct access to the transmissibility between two elements.
     */
    Scalar transmissibility(unsigned globalCenterElemIdx, unsigned globalElemIdx) const
    {
        return transmissibilities_.transmissibility(globalCenterElemIdx, globalElemIdx);
    }

    /*!
     * \copydoc EclTransmissiblity::diffusivity
     */
    template <class Context>
    Scalar diffusivity(const Context& context,
                       [[maybe_unused]] unsigned fromDofLocalIdx,
                       unsigned toDofLocalIdx) const
    {
        assert(fromDofLocalIdx == 0);
        return *pffDofData_.get(context.element(), toDofLocalIdx).diffusivity;
    }

    /*!
     * give the transmissibility for a face i.e. pair. should be symmetric?
     */
    Scalar diffusivity(const unsigned globalCellIn, const unsigned globalCellOut) const{
        return transmissibilities_.diffusivity(globalCellIn, globalCellOut);
    }

    /*!
     * give the dispersivity for a face i.e. pair.
     */
    Scalar dispersivity(const unsigned globalCellIn, const unsigned globalCellOut) const{
        return transmissibilities_.dispersivity(globalCellIn, globalCellOut);
    }

    /*!
     * \brief Direct access to a boundary transmissibility.
     */
    Scalar thermalTransmissibilityBoundary(const unsigned globalSpaceIdx,
                                    const unsigned boundaryFaceIdx) const
    {
        return transmissibilities_.thermalTransmissibilityBoundary(globalSpaceIdx, boundaryFaceIdx);
    }




    /*!
     * \copydoc EclTransmissiblity::transmissibilityBoundary
     */
    template <class Context>
    Scalar transmissibilityBoundary(const Context& elemCtx,
                                    unsigned boundaryFaceIdx) const
    {
        unsigned elemIdx = elemCtx.globalSpaceIndex(/*dofIdx=*/0, /*timeIdx=*/0);
        return transmissibilities_.transmissibilityBoundary(elemIdx, boundaryFaceIdx);
    }

    /*!
     * \brief Direct access to a boundary transmissibility.
     */
    Scalar transmissibilityBoundary(const unsigned globalSpaceIdx,
                                    const unsigned boundaryFaceIdx) const
    {
        return transmissibilities_.transmissibilityBoundary(globalSpaceIdx, boundaryFaceIdx);
    }


    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    Scalar thermalHalfTransmissibility(const unsigned globalSpaceIdxIn,
                                       const unsigned globalSpaceIdxOut) const
    {
        return transmissibilities_.thermalHalfTrans(globalSpaceIdxIn,globalSpaceIdxOut);
    }

    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    template <class Context>
    Scalar thermalHalfTransmissibilityIn(const Context& context,
                                         unsigned faceIdx,
                                         unsigned timeIdx) const
    {
        const auto& face = context.stencil(timeIdx).interiorFace(faceIdx);
        unsigned toDofLocalIdx = face.exteriorIndex();
        return *pffDofData_.get(context.element(), toDofLocalIdx).thermalHalfTransIn;
    }

    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    template <class Context>
    Scalar thermalHalfTransmissibilityOut(const Context& context,
                                          unsigned faceIdx,
                                          unsigned timeIdx) const
    {
        const auto& face = context.stencil(timeIdx).interiorFace(faceIdx);
        unsigned toDofLocalIdx = face.exteriorIndex();
        return *pffDofData_.get(context.element(), toDofLocalIdx).thermalHalfTransOut;
    }

    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    template <class Context>
    Scalar thermalHalfTransmissibilityBoundary(const Context& elemCtx,
                                               unsigned boundaryFaceIdx) const
    {
        unsigned elemIdx = elemCtx.globalSpaceIndex(/*dofIdx=*/0, /*timeIdx=*/0);
        return transmissibilities_.thermalHalfTransBoundary(elemIdx, boundaryFaceIdx);
    }

    /*!
     * \brief Return a reference to the object that handles the "raw" transmissibilities.
     */
    const typename Vanguard::TransmissibilityType& eclTransmissibilities() const
    { return transmissibilities_; }

    /*!
     * \copydoc BlackOilBaseProblem::thresholdPressure
     */
    Scalar thresholdPressure(unsigned elem1Idx, unsigned elem2Idx) const
    { return thresholdPressures_.thresholdPressure(elem1Idx, elem2Idx); }

    const EclThresholdPressure<TypeTag>& thresholdPressure() const
    { return thresholdPressures_; }

    EclThresholdPressure<TypeTag>& thresholdPressure()
    { return thresholdPressures_; }

    const EclTracerModel<TypeTag>& tracerModel() const
    { return tracerModel_; }

    EclTracerModel<TypeTag>& tracerModel()
    { return tracerModel_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     *
     * For the EclProblem, this method is identical to referencePorosity(). The intensive
     * quantities object may apply various multipliers (e.g. ones which model rock
     * compressibility and water induced rock compaction) to it which depend on the
     * current physical conditions.
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->porosity(globalSpaceIdx, timeIdx);
    }

    /*!
     * \brief Returns the depth of an degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depth of an element and is
     * thus slightly different from the depth of an element's centroid.
     */
    template <class Context>
    Scalar dofCenterDepth(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->dofCenterDepth(globalSpaceIdx);
    }

    /*!
     * \brief Direct indexed acces to the depth of an degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depth of an element and is
     * thus slightly different from the depth of an element's centroid.
     */
    Scalar dofCenterDepth(unsigned globalSpaceIdx) const
    {
        return this->simulator().vanguard().cellCenterDepth(globalSpaceIdx);
    }

    /*!
     * \copydoc BlackoilProblem::rockCompressibility
     */
    template <class Context>
    Scalar rockCompressibility(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->rockCompressibility(globalSpaceIdx);
    }

    /*!
     * \copydoc BlackoilProblem::rockReferencePressure
     */
    template <class Context>
    Scalar rockReferencePressure(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->rockReferencePressure(globalSpaceIdx);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return this->materialLawParams(globalSpaceIdx);
    }

    const MaterialLawParams& materialLawParams(unsigned globalDofIdx) const
    {
        return materialLawManager_->materialLawParams(globalDofIdx);
    }

    const MaterialLawParams& materialLawParams(unsigned globalDofIdx, FaceDir::DirEnum facedir) const
    {
        return materialLawManager_->materialLawParams(globalDofIdx, facedir);
    }

    /*!
     * \brief Return the parameters for the energy storage law of the rock
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context,
                         unsigned spaceIdx,
                         unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return thermalLawManager_->solidEnergyLawParams(globalSpaceIdx);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams &
    thermalConductionLawParams(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return thermalLawManager_->thermalConductionLawParams(globalSpaceIdx);
    }

    /*!
     * \brief Returns the ECL material law manager
     *
     * Note that this method is *not* part of the generic eWoms problem API because it
     * would force all problens use the ECL material laws.
     */
    std::shared_ptr<const EclMaterialLawManager> materialLawManager() const
    { return materialLawManager_; }

    template <class FluidState>
    void updateRelperms(
        std::array<Evaluation,numPhases> &mobility,
        DirectionalMobilityPtr &dirMob,
        FluidState &fluidState,
        unsigned globalSpaceIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(updateRelperms);
        {
            // calculate relative permeabilities. note that we store the result into the
            // mobility_ class attribute. the division by the phase viscosity happens later.
            const auto& materialParams = materialLawParams(globalSpaceIdx);
            MaterialLaw::relativePermeabilities(mobility, materialParams, fluidState);
            Valgrind::CheckDefined(mobility);
        }
        if (materialLawManager_->hasDirectionalRelperms()
               || materialLawManager_->hasDirectionalImbnum())
        {
            using Dir = FaceDir::DirEnum;
            constexpr int ndim = 3;
            dirMob = std::make_unique<DirectionalMobility<TypeTag, Evaluation>>();
            Dir facedirs[ndim] = {Dir::XPlus, Dir::YPlus, Dir::ZPlus};
            for (int i = 0; i<ndim; i++) {
                const auto& materialParams = materialLawParams(globalSpaceIdx, facedirs[i]);
                auto& mob_array = dirMob->getArray(i);
                MaterialLaw::relativePermeabilities(mob_array, materialParams, fluidState);
            }
        }
    }

    /*!
     * \copydoc materialLawManager()
     */
    std::shared_ptr<EclMaterialLawManager> materialLawManager()
    { return materialLawManager_; }

    using EclGenericProblem<GridView,FluidSystem,Scalar>::pvtRegionIndex;
    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned pvtRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return pvtRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    using EclGenericProblem<GridView,FluidSystem,Scalar>::satnumRegionIndex;
    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned satnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return this->satnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    using EclGenericProblem<GridView,FluidSystem,Scalar>::miscnumRegionIndex;
    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned miscnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return this->miscnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    using EclGenericProblem<GridView,FluidSystem,Scalar>::plmixnumRegionIndex;
    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned plmixnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return this->plmixnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    using EclGenericProblem<GridView,FluidSystem,Scalar>::maxPolymerAdsorption;
    /*!
     * \brief Returns the max polymer adsorption value
     */
    template <class Context>
    Scalar maxPolymerAdsorption(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return this->maxPolymerAdsorption(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return this->simulator().vanguard().caseName(); }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        // use the initial temperature of the DOF if temperature is not a primary
        // variable
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return initialFluidStates_[globalDofIdx].temperature(/*phaseIdx=*/0);
    }


    Scalar temperature(unsigned globalDofIdx, unsigned /*timeIdx*/) const
    {
        // use the initial temperature of the DOF if temperature is not a primary
        // variable
         return initialFluidStates_[globalDofIdx].temperature(/*phaseIdx=*/0);
    }

    const SolidEnergyLawParams&
    solidEnergyLawParams(unsigned globalSpaceIdx,
                         unsigned /*timeIdx*/) const
    {
        return this->thermalLawManager_->solidEnergyLawParams(globalSpaceIdx);
    }
    const ThermalConductionLawParams &
    thermalConductionLawParams(unsigned globalSpaceIdx,
                               unsigned /*timeIdx*/)const
    {
        return this->thermalLawManager_->thermalConductionLawParams(globalSpaceIdx);
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
        OPM_TIMEBLOCK_LOCAL(eclProblemBoundary);
        if (!context.intersection(spaceIdx).boundary())
            return;

        if constexpr (!enableEnergy || !enableThermalFluxBoundaries)
            values.setNoFlow();
        else {
            // in the energy case we need to specify a non-trivial boundary condition
            // because the geothermal gradient needs to be maintained. for this, we
            // simply assume the initial temperature at the boundary and specify the
            // thermal flow accordingly. in this context, "thermal flow" means energy
            // flow due to a temerature gradient while assuming no-flow for mass
            unsigned interiorDofIdx = context.interiorScvIndex(spaceIdx, timeIdx);
            unsigned globalDofIdx = context.globalSpaceIndex(interiorDofIdx, timeIdx);
            values.setThermalFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
        }

        if (nonTrivialBoundaryConditions()) {
            unsigned indexInInside  = context.intersection(spaceIdx).indexInInside();
            unsigned interiorDofIdx = context.interiorScvIndex(spaceIdx, timeIdx);
            unsigned globalDofIdx = context.globalSpaceIndex(interiorDofIdx, timeIdx);
            unsigned pvtRegionIdx = pvtRegionIndex(context, spaceIdx, timeIdx);
            const auto [type, massrate] = boundaryCondition(globalDofIdx, indexInInside);
            if (type == BCType::THERMAL)
                values.setThermalFlow(context, spaceIdx, timeIdx, boundaryFluidState(globalDofIdx, indexInInside));
            else if (type == BCType::FREE || type == BCType::DIRICHLET)
                values.setFreeFlow(context, spaceIdx, timeIdx, boundaryFluidState(globalDofIdx, indexInInside));
            else if (type == BCType::RATE)
                values.setMassRate(massrate, pvtRegionIdx);
        }
    }

    /*!
     * \brief Returns an element's historic maximum oil phase saturation that was
     *        observed during the simulation.
     *
     * In this context, "historic" means the the time before the current timestep began.
     *
     * This is a bit of a hack from the conceptional point of view, but it is required to
     * match the results of the 'flow' and ECLIPSE 100 simulators.
     */
    Scalar maxOilSaturation(unsigned globalDofIdx) const
    {
        if (!this->vapparsActive(this->episodeIndex()))
            return 0.0;

        return this->maxOilSaturation_[globalDofIdx];
    }

    /*!
     * \brief Sets an element's maximum oil phase saturation observed during the
     *        simulation.
     *
     * In this context, "historic" means the the time before the current timestep began.
     *
     * This a hack on top of the maxOilSaturation() hack but it is currently required to
     * do restart externally. i.e. from the flow code.
     */
    void setMaxOilSaturation(unsigned globalDofIdx, Scalar value)
    {
        if (!this->vapparsActive(this->episodeIndex()))
            return;

        this->maxOilSaturation_[globalDofIdx] = value;
    }

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
     * For quite technical reasons, the storage term cannot be recycled if either DRSDT
     * or DRVDT are active in ebos. Nor if the porosity is changes between timesteps
     * using a pore volume multiplier (i.e., poreVolumeMultiplier() != 1.0)
     */
    bool recycleFirstIterationStorage() const
    {
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

        SolventModule::assignPrimaryVars(values,
                                         enableSolvent ? this->solventSaturation_[globalDofIdx] : 0.0,
                                         enableSolvent ? this->solventRsw_[globalDofIdx] : 0.0);

        if constexpr (enablePolymer)
            values[Indices::polymerConcentrationIdx] = this->polymer_.concentration[globalDofIdx];

        if constexpr (enablePolymerMolarWeight)
            values[Indices::polymerMoleWeightIdx]= this->polymer_.moleWeight[globalDofIdx];

        if constexpr (enableBrine) {
            if (enableSaltPrecipitation && values.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp) {
                values[Indices::saltConcentrationIdx] = initialFluidStates_[globalDofIdx].saltSaturation();
            }
            else {
                values[Indices::saltConcentrationIdx] = initialFluidStates_[globalDofIdx].saltConcentration();
            }
        }

        if constexpr (enableMICP){
            values[Indices::microbialConcentrationIdx] = this->micp_.microbialConcentration[globalDofIdx];
            values[Indices::oxygenConcentrationIdx]= this->micp_.oxygenConcentration[globalDofIdx];
            values[Indices::ureaConcentrationIdx]= this->micp_.ureaConcentration[globalDofIdx];
            values[Indices::calciteConcentrationIdx]= this->micp_.calciteConcentration[globalDofIdx];
            values[Indices::biofilmConcentrationIdx]= this->micp_.biofilmConcentration[globalDofIdx];
        }

        values.checkDefined();
    }

    /*!
     * \copydoc FvBaseProblem::initialSolutionApplied()
     */
    void initialSolutionApplied()
    {
        // Calculate all intensive quantities.
        this->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx*/0);

        // We also need the intensive quantities for timeIdx == 1
        // corresponding to the start of the current timestep, if we
        // do not use the storage cache, or if we cannot recycle the
        // first iteration storage.
        if (!this->model().enableStorageCache() || !this->recycleFirstIterationStorage()) {
            this->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx*/1);
        }

        // initialize the wells. Note that this needs to be done after initializing the
        // intrinsic permeabilities and the after applying the initial solution because
        // the well model uses these...
        wellModel_.init();

        // let the object for threshold pressures initialize itself. this is done only at
        // this point, because determining the threshold pressures may require to access
        // the initial solution.
        thresholdPressures_.finishInit();

        updateCompositionChangeLimits_();

        if (enableAquifers_)
            aquiferModel_.initialSolutionApplied();

        if (this->simulator().episodeIndex() == 0) {
            eclWriter_->writeInitialFIPReport();
        }
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0 everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context,
                unsigned spaceIdx,
                unsigned timeIdx) const
    {
        const unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        source(rate, globalDofIdx, timeIdx);
    }

    void source(RateVector& rate,
                unsigned globalDofIdx,
                unsigned timeIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(eclProblemSource);
        rate = 0.0;

        // Add well contribution to source here.
        wellModel_.computeTotalRatesForDof(rate, globalDofIdx);

        // convert the source term from the total mass rate of the
        // cell to the one per unit of volume as used by the model.
        for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
            rate[eqIdx] /= this->model().dofTotalVolume(globalDofIdx);

            Valgrind::CheckDefined(rate[eqIdx]);
            assert(isfinite(rate[eqIdx]));
        }

        // Add non-well sources.
        addToSourceDense(rate, globalDofIdx, timeIdx);
    }

    void addToSourceDense(RateVector& rate,
                          unsigned globalDofIdx,
                          unsigned timeIdx) const
    {
        if (enableAquifers_)
            aquiferModel_.addToSource(rate, globalDofIdx, timeIdx);

        // Add source term from deck
        const auto& source = this->simulator().vanguard().schedule()[this->episodeIndex()].source();
        if (source.size() > 0) {
            std::array<int,3> ijk;
                this->simulator().vanguard().cartesianCoordinate(globalDofIdx, ijk);
            RateVector massRate(0.0);
            if ( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                massRate[Indices::canonicalToActiveComponentIndex(oilCompIdx)] = source.rate({ijk, SourceComponent::OIL}) / this->model().dofTotalVolume(globalDofIdx);
            }
            if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                massRate[Indices::canonicalToActiveComponentIndex(gasCompIdx)] = source.rate({ijk, SourceComponent::GAS}) / this->model().dofTotalVolume(globalDofIdx);
            }
            if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                massRate[Indices::canonicalToActiveComponentIndex(waterCompIdx)] = source.rate({ijk, SourceComponent::WATER}) / this->model().dofTotalVolume(globalDofIdx);
            }
            if constexpr (enableSolvent) {
                massRate[Indices::solventSaturationIdx] = source.rate({ijk, SourceComponent::SOLVENT}) / this->model().dofTotalVolume(globalDofIdx);
            }
            const int pvtRegionIdx = this->pvtRegionIndex(globalDofIdx);
            rate.setMassRate(massRate, pvtRegionIdx);
            if constexpr (enablePolymer) {
                rate[Indices::polymerConcentrationIdx] = source.rate({ijk, SourceComponent::POLYMER}) / this->model().dofTotalVolume(globalDofIdx);
            }
            if constexpr (enableEnergy) {
                rate[Indices::contiEnergyEqIdx] = source.hrate(ijk) / this->model().dofTotalVolume(globalDofIdx);
            }
        }

        // if requested, compensate systematic mass loss for cells which were "well
        // behaved" in the last time step
        // Note that we don't allow for drift compensation if there are no active wells.
        const bool compensateDrift = wellModel_.wellsActive();
        if (enableDriftCompensation_ && compensateDrift) {
            const auto& simulator = this->simulator();
            const auto& model = this->model();

            // we use a lower tolerance for the compensation too
            // assure the added drift from the last step does not
            // cause convergence issues on the current step
            Scalar maxCompensation = model.newtonMethod().tolerance()/10;
            Scalar poro = this->porosity(globalDofIdx, timeIdx);
            Scalar dt = simulator.timeStepSize();
            EqVector dofDriftRate = drift_[globalDofIdx];
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
     * \brief Returns a reference to the ECL well manager used by the problem.
     *
     * This can be used for inspecting wells outside of the problem.
     */
    const EclWellModel& wellModel() const
    { return wellModel_; }

    EclWellModel& wellModel()
    { return wellModel_; }

    const EclAquiferModel& aquiferModel() const
    { return aquiferModel_; }

    EclAquiferModel& mutableAquiferModel()
    { return aquiferModel_; }

    // temporary solution to facilitate output of initial state from flow
    const InitialFluidState& initialFluidState(unsigned globalDofIdx) const
    { return initialFluidStates_[globalDofIdx]; }

    const EclipseIO& eclIO() const
    { return eclWriter_->eclIO(); }

    void setSubStepReport(const SimulatorReportSingle& report)
    { return eclWriter_->setSubStepReport(report); }

    void setSimulationReport(const SimulatorReport& report)
    { return eclWriter_->setSimulationReport(report); }

    bool nonTrivialBoundaryConditions() const
    { return nonTrivialBoundaryConditions_; }

    const InitialFluidState boundaryFluidState(unsigned globalDofIdx, const int directionId) const
    {
        OPM_TIMEBLOCK_LOCAL(boundaryFluidState);
        const auto& bcprop = this->simulator().vanguard().schedule()[this->episodeIndex()].bcprop;
        if (bcprop.size() > 0) {
            FaceDir::DirEnum dir = FaceDir::FromIntersectionIndex(directionId);

            // index == 0: no boundary conditions for this
            // global cell and direction
            if (bcindex_(dir)[globalDofIdx] == 0)
                return initialFluidStates_[globalDofIdx];

            const auto& bc = bcprop[bcindex_(dir)[globalDofIdx]];
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
                    case BCComponent::NONE:
                        throw std::logic_error("you need to specify a valid component (OIL, WATER or GAS) when DIRICHLET type is set in BC");
                        break;
                }
                double pressure = initialFluidStates_[globalDofIdx].pressure(refPressurePhaseIdx_());
                const auto pressure_input = bc.pressure;
                if (pressure_input) {
                    pressure = *pressure_input;
                }

                std::array<Scalar, numPhases> pc = {0};
                const auto& matParams = materialLawParams(globalDofIdx);
                MaterialLaw::capillaryPressures(pc, matParams, fluidState);
                Valgrind::CheckDefined(pressure);
                Valgrind::CheckDefined(pc);
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx))
                        continue;

                    if (Indices::oilEnabled)
                        fluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
                    else if (Indices::gasEnabled)
                        fluidState.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[gasPhaseIdx]));
                    else if (Indices::waterEnabled)
                        //single (water) phase
                        fluidState.setPressure(phaseIdx, pressure);
                }
                
                double temperature = initialFluidStates_[globalDofIdx].temperature(0); // we only have one temperature
                const auto temperature_input = bc.temperature;
                if(temperature_input)
                    temperature = *temperature_input;
                fluidState.setTemperature(temperature);

                if (FluidSystem::enableDissolvedGas()) {
                    fluidState.setRs(0.0);
                    fluidState.setRv(0.0);
                }
                if (FluidSystem::enableDissolvedGasInWater()) {
                    fluidState.setRsw(0.0);
                }
                if (FluidSystem::enableVaporizedWater())
                    fluidState.setRvw(0.0);

                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx))
                        continue;

                    const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState, phaseIdx, pvtRegionIdx);
                    fluidState.setInvB(phaseIdx, b);

                    const auto& rho = FluidSystem::density(fluidState, phaseIdx, pvtRegionIdx);
                    fluidState.setDensity(phaseIdx, rho);

                }
                fluidState.checkDefined();
                return fluidState;
            }
        }
        return initialFluidStates_[globalDofIdx];
    }

    /*!
     * \brief Propose the size of the next time step to the simulator.
     *
     * This method is only called if the Newton solver does converge, the simulator
     * automatically cuts the time step in half without consultating this method again.
     */
    Scalar nextTimeStepSize() const
    {
        OPM_TIMEBLOCK(nexTimeStepSize);
        // allow external code to do the timestepping
        if (this->nextTimeStepSize_ > 0.0)
            return this->nextTimeStepSize_;

        const auto& simulator = this->simulator();
        int episodeIdx = simulator.episodeIndex();

        // for the initial episode, we use a fixed time step size
        if (episodeIdx < 0)
            return this->initialTimeStepSize_;

        // ask the newton method for a suggestion. This suggestion will be based on how
        // well the previous time step converged. After that, apply the runtime time
        // stepping constraints.
        const auto& newtonMethod = this->model().newtonMethod();
        return limitNextTimeStepSize_(newtonMethod.suggestTimeStepSize(simulator.timeStepSize()));
    }

    /*!
     * \brief Calculate the porosity multiplier due to water induced rock compaction.
     *
     * TODO: The API of this is a bit ad-hoc, it would be better to use context objects.
     */
    template <class LhsEval>
    LhsEval rockCompPoroMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(rockCompPoroMultiplier);
        if (this->rockCompPoroMult_.empty() && this->rockCompPoroMultWc_.empty())
            return 1.0;

        unsigned tableIdx = 0;
        if (!this->rockTableIdx_.empty())
            tableIdx = this->rockTableIdx_[elementIdx];

        const auto& fs = intQuants.fluidState();
        LhsEval effectivePressure = decay<LhsEval>(fs.pressure(refPressurePhaseIdx_()));
        if (!this->minRefPressure_.empty())
            // The pore space change is irreversible
            effectivePressure =
                min(decay<LhsEval>(fs.pressure(refPressurePhaseIdx_())),
                                   this->minRefPressure_[elementIdx]);

        if (!this->overburdenPressure_.empty())
            effectivePressure -= this->overburdenPressure_[elementIdx];


        if (!this->rockCompPoroMult_.empty()) {
            return this->rockCompPoroMult_[tableIdx].eval(effectivePressure, /*extrapolation=*/true);
        }

        // water compaction
        assert(!this->rockCompPoroMultWc_.empty());
        LhsEval SwMax = max(decay<LhsEval>(fs.saturation(waterPhaseIdx)), this->maxWaterSaturation_[elementIdx]);
        LhsEval SwDeltaMax = SwMax - initialFluidStates_[elementIdx].saturation(waterPhaseIdx);

        return this->rockCompPoroMultWc_[tableIdx].eval(effectivePressure, SwDeltaMax, /*extrapolation=*/true);
    }

    /*!
     * \brief Calculate the transmissibility multiplier due to water induced rock compaction.
     *
     * TODO: The API of this is a bit ad-hoc, it would be better to use context objects.
     */
    template <class LhsEval>
    LhsEval rockCompTransMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(rockCompTransMultiplier);
        if (this->rockCompTransMult_.empty() && this->rockCompTransMultWc_.empty())
            return 1.0;

        unsigned tableIdx = 0;
        if (!this->rockTableIdx_.empty())
            tableIdx = this->rockTableIdx_[elementIdx];

        const auto& fs = intQuants.fluidState();
        LhsEval effectivePressure = decay<LhsEval>(fs.pressure(refPressurePhaseIdx_()));

        if (!this->minRefPressure_.empty())
            // The pore space change is irreversible
            effectivePressure =
                min(decay<LhsEval>(fs.pressure(refPressurePhaseIdx_())),
                    this->minRefPressure_[elementIdx]);

        if (!this->overburdenPressure_.empty())
            effectivePressure -= this->overburdenPressure_[elementIdx];

        if (!this->rockCompTransMult_.empty())
            return this->rockCompTransMult_[tableIdx].eval(effectivePressure, /*extrapolation=*/true);

        // water compaction
        assert(!this->rockCompTransMultWc_.empty());
        LhsEval SwMax = max(decay<LhsEval>(fs.saturation(waterPhaseIdx)), this->maxWaterSaturation_[elementIdx]);
        LhsEval SwDeltaMax = SwMax - initialFluidStates_[elementIdx].saturation(waterPhaseIdx);

        return this->rockCompTransMultWc_[tableIdx].eval(effectivePressure, SwDeltaMax, /*extrapolation=*/true);
    }

    /*!
     * \brief Calculate the transmissibility multiplier due to porosity reduction.
     *
     * TODO: The API of this is a bit ad-hoc, it would be better to use context objects.
     */
    template <class LhsEval>
    LhsEval permFactTransMultiplier(const IntensiveQuantities& intQuants) const
    {
        OPM_TIMEBLOCK_LOCAL(permFactTransMultiplier);
        if (!enableSaltPrecipitation)
            return 1.0;
        
        const auto& fs = intQuants.fluidState();
        unsigned tableIdx = fs.pvtRegionIndex();
        LhsEval porosityFactor = decay<LhsEval>(1. - fs.saltSaturation());
        porosityFactor = min(porosityFactor, 1.0);
        const auto& permfactTable = BrineModule::permfactTable(tableIdx);
        return permfactTable.eval(porosityFactor, /*extrapolation=*/true);
    }

    /*!
     * \brief Return the well transmissibility multiplier due to rock changues.
     */
    template <class LhsEval>
    LhsEval wellTransMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(wellTransMultiplier);
        
        double trans_mult = this->simulator().problem().template rockCompTransMultiplier<double>(intQuants, elementIdx);
        trans_mult *= this->simulator().problem().template permFactTransMultiplier<double>(intQuants);
    
        return trans_mult;
    }

    std::pair<BCType, RateVector> boundaryCondition(const unsigned int globalSpaceIdx, const int directionId) const
    {
        OPM_TIMEBLOCK_LOCAL(boundaryCondition);
        if (!nonTrivialBoundaryConditions_) {
            return { BCType::NONE, RateVector(0.0) };
        }
        FaceDir::DirEnum dir = FaceDir::FromIntersectionIndex(directionId);
        const auto& schedule = this->simulator().vanguard().schedule();
        if (bcindex_(dir)[globalSpaceIdx] == 0) {
            return { BCType::NONE, RateVector(0.0) };
        }
        if (schedule[this->episodeIndex()].bcprop.size() == 0) {
            return { BCType::NONE, RateVector(0.0) };
        }
        const auto& bc = schedule[this->episodeIndex()].bcprop[bcindex_(dir)[globalSpaceIdx]];
        if (bc.bctype!=BCType::RATE) {
            return { bc.bctype, RateVector(0.0) };
        }

        RateVector rate = 0.0;
        switch (bc.component) {
        case BCComponent::OIL:
            rate[Indices::canonicalToActiveComponentIndex(oilCompIdx)] = bc.rate;
            break;
        case BCComponent::GAS:
            rate[Indices::canonicalToActiveComponentIndex(gasCompIdx)] = bc.rate;
            break;
        case BCComponent::WATER:
            rate[Indices::canonicalToActiveComponentIndex(waterCompIdx)] = bc.rate;
            break;
        case BCComponent::SOLVENT:
            if constexpr (!enableSolvent)
                throw std::logic_error("solvent is disabled and you're trying to add solvent to BC");

            rate[Indices::solventSaturationIdx] = bc.rate;
            break;
        case BCComponent::POLYMER:
            if constexpr (!enablePolymer)
                throw std::logic_error("polymer is disabled and you're trying to add polymer to BC");

            rate[Indices::polymerConcentrationIdx] = bc.rate;
            break;
        case BCComponent::NONE:
            throw std::logic_error("you need to specify the component when RATE type is set in BC");
            break;
        }
        //TODO add support for enthalpy rate
        return {bc.bctype, rate};
    }

    const std::unique_ptr<EclWriterType>& eclWriter() const
    {
        return eclWriter_;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<BaseType&>(*this));
        serializer(drift_);
        serializer(wellModel_);
        serializer(aquiferModel_);
        serializer(tracerModel_);
        serializer(*materialLawManager_);
        serializer(*eclWriter_);
    }
private:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }
protected:
    void updateExplicitQuantities_()
    {
        OPM_TIMEBLOCK(updateExplicitQuantities);
        const bool invalidateFromMaxWaterSat = updateMaxWaterSaturation_();
        const bool invalidateFromMinPressure = updateMinPressure_();

        // update hysteresis and max oil saturation used in vappars
        const bool invalidateFromHyst = updateHysteresis_();
        const bool invalidateFromMaxOilSat = updateMaxOilSaturation_();

        // the derivatives may have change
        bool invalidateIntensiveQuantities
            = invalidateFromMaxWaterSat || invalidateFromMinPressure || invalidateFromHyst || invalidateFromMaxOilSat;
        if (invalidateIntensiveQuantities) {
            OPM_TIMEBLOCK(beginTimeStepInvalidateIntensiveQuantities);
            this->model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        }

        if constexpr (getPropValue<TypeTag, Properties::EnablePolymer>())
            updateMaxPolymerAdsorption_();
    }

    template<class UpdateFunc>
    void updateProperty_(const std::string& failureMsg,
                         UpdateFunc func)
    {
        OPM_TIMEBLOCK(updateProperty);
        const auto& model = this->simulator().model();
        const auto& primaryVars = model.solution(/*timeIdx*/0);
        const auto& vanguard = this->simulator().vanguard();
        std::size_t numGridDof = primaryVars.size();
        OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
                const auto& iq = *model.cachedIntensiveQuantities(dofIdx, /*timeIdx=*/ 0);
                func(dofIdx, iq);
        }
        OPM_END_PARALLEL_TRY_CATCH(failureMsg, vanguard.grid().comm());
    }

    // update the parameters needed for DRSDT and DRVDT
    void updateCompositionChangeLimits_()
    {
        OPM_TIMEBLOCK(updateCompositionChangeLimits);
        // update the "last Rs" values for all elements, including the ones in the ghost
        // and overlap regions
        int episodeIdx = this->episodeIndex();
        std::array<bool,3> active{this->mixControls_.drsdtConvective(episodeIdx),
                                  this->mixControls_.drsdtActive(episodeIdx),
                                  this->mixControls_.drvdtActive(episodeIdx)};
        if (!active[0] && !active[1] && !active[2]) {
            return;
        }

        this->updateProperty_("EclProblem::updateCompositionChangeLimits_()) failed:",
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
                                                            pvtRegionIdx,
                                                            active);
                              }
            );
    }

    bool updateMaxOilSaturation_()
    {
        OPM_TIMEBLOCK(updateMaxOilSaturation);
        int episodeIdx = this->episodeIndex();

        // we use VAPPARS
        if (this->vapparsActive(episodeIdx)) {
            this->updateProperty_("EclProblem::updateMaxOilSaturation_() failed:",
                                  [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                                  {
                                      this->updateMaxOilSaturation_(compressedDofIdx,iq);
                                  });
            return true;
        }

        return false;
    }

    bool updateMaxOilSaturation_(unsigned compressedDofIdx, const IntensiveQuantities& iq)
    {
        OPM_TIMEBLOCK_LOCAL(updateMaxOilSaturation);
        const auto& fs = iq.fluidState();
        const Scalar So = decay<Scalar>(fs.saturation(refPressurePhaseIdx_()));
        auto& mos = this->maxOilSaturation_;
        if(mos[compressedDofIdx] < So){
            mos[compressedDofIdx] = So;
            return true;
        }else{
            return false;
        }
    }

    bool updateMaxWaterSaturation_()
    {
        OPM_TIMEBLOCK(updateMaxWaterSaturation);
        // water compaction is activated in ROCKCOMP
        if (this->maxWaterSaturation_.empty())
            return false;

        this->maxWaterSaturation_[/*timeIdx=*/1] = this->maxWaterSaturation_[/*timeIdx=*/0];
        this->updateProperty_("EclProblem::updateMaxWaterSaturation_() failed:",
                              [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                              {
                                  this->updateMaxWaterSaturation_(compressedDofIdx,iq);
                               });
        return true;
    }


    bool updateMaxWaterSaturation_(unsigned compressedDofIdx, const IntensiveQuantities& iq)
    {
        OPM_TIMEBLOCK_LOCAL(updateMaxWaterSaturation);
        const auto& fs = iq.fluidState();
        const Scalar Sw = decay<Scalar>(fs.saturation(waterPhaseIdx));
        auto& mow = this->maxWaterSaturation_;
        if(mow[compressedDofIdx]< Sw){
            mow[compressedDofIdx] = Sw;
            return true;
        }else{
            return false;
        }
    }

    bool updateMinPressure_()
    {
        OPM_TIMEBLOCK(updateMinPressure);
        // IRREVERS option is used in ROCKCOMP
        if (this->minRefPressure_.empty())
            return false;

        this->updateProperty_("EclProblem::updateMinPressure_() failed:",
                              [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                              {
                                  this->updateMinPressure_(compressedDofIdx,iq);
                              });
        return true;
    }

    bool updateMinPressure_(unsigned compressedDofIdx, const IntensiveQuantities& iq){
        OPM_TIMEBLOCK_LOCAL(updateMinPressure);
        const auto& fs = iq.fluidState();
        const Scalar min_pressure = getValue(fs.pressure(refPressurePhaseIdx_()));
        auto& min_pressures = this->minRefPressure_;
        if(min_pressures[compressedDofIdx]> min_pressure){
            min_pressures[compressedDofIdx] = min_pressure;
            return true;
        }else{
            return false;
        }
    }

    // \brief Function to assign field properties of type double, on the leaf grid view.
    //
    // For CpGrid with local grid refinement, the field property of a cell on the leaf
    // is inherited from its parent or equivalent (when has no parent) cell on level zero.
    std::function<std::vector<double>(const FieldPropsManager&, const std::string&, const unsigned int&)>
    fieldPropDoubleOnLeafAssigner_()
    {
        const auto& lookup = this->lookUpData_;
        return [&lookup](const FieldPropsManager& fieldPropManager, const std::string& propString,
                         const unsigned int& numElems)
        {
            return lookup.assignFieldPropsDoubleOnLeaf(fieldPropManager, propString, numElems);
        };
    }

    // \brief Function to assign field properties of type int, unsigned int, ..., on the leaf grid view.
    //
    // For CpGrid with local grid refinement, the field property of a cell on the leaf
    // is inherited from its parent or equivalent (when has no parent) cell on level zero.
    template<typename IntType>
    std::function<std::vector<IntType>(const FieldPropsManager&, const std::string&, const unsigned int&, bool)>
    fieldPropIntTypeOnLeafAssigner_()
    {
        const auto& lookup = this->lookUpData_;
        return [&lookup](const FieldPropsManager& fieldPropManager, const std::string& propString,
                         const unsigned int& numElems, bool needsTranslation)
        {
            return lookup.template assignFieldPropsIntOnLeaf<IntType>(fieldPropManager, propString,
                                                                      numElems, needsTranslation);
        };
    }

    void readMaterialParameters_()
    {
        OPM_TIMEBLOCK(readMaterialParameters);
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        // the PVT and saturation region numbers
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        this->updatePvtnum_();
        this->updateSatnum_();

        // the MISC region numbers (solvent model)
        this->updateMiscnum_();
        // the PLMIX region numbers (polymer model)
        this->updatePlmixnum_();

        // directional relative permeabilities
        this->updateKrnum_();
        OPM_END_PARALLEL_TRY_CATCH("Invalid region numbers: ", vanguard.gridView().comm());
        ////////////////////////////////
        // porosity
        updateReferencePorosity_();
        this->referencePorosity_[1] = this->referencePorosity_[0];
        ////////////////////////////////

        ////////////////////////////////
        // fluid-matrix interactions (saturation functions; relperm/capillary pressure)
        materialLawManager_ = std::make_shared<EclMaterialLawManager>();
        materialLawManager_->initFromState(eclState);
        materialLawManager_->initParamsForElements(eclState, this->model().numGridDof(),
                                                   this-> template fieldPropIntTypeOnLeafAssigner_<int>(),
                                                   this-> lookupIdxOnLevelZeroAssigner_());
        ////////////////////////////////
    }

    void readThermalParameters_()
    {
        if constexpr (enableEnergy)
        {
            const auto& simulator = this->simulator();
            const auto& vanguard = simulator.vanguard();
            const auto& eclState = vanguard.eclState();

            // fluid-matrix interactions (saturation functions; relperm/capillary pressure)
            thermalLawManager_ = std::make_shared<EclThermalLawManager>();
            thermalLawManager_->initParamsForElements(eclState, this->model().numGridDof(),
                                                      this-> fieldPropDoubleOnLeafAssigner_(),
                                                      this-> template fieldPropIntTypeOnLeafAssigner_<unsigned int>());
        }
    }

    void updateReferencePorosity_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        std::size_t numDof = this->model().numGridDof();

        this->referencePorosity_[/*timeIdx=*/0].resize(numDof);

        const auto& fp = eclState.fieldProps();
        const std::vector<double> porvData = fp.porv(false);
        const std::vector<int> actnumData = fp.actnum();
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            Scalar poreVolume = porvData[dofIdx];

            // we define the porosity as the accumulated pore volume divided by the
            // geometric volume of the element. Note that -- in pathetic cases -- it can
            // be larger than 1.0!
            Scalar dofVolume = simulator.model().dofTotalVolume(dofIdx);
            assert(dofVolume > 0.0);
            this->referencePorosity_[/*timeIdx=*/0][dofIdx] = poreVolume/dofVolume;
        }
    }

    void readInitialCondition_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        if (eclState.getInitConfig().hasEquil())
            readEquilInitialCondition_();
        else
            readExplicitInitialCondition_();

        if constexpr (enableSolvent || enablePolymer || enablePolymerMolarWeight || enableMICP)
            this->readBlackoilExtentionsInitialConditions_(this->model().numGridDof(),
                                                           enableSolvent,
                                                           enablePolymer,
                                                           enablePolymerMolarWeight,
                                                           enableMICP);

        //initialize min/max values
        std::size_t numElems = this->model().numGridDof();
        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            const auto& fs = initialFluidStates_[elemIdx];
            if (!this->maxWaterSaturation_.empty())
                this->maxWaterSaturation_[elemIdx] = std::max(this->maxWaterSaturation_[elemIdx], fs.saturation(waterPhaseIdx));
            if (!this->maxOilSaturation_.empty())
                this->maxOilSaturation_[elemIdx] = std::max(this->maxOilSaturation_[elemIdx], fs.saturation(oilPhaseIdx));
            if (!this->minRefPressure_.empty())
                this->minRefPressure_[elemIdx] = std::min(this->minRefPressure_[elemIdx], fs.pressure(refPressurePhaseIdx_()));
        }


    }

    void readEquilInitialCondition_()
    {
        const auto& simulator = this->simulator();

        // initial condition corresponds to hydrostatic conditions.
        using EquilInitializer = EclEquilInitializer<TypeTag>;
        EquilInitializer equilInitializer(simulator, *materialLawManager_);

        std::size_t numElems = this->model().numGridDof();
        initialFluidStates_.resize(numElems);
        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = initialFluidStates_[elemIdx];
            elemFluidState.assign(equilInitializer.initialFluidState(elemIdx));
        }
    }

    void readEclRestartSolution_()
    {
        // Throw an exception if the grid has LGRs. Refined grid are not supported for restart.
        if(this->simulator().vanguard().grid().maxLevel() > 0) {
            throw std::invalid_argument("Refined grids are not yet supported for restart ");
        }

        // Set the start time of the simulation
        auto& simulator = this->simulator();
        const auto& schedule = simulator.vanguard().schedule();
        const auto& eclState = simulator.vanguard().eclState();
        const auto& initconfig = eclState.getInitConfig();
        {
            int restart_step = initconfig.getRestartStep();

            simulator.setTime(schedule.seconds(restart_step));

            simulator.startNextEpisode(simulator.startTime() + simulator.time(),
                                       schedule.stepLength(restart_step));
            simulator.setEpisodeIndex(restart_step);
        }
        eclWriter_->beginRestart();

        Scalar dt = std::min(eclWriter_->restartTimeStepSize(), simulator.episodeLength());
        simulator.setTimeStepSize(dt);

        std::size_t numElems = this->model().numGridDof();
        initialFluidStates_.resize(numElems);
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

        if constexpr (enableMICP) {
            this->micp_.resize(numElems);
        }

        for (std::size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = initialFluidStates_[elemIdx];
            elemFluidState.setPvtRegionIndex(pvtRegionIndex(elemIdx));
            eclWriter_->eclOutputModule().initHysteresisParams(simulator, elemIdx);
            eclWriter_->eclOutputModule().assignToFluidState(elemFluidState, elemIdx);

            // Note: Function processRestartSaturations_() mutates the
            // 'ssol' argument--the value from the restart file--if solvent
            // is enabled.  Then, store the updated solvent saturation into
            // 'solventSaturation_'.  Otherwise, just pass a dummy value to
            // the function and discard the unchanged result.  Do not index
            // into 'solventSaturation_' unless solvent is enabled.
            {
                auto ssol = enableSolvent
                    ? eclWriter_->eclOutputModule().getSolventSaturation(elemIdx)
                    : Scalar(0);

                processRestartSaturations_(elemFluidState, ssol);

                if constexpr (enableSolvent) {
                    this->solventSaturation_[elemIdx] = ssol;
                    this->solventRsw_[elemIdx] = eclWriter_->eclOutputModule().getSolventRsw(elemIdx);
                }
            }

            this->mixControls_.updateLastValues(elemIdx, elemFluidState.Rs(), elemFluidState.Rv());

            if constexpr (enablePolymer)
                 this->polymer_.concentration[elemIdx] = eclWriter_->eclOutputModule().getPolymerConcentration(elemIdx);
            if constexpr (enableMICP){
                 this->micp_.microbialConcentration[elemIdx] = eclWriter_->eclOutputModule().getMicrobialConcentration(elemIdx);
                 this->micp_.oxygenConcentration[elemIdx] = eclWriter_->eclOutputModule().getOxygenConcentration(elemIdx);
                 this->micp_.ureaConcentration[elemIdx] = eclWriter_->eclOutputModule().getUreaConcentration(elemIdx);
                 this->micp_.biofilmConcentration[elemIdx] = eclWriter_->eclOutputModule().getBiofilmConcentration(elemIdx);
                 this->micp_.calciteConcentration[elemIdx] = eclWriter_->eclOutputModule().getCalciteConcentration(elemIdx);
            }
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
            initial(sol[elemIdx], elemCtx, /*spaceIdx=*/0, /*timeIdx=*/0);
        }

        // make sure that the ghost and overlap entities exhibit the correct
        // solution. alternatively, this could be done in the loop above by also
        // considering non-interior elements. Since the initial() method might not work
        // 100% correctly for such elements, let's play safe and explicitly synchronize
        // using message passing.
        this->model().syncOverlap();

        eclWriter_->endRestart();
    }

    void processRestartSaturations_(InitialFluidState& elemFluidState, Scalar& solventSaturation)
    {
        // each phase needs to be above certain value to be claimed to be existing
        // this is used to recover some RESTART running with the defaulted single-precision format
        const Scalar smallSaturationTolerance = 1.e-6;
        Scalar sumSaturation = 0.0;
        for (std::size_t phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                if (elemFluidState.saturation(phaseIdx) < smallSaturationTolerance)
                    elemFluidState.setSaturation(phaseIdx, 0.0);

                sumSaturation += elemFluidState.saturation(phaseIdx);
            }

        }
        if constexpr (enableSolvent) {
            if (solventSaturation < smallSaturationTolerance)
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

    void readExplicitInitialCondition_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& fp = eclState.fieldProps();
        bool has_swat     = fp.has_double("SWAT");
        bool has_sgas     = fp.has_double("SGAS");
        bool has_rs       = fp.has_double("RS");
        bool has_rv       = fp.has_double("RV");
        bool has_rvw       = fp.has_double("RVW");
        bool has_pressure = fp.has_double("PRESSURE");
        bool has_salt = fp.has_double("SALT");
        bool has_saltp = fp.has_double("SALTP");

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
                                     " dissolved gas is enabled");
        if (FluidSystem::enableVaporizedOil() && !has_rv)
            throw std::runtime_error("The ECL input file requires the RV keyword to be present if"
                                     " vaporized oil is enabled");
        if (FluidSystem::enableVaporizedWater() && !has_rvw)
            throw std::runtime_error("The ECL input file requires the RVW keyword to be present if"
                                     " vaporized water is enabled");
        if (enableBrine && !has_salt)
            throw std::runtime_error("The ECL input file requires the SALT keyword to be present if"
                                     " brine is enabled and the model is initialized explicitly");
        if (enableSaltPrecipitation && !has_saltp)
            throw std::runtime_error("The ECL input file requires the SALTP keyword to be present if"
                                     " salt precipitation is enabled and the model is initialized explicitly");

        std::size_t numDof = this->model().numGridDof();

        initialFluidStates_.resize(numDof);

        std::vector<double> waterSaturationData;
        std::vector<double> gasSaturationData;
        std::vector<double> pressureData;
        std::vector<double> rsData;
        std::vector<double> rvData;
        std::vector<double> rvwData;
        std::vector<double> tempiData;
        std::vector<double> saltData;
        std::vector<double> saltpData;

        if (FluidSystem::phaseIsActive(waterPhaseIdx) && Indices::numPhases > 1)
            waterSaturationData = fp.get_double("SWAT");
        else
            waterSaturationData.resize(numDof);

        if (FluidSystem::phaseIsActive(gasPhaseIdx) && FluidSystem::phaseIsActive(oilPhaseIdx))
            gasSaturationData = fp.get_double("SGAS");
        else
            gasSaturationData.resize(numDof);

        pressureData = fp.get_double("PRESSURE");
        if (FluidSystem::enableDissolvedGas())
            rsData = fp.get_double("RS");

        if (FluidSystem::enableVaporizedOil())
            rvData = fp.get_double("RV");

        if (FluidSystem::enableVaporizedWater())
            rvwData = fp.get_double("RVW");

        // initial reservoir temperature
        tempiData = fp.get_double("TEMPI");

        // initial salt concentration data
        if constexpr (enableBrine)
            saltData = fp.get_double("SALT");

         // initial precipitated salt saturation data
         if constexpr (enableSaltPrecipitation)
            saltpData = fp.get_double("SALTP");

        // calculate the initial fluid states
        for (std::size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofFluidState = initialFluidStates_[dofIdx];

            dofFluidState.setPvtRegionIndex(pvtRegionIndex(dofIdx));

            //////
            // set temperature
            //////
            Scalar temperatureLoc = tempiData[dofIdx];
            if (!std::isfinite(temperatureLoc) || temperatureLoc <= 0)
                temperatureLoc = FluidSystem::surfaceTemperature;
            dofFluidState.setTemperature(temperatureLoc);

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
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
                dofFluidState.setSaturation(FluidSystem::oilPhaseIdx,
                                            1.0
                                            - waterSaturationData[dofIdx]
                                            - gasSaturationData[dofIdx]);

            //////
            // set phase pressures
            //////
            Scalar pressure = pressureData[dofIdx]; // oil pressure (or gas pressure for water-gas system or water pressure for single phase)

            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            std::array<Scalar, numPhases> pc = {0};
            const auto& matParams = materialLawParams(dofIdx);
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

            if (FluidSystem::enableDissolvedGas())
                dofFluidState.setRs(rsData[dofIdx]);
            else if (Indices::gasEnabled && Indices::oilEnabled)
                dofFluidState.setRs(0.0);

            if (FluidSystem::enableVaporizedOil())
                dofFluidState.setRv(rvData[dofIdx]);
            else if (Indices::gasEnabled && Indices::oilEnabled)
                dofFluidState.setRv(0.0);

            if (FluidSystem::enableVaporizedWater())
                dofFluidState.setRvw(rvwData[dofIdx]);

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

    // update the hysteresis parameters of the material laws for the whole grid
    bool updateHysteresis_()
    {
        if (!materialLawManager_->enableHysteresis())
            return false;

        // we need to update the hysteresis data for _all_ elements (i.e., not just the
        // interior ones) to avoid desynchronization of the processes in the parallel case!
        this->updateProperty_("EclProblem::updateHysteresis_() failed:",
                              [this](unsigned compressedDofIdx, const IntensiveQuantities& iq)
                              {
                                  materialLawManager_->updateHysteresis(iq.fluidState(), compressedDofIdx);
                              });
        return true;
    }


    bool updateHysteresis_(unsigned compressedDofIdx, const IntensiveQuantities& iq)
    {
        OPM_TIMEBLOCK_LOCAL(updateHysteresis_);
        materialLawManager_->updateHysteresis(iq.fluidState(), compressedDofIdx);
        //TODO change materials to give a bool
        return true;
    }

    void updateMaxPolymerAdsorption_()
    {
        // we need to update the max polymer adsoption data for all elements
        this->updateProperty_("EclProblem::updateMaxPolymerAdsorption_() failed:",
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

private:
    struct PffDofData_
    {
        ConditionalStorage<enableEnergy, Scalar> thermalHalfTransIn;
        ConditionalStorage<enableEnergy, Scalar> thermalHalfTransOut;
        ConditionalStorage<enableDiffusion, Scalar> diffusivity;
        ConditionalStorage<enableDispersion, Scalar> dispersivity;
        Scalar transmissibility;
    };

    // update the prefetch friendly data object
    void updatePffDofData_()
    {
        const auto& distFn =
            [this](PffDofData_& dofData,
                   const Stencil& stencil,
                   unsigned localDofIdx)
            -> void
        {
            const auto& elementMapper = this->model().elementMapper();

            unsigned globalElemIdx = elementMapper.index(stencil.entity(localDofIdx));
            if (localDofIdx != 0) {
                unsigned globalCenterElemIdx = elementMapper.index(stencil.entity(/*dofIdx=*/0));
                dofData.transmissibility = transmissibilities_.transmissibility(globalCenterElemIdx, globalElemIdx);

                if constexpr (enableEnergy) {
                    *dofData.thermalHalfTransIn = transmissibilities_.thermalHalfTrans(globalCenterElemIdx, globalElemIdx);
                    *dofData.thermalHalfTransOut = transmissibilities_.thermalHalfTrans(globalElemIdx, globalCenterElemIdx);
                }
                if constexpr (enableDiffusion)
                    *dofData.diffusivity = transmissibilities_.diffusivity(globalCenterElemIdx, globalElemIdx);
                if (enableDispersion)
                    dofData.dispersivity = transmissibilities_.dispersivity(globalCenterElemIdx, globalElemIdx);
            }
        };

        pffDofData_.update(distFn);
    }

    void readBoundaryConditions_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& bcconfig = vanguard.eclState().getSimulationConfig().bcconfig();
        if (bcconfig.size() > 0) {
            nonTrivialBoundaryConditions_ = true;

            std::size_t numCartDof = vanguard.cartesianSize();
            unsigned numElems = vanguard.gridView().size(/*codim=*/0);
            std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);

            for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx)
                cartesianToCompressedElemIdx[vanguard.cartesianIndex(elemIdx)] = elemIdx;

            bcindex_.resize(numElems, 0);
            auto loopAndApply = [&cartesianToCompressedElemIdx,
                                 &vanguard](const auto& bcface,
                                            auto apply)
            {
                for (int i = bcface.i1; i <= bcface.i2; ++i) {
                    for (int j = bcface.j1; j <= bcface.j2; ++j) {
                        for (int k = bcface.k1; k <= bcface.k2; ++k) {
                            std::array<int, 3> tmp = {i,j,k};
                            auto elemIdx = cartesianToCompressedElemIdx[vanguard.cartesianIndex(tmp)];
                            if (elemIdx >= 0)
                                apply(elemIdx);
                        }
                    }
                }
            };
            for (const auto& bcface : bcconfig) {
                std::vector<int>& data = bcindex_(bcface.dir);
                const int index = bcface.index;
                    loopAndApply(bcface,
                                 [&data,index](int elemIdx)
                                 { data[elemIdx] = index; });
            }
        }
    }

    // this method applies the runtime constraints specified via the deck and/or command
    // line parameters for the size of the next time step.
    Scalar limitNextTimeStepSize_(Scalar dtNext) const
    {
        if constexpr (enableExperiments) {
            const auto& simulator = this->simulator();
            const auto& schedule = simulator.vanguard().schedule();
            int episodeIdx = simulator.episodeIndex();

            // first thing in the morning, limit the time step size to the maximum size
            Scalar maxTimeStepSize = EWOMS_GET_PARAM(TypeTag, double, SolverMaxTimeStepInDays)*24*60*60;
            int reportStepIdx = std::max(episodeIdx, 0);
            if (this->enableTuning_) {
                const auto& tuning = schedule[reportStepIdx].tuning();
                maxTimeStepSize = tuning.TSMAXZ;
            }

            dtNext = std::min(dtNext, maxTimeStepSize);

            Scalar remainingEpisodeTime =
                simulator.episodeStartTime() + simulator.episodeLength()
                - (simulator.startTime() + simulator.time());
            assert(remainingEpisodeTime >= 0.0);

            // if we would have a small amount of time left over in the current episode, make
            // two equal time steps instead of a big and a small one
            if (remainingEpisodeTime/2.0 < dtNext && dtNext < remainingEpisodeTime*(1.0 - 1e-5))
                // note: limiting to the maximum time step size here is probably not strictly
                // necessary, but it should not hurt and is more fool-proof
                dtNext = std::min(maxTimeStepSize, remainingEpisodeTime/2.0);

            if (simulator.episodeStarts()) {
                // if a well event occurred, respect the limit for the maximum time step after
                // that, too
                const auto& events = simulator.vanguard().schedule()[reportStepIdx].events();
                bool wellEventOccured =
                        events.hasEvent(ScheduleEvents::NEW_WELL)
                        || events.hasEvent(ScheduleEvents::PRODUCTION_UPDATE)
                        || events.hasEvent(ScheduleEvents::INJECTION_UPDATE)
                        || events.hasEvent(ScheduleEvents::WELL_STATUS_CHANGE);
                if (episodeIdx >= 0 && wellEventOccured && this->maxTimeStepAfterWellEvent_ > 0)
                    dtNext = std::min(dtNext, this->maxTimeStepAfterWellEvent_);
            }
        }

        return dtNext;
    }

    void computeAndSetEqWeights_()
    {
        std::vector<Scalar> sumInvB(numPhases, 0.0);
        const auto& gridView = this->gridView();
        ElementContext elemCtx(this->simulator());
        for(const auto& elem: elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updatePrimaryStencil(elem);
            int elemIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& dofFluidState = initialFluidStates_[elemIdx];
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
            unsigned solventCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
            unsigned activeSolventCompIdx = Indices::canonicalToActiveComponentIndex(solventCompIdx);
            this->model().setEqWeight(activeSolventCompIdx, avgB);
        }
    }

    int refPressurePhaseIdx_() const {
        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            return oilPhaseIdx;
        }
        else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            return gasPhaseIdx;
        }
        else {
            return waterPhaseIdx;
        }
    }

    typename Vanguard::TransmissibilityType transmissibilities_;

    std::shared_ptr<EclMaterialLawManager> materialLawManager_;
    std::shared_ptr<EclThermalLawManager> thermalLawManager_;

    EclThresholdPressure<TypeTag> thresholdPressures_;

    std::vector<InitialFluidState> initialFluidStates_;

    bool enableDriftCompensation_;
    GlobalEqVector drift_;

    EclWellModel wellModel_;
    bool enableAquifers_;
    EclAquiferModel aquiferModel_;

    bool enableEclOutput_;
    std::unique_ptr<EclWriterType> eclWriter_;
    
#if HAVE_DAMARIS
    bool enableDamarisOutput_ = false ;
    std::unique_ptr<DamarisWriterType> damarisWriter_;
#endif 

    PffGridVector<GridView, Stencil, PffDofData_, DofMapper> pffDofData_;
    TracerModel tracerModel_;

    EclActionHandler actionHandler_;

    template<class T>
    struct BCData
    {
        std::array<std::vector<T>,6> data;

        void resize(std::size_t size, T defVal)
        {
            for (auto& d : data)
                d.resize(size, defVal);
        }

        const std::vector<T>& operator()(FaceDir::DirEnum dir) const
        {
            if (dir == FaceDir::DirEnum::Unknown)
                throw std::runtime_error("Tried to access BC data for the 'Unknown' direction");
            int idx = 0;
            int div = static_cast<int>(dir);
            while ((div /= 2) >= 1)
              ++idx;
            assert(idx >= 0 && idx <= 5);
            return data[idx];
        }

        std::vector<T>& operator()(FaceDir::DirEnum dir)
        {
            return const_cast<std::vector<T>&>(std::as_const(*this)(dir));
        }
    };

    BCData<int> bcindex_;
    bool nonTrivialBoundaryConditions_ = false;
};

} // namespace Opm

#endif
