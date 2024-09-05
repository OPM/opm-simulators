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

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WetGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DeadOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityWaterPvt.hpp>

#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>

//#include <opm/output/eclipse/EclipseIO.hpp>

#include <opm/simulators/flow/FlowProblem.hpp>

#include <opm/simulators/flow/MixingRateControls.hpp>

#include <opm/simulators/flow/ActionHandler.hpp>
#if HAVE_DAMARIS
#include <opm/simulators/flow/DamarisWriter.hpp>
#endif

#include <algorithm>
#include <functional>
#include <set>
#include <string>
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
    // TODO: the naming of the Types will be adjusted
    using FlowProblemType = FlowProblem<TypeTag>;
    friend FlowProblemType;

    using typename FlowProblemType::Scalar;
    using typename FlowProblemType::Simulator;
    using typename FlowProblemType::GridView;
    using typename FlowProblemType::FluidSystem;
    using typename FlowProblemType::Vanguard;

    using FlowProblemType::dim;
    using FlowProblemType::dimWorld;
    using FlowProblemType::numEq;
    using FlowProblemType::numPhases;
    using FlowProblemType::numComponents;

    // TODO: potentially some cleaning up depending on the usage later here
    using FlowProblemType::enableConvectiveMixing;
    using FlowProblemType::enableBrine;
    using FlowProblemType::enableDiffusion;
    using FlowProblemType::enableDispersion;
    using FlowProblemType::enableEnergy;
    using FlowProblemType::enableExperiments;
    using FlowProblemType::enableExtbo;
    using FlowProblemType::enableFoam;
    using FlowProblemType::enableMICP;
    using FlowProblemType::enablePolymer;
    using FlowProblemType::enablePolymerMolarWeight;
    using FlowProblemType::enableSaltPrecipitation;
    using FlowProblemType::enableSolvent;
    using FlowProblemType::enableTemperature;
    using FlowProblemType::enableThermalFluxBoundaries;

    using FlowProblemType::gasPhaseIdx;
    using FlowProblemType::oilPhaseIdx;
    using FlowProblemType::waterPhaseIdx;

    using typename FlowProblemType::RateVector;
    using typename FlowProblemType::PrimaryVariables;
    using typename FlowProblemType::Indices;
    using typename FlowProblemType::IntensiveQuantities;
    using typename FlowProblemType::ElementContext;

    using typename FlowProblemType::MaterialLaw;
    using typename FlowProblemType::DimMatrix;

    using SolventModule = BlackOilSolventModule<TypeTag>;
    using PolymerModule = BlackOilPolymerModule<TypeTag>;
    using FoamModule = BlackOilFoamModule<TypeTag>;
    using BrineModule = BlackOilBrineModule<TypeTag>;
    using ExtboModule = BlackOilExtboModule<TypeTag>;
    using MICPModule = BlackOilMICPModule<TypeTag>;
    using DispersionModule = BlackOilDispersionModule<TypeTag, enableDispersion>;
    using DiffusionModule = BlackOilDiffusionModule<TypeTag, enableDiffusion>;
    using ConvectiveMixingModule = BlackOilConvectiveMixingModule<TypeTag, enableConvectiveMixing>;
    using ModuleParams = typename BlackOilLocalResidualTPFA<TypeTag>::ModuleParams;

    using InitialFluidState = typename EquilInitializer<TypeTag>::ScalarFluidState;
    using EclWriterType = EclWriter<TypeTag>;
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
    }

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    explicit FlowProblemBlackoil(Simulator& simulator)
        : FlowProblemType(simulator)
        , mixControls_(simulator.vanguard().schedule())
        , actionHandler_(simulator.vanguard().eclState(),
                         simulator.vanguard().schedule(),
                         simulator.vanguard().actionState(),
                         simulator.vanguard().summaryState(),
                         this->wellModel_,
                         simulator.vanguard().grid().comm())
    {
        // Tell the black-oil extensions to initialize their internal data structures
        const auto& vanguard = simulator.vanguard();

        BlackOilBrineParams<Scalar> brineParams;
        brineParams.template initFromState<enableBrine,
                                           enableSaltPrecipitation>(vanguard.eclState());
        BrineModule::setParams(std::move(brineParams));

        DiffusionModule::initFromState(vanguard.eclState());
        DispersionModule::initFromState(vanguard.eclState());

        BlackOilExtboParams<Scalar> extboParams;
        extboParams.template initFromState<enableExtbo>(vanguard.eclState());
        ExtboModule::setParams(std::move(extboParams));

        BlackOilFoamParams<Scalar> foamParams;
        foamParams.template initFromState<enableFoam>(vanguard.eclState());
        FoamModule::setParams(std::move(foamParams));

        BlackOilMICPParams<Scalar> micpParams;
        micpParams.template initFromState<enableMICP>(vanguard.eclState());
        MICPModule::setParams(std::move(micpParams));

        BlackOilPolymerParams<Scalar> polymerParams;
        polymerParams.template initFromState<enablePolymer, enablePolymerMolarWeight>(vanguard.eclState());
        PolymerModule::setParams(std::move(polymerParams));

        BlackOilSolventParams<Scalar> solventParams;
        solventParams.template initFromState<enableSolvent>(vanguard.eclState(), vanguard.schedule());
        SolventModule::setParams(std::move(solventParams));

        // create the ECL writer
        eclWriter_ = std::make_unique<EclWriterType>(simulator);
        enableEclOutput_ = Parameters::Get<Parameters::EnableEclOutput>();
#if HAVE_DAMARIS
        // create Damaris writer
        damarisWriter_ = std::make_unique<DamarisWriterType>(simulator);
        enableDamarisOutput_ = Parameters::Get<Parameters::EnableDamarisOutput>();
#endif
    }

    /*!
     * \brief Called by the simulator before an episode begins.
     */
    void beginEpisode()
    {
        FlowProblemType::beginEpisode();

        auto& simulator = this->simulator();
        int episodeIdx = simulator.episodeIndex();
        const auto& schedule = simulator.vanguard().schedule();

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

        ConvectiveMixingModule::beginEpisode(simulator.vanguard().eclState(), simulator.vanguard().schedule(), episodeIdx, moduleParams_.convectiveMixingModuleParam);
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        // TODO: some might be able to move back to the base class
        FlowProblemType::finishInit();
        auto& simulator = this->simulator();

        auto finishTransmissibilities = [updated = false, this]() mutable
        {
            if (updated) { return; }
            this->transmissibilities_.finishInit([&vg = this->simulator().vanguard()](const unsigned int it) {
                return vg.gridIdxToEquilGridIdx(it);
            });
            updated = true;
        };

        // calculating the TRANX, TRANY, TRANZ and NNC for output purpose
        // for parallel running, it is based on global trans_
        // for serial running, it is based on the transmissibilities_
        // we try to avoid for the parallel running, has both global trans_ and transmissibilities_ allocated at the same time
        if (enableEclOutput_) {
            if (simulator.vanguard().grid().comm().size() > 1) {
                if (simulator.vanguard().grid().comm().rank() == 0)
                    eclWriter_->setTransmissibilities(&simulator.vanguard().globalTransmissibility());
            } else {
                finishTransmissibilities();
                eclWriter_->setTransmissibilities(&simulator.problem().eclTransmissibilities());
            }

            std::function<unsigned int(unsigned int)> equilGridToGrid = [&simulator](unsigned int i) {
                return simulator.vanguard().gridEquilIdxToGridIdx(i);
            };
            eclWriter_->extractOutputTransAndNNC(equilGridToGrid);
        }
        simulator.vanguard().releaseGlobalTransmissibilities();

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
        if (Parameters::Get<Parameters::EnableGravity>())
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
        this->readMaterialParameters_();
        this->readThermalParameters_();

        // write the static output files (EGRID, INIT)
        if (enableEclOutput_) {
            eclWriter_->writeInit();
        }

        finishTransmissibilities();

        const auto& initconfig = eclState.getInitConfig();
        this->tracerModel_.init(initconfig.restartRequested());
        if (initconfig.restartRequested())
            this->readEclRestartSolution_();
        else
            this->readInitialCondition_();

        this->tracerModel_.prepareTracerBatches();

        this->updatePffDofData_();

        if constexpr (getPropValue<TypeTag, Properties::EnablePolymer>()) {
            const auto& vanguard = this->simulator().vanguard();
            const auto& gridView = vanguard.gridView();
            int numElements = gridView.size(/*codim=*/0);
            this->polymer_.maxAdsorption.resize(numElements, 0.0);
        }

        this->readBoundaryConditions_();

        // compute and set eq weights based on initial b values
        this->computeAndSetEqWeights_();

        if (this->enableDriftCompensation_) {
            this->drift_.resize(this->model().numGridDof());
            this->drift_ = 0.0;
        }

        if (this->enableVtkOutput_ && eclState.getIOConfig().initOnly()) {
            simulator.setTimeStepSize(0.0);
            // TODO: this might need to move back to the base class depending how we use it
            FlowProblemType::ParentType::writeOutput(true);
        }

        // after finishing the initialization and writing the initial solution, we move
        // to the first "real" episode/report step
        // for restart the episode index and start is already set
        if (!initconfig.restartRequested()) {
            simulator.startNextEpisode(schedule.seconds(1));
            simulator.setEpisodeIndex(0);
            simulator.setTimeStepIndex(0);
        }
        // const auto& eclState = this->simulator().vanguard().eclState();
        this->mixControls_.init(this->model().numGridDof(),
                                this->episodeIndex(),
                                eclState.runspec().tabdims().getNumPVTTables());



    }

    /*!
     * \brief Called by the simulator after each time integration.
     */
    void endTimeStep()
    {
        FlowProblemType::endTimeStep();

        const bool isSubStep = !this->simulator().episodeWillBeOver();

        // For CpGrid with LGRs, ecl/vtk output is not supported yet.
        const auto& grid = this->simulator().vanguard().gridView().grid();

        using GridType = std::remove_cv_t<std::remove_reference_t<decltype(grid)>>;
        constexpr bool isCpGrid = std::is_same_v<GridType, Dune::CpGrid>;
        if (!isCpGrid || (grid.maxLevel() == 0)) {
            this->eclWriter_->evalSummaryState(isSubStep);
        }

        {
            OPM_TIMEBLOCK(applyActions);

            const int episodeIdx = this->episodeIndex();
            auto& simulator = this->simulator();

            // Re-ordering in case of Alugrid
            this->actionHandler_
                .applyActions(episodeIdx, simulator.time() + simulator.timeStepSize(),
                              [this](const bool global)
            {
                using TransUpdateQuantities = typename Vanguard::TransmissibilityType::TransUpdateQuantities;
                this->transmissibilities_
                    .update(global,  TransUpdateQuantities::All, [&vg = this->simulator().vanguard()]
                            (const unsigned int i)
                    {
                        return vg.gridIdxToEquilGridIdx(i);
                    });
            });
        }

        // Deal with "clogging" for the MICP model
        if constexpr (enableMICP) {
            auto& model = this->model();
            const auto& residual = model.linearizer().residual();

            for (unsigned globalDofIdx = 0; globalDofIdx < residual.size(); ++globalDofIdx) {
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
        const int episodeIdx = this->episodeIndex();
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
                .evalUDQAssignments(episodeIdx, this->simulator().vanguard().udqState());

        FlowProblemType::endEpisode();
    }

    /*!
     * \brief Write the requested quantities of the current solution into the output
     *        files.
     */
    void writeOutput(const SimulatorTimer& timer, bool verbose = true)
    {
        FlowProblemType::writeOutput(verbose);

        bool isSubStep = !this->simulator().episodeWillBeOver();

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
     * \copydoc FvBaseProblem::initialSolutionApplied()
     */
    void initialSolutionApplied() {
        FlowProblemType::initialSolutionApplied();

        if (this->simulator().episodeIndex() == 0) {
            eclWriter_->writeInitialFIPReport();
        }
    }

    void addToSourceDense(RateVector& rate,
                          unsigned globalDofIdx,
                          unsigned timeIdx) const
    {
        FlowProblemType::addToSourceDense(rate, globalDofIdx, timeIdx);

        const auto& source = this->simulator().vanguard().schedule()[this->episodeIndex()].source();
        std::array<int,3> ijk;
        this->simulator().vanguard().cartesianCoordinate(globalDofIdx, ijk);

        if (source.hasSource(ijk)) {
            if constexpr (enableSolvent) {
                const int pvtRegionIdx = this->pvtRegionIndex(globalDofIdx);
                Scalar mass_rate =
                        source.rate({ijk, SourceComponent::SOLVENT}) / this->model().dofTotalVolume(globalDofIdx);
                if constexpr (getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>()) {
                    const auto& solventPvt = SolventModule::solventPvt();
                    mass_rate /= solventPvt.referenceDensity(pvtRegionIdx);
                }
                rate[Indices::contiSolventEqIdx] += mass_rate;
            }
            if constexpr (enablePolymer) {
                rate[Indices::polymerConcentrationIdx] +=
                        source.rate({ijk, SourceComponent::POLYMER}) / this->model().dofTotalVolume(globalDofIdx);
            }
        }
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

    // temporary solution to facilitate output of initial state from flow
    const InitialFluidState& initialFluidState(unsigned globalDofIdx) const
    { return initialFluidStates_[globalDofIdx]; }

    const EclipseIO& eclIO() const
    { return eclWriter_->eclIO(); }

    void setSubStepReport(const SimulatorReportSingle& report)
    { return eclWriter_->setSubStepReport(report); }

    void setSimulationReport(const SimulatorReport& report)
    { return eclWriter_->setSimulationReport(report); }

    InitialFluidState boundaryFluidState(unsigned globalDofIdx, const int directionId) const
    {
        OPM_TIMEBLOCK_LOCAL(boundaryFluidState);
        const auto& bcprop = this->simulator().vanguard().schedule()[this->episodeIndex()].bcprop;
        if (bcprop.size() > 0) {
            FaceDir::DirEnum dir = FaceDir::FromIntersectionIndex(directionId);

            // index == 0: no boundary conditions for this
            // global cell and direction
            if (this->bcindex_(dir)[globalDofIdx] == 0)
                return initialFluidStates_[globalDofIdx];

            const auto& bc = bcprop[this->bcindex_(dir)[globalDofIdx]];
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

                    fluidState.setPc(phaseIdx, pc[phaseIdx]);
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

                for (unsigned activePhaseIdx = 0; activePhaseIdx < FluidSystem::numActivePhases(); ++activePhaseIdx) {
                    const auto phaseIdx = FluidSystem::activeToCanonicalPhaseIdx(activePhaseIdx);

                    const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState, phaseIdx, pvtRegionIdx);
                    fluidState.setInvB(phaseIdx, b);

                    const auto& rho = FluidSystem::density(fluidState, phaseIdx, pvtRegionIdx);
                    fluidState.setDensity(phaseIdx, rho);
                    if (enableEnergy) {
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


    const std::unique_ptr<EclWriterType>& eclWriter() const
    {
        return eclWriter_;
    }

    void setConvData(const std::vector<std::vector<int>>& data)
    {
        eclWriter_->mutableOutputModule().setCnvData(data);
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
     * or DRVDT are active. Nor if the porosity is changes between timesteps
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


    Scalar drsdtcon(unsigned elemIdx, int episodeIdx) const
    {
        return this->mixControls_.drsdtcon(elemIdx, episodeIdx,
                                           this->pvtRegionIndex(elemIdx));
    }

    const ModuleParams& moduleParams() const
    {
        return moduleParams_;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<FlowProblemType&>(*this));
        serializer(mixControls_);
        serializer(*eclWriter_);
    }

protected:
    void updateExplicitQuantities_(int episodeIdx, int timeStepSize, const bool first_step_after_restart = false) {
        FlowProblemType::updateExplicitQuantities_(first_step_after_restart);

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
            unsigned solventCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
            unsigned activeSolventCompIdx = Indices::canonicalToActiveComponentIndex(solventCompIdx);
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
        // Throw an exception if the grid has LGRs. Refined grid are not supported for restart.
        if(this->simulator().vanguard().grid().maxLevel() > 0) {
            throw std::invalid_argument("Refined grids are not yet supported for restart ");
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

        if constexpr (enableMICP) {
            this->micp_.resize(numElems);
        }

        // Initialize mixing controls before trying to set any lastRx valuesx
        this->mixControls_.init(numElems, restart_step, eclState.runspec().tabdims().getNumPVTTables());

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

            // For CO2STORE and H2STORE we need to set the initial temperature for isothermal simulations
            bool isThermal = eclState.getSimulationConfig().isThermal();
            bool needTemperature = (eclState.runspec().co2Storage() || eclState.runspec().h2Storage());
            if (!isThermal && needTemperature) {
                const auto& fp = simulator.vanguard().eclState().fieldProps();
                elemFluidState.setTemperature(fp.get_double("TEMPI")[elemIdx]);
            }

            this->mixControls_.updateLastValues(elemIdx, elemFluidState.Rs(), elemFluidState.Rv());

            if constexpr (enablePolymer)
                this->polymer_.concentration[elemIdx] = this->eclWriter_->outputModule().getPolymerConcentration(elemIdx);
            if constexpr (enableMICP){
                this->micp_.microbialConcentration[elemIdx] = this->eclWriter_->outputModule().getMicrobialConcentration(elemIdx);
                this->micp_.oxygenConcentration[elemIdx] = this->eclWriter_->outputModule().getOxygenConcentration(elemIdx);
                this->micp_.ureaConcentration[elemIdx] = this->eclWriter_->outputModule().getUreaConcentration(elemIdx);
                this->micp_.biofilmConcentration[elemIdx] = this->eclWriter_->outputModule().getBiofilmConcentration(elemIdx);
                this->micp_.calciteConcentration[elemIdx] = this->eclWriter_->outputModule().getCalciteConcentration(elemIdx);
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
            this->initial(sol[elemIdx], elemCtx, /*spaceIdx=*/0, /*timeIdx=*/0);
        }

        // make sure that the ghost and overlap entities exhibit the correct
        // solution. alternatively, this could be done in the loop above by also
        // considering non-interior elements. Since the initial() method might not work
        // 100% correctly for such elements, let's play safe and explicitly synchronize
        // using message passing.
        this->model().syncOverlap();

        this->eclWriter_->endRestart();
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

    void readInitialCondition_()
    {
        FlowProblemType::readInitialCondition_();

        if constexpr (enableSolvent || enablePolymer || enablePolymerMolarWeight || enableMICP)
            this->readBlackoilExtentionsInitialConditions_(this->model().numGridDof(),
                                                           enableSolvent,
                                                           enablePolymer,
                                                           enablePolymerMolarWeight,
                                                           enableMICP);

    }

    virtual void handleSolventBC(const BCProp::BCFace& bc, RateVector& rate) const {
        if constexpr (!enableSolvent)
            throw std::logic_error("solvent is disabled and you're trying to add solvent to BC");

        rate[Indices::solventSaturationIdx] = bc.rate;
    }

    virtual void handlePolymerBC(const BCProp::BCFace& bc, RateVector& rate) const {
        if constexpr (!enablePolymer)
            throw std::logic_error("polymer is disabled and you're trying to add polymer to BC");

        rate[Indices::polymerConcentrationIdx] = bc.rate;
    }

    std::vector<InitialFluidState> initialFluidStates_;

    bool enableEclOutput_;
    std::unique_ptr<EclWriterType> eclWriter_;
#if HAVE_DAMARIS
    bool enableDamarisOutput_ = false ;
    std::unique_ptr<DamarisWriterType> damarisWriter_;
#endif
    MixingRateControls<FluidSystem> mixControls_;

    ActionHandler<Scalar> actionHandler_;

    ModuleParams moduleParams_;
};

} // namespace Opm

#endif // OPM_FLOW_PROBLEM_BLACK_HPP
