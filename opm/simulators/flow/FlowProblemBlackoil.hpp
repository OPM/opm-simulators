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

#include <opm/output/eclipse/EclipseIO.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/flow/ActionHandler.hpp>
#include <opm/simulators/flow/FlowProblem.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/flow/FlowThresholdPressure.hpp>
#include <opm/simulators/flow/MixingRateControls.hpp>
#include <opm/simulators/flow/OutputBlackoilModule.hpp>
#include <opm/simulators/flow/VtkTracerModule.hpp>

#include <opm/simulators/utils/satfunc/SatfuncConsistencyCheckManager.hpp>

#if HAVE_DAMARIS
#include <opm/simulators/flow/DamarisWriter.hpp>
#endif

#include <algorithm>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
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

    enum { enableDissolvedGas = Indices::compositionSwitchIdx >= 0 };
    enum { enableVapwat = getPropValue<TypeTag, Properties::EnableVapwat>() };
    enum { enableDisgasInWater = getPropValue<TypeTag, Properties::EnableDisgasInWater>() };

    using SolventModule = BlackOilSolventModule<TypeTag>;
    using PolymerModule = BlackOilPolymerModule<TypeTag>;
    using FoamModule = BlackOilFoamModule<TypeTag>;
    using BrineModule = BlackOilBrineModule<TypeTag>;
    using ExtboModule = BlackOilExtboModule<TypeTag>;
    using BioeffectsModule = BlackOilBioeffectsModule<TypeTag>;
    using DispersionModule = BlackOilDispersionModule<TypeTag, enableDispersion>;
    using DiffusionModule = BlackOilDiffusionModule<TypeTag, enableDiffusion>;
    using ConvectiveMixingModule = BlackOilConvectiveMixingModule<TypeTag, enableConvectiveMixing>;
    using ModuleParams = typename BlackOilLocalResidualTPFA<TypeTag>::ModuleParams;

    using InitialFluidState = typename EquilInitializer<TypeTag>::ScalarFluidState;
    using EclWriterType = EclWriter<TypeTag, OutputBlackOilModule<TypeTag> >;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
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
        , thresholdPressures_(simulator)
        , mixControls_(simulator.vanguard().schedule())
        , actionHandler_(simulator.vanguard().eclState(),
                         simulator.vanguard().schedule(),
                         simulator.vanguard().actionState(),
                         simulator.vanguard().summaryState(),
                         this->wellModel_,
                         simulator.vanguard().grid().comm())
    {
        this->model().addOutputModule(std::make_unique<VtkTracerModule<TypeTag>>(simulator));

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

        BlackOilBioeffectsParams<Scalar> bioeffectsParams;
        bioeffectsParams.template initFromState<enableBioeffects, enableMICP>(vanguard.eclState());
        BioeffectsModule::setParams(std::move(bioeffectsParams));

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
    void beginEpisode() override
    {
        FlowProblemType::beginEpisode();

        auto& simulator = this->simulator();

        const int episodeIdx = simulator.episodeIndex();
        const auto& schedule = simulator.vanguard().schedule();

        // Evaluate UDQ assign statements to make sure the settings are
        // available as UDA controls for the current report step.
        this->actionHandler_
            .evalUDQAssignments(episodeIdx, simulator.vanguard().udqState());

        if (episodeIdx >= 0) {
            const auto& oilVap = schedule[episodeIdx].oilvap();
            if (oilVap.getType() == OilVaporizationProperties::OilVaporization::VAPPARS) {
                FluidSystem::setVapPars(oilVap.vap1(), oilVap.vap2());
            }
            else {
                FluidSystem::setVapPars(0.0, 0.0);
            }
        }

        ConvectiveMixingModule::beginEpisode(simulator.vanguard().eclState(), schedule, episodeIdx,
                                             this->moduleParams_.convectiveMixingModuleParam);
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        // TODO: there should be room to remove duplication for this
        // function, but there is relatively complicated logic in the
        // function calls here.  Some refactoring is needed.
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

            this->eclWriter_->extractOutputTransAndNNC(equilGridToGrid);
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
        if (Parameters::Get<Parameters::EnableGravity>() &&
            eclState.getInitConfig().hasGravity())
        {
            // unit::gravity is 9.80665 m^2/s--i.e., standard measure at Tellus equator.
            this->gravity_[dim - 1] = unit::gravity;
        }

        if (this->enableTuning_) {
            // if support for the TUNING keyword is enabled, we get the initial time
            // steping parameters from it instead of from command line parameters
            const auto& tuning = schedule[0].tuning();
            this->initialTimeStepSize_ = tuning.TSINIT.has_value() ? tuning.TSINIT.value() : -1.0;
            this->maxTimeStepAfterWellEvent_ = tuning.TMAXWC;
        }

        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
            FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            this->maxOilSaturation_.resize(this->model().numGridDof(), 0.0);
        }

        this->readRockParameters_(simulator.vanguard().cellCenterDepths(),
                                  [&simulator](const unsigned idx)
                                  {
                                      std::array<int,dim> coords;
                                      simulator.vanguard().cartesianCoordinate(idx, coords);
                                      std::transform(coords.begin(), coords.end(), coords.begin(),
                                                     [](const auto c) { return c + 1; });
                                      return coords;
                                  });

        this->readMaterialParameters_();
        this->readThermalParameters_();

        // write the static output files (EGRID, INIT)
        if (enableEclOutput_) {
            this->eclWriter_->writeInit();
        }

        finishTransmissibilities();

        const auto& initconfig = eclState.getInitConfig();
        this->tracerModel_.init(initconfig.restartRequested());
        if (initconfig.restartRequested()) {
            this->readEclRestartSolution_();
        }
        else {
            this->readInitialCondition_();
        }

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

        if (this->enableDriftCompensation_) {
            this->drift_.resize(this->model().numGridDof());
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
        this->mixControls_.init(this->model().numGridDof(),
                                this->episodeIndex(),
                                eclState.runspec().tabdims().getNumPVTTables());

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
        // After the solution is updated, the values in output module needs
        // also updated.
        this->eclWriter().mutableOutputModule().invalidateLocalData();

        // For CpGrid with LGRs, ecl/vtk output is not supported yet.
        const auto& grid = this->simulator().vanguard().gridView().grid();

        using GridType = std::remove_cv_t<std::remove_reference_t<decltype(grid)>>;
        constexpr bool isCpGrid = std::is_same_v<GridType, Dune::CpGrid>;
        if (!isCpGrid || (grid.maxLevel() == 0)) {
            this->eclWriter_->evalSummaryState(!this->episodeWillBeOver());
        }

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
            this->eclWriter_->writeOutput(std::move(localCellData), isSubStep);
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

        // For CpGrid with LGRs, ecl-output is not supported yet.
        const auto& grid = this->simulator().vanguard().gridView().grid();

        using GridType = std::remove_cv_t<std::remove_reference_t<decltype(grid)>>;
        constexpr bool isCpGrid = std::is_same_v<GridType, Dune::CpGrid>;
        // Skip - for now -  calculate the initial fip values for CpGrid with LGRs.
        if (!isCpGrid || (grid.maxLevel() == 0)) {
            if (this->simulator().episodeIndex() == 0) {
                eclWriter_->writeInitialFIPReport();
            }
        }
    }

    void addToSourceDense(RateVector& rate,
                          unsigned globalDofIdx,
                          unsigned timeIdx) const override
    {
        this->aquiferModel_.addToSource(rate, globalDofIdx, timeIdx);

        // Add source term from deck
        const auto& source = this->simulator().vanguard().schedule()[this->episodeIndex()].source();
        std::array<int,3> ijk;
        this->simulator().vanguard().cartesianCoordinate(globalDofIdx, ijk);

        if (source.hasSource(ijk)) {
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
            if constexpr (enableEnergy) {
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
     *
     * TODO: The API of this is a bit ad-hoc, it would be better to use context objects.
     */
    template <class LhsEval>
    LhsEval permFactTransMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(permFactTransMultiplier, Subsystem::PvtProps);
        if constexpr (enableSaltPrecipitation) {
            const auto& fs = intQuants.fluidState();
            unsigned tableIdx = this->simulator().problem().satnumRegionIndex(elementIdx);
            LhsEval porosityFactor = decay<LhsEval>(1. - fs.saltSaturation());
            porosityFactor = min(porosityFactor, 1.0);
            const auto& permfactTable = BrineModule::permfactTable(tableIdx);
            return permfactTable.eval(porosityFactor, /*extrapolation=*/true);
        }
        else if constexpr (enableBioeffects) {
            return intQuants.permFactor().value();
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

    InitialFluidState boundaryFluidState(unsigned globalDofIdx, const int directionId) const
    {
        OPM_TIMEBLOCK_LOCAL(boundaryFluidState, Subsystem::Assembly);
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

                double temperature = initialFluidStates_[globalDofIdx].temperature(0); // we only have one temperature
                const auto temperature_input = bc.temperature;
                if(temperature_input)
                    temperature = *temperature_input;
                fluidState.setTemperature(temperature);

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
                    if constexpr (enableEnergy) {
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

        if constexpr (enableBioeffects) {
            values[Indices::microbialConcentrationIdx] = this->bioeffects_.microbialConcentration[globalDofIdx];
            values[Indices::biofilmVolumeFractionIdx]= this->bioeffects_.biofilmVolumeFraction[globalDofIdx];
            if constexpr (enableMICP) {
                values[Indices::oxygenConcentrationIdx]= this->bioeffects_.oxygenConcentration[globalDofIdx];
                values[Indices::ureaConcentrationIdx]= this->bioeffects_.ureaConcentration[globalDofIdx];
                values[Indices::calciteVolumeFractionIdx]= this->bioeffects_.calciteVolumeFraction[globalDofIdx];
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
            values.setThermalFlow(context, spaceIdx, timeIdx, this->initialFluidStates_[globalDofIdx] );
        }

        if (this->nonTrivialBoundaryConditions()) {
            unsigned indexInInside  = context.intersection(spaceIdx).indexInInside();
            unsigned interiorDofIdx = context.interiorScvIndex(spaceIdx, timeIdx);
            unsigned globalDofIdx = context.globalSpaceIndex(interiorDofIdx, timeIdx);
            unsigned pvtRegionIdx = pvtRegionIndex(context, spaceIdx, timeIdx);
            const auto [type, massrate] = this->boundaryCondition(globalDofIdx, indexInInside);
            if (type == BCType::THERMAL)
                values.setThermalFlow(context, spaceIdx, timeIdx, this->boundaryFluidState(globalDofIdx, indexInInside));
            else if (type == BCType::FREE || type == BCType::DIRICHLET)
                values.setFreeFlow(context, spaceIdx, timeIdx, this->boundaryFluidState(globalDofIdx, indexInInside));
            else if (type == BCType::RATE)
                values.setMassRate(massrate, pvtRegionIdx);
        }
    }

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
            this->mixControls_.init(this->model().numGridDof(),
                                    this->episodeIndex(),
                                    eclState.runspec().tabdims().getNumPVTTables());
        }
    }

    /*!
     * \copydoc BlackOilBaseProblem::thresholdPressure
     */
    Scalar thresholdPressure(unsigned elem1Idx, unsigned elem2Idx) const
    { return thresholdPressures_.thresholdPressure(elem1Idx, elem2Idx); }

    const FlowThresholdPressure<TypeTag>& thresholdPressure() const
    { return thresholdPressures_; }

    FlowThresholdPressure<TypeTag>& thresholdPressure()
    { return thresholdPressures_; }

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

        if (FluidSystem::enableDissolvedGasInWater() && has_rsw)
            rswData = fp.get_double("RSW");

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

    void handleSolventBC(const BCProp::BCFace& bc, RateVector& rate) const override
    {
        if constexpr (!enableSolvent)
            throw std::logic_error("solvent is disabled and you're trying to add solvent to BC");

        rate[Indices::solventSaturationIdx] = bc.rate;
    }

    void handlePolymerBC(const BCProp::BCFace& bc, RateVector& rate) const override
    {
        if constexpr (!enablePolymer)
            throw std::logic_error("polymer is disabled and you're trying to add polymer to BC");

        rate[Indices::polymerConcentrationIdx] = bc.rate;
    }

    void handleMicrBC(const BCProp::BCFace& bc, RateVector& rate) const override
    {
        if constexpr (!enableMICP)
            throw std::logic_error("MICP is disabled and you're trying to add microbes to BC");

        rate[Indices::microbialConcentrationIdx] = bc.rate;
    }

    void handleOxygBC(const BCProp::BCFace& bc, RateVector& rate) const override
    {
        if constexpr (!enableMICP)
            throw std::logic_error("MICP is disabled and you're trying to add oxygen to BC");

        rate[Indices::oxygenConcentrationIdx] = bc.rate;
    }

    void handleUreaBC(const BCProp::BCFace& bc, RateVector& rate) const override
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

    FlowThresholdPressure<TypeTag> thresholdPressures_;

    std::vector<InitialFluidState> initialFluidStates_;

    bool enableEclOutput_;
    std::unique_ptr<EclWriterType> eclWriter_;

    const Scalar smallSaturationTolerance_ = 1.e-6;
#if HAVE_DAMARIS
    bool enableDamarisOutput_ = false ;
    std::unique_ptr<DamarisWriterType> damarisWriter_;
#endif
    MixingRateControls<FluidSystem> mixControls_;

    ActionHandler<Scalar, IndexTraits> actionHandler_;

    ModuleParams moduleParams_;

private:
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
