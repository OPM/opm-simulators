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
/*!
 * \file
 *
 * \copydoc Opm::EclWriter
 */
#ifndef OPM_ECL_WRITER_HPP
#define OPM_ECL_WRITER_HPP

#include <dune/grid/common/partitionset.hh>

#include <opm/common/TimingMacros.hpp> // OPM_TIMEBLOCK
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/output/eclipse/Inplace.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/models/blackoil/blackoilproperties.hh> // Properties::EnableMech, EnableTemperature, EnableSolvent
#include <opm/models/common/multiphasebaseproperties.hh> // Properties::FluidSystem

#include <opm/simulators/flow/CollectDataOnIORank.hpp>
#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/flow/EclGenericWriter.hpp>
#include <opm/simulators/flow/FlowBaseVanguard.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelRestart.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>

#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Opm::Parameters {

// If available, write the ECL output in a non-blocking manner
struct EnableAsyncEclOutput { static constexpr bool value = true; };

// By default, use single precision for the ECL formated results
struct EclOutputDoublePrecision { static constexpr bool value = false; };

// Write all solutions for visualization, not just the ones for the
// report steps...
struct EnableWriteAllSolutions { static constexpr bool value = false; };

// Write ESMRY file for fast loading of summary data
struct EnableEsmry { static constexpr bool value = false; };

} // namespace Opm::Parameters

namespace Opm::Action {
    class State;
} // namespace Opm::Action

namespace Opm {
    class EclipseIO;
    class UDQState;
} // namespace Opm

namespace Opm {
/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Collects necessary output values and pass it to opm-common's ECL output.
 *
 * Caveats:
 * - For this class to do do anything meaningful, you will have to
 *   have the OPM module opm-common with ECL writing enabled.
 * - The only DUNE grid which is currently supported is Dune::CpGrid
 *   from the OPM module "opm-grid". Using another grid won't
 *   fail at compile time but you will provoke a fatal exception as
 *   soon as you try to write an ECL output file.
 * - This class requires to use the black oil model with the element
 *   centered finite volume discretization.
 */
template <class TypeTag, class OutputModule>
class EclWriter : public EclGenericWriter<GetPropType<TypeTag, Properties::Grid>,
                                          GetPropType<TypeTag, Properties::EquilGrid>,
                                          GetPropType<TypeTag, Properties::GridView>,
                                          GetPropType<TypeTag, Properties::ElementMapper>,
                                          GetPropType<TypeTag, Properties::Scalar>>
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    using BaseType = EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>;

    typedef Dune::MultipleCodimMultipleGeomTypeMapper< GridView > VertexMapper;

    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableMech = getPropValue<TypeTag, Properties::EnableMech>() };
    enum { enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };

public:

    static void registerParameters()
    {
        OutputModule::registerParameters();

        Parameters::Register<Parameters::EnableAsyncEclOutput>
            ("Write the ECL-formated results in a non-blocking way "
             "(i.e., using a separate thread).");
        Parameters::Register<Parameters::EnableEsmry>
            ("Write ESMRY file for fast loading of summary data.");
    }

    // The Simulator object should preferably have been const - the
    // only reason that is not the case is due to the SummaryState
    // object owned deep down by the vanguard.
    explicit EclWriter(Simulator& simulator)
        : BaseType(simulator.vanguard().schedule(),
                   simulator.vanguard().eclState(),
                   simulator.vanguard().summaryConfig(),
                   simulator.vanguard().grid(),
                   ((simulator.vanguard().grid().comm().rank() == 0)
                    ? &simulator.vanguard().equilGrid()
                    : nullptr),
                   simulator.vanguard().gridView(),
                   simulator.vanguard().cartesianIndexMapper(),
                   ((simulator.vanguard().grid().comm().rank() == 0)
                    ? &simulator.vanguard().equilCartesianIndexMapper()
                    : nullptr),
                   Parameters::Get<Parameters::EnableAsyncEclOutput>(),
                   Parameters::Get<Parameters::EnableEsmry>())
        , simulator_(simulator)
    {
#if HAVE_MPI
        if (this->simulator_.vanguard().grid().comm().size() > 1) {
            auto smryCfg = (this->simulator_.vanguard().grid().comm().rank() == 0)
                ? this->eclIO_->finalSummaryConfig()
                : SummaryConfig{};

            eclBroadcast(this->simulator_.vanguard().grid().comm(), smryCfg);

            this->outputModule_ = std::make_unique<OutputModule>
                (simulator, smryCfg, this->collectOnIORank_);
        }
        else
#endif
        {
            this->outputModule_ = std::make_unique<OutputModule>
                (simulator, this->eclIO_->finalSummaryConfig(), this->collectOnIORank_);
        }

        this->rank_ = this->simulator_.vanguard().grid().comm().rank();

        this->simulator_.vanguard().eclState().computeFipRegionStatistics();
    }

    ~EclWriter()
    {}

    const EquilGrid& globalGrid() const
    {
        return simulator_.vanguard().equilGrid();
    }

    /*!
     * \brief collect and pass data and pass it to eclIO writer
     */
    void evalSummaryState(bool isSubStep)
    {
        OPM_TIMEBLOCK(evalSummaryState);
        const int reportStepNum = simulator_.episodeIndex() + 1;

        /*
          The summary data is not evaluated for timestep 0, that is
          implemented with a:

             if (time_step == 0)
                 return;

          check somewhere in the summary code. When the summary code was
          split in separate methods Summary::eval() and
          Summary::add_timestep() it was necessary to pull this test out
          here to ensure that the well and group related keywords in the
          restart file, like XWEL and XGRP were "correct" also in the
          initial report step.

          "Correct" in this context means unchanged behavior, might very
          well be more correct to actually remove this if test.
        */

        if (reportStepNum == 0)
            return;

        const Scalar curTime = simulator_.time() + simulator_.timeStepSize();
        const Scalar totalCpuTime =
            simulator_.executionTimer().realTimeElapsed() +
            simulator_.setupTimer().realTimeElapsed() +
            simulator_.vanguard().setupTime();

        const auto localWellData            = simulator_.problem().wellModel().wellData();
        const auto localWBP                 = simulator_.problem().wellModel().wellBlockAveragePressures();
        const auto localGroupAndNetworkData = simulator_.problem().wellModel()
            .groupAndNetworkData(reportStepNum);

        const auto localAquiferData = simulator_.problem().aquiferModel().aquiferData();
        const auto localWellTestState = simulator_.problem().wellModel().wellTestState();
        this->prepareLocalCellData(isSubStep, reportStepNum);

        if (this->outputModule_->needInterfaceFluxes(isSubStep)) {
            this->captureLocalFluxData();
        }

        if (this->collectOnIORank_.isParallel()) {
            OPM_BEGIN_PARALLEL_TRY_CATCH()

            std::map<std::pair<std::string,int>,double> dummy;
            this->collectOnIORank_.collect({},
                                           outputModule_->getBlockData(),
                                           dummy,
                                           localWellData,
                                           localWBP,
                                           localGroupAndNetworkData,
                                           localAquiferData,
                                           localWellTestState,
                                           this->outputModule_->getInterRegFlows(),
                                           {},
                                           {});

            if (this->collectOnIORank_.isIORank()) {
                auto& iregFlows = this->collectOnIORank_.globalInterRegFlows();

                if (! iregFlows.readIsConsistent()) {
                    throw std::runtime_error {
                        "Inconsistent inter-region flow "
                        "region set names in parallel"
                    };
                }

                iregFlows.compress();
            }

            OPM_END_PARALLEL_TRY_CATCH("Collect to I/O rank: ",
                                       this->simulator_.vanguard().grid().comm());
        }


        std::map<std::string, double> miscSummaryData;
        std::map<std::string, std::vector<double>> regionData;
        Inplace inplace;

        {
            OPM_TIMEBLOCK(outputFipLogAndFipresvLog);

            inplace = outputModule_->calc_inplace(miscSummaryData, regionData, simulator_.gridView().comm());

            if (this->collectOnIORank_.isIORank()){
                inplace_ = inplace;
            }
        }

        // Add TCPU
        if (totalCpuTime != 0.0) {
            miscSummaryData["TCPU"] = totalCpuTime;
        }
        if (this->sub_step_report_.total_newton_iterations != 0) {
            miscSummaryData["NEWTON"] = this->sub_step_report_.total_newton_iterations;
        }
        if (this->sub_step_report_.total_linear_iterations != 0) {
            miscSummaryData["MLINEARS"] = this->sub_step_report_.total_linear_iterations;
        }
        if (this->sub_step_report_.total_newton_iterations != 0) {
            miscSummaryData["NLINEARS"] =  static_cast<float>(this->sub_step_report_.total_linear_iterations) / this->sub_step_report_.total_newton_iterations;
        }
        if (this->sub_step_report_.min_linear_iterations != std::numeric_limits<unsigned int>::max()) {
            miscSummaryData["NLINSMIN"] = this->sub_step_report_.min_linear_iterations;
        }
        if (this->sub_step_report_.max_linear_iterations != 0) {
            miscSummaryData["NLINSMAX"] = this->sub_step_report_.max_linear_iterations;
        }
        if (this->simulation_report_.success.total_linear_iterations != 0) {
            miscSummaryData["MSUMLINS"] = this->simulation_report_.success.total_linear_iterations;
        }
        if (this->simulation_report_.success.total_newton_iterations != 0) {
            miscSummaryData["MSUMNEWT"] = this->simulation_report_.success.total_newton_iterations;
        }

        {
            OPM_TIMEBLOCK(evalSummary);

            const auto& blockData = this->collectOnIORank_.isParallel()
                ? this->collectOnIORank_.globalBlockData()
                : this->outputModule_->getBlockData();

            const auto& interRegFlows = this->collectOnIORank_.isParallel()
                ? this->collectOnIORank_.globalInterRegFlows()
                : this->outputModule_->getInterRegFlows();

            this->evalSummary(reportStepNum,
                              curTime,
                              localWellData,
                              localWBP,
                              localGroupAndNetworkData,
                              localAquiferData,
                              blockData,
                              miscSummaryData,
                              regionData,
                              inplace,
                              this->outputModule_->initialInplace(),
                              interRegFlows,
                              this->summaryState(),
                              this->udqState());
        }
    }

    //! \brief Writes the initial FIP report as configured in RPTSOL.
    void writeInitialFIPReport()
    {
        const auto& gridView = simulator_.vanguard().gridView();
        const int num_interior = detail::
            countLocalInteriorCellsGridView(gridView);

        this->outputModule_->
            allocBuffers(num_interior, 0, false, false, /*isRestart*/ false);

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int dofIdx = 0; dofIdx < num_interior; ++dofIdx) {
            const auto& intQuants = *simulator_.model().cachedIntensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto totVolume = simulator_.model().dofTotalVolume(dofIdx);

            this->outputModule_->updateFluidInPlace(dofIdx, intQuants, totVolume);
        }

        // We always calculate the initial fip values as it may be used by various
        // keywords in the Schedule, e.g. FIP=2 in RPTSCHED but no FIP in RPTSOL
        outputModule_->calc_initial_inplace(simulator_.gridView().comm());

        // check if RPTSOL entry has FIP output
        const auto& fip = simulator_.vanguard().eclState().getEclipseConfig().fip();
        if (fip.output(FIPConfig::OutputField::FIELD) ||
            fip.output(FIPConfig::OutputField::RESV))
        {
            OPM_TIMEBLOCK(outputFipLogAndFipresvLog);
            boost::posix_time::ptime start_time =
                boost::posix_time::from_time_t(simulator_.vanguard().schedule().getStartTime());

            if (this->collectOnIORank_.isIORank()) {
                inplace_ = outputModule_->initialInplace().value();
                outputModule_->outputFipAndResvLog(inplace_, 0, 0.0, start_time,
                                                  false, simulator_.gridView().comm());
            }
        }

        outputModule_->outputFipAndResvLogToCSV(0, false, simulator_.gridView().comm());
    }

    void writeReports(const SimulatorTimer& timer)
    {
        if (! this->collectOnIORank_.isIORank()) {
            return;
        }

        // SimulatorTimer::reportStepNum() is the simulator's zero-based
        // "episode index".  This is generally the index value needed to
        // look up objects in the Schedule container.  That said, function
        // writeReports() is invoked at the *beginning* of a report
        // step/episode which means we typically need the objects from the
        // *previous* report step/episode.  We therefore need special case
        // handling for reportStepNum() == 0 in base runs and
        // reportStepNum() <= restart step in restarted runs.
        const auto firstStep = this->initialStep();
        const auto simStep =
            std::max(timer.reportStepNum() - 1, firstStep);

        const auto& rpt = this->schedule_[simStep].rpt_config();

        if (rpt.contains("WELSPECS") && (rpt.at("WELSPECS") > 0)) {
            // Requesting a well specification report is valid at all times,
            // including reportStepNum() == initialStep().
            this->writeWellspecReport(timer);
        }

        if (timer.reportStepNum() == firstStep) {
            // No dynamic flows at the beginning of the initialStep().
            return;
        }

        if (rpt.contains("WELLS") && rpt.at("WELLS") > 0) {
            this->writeWellflowReport(timer, simStep, rpt.at("WELLS"));
        }

        this->outputModule_->outputFipAndResvLog(this->inplace_,
                                                 timer.reportStepNum(),
                                                 timer.simulationTimeElapsed(),
                                                 timer.currentDateTime(),
                                                 /* isSubstep = */ false,
                                                 simulator_.gridView().comm());

        OpmLog::note("");   // Blank line after all reports.
    }

    void writeOutput(data::Solution&& localCellData, bool isSubStep)
    {
        OPM_TIMEBLOCK(writeOutput);

        const int reportStepNum = simulator_.episodeIndex() + 1;
        this->prepareLocalCellData(isSubStep, reportStepNum);
        this->outputModule_->outputErrorLog(simulator_.gridView().comm());

        // output using eclWriter if enabled
        auto localWellData = simulator_.problem().wellModel().wellData();
        auto localGroupAndNetworkData = simulator_.problem().wellModel()
            .groupAndNetworkData(reportStepNum);

        auto localAquiferData = simulator_.problem().aquiferModel().aquiferData();
        auto localWellTestState = simulator_.problem().wellModel().wellTestState();

        const bool isFlowsn = this->outputModule_->getFlows().hasFlowsn();
        auto flowsn = this->outputModule_->getFlows().getFlowsn();

        const bool isFloresn = this->outputModule_->getFlows().hasFloresn();
        auto floresn = this->outputModule_->getFlows().getFloresn();

        if (! isSubStep || Parameters::Get<Parameters::EnableWriteAllSolutions>()) {

            if (localCellData.empty()) {
                this->outputModule_->assignToSolution(localCellData);
            }

            // Add cell data to perforations for RFT output
            this->outputModule_->addRftDataToWells(localWellData,
                                                   reportStepNum,
                                                   simulator_.gridView().comm());
        }

        if (this->collectOnIORank_.isParallel() ||
            this->collectOnIORank_.doesNeedReordering())
        {
            // Note: We don't need WBP (well-block averaged pressures) or
            // inter-region flow rate values in order to create restart file
            // output.  There's consequently no need to collect those
            // properties on the I/O rank.

            this->collectOnIORank_.collect(localCellData,
                                           this->outputModule_->getBlockData(),
                                           this->outputModule_->getExtraBlockData(),
                                           localWellData,
                                           /* wbpData = */ {},
                                           localGroupAndNetworkData,
                                           localAquiferData,
                                           localWellTestState,
                                           /* interRegFlows = */ {},
                                           flowsn,
                                           floresn);
            if (this->collectOnIORank_.isIORank()) {
                this->outputModule_->assignGlobalFieldsToSolution(this->collectOnIORank_.globalCellData());
            }
        } else {
            this->outputModule_->assignGlobalFieldsToSolution(localCellData);
        }

        if (this->collectOnIORank_.isIORank()) {
            const Scalar curTime = simulator_.time() + simulator_.timeStepSize();
            const Scalar nextStepSize = simulator_.problem().nextTimeStepSize();
            std::optional<int> timeStepIdx;
            if (Parameters::Get<Parameters::EnableWriteAllSolutions>()) {
                timeStepIdx = simulator_.timeStepIndex();
            }
            this->doWriteOutput(reportStepNum, timeStepIdx, isSubStep,
                                std::move(localCellData),
                                std::move(localWellData),
                                std::move(localGroupAndNetworkData),
                                std::move(localAquiferData),
                                std::move(localWellTestState),
                                this->actionState(),
                                this->udqState(),
                                this->summaryState(),
                                this->simulator_.problem().thresholdPressure().getRestartVector(),
                                curTime, nextStepSize,
                                Parameters::Get<Parameters::EclOutputDoublePrecision>(),
                                isFlowsn, std::move(flowsn),
                                isFloresn, std::move(floresn));
        }
    }

    void beginRestart()
    {
        const auto enablePCHysteresis = simulator_.problem().materialLawManager()->enablePCHysteresis();
        const auto enableNonWettingHysteresis = simulator_.problem().materialLawManager()->enableNonWettingHysteresis();
        const auto enableWettingHysteresis = simulator_.problem().materialLawManager()->enableWettingHysteresis();
        const auto oilActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
        const auto gasActive = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
        const auto waterActive = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
        const auto enableSwatinit = simulator_.vanguard().eclState().fieldProps().has_double("SWATINIT");

        std::vector<RestartKey> solutionKeys {
            {"PRESSURE", UnitSystem::measure::pressure},
            {"SWAT",     UnitSystem::measure::identity,    waterActive},
            {"SGAS",     UnitSystem::measure::identity,    gasActive},
            {"TEMP",     UnitSystem::measure::temperature, enableEnergy},
            {"SSOLVENT", UnitSystem::measure::identity,    enableSolvent},

            {"RS",  UnitSystem::measure::gas_oil_ratio, FluidSystem::enableDissolvedGas()},
            {"RV",  UnitSystem::measure::oil_gas_ratio, FluidSystem::enableVaporizedOil()},
            {"RVW", UnitSystem::measure::oil_gas_ratio, FluidSystem::enableVaporizedWater()},
            {"RSW", UnitSystem::measure::gas_oil_ratio, FluidSystem::enableDissolvedGasInWater()},

            {"SGMAX", UnitSystem::measure::identity, enableNonWettingHysteresis && oilActive && gasActive},
            {"SHMAX", UnitSystem::measure::identity, enableWettingHysteresis && oilActive && gasActive},

            {"SOMAX", UnitSystem::measure::identity,
             (enableNonWettingHysteresis && oilActive && waterActive)
             || simulator_.problem().vapparsActive(simulator_.episodeIndex())},

            {"SOMIN", UnitSystem::measure::identity, enablePCHysteresis && oilActive && gasActive},
            {"SWHY1", UnitSystem::measure::identity, enablePCHysteresis && oilActive && waterActive},
            {"SWMAX", UnitSystem::measure::identity, enableWettingHysteresis && oilActive && waterActive},

            {"PPCW", UnitSystem::measure::pressure, enableSwatinit},
        };

        {
            const auto& tracers = simulator_.vanguard().eclState().tracer();

            for (const auto& tracer : tracers) {
                const auto enableSolTracer =
                    ((tracer.phase == Phase::GAS) && FluidSystem::enableDissolvedGas()) ||
                    ((tracer.phase == Phase::OIL) && FluidSystem::enableVaporizedOil());

                solutionKeys.emplace_back(tracer.fname(), UnitSystem::measure::identity, true);
                solutionKeys.emplace_back(tracer.sname(), UnitSystem::measure::identity, enableSolTracer);
            }
        }

        const auto& inputThpres = eclState().getSimulationConfig().getThresholdPressure();
        const std::vector<RestartKey> extraKeys {
            {"OPMEXTRA", UnitSystem::measure::identity, false},
            {"THRESHPR", UnitSystem::measure::pressure, inputThpres.active()},
        };

        const auto& gridView = this->simulator_.vanguard().gridView();
        const auto numElements = gridView.size(/*codim=*/0);

        // Try to load restart step 0 to calculate initial FIP
        {
            this->outputModule_->allocBuffers(numElements,
                                              0,
                                              /*isSubStep = */false,
                                              /*log = */      false,
                                              /*isRestart = */true);

            const auto restartSolution =
                loadParallelRestartSolution(this->eclIO_.get(),
                                            solutionKeys, gridView.comm(), 0);

            if (!restartSolution.empty()) {
                for (auto elemIdx = 0*numElements; elemIdx < numElements; ++elemIdx) {
                    const auto globalIdx = this->collectOnIORank_.localIdxToGlobalIdx(elemIdx);
                    this->outputModule_->setRestart(restartSolution, elemIdx, globalIdx);
                }

                this->simulator_.problem().readSolutionFromOutputModule(0, true);
                ElementContext elemCtx(this->simulator_);
                for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                    elemCtx.updatePrimaryStencil(elem);
                    elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                    this->outputModule_->updateFluidInPlace(elemCtx);
                }

                this->outputModule_->calc_initial_inplace(this->simulator_.gridView().comm());
            }
        }

        {
            // The episodeIndex is rewound one step back before calling
            // beginRestart() and cannot be used here.  We just ask the
            // initconfig directly to be sure that we use the correct index.
            const auto restartStepIdx = this->simulator_.vanguard()
                .eclState().getInitConfig().getRestartStep();

            this->outputModule_->allocBuffers(numElements,
                                              restartStepIdx,
                                              /*isSubStep = */false,
                                              /*log = */      false,
                                              /*isRestart = */true);
        }

        {
            const auto restartValues =
                loadParallelRestart(this->eclIO_.get(),
                                    this->actionState(),
                                    this->summaryState(),
                                    solutionKeys, extraKeys, gridView.comm());

            for (auto elemIdx = 0*numElements; elemIdx < numElements; ++elemIdx) {
                const auto globalIdx = this->collectOnIORank_.localIdxToGlobalIdx(elemIdx);
                this->outputModule_->setRestart(restartValues.solution, elemIdx, globalIdx);
            }

            auto& tracer_model = simulator_.problem().tracerModel();
            for (int tracer_index = 0; tracer_index < tracer_model.numTracers(); ++tracer_index) {
                // Free tracers
                {
                    const auto& free_tracer_name = tracer_model.fname(tracer_index);
                    const auto& free_tracer_solution = restartValues.solution
                        .template data<double>(free_tracer_name);

                    for (auto elemIdx = 0*numElements; elemIdx < numElements; ++elemIdx) {
                        const auto globalIdx = this->collectOnIORank_.localIdxToGlobalIdx(elemIdx);
                        tracer_model.setFreeTracerConcentration
                            (tracer_index, elemIdx, free_tracer_solution[globalIdx]);
                    }
                }

                // Solution tracer (only if DISGAS/VAPOIL are active for gas/oil tracers)
                if ((tracer_model.phase(tracer_index) == Phase::GAS && FluidSystem::enableDissolvedGas()) ||
                    (tracer_model.phase(tracer_index) == Phase::OIL && FluidSystem::enableVaporizedOil()))
                {
                    tracer_model.setEnableSolTracers(tracer_index, true);

                    const auto& sol_tracer_name = tracer_model.sname(tracer_index);
                    const auto& sol_tracer_solution = restartValues.solution
                        .template data<double>(sol_tracer_name);

                    for (auto elemIdx = 0*numElements; elemIdx < numElements; ++elemIdx) {
                        const auto globalIdx = this->collectOnIORank_.localIdxToGlobalIdx(elemIdx);
                        tracer_model.setSolTracerConcentration
                            (tracer_index, elemIdx, sol_tracer_solution[globalIdx]);
                    }
                }
                else {
                    tracer_model.setEnableSolTracers(tracer_index, false);

                    for (auto elemIdx = 0*numElements; elemIdx < numElements; ++elemIdx) {
                        tracer_model.setSolTracerConcentration(tracer_index, elemIdx, 0.0);
                    }
                }
            }

            if (inputThpres.active()) {
                const_cast<Simulator&>(this->simulator_)
                    .problem().thresholdPressure()
                    .setFromRestart(restartValues.getExtra("THRESHPR"));
            }

            restartTimeStepSize_ = restartValues.getExtra("OPMEXTRA")[0];
            if (restartTimeStepSize_ <= 0) {
                restartTimeStepSize_ = std::numeric_limits<double>::max();
            }

            // Initialize the well model from restart values
            this->simulator_.problem().wellModel()
                .initFromRestartFile(restartValues);

            if (!restartValues.aquifer.empty()) {
                this->simulator_.problem().mutableAquiferModel()
                    .initFromRestart(restartValues.aquifer);
            }
        }
    }

    void endRestart()
    {
        // Calculate initial in-place volumes.
        // Does nothing if they have already been calculated,
        // e.g. from restart data at T=0.
        this->outputModule_->calc_initial_inplace(this->simulator_.gridView().comm());

        if (this->collectOnIORank_.isIORank()) {
            if (this->outputModule_->initialInplace().has_value()) {
                this->inplace_ = this->outputModule_->initialInplace().value();
            }
        }
    }

    const OutputModule& outputModule() const
    { return *outputModule_; }

    OutputModule& mutableOutputModule() const
    { return *outputModule_; }

    Scalar restartTimeStepSize() const
    { return restartTimeStepSize_; }

    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(*outputModule_);
    }

private:
    static bool enableEclOutput_()
    {
        static bool enable = Parameters::Get<Parameters::EnableEclOutput>();
        return enable;
    }

    const EclipseState& eclState() const
    { return simulator_.vanguard().eclState(); }

    SummaryState& summaryState()
    { return simulator_.vanguard().summaryState(); }

    Action::State& actionState()
    { return simulator_.vanguard().actionState(); }

    UDQState& udqState()
    { return simulator_.vanguard().udqState(); }

    const Schedule& schedule() const
    { return simulator_.vanguard().schedule(); }

    void prepareLocalCellData(const bool isSubStep,
                              const int  reportStepNum)
    {
        OPM_TIMEBLOCK(prepareLocalCellData);

        if (this->outputModule_->localDataValid()) {
            return;
        }

        const auto& gridView = simulator_.vanguard().gridView();
        const bool log = this->collectOnIORank_.isIORank();

        const int num_interior = detail::
            countLocalInteriorCellsGridView(gridView);
        this->outputModule_->
            allocBuffers(num_interior, reportStepNum,
                         isSubStep && !Parameters::Get<Parameters::EnableWriteAllSolutions>(),
                         log, /*isRestart*/ false);

        ElementContext elemCtx(simulator_);

        OPM_BEGIN_PARALLEL_TRY_CATCH();

        {
            OPM_TIMEBLOCK(prepareCellBasedData);

            this->outputModule_->prepareDensityAccumulation();
            this->outputModule_->setupExtractors(isSubStep, reportStepNum);
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                this->outputModule_->processElement(elemCtx);
                this->outputModule_->processElementBlockData(elemCtx);
            }
            this->outputModule_->clearExtractors();

            this->outputModule_->accumulateDensityParallel();
        }

        {
            OPM_TIMEBLOCK(prepareFluidInPlace);

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int dofIdx = 0; dofIdx < num_interior; ++dofIdx) {
                const auto& intQuants = *simulator_.model().cachedIntensiveQuantities(dofIdx, /*timeIdx=*/0);
                const auto totVolume = simulator_.model().dofTotalVolume(dofIdx);

                this->outputModule_->updateFluidInPlace(dofIdx, intQuants, totVolume);
            }
        }

        this->outputModule_->validateLocalData();

        OPM_END_PARALLEL_TRY_CATCH("EclWriter::prepareLocalCellData() failed: ",
                                   this->simulator_.vanguard().grid().comm());
    }

    void captureLocalFluxData()
    {
        OPM_TIMEBLOCK(captureLocalData);

        const auto& gridView = this->simulator_.vanguard().gridView();
        const auto timeIdx = 0u;

        auto elemCtx = ElementContext { this->simulator_ };

        const auto elemMapper = ElementMapper { gridView, Dune::mcmgElementLayout() };
        const auto activeIndex = [&elemMapper](const Element& e)
        {
            return elemMapper.index(e);
        };

        const auto cartesianIndex = [this](const int elemIndex)
        {
            return this->cartMapper_.cartesianIndex(elemIndex);
        };

        this->outputModule_->initializeFluxData();

        OPM_BEGIN_PARALLEL_TRY_CATCH();

        for (const auto& elem : elements(gridView, Dune::Partitions::interiorBorder)) {
            elemCtx.updateStencil(elem);
            elemCtx.updateIntensiveQuantities(timeIdx);
            elemCtx.updateExtensiveQuantities(timeIdx);

            this->outputModule_->processFluxes(elemCtx, activeIndex, cartesianIndex);
        }

        OPM_END_PARALLEL_TRY_CATCH("EclWriter::captureLocalFluxData() failed: ",
                                   this->simulator_.vanguard().grid().comm())

        this->outputModule_->finalizeFluxData();
    }

    void writeWellspecReport(const SimulatorTimer& timer) const
    {
        const auto changedWells = this->schedule_
            .changed_wells(timer.reportStepNum(), this->initialStep());

        if (changedWells.empty()) {
            return;
        }

        this->outputModule_->outputWellspecReport(changedWells,
                                                  timer.reportStepNum(),
                                                  timer.simulationTimeElapsed(),
                                                  timer.currentDateTime());
    }

    void writeWellflowReport(const SimulatorTimer& timer,
                             const int             simStep,
                             const int             wellsRequest) const
    {
        this->outputModule_->outputTimeStamp("WELLS",
                                             timer.simulationTimeElapsed(),
                                             timer.reportStepNum(),
                                             timer.currentDateTime());

        const auto wantConnData = wellsRequest > 1;

        this->outputModule_->outputProdLog(simStep, wantConnData);
        this->outputModule_->outputInjLog(simStep, wantConnData);
        this->outputModule_->outputCumLog(simStep, wantConnData);
        this->outputModule_->outputMSWLog(simStep);
    }

    int initialStep() const
    {
        const auto& initConfig = this->eclState().cfg().init();

        return initConfig.restartRequested()
            ? initConfig.getRestartStep()
            : 0;
    }

    Simulator& simulator_;
    std::unique_ptr<OutputModule> outputModule_;
    Scalar restartTimeStepSize_;
    int rank_ ;
    Inplace inplace_;
};

} // namespace Opm

#endif // OPM_ECL_WRITER_HPP
