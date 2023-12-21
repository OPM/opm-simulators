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
#ifndef EWOMS_ECL_WRITER_HH
#define EWOMS_ECL_WRITER_HH

#include <dune/grid/common/partitionset.hh>

#include <ebos/collecttoiorank.hh>
#include <ebos/eclbasevanguard.hh>
#include <ebos/eclgenericwriter.hh>
#include <ebos/ecloutputblackoilmodule.hh>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelRestart.hpp>
#include <opm/simulators/flow/countGlobalCells.hpp>


#include <opm/common/OpmLog/OpmLog.hpp>

#include <limits>
#include <stdexcept>
#include <string>


// #include <opm/simulators/utils/GridDataOutput.hpp>


namespace Opm::Properties {

template<class TypeTag, class MyTypeTag>
struct EnableEclOutput {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableAsyncEclOutput {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EclOutputDoublePrecision {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableEsmry {
    using type = UndefinedProperty;
};
} // namespace Opm::Properties

namespace Opm {

namespace Action { class State; }
class EclipseIO;
class UDQState;

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Collects necessary output values and pass it to opm-output.
 *
 * Caveats:
 * - For this class to do do anything meaningful, you will have to
 *   have the OPM module opm-output.
 * - The only DUNE grid which is currently supported is Dune::CpGrid
 *   from the OPM module "opm-grid". Using another grid won't
 *   fail at compile time but you will provoke a fatal exception as
 *   soon as you try to write an ECL output file.
 * - This class requires to use the black oil model with the element
 *   centered finite volume discretization.
 */
template <class TypeTag>
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
        EclOutputBlackOilModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableAsyncEclOutput,
                             "Write the ECL-formated results in a non-blocking way (i.e., using a separate thread).");

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableEsmry,
                             "Write ESMRY file for fast loading of summary data.");
                                                        
    }

    // The Simulator object should preferably have been const - the
    // only reason that is not the case is due to the SummaryState
    // object owned deep down by the vanguard.
    EclWriter(Simulator& simulator)
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
                   EWOMS_GET_PARAM(TypeTag, bool, EnableAsyncEclOutput),
                   EWOMS_GET_PARAM(TypeTag, bool, EnableEsmry))
        , simulator_(simulator)
    {
        this->eclOutputModule_ = std::make_unique<EclOutputBlackOilModule<TypeTag>>
            (simulator, this->collectToIORank_);
            
        rank_ = simulator_.vanguard().grid().comm().rank() ;
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

        if (this->eclOutputModule_->needInterfaceFluxes(isSubStep)) {
            this->captureLocalFluxData();
        }

        if (this->collectToIORank_.isParallel()) {
            OPM_BEGIN_PARALLEL_TRY_CATCH()

            this->collectToIORank_.collect({},
                                           eclOutputModule_->getBlockData(),
                                           localWellData,
                                           localWBP,
                                           localGroupAndNetworkData,
                                           localAquiferData,
                                           localWellTestState,
                                           this->eclOutputModule_->getInterRegFlows(),
                                           {},
                                           {});

            if (this->collectToIORank_.isIORank()) {
                auto& iregFlows = this->collectToIORank_.globalInterRegFlows();

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
            inplace = eclOutputModule_->outputFipLog(miscSummaryData, regionData, reportStepNum,
                                                     isSubStep, simulator_.gridView().comm());
            eclOutputModule_->outputFipresvLog(miscSummaryData, regionData, reportStepNum,
                                               isSubStep, simulator_.gridView().comm());
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
        if (this->simulation_report_.success.total_newton_iterations != 0) {
            miscSummaryData["MSUMLINS"] = this->simulation_report_.success.total_linear_iterations;
        }
        if (this->simulation_report_.success.total_newton_iterations != 0) {
            miscSummaryData["MSUMNEWT"] = this->simulation_report_.success.total_newton_iterations;
        }

        {
            OPM_TIMEBLOCK(evalSummary);

            const auto& blockData = this->collectToIORank_.isParallel()
                ? this->collectToIORank_.globalBlockData()
                : this->eclOutputModule_->getBlockData();

            const auto& interRegFlows = this->collectToIORank_.isParallel()
                ? this->collectToIORank_.globalInterRegFlows()
                : this->eclOutputModule_->getInterRegFlows();

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
                              this->eclOutputModule_->initialInplace(),
                              interRegFlows,
                              this->summaryState(),
                              this->udqState());
        }

        if (! isSubStep) {
            OPM_TIMEBLOCK(outputProdInjLogs);

            eclOutputModule_->outputProdLog(reportStepNum);
            eclOutputModule_->outputInjLog(reportStepNum);
            eclOutputModule_->outputCumLog(reportStepNum);

            OpmLog::note("");   // Blank line after all reports.
        }
    }

    //! \brief Writes the initial FIP report as configured in RPTSOL.
    void writeInitialFIPReport()
    {
        const auto& fip = simulator_.vanguard().eclState().getEclipseConfig().fip();
        if (!fip.output(FIPConfig::OutputField::FIELD) &&
            !fip.output(FIPConfig::OutputField::RESV)) {
            return;
        }

        const auto& gridView = simulator_.vanguard().gridView();
        const int num_interior = detail::
            countLocalInteriorCellsGridView(gridView);

        this->eclOutputModule_->
            allocBuffers(num_interior, 0, false, false, /*isRestart*/ false);

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int dofIdx = 0; dofIdx < num_interior; ++dofIdx) {
            const auto& intQuants = *simulator_.model().cachedIntensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto totVolume = simulator_.model().dofTotalVolume(dofIdx);

            this->eclOutputModule_->updateFluidInPlace(dofIdx, intQuants, totVolume);
        }

        std::map<std::string, double> miscSummaryData;
        std::map<std::string, std::vector<double>> regionData;
        Inplace inplace;
        {
            OPM_TIMEBLOCK(outputFipLogAndFipresvLog);
            inplace = eclOutputModule_->outputFipLog(miscSummaryData, regionData, 0,
                                                     false, simulator_.gridView().comm());
            eclOutputModule_->outputFipresvLog(miscSummaryData, regionData, 0,
                                               false, simulator_.gridView().comm());
        }
    }

    void writeOutput(data::Solution&& localCellData, bool isSubStep)
    {
        OPM_TIMEBLOCK(writeOutput);

        const int reportStepNum = simulator_.episodeIndex() + 1;
        this->prepareLocalCellData(isSubStep, reportStepNum);
        this->eclOutputModule_->outputErrorLog(simulator_.gridView().comm());

        // output using eclWriter if enabled
        auto localWellData = simulator_.problem().wellModel().wellData();
        auto localGroupAndNetworkData = simulator_.problem().wellModel()
            .groupAndNetworkData(reportStepNum);

        auto localAquiferData = simulator_.problem().aquiferModel().aquiferData();
        auto localWellTestState = simulator_.problem().wellModel().wellTestState();

        const bool isFlowsn = this->eclOutputModule_->hasFlowsn();
        auto flowsn = this->eclOutputModule_->getFlowsn();

        const bool isFloresn = this->eclOutputModule_->hasFloresn();
        auto floresn = this->eclOutputModule_->getFloresn();

        // data::Solution localCellData = {};
        if (! isSubStep) {
            
            if (localCellData.empty()) {
                this->eclOutputModule_->assignToSolution(localCellData);
            }

            // Add cell data to perforations for RFT output
            this->eclOutputModule_->addRftDataToWells(localWellData, reportStepNum);
        }

        if (this->collectToIORank_.isParallel() ||
            this->collectToIORank_.doesNeedReordering())
        {
            // Note: We don't need WBP (well-block averaged pressures) or
            // inter-region flow rate values in order to create restart file
            // output.  There's consequently no need to collect those
            // properties on the I/O rank.

            this->collectToIORank_.collect(localCellData,
                                           this->eclOutputModule_->getBlockData(),
                                           localWellData,
                                           /* wbpData = */ {},
                                           localGroupAndNetworkData,
                                           localAquiferData,
                                           localWellTestState,
                                           /* interRegFlows = */ {},
                                           flowsn,
                                           floresn);
        }

        if (this->collectToIORank_.isIORank()) {
            const Scalar curTime = simulator_.time() + simulator_.timeStepSize();
            const Scalar nextStepSize = simulator_.problem().nextTimeStepSize();

            this->doWriteOutput(reportStepNum, isSubStep,
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
                                EWOMS_GET_PARAM(TypeTag, bool, EclOutputDoublePrecision),
                                isFlowsn, std::move(flowsn),
                                isFloresn, std::move(floresn));
        }
    }

    void beginRestart()
    {
        bool enableHysteresis = simulator_.problem().materialLawManager()->enableHysteresis();
        bool enableSwatinit = simulator_.vanguard().eclState().fieldProps().has_double("SWATINIT");
        bool opm_rst_file = EWOMS_GET_PARAM(TypeTag, bool, EnableOpmRstFile);
        bool read_temp = enableEnergy || (opm_rst_file && enableTemperature);
        std::vector<RestartKey> solutionKeys{
            {"PRESSURE", UnitSystem::measure::pressure},
            {"SWAT", UnitSystem::measure::identity, static_cast<bool>(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))},
            {"SGAS", UnitSystem::measure::identity, static_cast<bool>(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))},
            {"TEMP" , UnitSystem::measure::temperature, read_temp},
            {"SSOLVENT" , UnitSystem::measure::identity, enableSolvent},
            {"RS", UnitSystem::measure::gas_oil_ratio, FluidSystem::enableDissolvedGas()},
            {"RV", UnitSystem::measure::oil_gas_ratio, FluidSystem::enableVaporizedOil()},
            {"RVW", UnitSystem::measure::oil_gas_ratio, FluidSystem::enableVaporizedWater()},
            {"SOMAX", UnitSystem::measure::identity, simulator_.problem().vapparsActive(simulator_.episodeIndex())},
            {"PCSWM_OW", UnitSystem::measure::identity, enableHysteresis},
            {"KRNSW_OW", UnitSystem::measure::identity, enableHysteresis},
            {"PCSWM_GO", UnitSystem::measure::identity, enableHysteresis},
            {"KRNSW_GO", UnitSystem::measure::identity, enableHysteresis},
            {"PPCW", UnitSystem::measure::pressure, enableSwatinit}
        };

        const auto& inputThpres = eclState().getSimulationConfig().getThresholdPressure();
        std::vector<RestartKey> extraKeys = {{"OPMEXTRA", UnitSystem::measure::identity, false},
                                             {"THRESHPR", UnitSystem::measure::pressure, inputThpres.active()}};

        {
            const auto& tracers = simulator_.vanguard().eclState().tracer();
            for (const auto& tracer : tracers)
                solutionKeys.emplace_back(tracer.fname(), UnitSystem::measure::identity, true);
        }

        // The episodeIndex is rewined one back before beginRestart is called
        // and can not be used here.
        // We just ask the initconfig directly to be sure that we use the correct
        // index.
        const auto& initconfig = simulator_.vanguard().eclState().getInitConfig();
        int restartStepIdx = initconfig.getRestartStep();

        const auto& gridView = simulator_.vanguard().gridView();
        unsigned numElements = gridView.size(/*codim=*/0);
        eclOutputModule_->allocBuffers(numElements, restartStepIdx, /*isSubStep=*/false, /*log=*/false, /*isRestart*/ true);

        {
            SummaryState& summaryState = simulator_.vanguard().summaryState();
            Action::State& actionState = simulator_.vanguard().actionState();
            auto restartValues = loadParallelRestart(this->eclIO_.get(), actionState, summaryState, solutionKeys, extraKeys,
                                                     gridView.grid().comm());
            for (unsigned elemIdx = 0; elemIdx < numElements; ++elemIdx) {
                unsigned globalIdx = this->collectToIORank_.localIdxToGlobalIdx(elemIdx);
                eclOutputModule_->setRestart(restartValues.solution, elemIdx, globalIdx);
            }

            auto& tracer_model = simulator_.problem().tracerModel();
            for (int tracer_index = 0; tracer_index < tracer_model.numTracers(); tracer_index++) {
                const auto& tracer_name = tracer_model.fname(tracer_index);
                const auto& tracer_solution = restartValues.solution.template data<double>(tracer_name);
                for (unsigned elemIdx = 0; elemIdx < numElements; ++elemIdx) {
                    unsigned globalIdx = this->collectToIORank_.localIdxToGlobalIdx(elemIdx);
                    tracer_model.setTracerConcentration(tracer_index, globalIdx, tracer_solution[globalIdx]);
                }
            }

            if (inputThpres.active()) {
                Simulator& mutableSimulator = const_cast<Simulator&>(simulator_);
                auto& thpres = mutableSimulator.problem().thresholdPressure();
                const auto& thpresValues = restartValues.getExtra("THRESHPR");
                thpres.setFromRestart(thpresValues);
            }
            restartTimeStepSize_ = restartValues.getExtra("OPMEXTRA")[0];

            // initialize the well model from restart values
            simulator_.problem().wellModel().initFromRestartFile(restartValues);

            if (!restartValues.aquifer.empty())
                simulator_.problem().mutableAquiferModel().initFromRestart(restartValues.aquifer);
        }
    }

    void endRestart()
    {}

    const EclOutputBlackOilModule<TypeTag>& eclOutputModule() const
    { return *eclOutputModule_; }

    EclOutputBlackOilModule<TypeTag>& mutableEclOutputModule() const
    { return *eclOutputModule_; }

    Scalar restartTimeStepSize() const
    { return restartTimeStepSize_; }

    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(*eclOutputModule_);
    }

private:
    static bool enableEclOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput); }

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

        if (this->eclOutputModule_->localDataValid()) {
            return;
        }

        const auto& gridView = simulator_.vanguard().gridView();
        const bool log = this->collectToIORank_.isIORank();

        const int num_interior = detail::
            countLocalInteriorCellsGridView(gridView);
        this->eclOutputModule_->
            allocBuffers(num_interior, reportStepNum,
                         isSubStep, log, /*isRestart*/ false);

        ElementContext elemCtx(simulator_);

        OPM_BEGIN_PARALLEL_TRY_CATCH();

        {
            OPM_TIMEBLOCK(prepareCellBasedData);
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                this->eclOutputModule_->processElement(elemCtx);
            }
        }

        if constexpr (enableMech) {
            if (simulator_.vanguard().eclState().runspec().mech()) {
                OPM_TIMEBLOCK(prepareMechData);
                for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                    elemCtx.updatePrimaryStencil(elem);
                    elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                    eclOutputModule_->processElementMech(elemCtx);
                }
            }
        }

        if (! this->simulator_.model().linearizer().getFlowsInfo().empty()) {
            OPM_TIMEBLOCK(prepareFlowsData);
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                this->eclOutputModule_->processElementFlows(elemCtx);
            }
        }

        {
            OPM_TIMEBLOCK(prepareBlockData);
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                this->eclOutputModule_->processElementBlockData(elemCtx);
            }
        }

        {
            OPM_TIMEBLOCK(prepareFluidInPlace);

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int dofIdx = 0; dofIdx < num_interior; ++dofIdx) {
                const auto& intQuants = *simulator_.model().cachedIntensiveQuantities(dofIdx, /*timeIdx=*/0);
                const auto totVolume = simulator_.model().dofTotalVolume(dofIdx);

                this->eclOutputModule_->updateFluidInPlace(dofIdx, intQuants, totVolume);
            }
        }

        this->eclOutputModule_->validateLocalData();

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

        this->eclOutputModule_->initializeFluxData();

        OPM_BEGIN_PARALLEL_TRY_CATCH();

        for (const auto& elem : elements(gridView, Dune::Partitions::interiorBorder)) {
            elemCtx.updateStencil(elem);
            elemCtx.updateIntensiveQuantities(timeIdx);
            elemCtx.updateExtensiveQuantities(timeIdx);

            this->eclOutputModule_->processFluxes(elemCtx, activeIndex, cartesianIndex);
        }

        OPM_END_PARALLEL_TRY_CATCH("EclWriter::captureLocalFluxData() failed: ",
                                   this->simulator_.vanguard().grid().comm())

        this->eclOutputModule_->finalizeFluxData();
    }

    Simulator& simulator_;
    std::unique_ptr<EclOutputBlackOilModule<TypeTag> > eclOutputModule_;
    Scalar restartTimeStepSize_;
    int rank_ ;
};

} // namespace Opm

#endif
