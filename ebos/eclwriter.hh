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
 * \copydoc Ewoms::EclWriter
 */
#ifndef EWOMS_ECL_WRITER_HH
#define EWOMS_ECL_WRITER_HH

#include "collecttoiorank.hh"
#include "ecloutputblackoilmodule.hh"

#include <ewoms/disc/ecfv/ecfvdiscretization.hh>
#include <ewoms/io/baseoutputwriter.hh>
#include <ewoms/parallel/tasklets.hh>

#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/parser/eclipse/Units/UnitSystem.hpp>

#include <opm/grid/GridHelpers.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <list>
#include <utility>
#include <string>

BEGIN_PROPERTIES

NEW_PROP_TAG(EnableEclOutput);
NEW_PROP_TAG(EnableAsyncEclOutput);
NEW_PROP_TAG(EclOutputDoublePrecision);

END_PROPERTIES

namespace Ewoms {

template <class TypeTag>
class EclWriter;

template <class TypeTag>
class EclOutputBlackOilModule;

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
class EclWriter
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Vanguard) Vanguard;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef CollectDataToIORank<Vanguard> CollectDataToIORankType;

    typedef std::vector<Scalar> ScalarBuffer;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

public:
    static void registerParameters()
    {
        EclOutputBlackOilModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableAsyncEclOutput,
                             "Write the ECL-formated results in a non-blocking way (i.e., using a separate thread).");
    }

    EclWriter(const Simulator& simulator)
        : simulator_(simulator)
        , collectToIORank_(simulator_.vanguard())
        , eclOutputModule_(simulator, collectToIORank_)
    {
        globalGrid_ = simulator_.vanguard().grid();
        globalGrid_.switchToGlobalView();
        eclIO_.reset(new Opm::EclipseIO(simulator_.vanguard().eclState(),
                                        Opm::UgGridHelpers::createEclipseGrid(globalGrid_, simulator_.vanguard().eclState().getInputGrid()),
                                        simulator_.vanguard().schedule(),
                                        simulator_.vanguard().summaryConfig()));

        // create output thread if enabled and rank is I/O rank
        // async output is enabled by default if pthread are enabled
        bool enableAsyncOutput = EWOMS_GET_PARAM(TypeTag, bool, EnableAsyncEclOutput);
        int numWorkerThreads = 0;
        if (enableAsyncOutput && collectToIORank_.isIORank())
            numWorkerThreads = 1;
        taskletRunner_.reset(new TaskletRunner(numWorkerThreads));
    }

    ~EclWriter()
    { }

    const Opm::EclipseIO& eclIO() const
    { return *eclIO_; }

    void writeInit()
    {
        if (collectToIORank_.isIORank()) {
            std::map<std::string, std::vector<int> > integerVectors;
            if (collectToIORank_.isParallel())
                integerVectors.emplace("MPI_RANK", collectToIORank_.globalRanks());
            eclIO_->writeInitial(computeTrans_(), integerVectors, exportNncStructure_());
        }
    }

    /*!
     * \brief collect and pass data and pass it to eclIO writer
     */
    void writeOutput(bool isSubStep)
    {
        Scalar curTime = simulator_.time() + simulator_.timeStepSize();
        Scalar totalCpuTime =
            simulator_.executionTimer().realTimeElapsed() +
            simulator_.setupTimer().realTimeElapsed() +
            simulator_.vanguard().externalSetupTime();
        Scalar nextStepSize = simulator_.problem().nextTimeStepSize();

        // output using eclWriter if enabled
        Opm::data::Wells localWellData = simulator_.problem().wellModel().wellData();

        int episodeIdx = simulator_.episodeIndex() + 1;
        const auto& gridView = simulator_.vanguard().gridView();
        int numElements = gridView.size(/*codim=*/0);
        bool log = collectToIORank_.isIORank();
        eclOutputModule_.allocBuffers(numElements, episodeIdx, isSubStep, log);

        ElementContext elemCtx(simulator_);
        ElementIterator elemIt = gridView.template begin</*codim=*/0>();
        const ElementIterator& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            eclOutputModule_.processElement(elemCtx);
        }
        eclOutputModule_.outputErrorLog();

        // collect all data to I/O rank and assign to sol
        Opm::data::Solution localCellData = {};
        if (!isSubStep)
            eclOutputModule_.assignToSolution(localCellData);

        // add cell data to perforations for Rft output
        if (!isSubStep)
            eclOutputModule_.addRftDataToWells(localWellData, episodeIdx);

        if (collectToIORank_.isParallel())
            collectToIORank_.collect(localCellData, eclOutputModule_.getBlockData(), localWellData);

        std::map<std::string, double> miscSummaryData;
        std::map<std::string, std::vector<double>> regionData;
        eclOutputModule_.outputFipLog(miscSummaryData, regionData, isSubStep);

        // write output on I/O rank
        if (collectToIORank_.isIORank()) {
            const auto& eclState = simulator_.vanguard().eclState();
            const auto& simConfig = eclState.getSimulationConfig();

            // Add TCPU
            if (totalCpuTime != 0.0)
                miscSummaryData["TCPU"] = totalCpuTime;

            bool enableDoublePrecisionOutput = EWOMS_GET_PARAM(TypeTag, bool, EclOutputDoublePrecision);
            const Opm::data::Solution& cellData = collectToIORank_.isParallel() ? collectToIORank_.globalCellData() : localCellData;
            const Opm::data::Wells& wellData = collectToIORank_.isParallel() ? collectToIORank_.globalWellData() : localWellData;
            Opm::RestartValue restartValue(cellData, wellData);

            const std::map<std::pair<std::string, int>, double>& blockData
                = collectToIORank_.isParallel()
                ? collectToIORank_.globalBlockData()
                : eclOutputModule_.getBlockData();

            // Add suggested next timestep to extra data.
            if (!isSubStep)
                restartValue.addExtra("OPMEXTRA", std::vector<double>(1, nextStepSize));

            if (simConfig.useThresholdPressure())
                restartValue.addExtra("THRESHPR", Opm::UnitSystem::measure::pressure, simulator_.problem().thresholdPressure().data());

            // first, create a tasklet to write the data for the current time step to disk
            auto eclWriteTasklet = std::make_shared<EclWriteTasklet>(*eclIO_,
                                                                     episodeIdx,
                                                                     isSubStep,
                                                                     curTime,
                                                                     restartValue,
                                                                     miscSummaryData,
                                                                     regionData,
                                                                     blockData,
                                                                     enableDoublePrecisionOutput);

            // then, make sure that the previous I/O request has been completed and the
            // number of incomplete tasklets does not increase between time steps
            taskletRunner_->barrier();

            // finally, start a new output writing job
            taskletRunner_->dispatch(eclWriteTasklet);
        }
    }

    void beginRestart()
    {
        bool enableHysteresis = simulator_.problem().materialLawManager()->enableHysteresis();
        bool enableSwatinit = simulator_.vanguard().eclState().get3DProperties().hasDeckDoubleGridProperty("SWATINIT");
        std::vector<Opm::RestartKey> solutionKeys{
            {"PRESSURE", Opm::UnitSystem::measure::pressure},
            {"SWAT", Opm::UnitSystem::measure::identity, static_cast<bool>(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))},
            {"SGAS", Opm::UnitSystem::measure::identity, static_cast<bool>(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))},
            {"TEMP" , Opm::UnitSystem::measure::temperature, enableEnergy},
            {"RS", Opm::UnitSystem::measure::gas_oil_ratio, FluidSystem::enableDissolvedGas()},
            {"RV", Opm::UnitSystem::measure::oil_gas_ratio, FluidSystem::enableVaporizedOil()},
            {"SOMAX", Opm::UnitSystem::measure::identity, simulator_.problem().vapparsActive()},
            {"PCSWM_OW", Opm::UnitSystem::measure::identity, enableHysteresis},
            {"KRNSW_OW", Opm::UnitSystem::measure::identity, enableHysteresis},
            {"PCSWM_GO", Opm::UnitSystem::measure::identity, enableHysteresis},
            {"KRNSW_GO", Opm::UnitSystem::measure::identity, enableHysteresis},
            {"PPCW", Opm::UnitSystem::measure::pressure, enableSwatinit}
        };

        const auto& inputThpres = eclState().getSimulationConfig().getThresholdPressure();
        std::vector<Opm::RestartKey> extraKeys = {{"OPMEXTRA", Opm::UnitSystem::measure::identity, false},
                                                  {"THRESHPR", Opm::UnitSystem::measure::pressure, inputThpres.active()}};

        unsigned episodeIdx = simulator_.episodeIndex();
        const auto& gridView = simulator_.vanguard().gridView();
        unsigned numElements = gridView.size(/*codim=*/0);
        eclOutputModule_.allocBuffers(numElements, episodeIdx, /*isSubStep=*/false, /*log=*/false);

        auto restartValues = eclIO_->loadRestart(solutionKeys, extraKeys);
        for (unsigned elemIdx = 0; elemIdx < numElements; ++elemIdx) {
            unsigned globalIdx = collectToIORank_.localIdxToGlobalIdx(elemIdx);
            eclOutputModule_.setRestart(restartValues.solution, elemIdx, globalIdx);
        }
        if (inputThpres.active()) {
            Simulator& mutableSimulator = const_cast<Simulator&>(simulator_);
            auto& thpres = mutableSimulator.problem().thresholdPressure();
            const auto& thpresValues = restartValues.getExtra("THRESHPR");
            thpres.setFromRestart(thpresValues);
        }
        restartTimeStepSize_ = restartValues.getExtra("OPMEXTRA")[0];
    }

    void endRestart()
    {}

    const EclOutputBlackOilModule<TypeTag>& eclOutputModule() const
    { return eclOutputModule_; }

    Scalar restartTimeStepSize() const
    { return restartTimeStepSize_; }


private:
    static bool enableEclOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput); }

    Opm::data::Solution computeTrans_() const
    {
        const auto& cartMapper = simulator_.vanguard().cartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();
        const int globalSize = cartDims[0]*cartDims[1]*cartDims[2];

        Opm::data::CellData tranx = {Opm::UnitSystem::measure::transmissibility, std::vector<double>(globalSize), Opm::data::TargetType::INIT};
        Opm::data::CellData trany = {Opm::UnitSystem::measure::transmissibility, std::vector<double>(globalSize), Opm::data::TargetType::INIT};
        Opm::data::CellData tranz = {Opm::UnitSystem::measure::transmissibility, std::vector<double>(globalSize), Opm::data::TargetType::INIT};

        for (size_t i = 0; i < tranx.data.size(); ++i) {
            tranx.data[0] = 0.0;
            trany.data[0] = 0.0;
            tranz.data[0] = 0.0;
        }

        const auto& globalGridView = globalGrid_.leafGridView();
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> ElementMapper;
        ElementMapper globalElemMapper(globalGridView, Dune::mcmgElementLayout());
#else
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout> ElementMapper;
        ElementMapper globalElemMapper(globalGridView);
#endif

        const auto& cartesianCellIdx = globalGrid_.globalCell();
        const auto* globalTrans = &(simulator_.vanguard().globalTransmissibility());

        if (!collectToIORank_.isParallel())
            // in the sequential case we must use the transmissibilites defined by
            // the problem. (because in the sequential case, the grid manager does
            // not compute "global" transmissibilities for performance reasons. in
            // the parallel case, the problem's transmissibilities can't be used
            // because this object refers to the distributed grid and we need the
            // sequential version here.)
            globalTrans = &simulator_.problem().eclTransmissibilities();

        auto elemIt = globalGridView.template begin</*codim=*/0>();
        const auto& elemEndIt = globalGridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            const auto& elem = *elemIt;

            auto isIt = globalGridView.ibegin(elem);
            const auto& isEndIt = globalGridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                const auto& is = *isIt;

                if (!is.neighbor())
                    continue; // intersection is on the domain boundary

                unsigned c1 = globalElemMapper.index(is.inside());
                unsigned c2 = globalElemMapper.index(is.outside());

                if (c1 > c2)
                    continue; // we only need to handle each connection once, thank you.


                int gc1 = std::min(cartesianCellIdx[c1], cartesianCellIdx[c2]);
                int gc2 = std::max(cartesianCellIdx[c1], cartesianCellIdx[c2]);
                if (gc2 - gc1 == 1)
                    tranx.data[gc1] = globalTrans->transmissibility(c1, c2);

                if (gc2 - gc1 == cartDims[0])
                    trany.data[gc1] = globalTrans->transmissibility(c1, c2);

                if (gc2 - gc1 == cartDims[0]*cartDims[1])
                    tranz.data[gc1] = globalTrans->transmissibility(c1, c2);
            }
        }

        return {{"TRANX", tranx},
                {"TRANY", trany},
                {"TRANZ", tranz}};
    }

    Opm::NNC exportNncStructure_() const
    {
        Opm::NNC nnc = eclState().getInputNNC();
        int nx = eclState().getInputGrid().getNX();
        int ny = eclState().getInputGrid().getNY();

        const auto& globalGridView = globalGrid_.leafGridView();
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView> ElementMapper;
        ElementMapper globalElemMapper(globalGridView, Dune::mcmgElementLayout());

#else
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout> ElementMapper;
        ElementMapper globalElemMapper(globalGridView);
#endif

        const auto* globalTrans = &(simulator_.vanguard().globalTransmissibility());
        if (!collectToIORank_.isParallel()) {
            // in the sequential case we must use the transmissibilites defined by
            // the problem. (because in the sequential case, the grid manager does
            // not compute "global" transmissibilities for performance reasons. in
            // the parallel case, the problem's transmissibilities can't be used
            // because this object refers to the distributed grid and we need the
            // sequential version here.)
            globalTrans = &simulator_.problem().eclTransmissibilities();
        }

        auto elemIt = globalGridView.template begin</*codim=*/0>();
        const auto& elemEndIt = globalGridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            const auto& elem = *elemIt;

            auto isIt = globalGridView.ibegin(elem);
            const auto& isEndIt = globalGridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                const auto& is = *isIt;

                if (!is.neighbor())
                    continue; // intersection is on the domain boundary

                unsigned c1 = globalElemMapper.index(is.inside());
                unsigned c2 = globalElemMapper.index(is.outside());

                if (c1 > c2)
                    continue; // we only need to handle each connection once, thank you.

                // TODO (?): use the cartesian index mapper to make this code work
                // with grids other than Dune::CpGrid. The problem is that we need
                // the a mapper for the sequential grid, not for the distributed one.
                int cc1 = globalGrid_.globalCell()[c1];
                int cc2 = globalGrid_.globalCell()[c2];

                if (std::abs(cc1 - cc2) != 1 &&
                    std::abs(cc1 - cc2) != nx &&
                    std::abs(cc1 - cc2) != nx*ny)
                {
                    nnc.addNNC(cc1, cc2, globalTrans->transmissibility(c1, c2));
                }
            }
        }
        return nnc;
    }

    struct EclWriteTasklet
        : public TaskletInterface
    {
        Opm::EclipseIO& eclIO_;
        int episodeIdx_;
        bool isSubStep_;
        double secondsElapsed_;
        Opm::RestartValue restartValue_;
        std::map<std::string, double> singleSummaryValues_;
        std::map<std::string, std::vector<double>> regionSummaryValues_;
        std::map<std::pair<std::string, int>, double> blockSummaryValues_;
        bool writeDoublePrecision_;

        explicit EclWriteTasklet(Opm::EclipseIO& eclIO,
                                 int episodeIdx,
                                 bool isSubStep,
                                 double secondsElapsed,
                                 Opm::RestartValue restartValue,
                                 const std::map<std::string, double>& singleSummaryValues,
                                 const std::map<std::string, std::vector<double>>& regionSummaryValues,
                                 const std::map<std::pair<std::string, int>, double>& blockSummaryValues,
                                 bool writeDoublePrecision)
            : eclIO_(eclIO)
            , episodeIdx_(episodeIdx)
            , isSubStep_(isSubStep)
            , secondsElapsed_(secondsElapsed)
            , restartValue_(restartValue)
            , singleSummaryValues_(singleSummaryValues)
            , regionSummaryValues_(regionSummaryValues)
            , blockSummaryValues_(blockSummaryValues)
            , writeDoublePrecision_(writeDoublePrecision)
        { }

        // callback to eclIO serial writeTimeStep method
        void run()
        {
            eclIO_.writeTimeStep(episodeIdx_,
                                 isSubStep_,
                                 secondsElapsed_,
                                 restartValue_,
                                 singleSummaryValues_,
                                 regionSummaryValues_,
                                 blockSummaryValues_,
                                 writeDoublePrecision_);
        }
    };

    const Opm::EclipseState& eclState() const
    { return simulator_.vanguard().eclState(); }

    const Simulator& simulator_;
    CollectDataToIORankType collectToIORank_;
    EclOutputBlackOilModule<TypeTag> eclOutputModule_;
    std::unique_ptr<Opm::EclipseIO> eclIO_;
    Grid globalGrid_;
    std::unique_ptr<TaskletRunner> taskletRunner_;
    Scalar restartTimeStepSize_;


};
} // namespace Ewoms

#endif
