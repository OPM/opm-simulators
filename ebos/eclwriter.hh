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

#include <ewoms/models/blackoil/blackoilmodel.hh>
#include <ewoms/disc/ecfv/ecfvdiscretization.hh>
#include <ewoms/io/baseoutputwriter.hh>
#include <ewoms/parallel/tasklets.hh>

#include <ebos/nncsorter.hpp>

#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/parser/eclipse/Units/UnitSystem.hpp>

#include <opm/grid/GridHelpers.hpp>
#include <opm/grid/utility/cartesianToCompressed.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

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
 * \brief Detect whether two cells are direct vertical neighbours.
 *
 * I.e. have the same i and j index and all cartesian cells between them
 * along the vertical column are inactive.
 *
 * \tparam CM The type of the cartesian index mapper.
 * \param cartMapper The mapper onto cartesian indices.
 * \param cartesianToActive The mapping of cartesian indices to active indices.
 * \param smallGlobalIndex The cartesian cell index of the cell with smaller index
 * \param largeGlobalIndex The cartesian cell index of the cell with larger index
 * \return True if the cells have the same i and j indices and all cartesian cells
 *         between them are inactive.
 */
inline
bool directVerticalNeighbors(const std::array<int, 3>& cartDims,
                             const std::unordered_map<int,int>& cartesianToActive,
                             int smallGlobalIndex, int largeGlobalIndex)
{
    assert(smallGlobalIndex <= largeGlobalIndex);
    std::array<int, 3> ijk1, ijk2;
    auto globalToIjk = [cartDims](int gc) {
                           std::array<int, 3> ijk;
                           ijk[0] = gc % cartDims[0];
                           gc /= cartDims[0];
                           ijk[1] = gc % cartDims[1];
                           ijk[2] = gc / cartDims[1];
                           return ijk;
                       };
    ijk1 = globalToIjk(smallGlobalIndex);
    ijk2 = globalToIjk(largeGlobalIndex);
    assert(ijk2[2]>=ijk1[2]);

    if ( ijk1[0] == ijk2[0] && ijk1[1] == ijk2[1] && (ijk2[2] - ijk1[2]) > 1)
    {
        assert((largeGlobalIndex-smallGlobalIndex)%(cartDims[0]*cartDims[1])==0);
        for ( int gi = smallGlobalIndex + cartDims[0] * cartDims[1]; gi < largeGlobalIndex;
              gi += cartDims[0] * cartDims[1] )
        {
            if ( cartesianToActive.find( gi ) != cartesianToActive.end() )
            {
                return false;
            }
        }
        return true;
    } else
        return false;
}

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

    // The Simulator object should preferably have been const - the
    // only reason that is not the case is due to the SummaryState
    // object owned deep down by the vanguard.
    EclWriter(Simulator& simulator)
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
            auto cartMap = Opm::cartesianToCompressed(globalGrid_.size(0),
                                                      Opm::UgGridHelpers::globalCell(globalGrid_));
            eclIO_->writeInitial(computeTrans_(cartMap), integerVectors, exportNncStructure_(cartMap));
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

        int reportStepNum = simulator_.episodeIndex() + 1;
        const auto& gridView = simulator_.vanguard().gridView();
        int numElements = gridView.size(/*codim=*/0);
        bool log = collectToIORank_.isIORank();
        eclOutputModule_.allocBuffers(numElements, reportStepNum, isSubStep, log);

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
            eclOutputModule_.addRftDataToWells(localWellData, reportStepNum);

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
            auto eclWriteTasklet = std::make_shared<EclWriteTasklet>(summaryState(),
                                                                     eclState,
                                                                     schedule(),
                                                                     *eclIO_,
                                                                     reportStepNum,
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

        // The episodeIndex is rewined one back before beginRestart is called
        // and can not be used here.
        // We just ask the initconfig directly to be sure that we use the correct
        // index.
        const auto& initconfig = simulator_.vanguard().eclState().getInitConfig();
        int restartStepIdx = initconfig.getRestartStep();

        const auto& gridView = simulator_.vanguard().gridView();
        unsigned numElements = gridView.size(/*codim=*/0);
        eclOutputModule_.allocBuffers(numElements, restartStepIdx, /*isSubStep=*/false, /*log=*/false);

        {
            /*
              When running a restarted simulation the restart file is loaded
              twice, first here as part of the state initialization and then
              subsequently in the Simulator::run() method. The global
              SummaryState instance is accumulates total variables like FOPT, if
              the same instance is used twice when loading the restart file, the
              cumulatives will be counted doubly, we therefor use a temporary
              SummaryState instance in this call to loadRestart().
            */
            Opm::SummaryState summaryState;
            auto restartValues = eclIO_->loadRestart(summaryState, solutionKeys, extraKeys);

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

    Opm::data::Solution computeTrans_(const std::unordered_map<int,int>& cartesianToActive) const
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

                // Ordering of compressed and uncompressed index should be the same
                assert(cartesianCellIdx[c1] <= cartesianCellIdx[c2]);
                int gc1 = std::min(cartesianCellIdx[c1], cartesianCellIdx[c2]);
                int gc2 = std::max(cartesianCellIdx[c1], cartesianCellIdx[c2]);

                if (gc2 - gc1 == 1) {
                    tranx.data[gc1] = globalTrans->transmissibility(c1, c2);
                    continue; // skip other if clauses as they are false, last one needs some computation
                }

                if (gc2 - gc1 == cartDims[0]) {
                    trany.data[gc1] = globalTrans->transmissibility(c1, c2);
                    continue; // skipt next if clause as it needs some computation
                }

                if ( gc2 - gc1 == cartDims[0]*cartDims[1] ||
                     directVerticalNeighbors(cartDims, cartesianToActive, gc1, gc2))
                    tranz.data[gc1] = globalTrans->transmissibility(c1, c2);
            }
        }

        return {{"TRANX", tranx},
                {"TRANY", trany},
                {"TRANZ", tranz}};
    }

    Opm::NNC exportNncStructure_(const std::unordered_map<int,int>& cartesianToActive) const
    {
        std::size_t nx = eclState().getInputGrid().getNX();
        std::size_t ny = eclState().getInputGrid().getNY();
        auto nncData = sortNncAndApplyEditnnc(eclState().getInputNNC().nncdata(),
                                              eclState().getInputEDITNNC().data());
        const auto& unitSystem = simulator_.vanguard().deck().getActiveUnitSystem();
        std::vector<Opm::NNCdata> outputNnc;
        std::size_t index = 0;

        for( const auto& entry : nncData ) {
            // test whether NNC is not a neighboring connection
            // cell2>=cell1 holds due to sortNncAndApplyEditnnc
            assert( entry.cell2 >= entry.cell1 );
            auto cellDiff = entry.cell2 - entry.cell1;

            if (cellDiff != 1 && cellDiff != nx && cellDiff != nx*ny) {
                auto tt = unitSystem.from_si(Opm::UnitSystem::measure::transmissibility, entry.trans);
                // Eclipse ignores NNCs (with EDITNNC applied) that are small. Seems like the threshold is 1.0e-6
                if ( tt >= 1.0e-6 )
                    outputNnc.emplace_back(entry.cell1, entry.cell2, entry.trans);
            }
            ++index;
        }

        auto nncCompare =  []( const Opm::NNCdata& nnc1, const Opm::NNCdata& nnc2){
                               return nnc1.cell1 < nnc2.cell1 ||
                                      ( nnc1.cell1 == nnc2.cell1 && nnc1.cell2 < nnc2.cell2);};
        // Sort the nncData values from the deck as they need to be
        // Checked when writing NNC transmissibilities from the simulation.
        std::sort(nncData.begin(), nncData.end(), nncCompare);

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

        auto cartDims = simulator_.vanguard().cartesianIndexMapper().cartesianDimensions();
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
                std::size_t cc1 = globalGrid_.globalCell()[c1];
                std::size_t cc2 = globalGrid_.globalCell()[c2];

                if ( cc2 < cc1 )
                    std::swap(cc1, cc2);

                auto cellDiff = cc2 - cc1;

                if (cellDiff != 1 &&
                    cellDiff != nx &&
                    cellDiff != nx*ny &&
                    ! directVerticalNeighbors(cartDims, cartesianToActive, cc1, cc2)) {
                    // We need to check whether an NNC for this face was also specified
                    // via the NNC keyword in the deck (i.e. in the first origNncSize entries.
                    auto t = globalTrans->transmissibility(c1, c2);
                    auto candidate = std::lower_bound(nncData.begin(), nncData.end(), Opm::NNCdata(cc1, cc2, 0.0), nncCompare);

                    while ( candidate != nncData.end() && candidate->cell1 == cc1
                         && candidate->cell2 == cc2) {
                        t -= candidate->trans;
                        ++candidate;
                    }
                    // eclipse ignores NNCs with zero transmissibility (different threshold than for NNC
                    // with corresponding EDITNNC above). In addition we do set small transmissibilties
                    // to zero when setting up the simulator. These will be ignored here, too.
                    auto tt = unitSystem.from_si(Opm::UnitSystem::measure::transmissibility, std::abs(t));
                    if ( tt > 1e-12 )
                        outputNnc.push_back({cc1, cc2, t});
                }
            }
        }
        Opm::NNC ret;
        for(const auto& nncItem: outputNnc)
            ret.addNNC(nncItem.cell1, nncItem.cell2, nncItem.trans);
        return ret;
    }

    struct EclWriteTasklet
        : public TaskletInterface
    {
        Opm::SummaryState& summaryState;
        const Opm::EclipseState& eclState;
        const Opm::Schedule& schedule;
        Opm::EclipseIO& eclIO_;
        int reportStepNum_;
        bool isSubStep_;
        double secondsElapsed_;
        Opm::RestartValue restartValue_;
        std::map<std::string, double> singleSummaryValues_;
        std::map<std::string, std::vector<double>> regionSummaryValues_;
        std::map<std::pair<std::string, int>, double> blockSummaryValues_;
        bool writeDoublePrecision_;

        explicit EclWriteTasklet(Opm::SummaryState& summaryState,
                                 const Opm::EclipseState& eclState,
                                 const Opm::Schedule& schedule,
                                 Opm::EclipseIO& eclIO,
                                 int reportStepNum,
                                 bool isSubStep,
                                 double secondsElapsed,
                                 Opm::RestartValue restartValue,
                                 const std::map<std::string, double>& singleSummaryValues,
                                 const std::map<std::string, std::vector<double>>& regionSummaryValues,
                                 const std::map<std::pair<std::string, int>, double>& blockSummaryValues,
                                 bool writeDoublePrecision)
            : summaryState(summaryState)
            , eclState(eclState)
            , schedule(schedule)
            , eclIO_(eclIO)
            , reportStepNum_(reportStepNum)
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
            const auto& summary = eclIO_.summary();

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
            if (reportStepNum_ > 0)
                summary.eval(summaryState,
                             reportStepNum_,
                             secondsElapsed_,
                             eclState,
                             schedule,
                             restartValue_.wells,
                             singleSummaryValues_,
                             regionSummaryValues_,
                             blockSummaryValues_);

            eclIO_.writeTimeStep(summaryState,
                                 reportStepNum_,
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

    Opm::SummaryState& summaryState()
    { return simulator_.vanguard().summaryState(); }

    const Opm::Schedule& schedule() const
    { return simulator_.vanguard().schedule(); }

    Simulator& simulator_;
    CollectDataToIORankType collectToIORank_;
    EclOutputBlackOilModule<TypeTag> eclOutputModule_;
    std::unique_ptr<Opm::EclipseIO> eclIO_;
    Grid globalGrid_;
    std::unique_ptr<TaskletRunner> taskletRunner_;
    Scalar restartTimeStepSize_;


};
} // namespace Ewoms

#endif
