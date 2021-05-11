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

#include "collecttoiorank.hh"
#include "ecloutputblackoilmodule.hh"

#include <opm/models/blackoil/blackoilmodel.hh>

#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/io/baseoutputwriter.hh>
#include <opm/models/parallel/tasklets.hh>

#include <opm/output/eclipse/EclipseIO.hpp>

#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/output/eclipse/Summary.hpp>
#include <opm/parser/eclipse/Units/UnitSystem.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/State.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/PAvgCalculatorCollection.hpp>

#include <opm/simulators/utils/ParallelRestart.hpp>
#include <opm/grid/GridHelpers.hpp>
#include <opm/grid/utility/cartesianToCompressed.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <list>
#include <utility>
#include <string>
#include <chrono>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

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

} // namespace Opm::Properties

namespace Opm {

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
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;

    using CollectDataToIORankType = CollectDataToIORank<Grid,EquilGrid,GridView>;

    typedef std::vector<Scalar> ScalarBuffer;

    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };


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
        , collectToIORank_(simulator_.vanguard().grid(),
                           simulator_.vanguard().grid().comm().rank() == 0 ?
                           &simulator_.vanguard().equilGrid() : nullptr,
                           simulator_.vanguard().gridView(),
                           simulator_.vanguard().cartesianIndexMapper(),
                           simulator_.vanguard().grid().comm().rank() == 0 ?
                           &simulator_.vanguard().equilCartesianIndexMapper() : nullptr)

    {
        std::vector<std::size_t> wbp_index_list;
        if (collectToIORank_.isIORank()) {
            const auto& schedule = simulator_.vanguard().schedule();
            eclIO_.reset(new EclipseIO(simulator_.vanguard().eclState(),
                                       UgGridHelpers::createEclipseGrid(globalGrid(), simulator_.vanguard().eclState().getInputGrid()),
                                       schedule,
                                       simulator_.vanguard().summaryConfig()));

            const auto& wbp_calculators = eclIO_->summary().wbp_calculators( schedule.size() - 1 );
            wbp_index_list = wbp_calculators.index_list();
        }
        if (collectToIORank_.isParallel()) {
            const auto& comm = simulator_.vanguard().grid().comm();
            unsigned long size = wbp_index_list.size();
            comm.broadcast(&size, 1, collectToIORank_.ioRank);
            if (!collectToIORank_.isIORank())
                wbp_index_list.resize( size );
            comm.broadcast(wbp_index_list.data(), size, collectToIORank_.ioRank);
        }
        // create output thread if enabled and rank is I/O rank
        // async output is enabled by default if pthread are enabled
        bool enableAsyncOutput = EWOMS_GET_PARAM(TypeTag, bool, EnableAsyncEclOutput);
        int numWorkerThreads = 0;
        if (enableAsyncOutput && collectToIORank_.isIORank())
            numWorkerThreads = 1;
        taskletRunner_.reset(new TaskletRunner(numWorkerThreads));

        this->eclOutputModule_ = std::make_unique<EclOutputBlackOilModule<TypeTag>>(simulator, wbp_index_list, this->collectToIORank_);
    }

    ~EclWriter()
    { }

    const EclipseIO& eclIO() const
    {
        assert(eclIO_);
        return *eclIO_;
    }

    const EquilGrid& globalGrid() const
    {
        return simulator_.vanguard().equilGrid();
    }

    void writeInit()
    {
        if (collectToIORank_.isIORank()) {
            std::map<std::string, std::vector<int> > integerVectors;
            if (collectToIORank_.isParallel())
                integerVectors.emplace("MPI_RANK", collectToIORank_.globalRanks());
            auto cartMap = cartesianToCompressed(globalGrid().size(0),
                                                 UgGridHelpers::globalCell(globalGrid()));
            eclIO_->writeInitial(computeTrans_(cartMap), integerVectors, exportNncStructure_(cartMap));
        }
    }

    /*!
     * \brief collect and pass data and pass it to eclIO writer
     */

    void evalSummaryState(bool isSubStep)
    {
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
            simulator_.vanguard().externalSetupTime();

        const auto localWellData            = simulator_.problem().wellModel().wellData();
        const auto localGroupAndNetworkData = simulator_.problem().wellModel()
            .groupAndNetworkData(reportStepNum, simulator_.vanguard().schedule());


        const auto localAquiferData = simulator_.problem().mutableAquiferModel().aquiferData();

        this->prepareLocalCellData(isSubStep, reportStepNum);

        if (collectToIORank_.isParallel())
            collectToIORank_.collect({},
                                     eclOutputModule_->getBlockData(),
                                     eclOutputModule_->getWBPData(),
                                     localWellData,
                                     localGroupAndNetworkData,
                                     localAquiferData);


        std::map<std::string, double> miscSummaryData;
        std::map<std::string, std::vector<double>> regionData;
        auto inplace = eclOutputModule_->outputFipLog(miscSummaryData, regionData, isSubStep);

        bool forceDisableProdOutput = false;
        bool forceDisableInjOutput = false;
        bool forceDisableCumOutput = false;
        eclOutputModule_->outputProdLog(reportStepNum, isSubStep, forceDisableProdOutput);
        eclOutputModule_->outputInjLog(reportStepNum, isSubStep, forceDisableInjOutput);
        eclOutputModule_->outputCumLog(reportStepNum, isSubStep, forceDisableCumOutput);


        std::vector<char> buffer;
        if (this->collectToIORank_.isIORank()) {
            const auto& summary = eclIO_->summary();
            auto wbp_calculators = summary.wbp_calculators(reportStepNum);
            const auto& wbpData
                = this->collectToIORank_.isParallel()
                ? this->collectToIORank_.globalWBPData()
                : this->eclOutputModule_->getWBPData();

            for (const auto& [global_index, pressure] : wbpData)
                wbp_calculators.add_pressure( global_index, pressure );

            // Add TCPU
            if (totalCpuTime != 0.0) {
                miscSummaryData["TCPU"] = totalCpuTime;
            }

            const auto& wellData = this->collectToIORank_.isParallel()
                ? this->collectToIORank_.globalWellData()
                : localWellData;

            const auto& groupAndNetworkData = this->collectToIORank_.isParallel()
                ? this->collectToIORank_.globalGroupAndNetworkData()
                : localGroupAndNetworkData;

            const auto& aquiferData = this->collectToIORank_.isParallel()
                ? this->collectToIORank_.globalAquiferData()
                : localAquiferData;

            const auto& blockData
                = this->collectToIORank_.isParallel()
                ? this->collectToIORank_.globalBlockData()
                : this->eclOutputModule_->getBlockData();

            summary.eval(summaryState(),
                         reportStepNum,
                         curTime,
                         wellData,
                         groupAndNetworkData,
                         miscSummaryData,
                         eclOutputModule_->initialInplace(),
                         inplace,
                         wbp_calculators,
                         regionData,
                         blockData,
                         aquiferData);

            /*
              Off-by-one-fun: The reportStepNum argument corresponds to the
              report step these results will be written to, whereas the argument
              to UDQ function evaluation corresponds to the report step we are
              currently on.
            */
            auto udq_step = reportStepNum - 1;
            const auto& udq_config = schedule().getUDQConfig(udq_step);
            udq_config.eval( udq_step, schedule().wellMatcher(udq_step), summaryState(), udqState() );

            buffer = summaryState().serialize();
        }

        if (collectToIORank_.isParallel()) {
#ifdef HAVE_MPI
            unsigned long buffer_size = buffer.size();
            MPI_Bcast(&buffer_size, 1, MPI_UNSIGNED_LONG, collectToIORank_.ioRank, MPI_COMM_WORLD);
            if (!collectToIORank_.isIORank())
                buffer.resize( buffer_size );

            MPI_Bcast(buffer.data(), buffer_size, MPI_CHAR, collectToIORank_.ioRank, MPI_COMM_WORLD);
            if (!collectToIORank_.isIORank()) {
                SummaryState& st = summaryState();
                st.deserialize(buffer);
            }
#endif
        }
    }


    void writeOutput(bool isSubStep)
    {
        const int reportStepNum = simulator_.episodeIndex() + 1;

        this->prepareLocalCellData(isSubStep, reportStepNum);
        this->eclOutputModule_->outputErrorLog();

        // output using eclWriter if enabled
        auto localWellData = simulator_.problem().wellModel().wellData();
        auto localGroupAndNetworkData = simulator_.problem().wellModel()
            .groupAndNetworkData(reportStepNum, simulator_.vanguard().schedule());

        data::Solution localCellData = {};
        if (! isSubStep) {
            this->eclOutputModule_->assignToSolution(localCellData);

            // add cell data to perforations for Rft output
            this->eclOutputModule_->addRftDataToWells(localWellData, reportStepNum);
        }

        if (this->collectToIORank_.isParallel()) {
            collectToIORank_.collect(localCellData,
                                     eclOutputModule_->getBlockData(),
                                     eclOutputModule_->getWBPData(),
                                     localWellData,
                                     localGroupAndNetworkData,
                                     {});
        }

        if (this->collectToIORank_.isIORank()) {
            this->writeOutput(reportStepNum, isSubStep,
                              std::move(localCellData),
                              std::move(localWellData),
                              std::move(localGroupAndNetworkData));
        }
    }

    void beginRestart()
    {
        bool enableHysteresis = simulator_.problem().materialLawManager()->enableHysteresis();
        bool enableSwatinit = simulator_.vanguard().eclState().fieldProps().has_double("SWATINIT");
        std::vector<RestartKey> solutionKeys{
            {"PRESSURE", UnitSystem::measure::pressure},
            {"SWAT", UnitSystem::measure::identity, static_cast<bool>(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))},
            {"SGAS", UnitSystem::measure::identity, static_cast<bool>(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))},
            {"TEMP" , UnitSystem::measure::temperature, enableEnergy},
            {"SSOLVENT" , UnitSystem::measure::identity, enableSolvent},
            {"RS", UnitSystem::measure::gas_oil_ratio, FluidSystem::enableDissolvedGas()},
            {"RV", UnitSystem::measure::oil_gas_ratio, FluidSystem::enableVaporizedOil()},
            {"SOMAX", UnitSystem::measure::identity, simulator_.problem().vapparsActive()},
            {"PCSWM_OW", UnitSystem::measure::identity, enableHysteresis},
            {"KRNSW_OW", UnitSystem::measure::identity, enableHysteresis},
            {"PCSWM_GO", UnitSystem::measure::identity, enableHysteresis},
            {"KRNSW_GO", UnitSystem::measure::identity, enableHysteresis},
            {"PPCW", UnitSystem::measure::pressure, enableSwatinit}
        };

        const auto& inputThpres = eclState().getSimulationConfig().getThresholdPressure();
        std::vector<RestartKey> extraKeys = {{"OPMEXTRA", UnitSystem::measure::identity, false},
                                             {"THRESHPR", UnitSystem::measure::pressure, inputThpres.active()}};

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
            auto restartValues = loadParallelRestart(eclIO_.get(), actionState, summaryState, solutionKeys, extraKeys,
                                                     gridView.grid().comm());
            for (unsigned elemIdx = 0; elemIdx < numElements; ++elemIdx) {
                unsigned globalIdx = collectToIORank_.localIdxToGlobalIdx(elemIdx);
                eclOutputModule_->setRestart(restartValues.solution, elemIdx, globalIdx);
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

    Scalar restartTimeStepSize() const
    { return restartTimeStepSize_; }


private:
    static bool enableEclOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput); }

    data::Solution computeTrans_(const std::unordered_map<int,int>& cartesianToActive) const
    {
        const auto& cartMapper = simulator_.vanguard().equilCartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();
        const int globalSize = cartDims[0]*cartDims[1]*cartDims[2];

        data::CellData tranx = {UnitSystem::measure::transmissibility, std::vector<double>(globalSize), data::TargetType::INIT};
        data::CellData trany = {UnitSystem::measure::transmissibility, std::vector<double>(globalSize), data::TargetType::INIT};
        data::CellData tranz = {UnitSystem::measure::transmissibility, std::vector<double>(globalSize), data::TargetType::INIT};

        for (size_t i = 0; i < tranx.data.size(); ++i) {
            tranx.data[0] = 0.0;
            trany.data[0] = 0.0;
            tranz.data[0] = 0.0;
        }

        typedef typename EquilGrid :: LeafGridView  GlobalGridView;
        const GlobalGridView& globalGridView = globalGrid().leafGridView();
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView> ElementMapper;
        ElementMapper globalElemMapper(globalGridView, Dune::mcmgElementLayout());

        const EclTransmissibility<TypeTag>* globalTrans;

        if (!collectToIORank_.isParallel())
        {
            // in the sequential case we must use the transmissibilites defined by
            // the problem. (because in the sequential case, the grid manager does
            // not compute "global" transmissibilities for performance reasons. in
            // the parallel case, the problem's transmissibilities can't be used
            // because this object refers to the distributed grid and we need the
            // sequential version here.)
            globalTrans = &simulator_.problem().eclTransmissibilities();
        }
        else
        {
            globalTrans = &(simulator_.vanguard().globalTransmissibility());
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

                // Ordering of compressed and uncompressed index should be the same
                const int cartIdx1 = cartMapper.cartesianIndex( c1 );
                const int cartIdx2 = cartMapper.cartesianIndex( c2 );
                // Ordering of compressed and uncompressed index should be the same
                assert(cartIdx1 <= cartIdx2);
                int gc1 = std::min(cartIdx1, cartIdx2);
                int gc2 = std::max(cartIdx1, cartIdx2);

                if (gc2 - gc1 == 1 && cartDims[0] > 1 ) {
                    tranx.data[gc1] = globalTrans->transmissibility(c1, c2);
                    continue; // skip other if clauses as they are false, last one needs some computation
                }

                if (gc2 - gc1 == cartDims[0] && cartDims[1] > 1) {
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

    std::vector<NNCdata> exportNncStructure_(const std::unordered_map<int,int>& cartesianToActive) const
    {
        std::size_t nx = eclState().getInputGrid().getNX();
        std::size_t ny = eclState().getInputGrid().getNY();
        auto nncData = eclState().getInputNNC().input();
        const auto& unitSystem = simulator_.vanguard().eclState().getDeckUnitSystem();
        std::vector<NNCdata> outputNnc;
        std::size_t index = 0;

        for( const auto& entry : nncData ) {
            // test whether NNC is not a neighboring connection
            // cell2>=cell1 holds due to sortNncAndApplyEditnnc
            assert( entry.cell2 >= entry.cell1 );
            auto cellDiff = entry.cell2 - entry.cell1;

            if (cellDiff != 1 && cellDiff != nx && cellDiff != nx*ny) {
                auto tt = unitSystem.from_si(UnitSystem::measure::transmissibility, entry.trans);
                // Eclipse ignores NNCs (with EDITNNC applied) that are small. Seems like the threshold is 1.0e-6
                if ( tt >= 1.0e-6 )
                    outputNnc.emplace_back(entry.cell1, entry.cell2, entry.trans);
            }
            ++index;
        }

        typedef typename EquilGrid :: LeafGridView  GlobalGridView;
        const GlobalGridView& globalGridView = globalGrid().leafGridView();
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView> ElementMapper;
        ElementMapper globalElemMapper(globalGridView, Dune::mcmgElementLayout());

        const EclTransmissibility<TypeTag>* globalTrans;
        if (!collectToIORank_.isParallel()) {
            // in the sequential case we must use the transmissibilites defined by
            // the problem. (because in the sequential case, the grid manager does
            // not compute "global" transmissibilities for performance reasons. in
            // the parallel case, the problem's transmissibilities can't be used
            // because this object refers to the distributed grid and we need the
            // sequential version here.)
            globalTrans = &simulator_.problem().eclTransmissibilities();
        }
        else
        {
            globalTrans = &(simulator_.vanguard().globalTransmissibility());
        }

        // Cartesian index mapper for the serial I/O grid
        const auto& equilCartMapper =  simulator_.vanguard().equilCartesianIndexMapper();
        const auto& cartDims = simulator_.vanguard().cartesianIndexMapper().cartesianDimensions();
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

                std::size_t cc1 = equilCartMapper.cartesianIndex( c1 ); //globalIOGrid_.globalCell()[c1];
                std::size_t cc2 = equilCartMapper.cartesianIndex( c2 ); //globalIOGrid_.globalCell()[c2];

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
                    auto candidate = std::lower_bound(nncData.begin(), nncData.end(), NNCdata(cc1, cc2, 0.0));

                    while ( candidate != nncData.end() && candidate->cell1 == cc1
                         && candidate->cell2 == cc2) {
                        t -= candidate->trans;
                        ++candidate;
                    }
                    // eclipse ignores NNCs with zero transmissibility (different threshold than for NNC
                    // with corresponding EDITNNC above). In addition we do set small transmissibilties
                    // to zero when setting up the simulator. These will be ignored here, too.
                    auto tt = unitSystem.from_si(UnitSystem::measure::transmissibility, std::abs(t));
                    if ( tt > 1e-12 )
                        outputNnc.push_back({cc1, cc2, t});
                }
            }
        }
        return outputNnc;
    }

    struct EclWriteTasklet
        : public TaskletInterface
    {
        Action::State actionState_;
        SummaryState summaryState_;
        UDQState udqState_;
        EclipseIO& eclIO_;
        int reportStepNum_;
        bool isSubStep_;
        double secondsElapsed_;
        RestartValue restartValue_;
        bool writeDoublePrecision_;

        explicit EclWriteTasklet(const Action::State& actionState,
                                 const SummaryState& summaryState,
                                 const UDQState& udqState,
                                 EclipseIO& eclIO,
                                 int reportStepNum,
                                 bool isSubStep,
                                 double secondsElapsed,
                                 RestartValue restartValue,
                                 bool writeDoublePrecision)
            : actionState_(actionState)
            , summaryState_(summaryState)
            , udqState_(udqState)
            , eclIO_(eclIO)
            , reportStepNum_(reportStepNum)
            , isSubStep_(isSubStep)
            , secondsElapsed_(secondsElapsed)
            , restartValue_(restartValue)
            , writeDoublePrecision_(writeDoublePrecision)
        { }

        // callback to eclIO serial writeTimeStep method
        void run()
        {
            eclIO_.writeTimeStep(actionState_,
                                 summaryState_,
                                 udqState_,
                                 reportStepNum_,
                                 isSubStep_,
                                 secondsElapsed_,
                                 restartValue_,
                                 writeDoublePrecision_);
        }
    };

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
        const auto& gridView = simulator_.vanguard().gridView();
        const int numElements = gridView.size(/*codim=*/0);
        const bool log = collectToIORank_.isIORank();

        eclOutputModule_->allocBuffers(numElements, reportStepNum,
                                      isSubStep, log, /*isRestart*/ false);

        ElementContext elemCtx(simulator_);
        ElementIterator elemIt = gridView.template begin</*codim=*/0>();

        const ElementIterator& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            eclOutputModule_->processElement(elemCtx);
        }
    }

    void writeOutput(const int                            reportStepNum,
                     const bool                           isSubStep,
                     ::Opm::data::Solution&&              localCellData,
                     ::Opm::data::Wells&&                 localWellData,
                     ::Opm::data::GroupAndNetworkValues&& localGroupAndNetworkData)
    {
        const Scalar curTime = simulator_.time() + simulator_.timeStepSize();
        const Scalar nextStepSize = simulator_.problem().nextTimeStepSize();
        const auto isParallel = this->collectToIORank_.isParallel();

        RestartValue restartValue {
            isParallel ? this->collectToIORank_.globalCellData()
                       : std::move(localCellData),

            isParallel ? this->collectToIORank_.globalWellData()
                       : std::move(localWellData),

            isParallel ? this->collectToIORank_.globalGroupAndNetworkData()
                       : std::move(localGroupAndNetworkData)
        };

        if (simulator_.vanguard().eclState().getSimulationConfig().useThresholdPressure()) {
            restartValue.addExtra("THRESHPR", UnitSystem::measure::pressure,
                                  simulator_.problem().thresholdPressure().data());
        }

        // Add suggested next timestep to extra data.
        if (! isSubStep) {
            restartValue.addExtra("OPMEXTRA", std::vector<double>(1, nextStepSize));
        }

        // first, create a tasklet to write the data for the current time
        // step to disk
        auto eclWriteTasklet = std::make_shared<EclWriteTasklet>(
            this->actionState(), this->summaryState(), this->udqState(), *this->eclIO_,
            reportStepNum, isSubStep, curTime, std::move(restartValue),
            EWOMS_GET_PARAM(TypeTag, bool, EclOutputDoublePrecision)
            );

        // then, make sure that the previous I/O request has been completed
        // and the number of incomplete tasklets does not increase between
        // time steps
        this->taskletRunner_->barrier();

        // finally, start a new output writing job
        this->taskletRunner_->dispatch(std::move(eclWriteTasklet));
    }

    Simulator& simulator_;
    CollectDataToIORankType collectToIORank_;
    std::unique_ptr<EclOutputBlackOilModule<TypeTag>> eclOutputModule_;
    std::unique_ptr<EclipseIO> eclIO_;
    std::unique_ptr<TaskletRunner> taskletRunner_;
    Scalar restartTimeStepSize_;
};
} // namespace Opm

#endif
