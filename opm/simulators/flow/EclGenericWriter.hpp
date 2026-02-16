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
#ifndef OPM_ECL_GENERIC_WRITER_HPP
#define OPM_ECL_GENERIC_WRITER_HPP

#include <opm/models/parallel/tasklets.hpp>

#include <opm/simulators/flow/CollectDataOnIORank.hpp>
#include <opm/simulators/flow/Transmissibility.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace Opm {

class EclipseIO;
class EclipseState;
class InterRegFlowMap;
class Inplace;
template <class Grid> class LevelCartesianIndexMapper;
struct NNCdata;
class Schedule;
class SummaryConfig;
class SummaryState;
class UDQState;

} // namespace Opm

namespace Opm { namespace Action {
class State;
}} // namespace Opm::Action

namespace Opm {

template <class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
class EclGenericWriter
{
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using EquilCartesianIndexMapper = Dune::CartesianIndexMapper<EquilGrid>;
    using CollectDataOnIORankType = CollectDataOnIORank<Grid,EquilGrid,GridView>;
    using TransmissibilityType = Transmissibility<Grid,GridView,ElementMapper,CartesianIndexMapper,Scalar>;

public:
    EclGenericWriter(const Schedule& schedule,
                     const EclipseState& eclState,
                     const SummaryConfig& summaryConfig,
                     const Grid& grid,
                     const EquilGrid* equilGrid,
                     const GridView& gridView,
                     const Dune::CartesianIndexMapper<Grid>& cartMapper,
                     const Dune::CartesianIndexMapper<EquilGrid>* equilCartMapper,
                     bool enableAsyncOutput,
                     bool enableEsmry);

    const EclipseIO& eclIO() const;

    void writeInit();

    void setTransmissibilities(const TransmissibilityType* globalTrans)
    {
        globalTrans_ = globalTrans;
    }

    void setSubStepReport(const SimulatorReportSingle& report)
    {
        sub_step_report_ = report;
    }
    void setSimulationReport(const SimulatorReport& report)
    {
        simulation_report_ = report;
    }

    const std::vector<NNCdata>& getOutputNnc() const
    {
        return outputNnc_;
    }

    const CollectDataOnIORankType& collectOnIORank() const
    {
        return collectOnIORank_;
    }

    void extractOutputTransAndNNC(const std::function<unsigned int(unsigned int)>& map);

protected:
    const TransmissibilityType& globalTrans() const;
    unsigned int gridEquilIdxToGridIdx(unsigned int elemIndex) const;

    void doWriteOutput(const int                          reportStepNum,
                       const std::optional<int>           timeStepNum,
                       const bool                         isSubStep,
                       data::Solution&&                   localCellData,
                       data::Wells&&                      localWellData,
                       data::GroupAndNetworkValues&&      localGroupAndNetworkData,
                       data::Aquifers&&                   localAquiferData,
                       WellTestState&&                    localWTestState,
                       const Action::State&               actionState,
                       const UDQState&                    udqState,
                       const SummaryState&                summaryState,
                       const std::vector<Scalar>&         thresholdPressure,
                       Scalar                             curTime,
                       Scalar                             nextStepSize,
                       bool                               doublePrecision,
                       bool                               isFlowsn,
                       std::array<FlowsData<double>,3>&&  flowsn,
                       bool                               isFloresn,
                       std::array<FlowsData<double>, 3>&& floresn);

    void evalSummary(int                                                  reportStepNum,
                     Scalar                                               curTime,
                     const data::Wells&                                   localWellData,
                     const data::WellBlockAveragePressures&               localWBPData,
                     const data::GroupAndNetworkValues&                   localGroupAndNetworkData,
                     const std::map<int,data::AquiferData>&               localAquiferData,
                     const std::map<std::pair<std::string, int>, double>& blockData,
                     const std::map<std::string, double>&                 miscSummaryData,
                     const std::map<std::string, std::vector<double>>&    regionData,
                     const Inplace&                                       inplace,
                     const Inplace*                                       initialInPlace,
                     const InterRegFlowMap&                               interRegFlows,
                     SummaryState&                                        summaryState,
                     UDQState&                                            udqState);

    CollectDataOnIORankType collectOnIORank_;
    const Grid& grid_;
    const GridView& gridView_;
    const Schedule& schedule_;
    const EclipseState& eclState_;
    std::unique_ptr<EclipseIO> eclIO_;
    std::unique_ptr<TaskletRunner> taskletRunner_;
    Scalar restartTimeStepSize_;
    const TransmissibilityType* globalTrans_ = nullptr;
    const Dune::CartesianIndexMapper<Grid>& cartMapper_;
    const Dune::CartesianIndexMapper<EquilGrid>* equilCartMapper_;
    const EquilGrid* equilGrid_;
    SimulatorReportSingle sub_step_report_;
    SimulatorReport simulation_report_;
    mutable std::vector<NNCdata> outputNnc_; 
    
    // Regular NNCs per grid: internal to a grid.
    // Both cells belong to the same level grid, either the main grid or a level/local grid. 
    // nnc.cell1 (NNC1 in *.EGRID) Level/Local Cartesian index of cell1 
    // nnc.cell2 (NNC2 in *.EGRID) Level/Local Cartesian index of cell2
    // Equivalent to TRANNNC in *.INIT
     mutable std::vector<std::vector<NNCdata>> outputInnerLevelNnc_;

    // NNCs between main (level zero) grid and local grids:
    // nnc.cell1 (NNCG in *.EGRID) Cartesian index of cell1 in the main grid where the cell belongs to. 
    // nnc.cell2 (NNCL in *.EGRID) Level/Local Cartesian index of cell2 in the refined level grid where the cell belongs to.
    // Equivalent to TRANGL in *.INIT
    mutable std::vector<std::vector<NNCdata>> outputNncGlobalLocal_; // here GLOBAL refers to level 0 grid, local to any LGR (refined grid)

    // Amalgamated NNCs: nncs between different LGRs. For example, nested refinement or neighboring LGRs.
    // The cells belong to different refined level grids. 
    // nnc.cell1 (NNA1 in *.EGRID) Level/Local Cartesian index of cell1 (in its level grid: level1)
    // nnc.cell2 (NNA2 in *.EGRID) Level/Local Cartesian index of cell2 (in its level grid: level2, with level2 > level1).
    // Equivalent to TRANLL in *.INIT
    mutable std::vector<std::vector<std::vector<NNCdata>>> outputAmalgamatedNnc_;
    
    mutable std::unique_ptr<std::vector<data::Solution>> outputTrans_;

private:
    template<typename LevelIndicesFunction>
    void computeTrans_(const Opm::LevelCartesianIndexMapper<EquilGrid>& levelCartMapp,
                       const std::vector<std::unordered_map<int,int>>& levelCartToLevelCompressed,
                       const std::function<unsigned int(unsigned int)>& map,
                       const LevelIndicesFunction& computeLevelIndices,
                       const std::function<int(int, int)>& computeLevelCartIdx) const;

    template<typename LevelIndicesFunction>
    std::vector<NNCdata> exportNncStructure_(const Opm::LevelCartesianIndexMapper<EquilGrid>& levelCartMapp,
                                             const std::unordered_map<int,int>& cartesianToActive,
                                             const std::function<unsigned int(unsigned int)>& map,
                                             const LevelIndicesFunction& computeLevelIndices,
                                             const std::function<int(int, int)>& computeLevelCartIdx) const;

    /// Returns true if the given Cartesian cell index belongs to a numerical aquifer.
    /// If no numerical aquifer exists, always returns false.
    bool isNumAquCell_(const std::size_t cartIdx) const;
    /// Returns true if either of the two connected cells belongs to a numerical aquifer.
    bool isNumAquConn_(const std::size_t cartIdx1, const std::size_t cartIdx2) const;

    /// Create LevelCartesianIndexMapper.
    ///
    /// For CpGrid, the LevelCartesianIndexMapper constructor takes the grid
    /// For other grid types,the CartesianIndexMapper.
    /// Note: the design of the LevelCartesianIndexMapper class and its relation/duplication
    ///       with CartesianIndexMapper will be improved "soon".
    /// @tparam equilGridIsCpGrid Compile-time flag indicating whether the equilGrid_ is a Dune::CpGrid
    ///                           (std::is_same_v<EquilGrid, Dune::CpGrid>).
    template<bool equilGridIsCpGrid>
    Opm::LevelCartesianIndexMapper<EquilGrid> createLevelCartMapp_() const;

    /// Create maps from [level Cartesian index] to [level active/compressed index]
    ///
    /// Note: Refinement is supported only for CpGrid for now. For other grid types,
    ///       there is only one map carteisanToActive since level zero and leaf grids
    ///       coincide.
    /// @tparam equilGridIsCpGrid Compile-time flag indicating whether the equilGrid_ is a Dune::CpGrid
    ///                           (std::is_same_v<EquilGrid, Dune::CpGrid>).
    template<bool equilGridIsCpGrid>
    std::vector<std::unordered_map<int,int>> createCartesianToActiveMaps_(const Opm::LevelCartesianIndexMapper<EquilGrid>& levelCartMapp) const;

    /// Return function to compute (level) Cartesian dimensions
    ///
    /// Note: Refinement is supported only for CpGrid for now.
    /// @tparam equilGridIsCpGrid Compile-time flag indicating whether the equilGrid_ is a Dune::CpGrid
    ///                           (std::is_same_v<EquilGrid, Dune::CpGrid>).
    template<bool equilGridIsCpGrid>
    std::function<std::array<int,3>(int)> computeLevelCartDimensions_(const Opm::LevelCartesianIndexMapper<EquilGrid>& levelCartMapp,
                                                                      const Dune::CartesianIndexMapper<EquilGrid>& equilCartMapp) const;

    /// Return function to compute (level) Cartesian indices of the cells adjacent to an intersection.
    ///
    /// @tparam equilGridIsCpGrid Compile-time flag indicating whether the equilGrid_ is a Dune::CpGrid
    ///                           (std::is_same_v<EquilGrid, Dune::CpGrid>).
    template<bool equilGridIsCpGrid>
    std::function<int(int, int)> computeLevelCartIdx_(const Opm::LevelCartesianIndexMapper<EquilGrid>& levelCartMapp,
                                                      const Dune::CartesianIndexMapper<EquilGrid>& equilCartMapp) const;

    /// Return function to compute level (local/compressed) indices of the cells adjacent to an intersection.
    ///
    /// @tparam equilGridIsCpGrid Compile-time flag indicating whether the equilGrid_ is a Dune::CpGrid
    ///                           (std::is_same_v<EquilGrid, Dune::CpGrid>).
    template <bool equilGridIsCpGrid>
    auto computeLevelIndices_() const;

    /// Creates/allocates CellData for TRANX/Y/Z for level grids.
    /// Only for CpGrid. Other grid types do not support refinement yet.
    void allocateLevelTrans_(const std::array<int,3>& levelCartDims,
                             data::Solution& levelTrans) const;
    /// Allocates NNCdata  vectors for
    /// - TRANNNC for level grids,
    /// - TRANGL for NNCs between level zero grid and an LGR (refined level grid)
    /// - TRANLL for NNCs between two different refined level grids. 
    /// Only for CpGrid. Other grid types do not support refinement yet.
    void allocateAllNncs_(int maxLevel) const;
};

} // namespace Opm

#endif // OPM_ECL_GENERIC_WRITER_HPP
