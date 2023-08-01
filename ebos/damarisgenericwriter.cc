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

#include <config.h>

#ifdef HAVE_DAMARIS

#include <ebos/damarisgenericwriter.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/utility/cartesianToCompressed.hpp>
#if HAVE_DUNE_ALUGRID
#include "eclalugridvanguard.hh"
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/3d/gridview.hh>
#endif // HAVE_DUNE_ALUGRID

// #include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/output/eclipse/Summary.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvgCalculatorCollection.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#if HAVE_MPI
#include <ebos/eclmpiserializer.hh>
#endif

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#include <ebos/femcpgridcompat.hh>
#endif // HAVE_DUNE_FEM

#if HAVE_MPI
#include <mpi.h>
#endif

#include <array>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>


namespace Opm {

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
DamarisGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
DamarisGenericWriter(const Schedule& schedule,
                 const EclipseState& eclState,
                 const SummaryConfig& summaryConfig,
                 const Grid& grid,
                 const EquilGrid* equilGrid,
                 const GridView& gridView,
                 const Dune::CartesianIndexMapper<Grid>& cartMapper,
                 const Dune::CartesianIndexMapper<EquilGrid>* equilCartMapper,
                 bool enableAsyncOutput,
                 bool enableEsmry )
    : collectToIORank_(grid,
                       equilGrid,
                       gridView,
                       cartMapper,
                       equilCartMapper,
                       summaryConfig.fip_regions_interreg_flow())
    , grid_(grid)
    , gridView_(gridView)
    , schedule_(schedule)
    , eclState_(eclState)
    , summaryConfig_(summaryConfig)
    , cartMapper_(cartMapper)
    , equilCartMapper_(equilCartMapper)
    , equilGrid_(equilGrid)
{
 
}


template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void DamarisGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
writeInit(const std::function<unsigned int(unsigned int)>& map)
{

}


template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void DamarisGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
doWriteOutput(const int                     reportStepNum,
              const bool                    isSubStep,
              data::Solution&&              localCellData,
              data::Wells&&                 localWellData,
              data::GroupAndNetworkValues&& localGroupAndNetworkData,
              data::Aquifers&&              localAquiferData,
              WellTestState&&               localWTestState,
              const Action::State& actionState,
              const UDQState& udqState,
              const SummaryState& summaryState,
              const std::vector<Scalar>& thresholdPressure,
              Scalar curTime,
              Scalar nextStepSize,
              bool doublePrecision,
              bool isFlowsn,
              std::array<std::pair<std::string, std::pair<std::vector<int>, std::vector<double>>>, 3>&& flowsn,
              bool isFloresn,
              std::array<std::pair<std::string, std::pair<std::vector<int>, std::vector<double>>>, 3>&& floresn)
{
    const auto isParallel = this->collectToIORank_.isParallel();
    const bool needsReordering = this->collectToIORank_.doesNeedReordering();

    RestartValue restartValue {
        (isParallel || needsReordering) ? this->collectToIORank_.globalCellData()
                   : std::move(localCellData),

        isParallel ? this->collectToIORank_.globalWellData()
                   : std::move(localWellData),

        isParallel ? this->collectToIORank_.globalGroupAndNetworkData()
                   : std::move(localGroupAndNetworkData),

        isParallel ? this->collectToIORank_.globalAquiferData()
                   : std::move(localAquiferData)
    };

    if (eclState_.getSimulationConfig().useThresholdPressure()) {
        restartValue.addExtra("THRESHPR", UnitSystem::measure::pressure,
                              thresholdPressure);
    }

    // Add suggested next timestep to extra data.
    if (! isSubStep) {
        restartValue.addExtra("OPMEXTRA", std::vector<double>(1, nextStepSize));
    }

    // Add nnc flows and flores.
    if (isFlowsn) {
        const auto flowsn_global = isParallel ? this->collectToIORank_.globalFlowsn() : std::move(flowsn);
        for (const auto& flows : flowsn_global) { 
            if (flows.first.empty())
                continue;
            if (flows.first == "FLOGASN+") {
                restartValue.addExtra(flows.first, UnitSystem::measure::gas_surface_rate, flows.second.second);
            }
            else {
                restartValue.addExtra(flows.first, UnitSystem::measure::liquid_surface_rate, flows.second.second);
            }
        }
    }
    if (isFloresn) {
        const auto floresn_global = isParallel ? this->collectToIORank_.globalFloresn() : std::move(floresn);
        for (const auto& flores : floresn_global) {
            if (flores.first.empty())
                continue;
            restartValue.addExtra(flores.first, UnitSystem::measure::rate, flores.second.second);
        }
    }

    /*
    // We do not use tasklet to write the data when offloading to Damaris
    */
}




template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
const typename DamarisGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::TransmissibilityType&
DamarisGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
globalTrans() const
{
    assert (globalTrans_);
    return *globalTrans_;
}

#if HAVE_DUNE_FEM
template class DamarisGenericWriter<Dune::CpGrid,
                                Dune::CpGrid,
                                Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>, Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>>,
                                double>;
template class DamarisGenericWriter<Dune::CpGrid,
                                Dune::CpGrid,
                                Dune::Fem::GridPart2GridViewImpl<
                                    Dune::Fem::AdaptiveLeafGridPart<
                                        Dune::CpGrid,
                                        Dune::PartitionIteratorType(4),
                                        false>>,
                                Dune::MultipleCodimMultipleGeomTypeMapper<
                                    Dune::Fem::GridPart2GridViewImpl<
                                        Dune::Fem::AdaptiveLeafGridPart<
                                            Dune::CpGrid,
                                            Dune::PartitionIteratorType(4),
                                            false>>>,
                                double>;


#if HAVE_DUNE_ALUGRID
#if HAVE_MPI
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>;
#else
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI
                                                               
template class DamarisGenericWriter<ALUGrid3CN,
                                Dune::CpGrid,
                                Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false>>>, Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false>>>>,
                                double>;
                                
template class DamarisGenericWriter<ALUGrid3CN,
                                Dune::CpGrid,
                                Dune::Fem::GridPart2GridViewImpl<
                                    Dune::Fem::AdaptiveLeafGridPart<
                                        ALUGrid3CN,
                                        Dune::PartitionIteratorType(4),
                                        false>>,
                                Dune::MultipleCodimMultipleGeomTypeMapper<
                                    Dune::Fem::GridPart2GridViewImpl<
                                        Dune::Fem::AdaptiveLeafGridPart<
                                            ALUGrid3CN,
                                            Dune::PartitionIteratorType(4),
                                            false>>>,
                                double>;                                
#endif // HAVE_DUNE_ALUGRID

#else // !HAVE_DUNE_FEM
template class DamarisGenericWriter<Dune::CpGrid,
                                Dune::CpGrid,
                                Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>,
                                double>;
#if HAVE_DUNE_ALUGRID
#if HAVE_MPI
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>;
#else
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI
template class DamarisGenericWriter<ALUGrid3CN,
                                Dune::CpGrid,
                                Dune::GridView<Dune::ALU3dLeafGridViewTraits<const ALUGrid3CN,Dune::PartitionIteratorType(4)>>,  Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::ALU3dLeafGridViewTraits<const ALUGrid3CN, Dune::PartitionIteratorType(4)>>>,
                                double>;

#endif // HAVE_DUNE_ALUGRID
#endif // !HAVE_DUNE_FEM

template class DamarisGenericWriter<Dune::PolyhedralGrid<3,3,double>,
                                Dune::PolyhedralGrid<3,3,double>,
                                Dune::GridView<Dune::PolyhedralGridViewTraits<3, 3, double, Dune::PartitionIteratorType(4)>>,                              Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>>,
                                double>;
                                                                  
} // namespace Opm

#endif  // #ifdef HAVE_DAMARIS
