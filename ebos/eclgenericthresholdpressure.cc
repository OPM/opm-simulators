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
#include <ebos/eclgenericthresholdpressure.hh>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Eqldims.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/SimulationConfig.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/ThresholdPressure.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/polyhedralgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/3d/gridview.hh>
#include "alucartesianindexmapper.hh"
#endif // HAVE_DUNE_ALUGRID

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#include <ebos/femcpgridcompat.hh>
#endif // HAVE_DUNE_FEM

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <stdexcept>

namespace Opm {

template<class Grid, class GridView, class ElementMapper, class Scalar>
EclGenericThresholdPressure<Grid,GridView,ElementMapper,Scalar>::
EclGenericThresholdPressure(const CartesianIndexMapper& cartMapper,
                            const GridView& gridView,
                            const ElementMapper& elementMapper,
                            const EclipseState& eclState)
    : cartMapper_(cartMapper)
    , gridView_(gridView)
    , elementMapper_(elementMapper)
    , eclState_(eclState)
{
}

template<class Grid, class GridView, class ElementMapper,class Scalar>
Scalar EclGenericThresholdPressure<Grid,GridView,ElementMapper,Scalar>::
thresholdPressure(int elem1Idx, int elem2Idx) const
{
    if (!enableThresholdPressure_)
        return 0.0;

    // threshold pressure accross faults
    if (!thpresftValues_.empty()) {
        int cartElem1Idx = cartMapper_.cartesianIndex(elem1Idx);
        int cartElem2Idx = cartMapper_.cartesianIndex(elem2Idx);

        assert(0 <= cartElem1Idx && static_cast<int>(cartElemFaultIdx_.size()) > cartElem1Idx);
        assert(0 <= cartElem2Idx && static_cast<int>(cartElemFaultIdx_.size()) > cartElem2Idx);

        int fault1Idx = cartElemFaultIdx_[cartElem1Idx];
        int fault2Idx = cartElemFaultIdx_[cartElem2Idx];
        if (fault1Idx != -1 && fault1Idx == fault2Idx)
            // inside a fault there's no threshold pressure, even accross EQUIL
            // regions.
            return 0.0;
        if (fault1Idx != fault2Idx) {
            // TODO: which value if a cell is part of multiple faults? we take
            // the maximum here.
            Scalar val1 = (fault1Idx >= 0) ? thpresftValues_[fault1Idx] : 0.0;
            Scalar val2 = (fault2Idx >= 0) ? thpresftValues_[fault2Idx] : 0.0;
            return std::max(val1, val2);
        }
    }

    // threshold pressure accross EQUIL regions
    unsigned short equilRegion1Idx = elemEquilRegion_[elem1Idx];
    unsigned short equilRegion2Idx = elemEquilRegion_[elem2Idx];

    if (equilRegion1Idx == equilRegion2Idx)
        return 0.0;

    return thpres_[equilRegion1Idx*numEquilRegions_ + equilRegion2Idx];
}

template<class Grid, class GridView, class ElementMapper, class Scalar>
void EclGenericThresholdPressure<Grid,GridView,ElementMapper,Scalar>::
finishInit()
{
    unsigned numElements = gridView_.size(/*codim=*/0);

    const auto& simConfig = eclState_.getSimulationConfig();

    enableThresholdPressure_ = simConfig.useThresholdPressure();
    if (!enableThresholdPressure_)
        return;

    numEquilRegions_ = eclState_.getTableManager().getEqldims().getNumEquilRegions();
    if (numEquilRegions_ > 0xff) {
        // make sure that the index of an equilibration region can be stored in a
        // single byte
        throw std::runtime_error("The maximum number of supported equilibration regions is 255!");
    }

    // internalize the data specified using the EQLNUM keyword
    const auto& fp = eclState_.fieldProps();
    const auto& equilRegionData = fp.get_int("EQLNUM");
    elemEquilRegion_.resize(numElements, 0);
    for (unsigned elemIdx = 0; elemIdx < numElements; ++elemIdx) {
        elemEquilRegion_[elemIdx] = equilRegionData[elemIdx] - 1;
    }

    /*
      If this is a restart run the ThresholdPressure object will be active,
      but it will *not* be properly initialized with numerical values. The
      values must instead come from the THPRES vector in the restart file.
    */
    if (simConfig.getThresholdPressure().restart())
        return;

    // allocate the array which specifies the threshold pressures
    thpres_.resize(numEquilRegions_*numEquilRegions_, 0.0);
    thpresDefault_.resize(numEquilRegions_*numEquilRegions_, 0.0);
}

template<class Grid, class GridView, class ElementMapper, class Scalar>
void EclGenericThresholdPressure<Grid,GridView,ElementMapper,Scalar>::
applyExplicitThresholdPressures_()
{
    const SimulationConfig& simConfig = eclState_.getSimulationConfig();
    const auto& thpres = simConfig.getThresholdPressure();

    // set the threshold pressures for all EQUIL region boundaries which have a
    // intersection in the grid
    for (const auto& elem : elements(gridView_, Dune::Partitions::interior)) {
        for (const auto& intersection : intersections(gridView_, elem)) {
            if (intersection.boundary())
                continue; // ignore boundary intersections for now (TODO?)
            else if (!intersection.neighbor()) //processor boundary but not domain boundary
                continue;

            const auto& inside = intersection.inside();
            const auto& outside = intersection.outside();

            unsigned insideElemIdx = elementMapper_.index(inside);
            unsigned outsideElemIdx = elementMapper_.index(outside);

            unsigned equilRegionInside = elemEquilRegion_[insideElemIdx];
            unsigned equilRegionOutside = elemEquilRegion_[outsideElemIdx];
            if (thpres.hasRegionBarrier(equilRegionInside + 1, equilRegionOutside + 1)) {
                Scalar pth = 0.0;
                if (thpres.hasThresholdPressure(equilRegionInside + 1, equilRegionOutside + 1)) {
                    // threshold pressure explicitly specified
                    pth = thpres.getThresholdPressure(equilRegionInside + 1, equilRegionOutside + 1);
                }
                else {
                    // take the threshold pressure from the initial condition
                    unsigned offset = equilRegionInside*numEquilRegions_ + equilRegionOutside;
                    pth = thpresDefault_[offset];
                }

                unsigned offset1 = equilRegionInside*numEquilRegions_ + equilRegionOutside;
                unsigned offset2 = equilRegionOutside*numEquilRegions_ + equilRegionInside;

                thpres_[offset1] = pth;
                thpres_[offset2] = pth;
            }
        }
    }

    // apply threshold pressures across faults
    if (thpres.ftSize() > 0)
        configureThpresft_();
}

template<class Grid, class GridView, class ElementMapper, class Scalar>
void EclGenericThresholdPressure<Grid,GridView,ElementMapper,Scalar>::
configureThpresft_()
{
    // retrieve the faults collection.
    const FaultCollection& faults = eclState_.getFaults();
    const SimulationConfig& simConfig = eclState_.getSimulationConfig();
    const auto& thpres = simConfig.getThresholdPressure();

    // extract the multipliers
    int numFaults = faults.size();
    int numCartesianElem = eclState_.getInputGrid().getCartesianSize();
    thpresftValues_.resize(numFaults, -1.0);
    cartElemFaultIdx_.resize(numCartesianElem, -1);
    for (size_t faultIdx = 0; faultIdx < faults.size(); faultIdx++) {
        auto& fault = faults.getFault(faultIdx);
        thpresftValues_[faultIdx] = thpres.getThresholdPressureFault(faultIdx);
        for (const FaultFace& face : fault)
            // "face" is a misnomer because the object describes a set of cell
            // indices, but we go with the conventions of the parser here...
            for (size_t cartElemIdx : face)
                cartElemFaultIdx_[cartElemIdx] = faultIdx;
    }
}

template<class Grid, class GridView, class ElementMapper, class Scalar>
std::vector<Scalar>
EclGenericThresholdPressure<Grid,GridView,ElementMapper,Scalar>::
getRestartVector() const
{
    if (!enableThresholdPressure_)
        return {};

    std::vector<Scalar> result(numEquilRegions_ * numEquilRegions_, 0.0);
    const auto& simConfig = eclState_.getSimulationConfig();
    const auto& thpres = simConfig.getThresholdPressure();

    std::size_t idx = 0;
    for (unsigned j = 1; j <= numEquilRegions_; ++j) {
        for (unsigned i = 1; i <= numEquilRegions_; ++i, ++idx) {
            if (thpres.hasRegionBarrier(i, j)) {
                if (thpres.hasThresholdPressure(i, j)) {
                    result[idx] = thpres.getThresholdPressure(i, j);
                } else {
                    result[idx] = this->thpresDefault_[idx];
                }
            }
        }
    }

    return result;
}

template<class Grid, class GridView, class ElementMapper, class Scalar>
void
EclGenericThresholdPressure<Grid,GridView,ElementMapper,Scalar>::
logPressures()
{
    if (!enableThresholdPressure_)
        return;

    auto lineFormat = [this](unsigned i, unsigned j, double val)
    {
        const auto& units = eclState_.getUnits();
        return fmt::format("{:4}{:>6}{:23}{:>6}{:24}{:>11.07}{:7}{}\n",
                           " ", i,
                           " ", j,
                           " ", units.from_si(UnitSystem::measure::pressure, val),
                           " ", units.name(UnitSystem::measure::pressure));
    };
    auto lineFormatS = [](auto s1, auto s2, auto s3)
    {
        return fmt::format("{:4}{:^16}{:13}{:^9}{:21}{:^18}\n",
                           " ", s1, " ", s2, " ", s3);
    };

    std::string str = "\nLIST OF ALL NON-ZERO THRESHOLD PRESSURES\n"
                      "----------------------------------------\n"
                      "\n";
    str += lineFormatS("FLOW FROM REGION", "TO REGION", "THRESHOLD PRESSURE");
    str += lineFormatS(std::string(16, '-'), std::string(9, '-'), std::string(18, '-'));
    const auto& simConfig = eclState_.getSimulationConfig();
    const auto& thpres = simConfig.getThresholdPressure();

    for (unsigned i = 1; i <= numEquilRegions_; ++i) {
        for (unsigned j = (thpres.irreversible() ? 1 : i); j <= numEquilRegions_; ++j) {
            if (thpres.hasRegionBarrier(i, j)) {
                if (thpres.hasThresholdPressure(i, j)) {
                    str += lineFormat(i, j, thpres.getThresholdPressure(j, i));
                } else {
                    std::size_t idx = (j - 1) * numEquilRegions_ + (i - 1);
                    str += lineFormat(i, j, this->thpresDefault_[idx]);
                }
            }
        }
    }
    str += lineFormatS(std::string(16, '-'), std::string(9, '-'), std::string(18, '-'));
    OpmLog::note(str);
}

#if HAVE_DUNE_FEM
template class EclGenericThresholdPressure<Dune::CpGrid,
                                           Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>,
                                           Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>>,
                                           double>;
template class EclGenericThresholdPressure<Dune::CpGrid,
                                            Dune::Fem::GridPart2GridViewImpl<
                                                Dune::Fem::AdaptiveLeafGridPart<
                                                    Dune::CpGrid,
                                                    Dune::PartitionIteratorType(4),
                                                    false> >,
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

template class EclGenericThresholdPressure<ALUGrid3CN,
                                           Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false>>>,
                                           Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false>>>>,
                                           double>;
template class EclGenericThresholdPressure<ALUGrid3CN,
                                            Dune::Fem::GridPart2GridViewImpl<
                                                Dune::Fem::AdaptiveLeafGridPart<
                                                    ALUGrid3CN,
                                                    Dune::PartitionIteratorType(4),
                                                    false> >,
                                            Dune::MultipleCodimMultipleGeomTypeMapper<
                                                Dune::Fem::GridPart2GridViewImpl<
                                                    Dune::Fem::AdaptiveLeafGridPart<
                                                        ALUGrid3CN,
                                                        Dune::PartitionIteratorType(4),
                                                        false>>>,
                                            double>;                                           
                                           
                                           

#endif //HAVE_DUNE_ALUGRID
#else
template class EclGenericThresholdPressure<Dune::CpGrid,
                                           Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                           Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>,
                                           double>;
#if HAVE_DUNE_ALUGRID

#if HAVE_MPI
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>;
#else
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI
template class EclGenericThresholdPressure<ALUGrid3CN,
                                           Dune::GridView<Dune::ALU3dLeafGridViewTraits<const ALUGrid3CN,Dune::PartitionIteratorType(4)>>,
                                           Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::ALU3dLeafGridViewTraits<const ALUGrid3CN,Dune::PartitionIteratorType(4)>>>,
                                           double>;
#endif //HAVE_DUNE_ALUGRID

#endif

template class EclGenericThresholdPressure<Dune::PolyhedralGrid<3,3,double>,
                                           Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>,
                                           Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>>,
                                           double>;

} // namespace Opm
