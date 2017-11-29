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
 * \copydoc Ewoms::EclThresholdPressure
 */
#ifndef EWOMS_ECL_THRESHOLD_PRESSURE_HH
#define EWOMS_ECL_THRESHOLD_PRESSURE_HH

#include <ewoms/common/propertysystem.hh>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/GridProperty.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Eqldims.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/SimulationConfig.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/ThresholdPressure.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <dune/grid/common/gridenums.hh>
#include <dune/common/version.hh>

#include <array>
#include <vector>
#include <unordered_map>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Evaluation);
NEW_PROP_TAG(ElementContext);
NEW_PROP_TAG(FluidSystem);
}

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This class calculates the threshold pressure for grid faces according to the
 *        Eclipse Reference Manual.
 *
 * If the difference of the pressure potential between two cells is below the threshold
 * pressure, the pressure potential difference is assumed to be zero, if it is larger
 * than the threshold pressure, it is reduced by the threshold pressure.
 */
template <class TypeTag>
class EclThresholdPressure
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { numPhases = FluidSystem::numPhases };

public:
    EclThresholdPressure(const Simulator& simulator)
        : simulator_(simulator)
    {
        enableThresholdPressure_ = false;
    }

    /*!
     * \brief Actually compute the threshold pressures over a face as a pre-compute step.
     */
    void finishInit()
    {
        const auto& gridView = simulator_.gridView();

        unsigned numElements = gridView.size(/*codim=*/0);

        // this code assumes that the DOFs are the elements. (i.e., an
        // ECFV spatial discretization with TPFA). if you try to use
        // it with something else, you're currently out of luck,
        // sorry!
        assert(simulator_.model().numGridDof() == numElements);

        const auto& gridManager = simulator_.gridManager();
        const auto& eclState = gridManager.eclState();
        const auto& simConfig = eclState.getSimulationConfig();

        enableThresholdPressure_ = simConfig.hasThresholdPressure();
        if (!enableThresholdPressure_)
            return;

        numEquilRegions_ = eclState.getTableManager().getEqldims().getNumEquilRegions();
        if (numEquilRegions_ > 0xff) {
            // make sure that the index of an equilibration region can be stored in a
            // single byte
            OPM_THROW(std::runtime_error,
                      "The maximum number of supported equilibration regions is 255!");
        }

        // allocate the array which specifies the threshold pressures
        thpres_.resize(numEquilRegions_*numEquilRegions_, 0.0);
        thpresDefault_.resize(numEquilRegions_*numEquilRegions_, 0.0);

        // internalize the data specified using the EQLNUM keyword
        const std::vector<int>& equilRegionData =
            eclState.get3DProperties().getIntGridProperty("EQLNUM").getData();
        elemEquilRegion_.resize(numElements, 0);
        for (unsigned elemIdx = 0; elemIdx < numElements; ++elemIdx) {
            int cartElemIdx = gridManager.cartesianIndex(elemIdx);

            // ECL uses Fortran-style indices but we want C-style ones!
            elemEquilRegion_[elemIdx] = equilRegionData[cartElemIdx] - 1;
        }

        computeDefaultThresholdPressures_();
        applyExplicitThresholdPressures_();
    }

    /*!
     * \brief Returns the theshold pressure [Pa] for the intersection between two elements.
     *
     * This is tailor made for the E100 threshold pressure mechanism and it is thus quite
     * a hack: First of all threshold pressures in general are unphysical, and second,
     * they should be different for the fluid phase but are not. Anyway, this seems to be
     * E100's way of doing things, so we do it the same way.
     */
    Scalar thresholdPressure(int elemIdx1, int elemIdx2) const
    {
        if (!enableThresholdPressure_)
            return 0.0;

        unsigned short equilRegion1Idx = elemEquilRegion_[elemIdx1];
        unsigned short equilRegion2Idx = elemEquilRegion_[elemIdx2];

        if (equilRegion1Idx == equilRegion2Idx)
            return 0.0;

        return thpres_[equilRegion1Idx*numEquilRegions_ + equilRegion2Idx];
    }

private:
    // compute the defaults of the threshold pressures using the initial condition
    void computeDefaultThresholdPressures_()
    {
        const auto& gridManager = simulator_.gridManager();
        const auto& gridView = gridManager.gridView();

        typedef Opm::MathToolbox<Evaluation> Toolbox;
        // loop over the whole grid and compute the maximum gravity adjusted pressure
        // difference between two EQUIL regions.
        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();
        ElementContext elemCtx(simulator_);
        for (; elemIt != elemEndIt; ++elemIt) {

            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updateAll(elem);
            const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);

            for (unsigned scvfIdx = 0; scvfIdx < stencil.numInteriorFaces(); ++ scvfIdx) {
                const auto& face = stencil.interiorFace(scvfIdx);

                unsigned i = face.interiorIndex();
                unsigned j = face.exteriorIndex();

                unsigned insideElemIdx = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
                unsigned outsideElemIdx = elemCtx.globalSpaceIndex(j, /*timeIdx=*/0);

                unsigned equilRegionInside = elemEquilRegion_[insideElemIdx];
                unsigned equilRegionOutside = elemEquilRegion_[outsideElemIdx];

                if (equilRegionInside == equilRegionOutside)
                    // the current face is not at the boundary between EQUIL regions!
                    continue;

                // don't include connections with negligible flow
                const Scalar& trans = simulator_.problem().transmissibility(elemCtx, i, j);
                const Scalar& faceArea = face.area();
                if ( std::abs(faceArea * trans) < 1e-18)
                    continue;

                // determine the maximum difference of the pressure of any phase over the
                // intersection
                Scalar pth = 0.0;
                const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, /*timeIdx=*/0);
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    unsigned upIdx = extQuants.upstreamIndex(phaseIdx);
                    const auto& up = elemCtx.intensiveQuantities(upIdx, /*timeIdx=*/0);

                    if (up.mobility(phaseIdx) > 0.0) {
                        Scalar phaseVal = Toolbox::value(extQuants.pressureDifference(phaseIdx));
                        pth = std::max(pth, std::abs(phaseVal));
                    }
                }

                int offset1 = equilRegionInside*numEquilRegions_ + equilRegionOutside;
                int offset2 = equilRegionOutside*numEquilRegions_ + equilRegionInside;

                thpresDefault_[offset1] = std::max(thpresDefault_[offset1], pth);
                thpresDefault_[offset2] = std::max(thpresDefault_[offset2], pth);
            }
        }

        // make sure that the threshold pressures is consistent for parallel
        // runs. (i.e. take the maximum of all processes)
        for (unsigned i = 0; i < thpresDefault_.size(); ++i)
            thpresDefault_[i] = gridView.comm().max(thpresDefault_[i]);
    }

    // internalize the threshold pressures which where explicitly specified via the
    // THPRES keyword.
    void applyExplicitThresholdPressures_()
    {
        const auto& gridManager = simulator_.gridManager();
        const auto& gridView = gridManager.gridView();
        const auto& elementMapper = simulator_.model().elementMapper();
        const auto& eclState = simulator_.gridManager().eclState();
        const Opm::SimulationConfig& simConfig = eclState.getSimulationConfig();
        const auto& thpres = simConfig.getThresholdPressure();

        // set the threshold pressures for all EQUIL region boundaries which have a
        // intersection in the grid
        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue;

            auto isIt = gridView.ibegin(elem);
            const auto& isEndIt = gridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                // store intersection, this might be costly
                const auto& intersection = *isIt;

                // ignore boundary intersections for now (TODO?)
                if (intersection.boundary())
                    continue;

                const auto& inside = intersection.inside();
                const auto& outside = intersection.outside();

                unsigned insideElemIdx = elementMapper.index(inside);
                unsigned outsideElemIdx = elementMapper.index(outside);

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
    }

    const Simulator& simulator_;

    std::vector<Scalar> thpresDefault_;
    std::vector<Scalar> thpres_;
    unsigned numEquilRegions_;
    std::vector<unsigned char> elemEquilRegion_;

    bool enableThresholdPressure_;
};

} // namespace Ewoms

#endif
