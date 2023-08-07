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
 * \copydoc Opm::EclThresholdPressure
 */
#ifndef EWOMS_ECL_THRESHOLD_PRESSURE_HH
#define EWOMS_ECL_THRESHOLD_PRESSURE_HH

#include <ebos/eclgenericthresholdpressure.hh>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <algorithm>
#include <vector>

namespace Opm {

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
class EclThresholdPressure : public EclGenericThresholdPressure<GetPropType<TypeTag, Properties::Grid>,
                                                                GetPropType<TypeTag, Properties::GridView>,
                                                                GetPropType<TypeTag, Properties::ElementMapper>,
                                                                GetPropType<TypeTag, Properties::Scalar>>
{
    using BaseType = EclGenericThresholdPressure<GetPropType<TypeTag, Properties::Grid>,
                                                 GetPropType<TypeTag, Properties::GridView>,
                                                 GetPropType<TypeTag, Properties::ElementMapper>,
                                                 GetPropType<TypeTag, Properties::Scalar>>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    enum { numPhases = FluidSystem::numPhases };

public:
    EclThresholdPressure(const Simulator& simulator)
        : BaseType(simulator.vanguard().cartesianIndexMapper(),
                   simulator.vanguard().gridView(),
                   simulator.model().elementMapper(),
                   simulator.vanguard().eclState())
        , simulator_(simulator)
    {
    }

    /*!
     * \brief Actually compute the threshold pressures over a face as a pre-compute step.
     */
    void finishInit()
    {
        this->BaseType::finishInit();
        if (this->enableThresholdPressure_ && !this->thpresDefault_.empty()) {
            this->computeDefaultThresholdPressures_();
            this->applyExplicitThresholdPressures_();
        }
    }

private:
    // compute the defaults of the threshold pressures using the initial condition
    void computeDefaultThresholdPressures_()
    {
        const auto& vanguard = simulator_.vanguard();
        const auto& gridView = vanguard.gridView();

        using Toolbox = MathToolbox<Evaluation>;
        // loop over the whole grid and compute the maximum gravity adjusted pressure
        // difference between two EQUIL regions.
        ElementContext elemCtx(simulator_);
        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            elemCtx.updateAll(elem);
            const auto& stencil = elemCtx.stencil(/*timeIdx=*/0);

            for (unsigned scvfIdx = 0; scvfIdx < stencil.numInteriorFaces(); ++ scvfIdx) {
                const auto& face = stencil.interiorFace(scvfIdx);

                unsigned i = face.interiorIndex();
                unsigned j = face.exteriorIndex();

                unsigned insideElemIdx = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
                unsigned outsideElemIdx = elemCtx.globalSpaceIndex(j, /*timeIdx=*/0);

                unsigned equilRegionInside = this->elemEquilRegion_[insideElemIdx];
                unsigned equilRegionOutside = this->elemEquilRegion_[outsideElemIdx];

                if (equilRegionInside == equilRegionOutside)
                    // the current face is not at the boundary between EQUIL regions!
                    continue;

                // don't include connections with negligible flow
                const Evaluation& trans = simulator_.problem().transmissibility(elemCtx, i, j);
                Scalar faceArea = face.area();
                if (std::abs(faceArea*getValue(trans)) < 1e-18)
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

                int offset1 = equilRegionInside*this->numEquilRegions_ + equilRegionOutside;
                int offset2 = equilRegionOutside*this->numEquilRegions_ + equilRegionInside;

                this->thpresDefault_[offset1] = std::max(this->thpresDefault_[offset1], pth);
                this->thpresDefault_[offset2] = std::max(this->thpresDefault_[offset2], pth);
            }
        }

        // make sure that the threshold pressures is consistent for parallel
        // runs. (i.e. take the maximum of all processes)
        for (unsigned i = 0; i < this->thpresDefault_.size(); ++i)
            this->thpresDefault_[i] = gridView.comm().max(this->thpresDefault_[i]);

        this->logPressures();
    }

    const Simulator& simulator_;
};

} // namespace Opm

#endif
