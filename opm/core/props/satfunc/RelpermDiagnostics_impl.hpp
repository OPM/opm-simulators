/*
  Copyright 2016 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_RELPERMDIAGNOSTICS_IMPL_HEADER_INCLUDED
#define OPM_RELPERMDIAGNOSTICS_IMPL_HEADER_INCLUDED

#include <vector>
#include <utility>

#include <opm/material/fluidmatrixinteractions/EclEpsGridProperties.hpp>
#include <opm/core/props/satfunc/RelpermDiagnostics.hpp>
#include <opm/grid/utility/compressedToCartesian.hpp>
#include <opm/grid/GridHelpers.hpp>


namespace Opm {

    template <class CartesianIndexMapper>
    void RelpermDiagnostics::diagnosis(const Opm::EclipseState& eclState,
                                       const CartesianIndexMapper& cartesianIndexMapper)
    {
        OpmLog::info("\n===============Saturation Functions Diagnostics===============\n");
        phaseCheck_(eclState);
        satFamilyCheck_(eclState);
        tableCheck_(eclState);
        unscaledEndPointsCheck_(eclState);
        scaledEndPointsCheck_(eclState, cartesianIndexMapper);
    }

    template <class CartesianIndexMapper>
    void RelpermDiagnostics::scaledEndPointsCheck_(const EclipseState& eclState,
                                                   const CartesianIndexMapper& cartesianIndexMapper)
    {
        // All end points are subject to round-off errors, checks should account for it
        const float tolerance = 1e-6;
        const int nc = cartesianIndexMapper.compressedSize();
        const bool threepoint = eclState.runspec().endpointScaling().threepoint();
        scaledEpsInfo_.resize(nc);
        EclEpsGridProperties epsGridProperties(eclState, false);
        const std::string tag = "Scaled endpoints";
        for (int c = 0; c < nc; ++c) {
            const std::string satnumIdx = std::to_string(epsGridProperties.satRegion(c));
            std::string cellIdx;
            {
                std::array<int, 3> ijk;
                cartesianIndexMapper.cartesianCoordinate(c, ijk);
                cellIdx = "(" + std::to_string(ijk[0]) + ", " +
                    std::to_string(ijk[1]) + ", " +
                    std::to_string(ijk[2]) + ")";
            }
            scaledEpsInfo_[c].extractScaled(eclState, epsGridProperties, c);

            // SGU <= 1.0 - SWL
            if (scaledEpsInfo_[c].Sgu > (1.0 - scaledEpsInfo_[c].Swl + tolerance)) {
                const std::string msg = "For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SGU exceed 1.0 - SWL";
                OpmLog::warning(tag, msg);
            }

            // SGL <= 1.0 - SWU
            if (scaledEpsInfo_[c].Sgl > (1.0 - scaledEpsInfo_[c].Swu + tolerance)) {
                const std::string msg = "For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SGL exceed 1.0 - SWU";
                OpmLog::warning(tag, msg);
            }

            if (threepoint && fluidSystem_ == FluidSystem::BlackOil) {
                // Mobilility check.
		    if ((scaledEpsInfo_[c].Sowcr + scaledEpsInfo_[c].Swcr) >= (1.0 + tolerance)) {
                    const std::string msg = "For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SOWCR + SWCR exceed 1.0";
                    OpmLog::warning(tag, msg);
                }

            if ((scaledEpsInfo_[c].Sogcr + scaledEpsInfo_[c].Sgcr + scaledEpsInfo_[c].Swl) >= (1.0 + tolerance)) {
                    const std::string msg = "For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SOGCR + SGCR + SWL exceed 1.0";
                    OpmLog::warning(tag, msg);
                }
            }
        }
    }

} //namespace Opm

#endif // OPM_RELPERMDIAGNOSTICS_IMPL_HEADER_INCLUDED
