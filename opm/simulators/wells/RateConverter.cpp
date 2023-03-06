/*
  Copyright 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2017, IRIS

  This file is part of the Open Porous Media Project (OPM).

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
#include <config.h>
#include <opm/simulators/wells/RateConverter.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <algorithm>

namespace Opm {
namespace RateConverter {

template <class FluidSystem, class Region>
void SurfaceToReservoirVoidage<FluidSystem,Region>::
sumRates(std::unordered_map<RegionId,Attributes>& attributes_hpv,
         std::unordered_map<RegionId,Attributes>& attributes_pv,
         Parallel::Communication comm)
{
    for (const auto& reg : rmap_.activeRegions()) {
        // Calculating first using the hydrocarbon pore volumes
        auto& ra = attr_.attributes(reg);
        const auto& attri_hpv = attributes_hpv[reg];
        ra.data = attri_hpv.data;
        comm.sum(ra.data.data(), ra.data.size());
        // TODO: should we have some epsilon here instead of zero?
        // No hydrocarbon pore volume, redo but this time using full pore volumes.
        if (ra.pv <= 0.) {
            // using the pore volume to do the averaging
            const auto& attri_pv = attributes_pv[reg];
            ra.data = attri_pv.data;
            comm.sum(ra.data.data(), ra.data.size());
            assert(ra.pv > 0.);
        }
        const double pv_sum = ra.pv;
        for (double& d : ra.data)
            d /= pv_sum;
        ra.pv = pv_sum;
    }
}

template <class FluidSystem, class Region>
template <class Coeff>
void SurfaceToReservoirVoidage<FluidSystem,Region>::
calcInjCoeff(const RegionId r, const int pvtRegionIdx, Coeff& coeff) const
{
    const auto& pu = phaseUsage_;
    const auto& ra = attr_.attributes(r);
    const double p = ra.pressure;
    const double T = ra.temperature;
    const double saltConcentration = ra.saltConcentration;

    const int   iw = RegionAttributeHelpers::PhasePos::water(pu);
    const int   io = RegionAttributeHelpers::PhasePos::oil  (pu);
    const int   ig = RegionAttributeHelpers::PhasePos::gas  (pu);

    std::fill(& coeff[0], & coeff[0] + phaseUsage_.num_phases, 0.0);

    if (RegionAttributeHelpers::PhaseUsed::water(pu)) {
        // q[w]_r = q[w]_s / bw

        const double bw = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, 0.0, saltConcentration);

        coeff[iw] = 1.0 / bw;
    }

    if (RegionAttributeHelpers::PhaseUsed::oil(pu)) {
        const double bo = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, 0.0);
        coeff[io] += 1.0 / bo;
    }

    if (RegionAttributeHelpers::PhaseUsed::gas(pu)) {
        const double bg = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, 0.0, 0.0);
        coeff[ig] += 1.0 / bg;
    }
}

using FS = BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>;
template void SurfaceToReservoirVoidage<FS,std::vector<int>>::
    sumRates(std::unordered_map<int,Attributes>&,
             std::unordered_map<int,Attributes>&,
             Parallel::Communication);

template void SurfaceToReservoirVoidage<FS,std::vector<int>>::
    calcInjCoeff<std::vector<double>>(const int, const int, std::vector<double>&) const;

} // namespace RateConverter
} // namespace Opm
