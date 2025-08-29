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

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

template <typename Scalar, typename Rates>
std::pair<Scalar, Scalar>
dissolvedVaporisedRatio(const int    io,
                        const int    ig,
                        const Scalar rs,
                        const Scalar rv,
                        const Rates& surface_rates)
{
    if ((io < 0) || (ig < 0)) {
        return { rs, rv };
    }
    auto eps = std::copysign(1.0e-15, surface_rates[io]);
    const auto Rs = surface_rates[ig] / (surface_rates[io] + eps);

    eps = std::copysign(1.0e-15, surface_rates[ig]);
    const auto Rv = surface_rates[io] / (surface_rates[ig] + eps);

    return {
        std::clamp(static_cast<Scalar>(Rs), Scalar{0.0}, rs),
        std::clamp(static_cast<Scalar>(Rv), Scalar{0.0}, rv)
    };
}

}

namespace Opm {
namespace RateConverter {

template <class FluidSystem, class Region>
void SurfaceToReservoirVoidage<FluidSystem, Region>::
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
        const Scalar pv_sum = ra.pv;
        std::transform(ra.data.begin(), ra.data.end(), ra.data.begin(),
                       [pv_sum](const auto d) { return d / pv_sum; });
        ra.pv = pv_sum;
    }
}

template <class FluidSystem, class Region>
template <class Coeff>
void SurfaceToReservoirVoidage<FluidSystem,  Region>::
calcInjCoeff(const RegionId r, const int pvtRegionIdx, Coeff& coeff) const
{
    const auto& ra = attr_.attributes(r);
    const Scalar p = ra.pressure;
    const Scalar T = ra.temperature;
    const Scalar saltConcentration = ra.saltConcentration;

    const int   iw = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
    const int   io = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
    const int   ig = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);

    std::fill(& coeff[0], & coeff[0] + FluidSystem::numActivePhases(), 0.0);

    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        // q[w]_r = q[w]_s / bw

        const Scalar bw = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx,
                                                                               T,
                                                                               p,
                                                                               Scalar{0.0},
                                                                               saltConcentration);

        coeff[iw] = 1.0 / bw;
    }

    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        const Scalar bo = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIdx,
                                                                             T,
                                                                             p,
                                                                             Scalar{0.0});
        coeff[io] += 1.0 / bo;
    }

    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        const Scalar bg = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx,
                                                                             T,
                                                                             p,
                                                                             Scalar{0.0},
                                                                             Scalar{0.0});
        coeff[ig] += 1.0 / bg;
    }
}

template <class FluidSystem, class Region>
template <class Coeff>
void SurfaceToReservoirVoidage<FluidSystem,  Region>::
calcCoeff(const RegionId r, const int pvtRegionIdx, Coeff& coeff) const
{
    const auto& ra = attr_.attributes(r);
    calcCoeff(pvtRegionIdx, ra.pressure, ra.rs, ra.rv, ra.rsw, ra.rvw, ra.temperature, ra.saltConcentration, coeff);
}

template <class FluidSystem, class Region>
template <class Coeff, class Rates>
void SurfaceToReservoirVoidage<FluidSystem,  Region>::
calcCoeff(const RegionId r, const int pvtRegionIdx, const Rates& surface_rates, Coeff& coeff) const
{
    const auto& ra = attr_.attributes(r);
    const int   iw = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
    const int   io = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
    const int   ig = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
    const auto [Rs, Rv] =
        dissolvedVaporisedRatio(io, ig, ra.rs, ra.rv, surface_rates);

    const auto [Rsw, Rvw] =
        dissolvedVaporisedRatio(iw, ig, ra.rsw, ra.rvw, surface_rates);

    calcCoeff(pvtRegionIdx, ra.pressure, Rs, Rv, Rsw, Rvw, ra.temperature, ra.saltConcentration, coeff);
}

template <class FluidSystem, class Region>
template <class Coeff>
void SurfaceToReservoirVoidage<FluidSystem,  Region>::
calcCoeff(const int pvtRegionIdx,
          const Scalar p,
          const Scalar Rs,
          const Scalar Rv,
          const Scalar Rsw,
          const Scalar Rvw,
          const Scalar T,
          const Scalar saltConcentration,
          Coeff& coeff) const
{
    const int   iw = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
    const int   io = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
    const int   ig = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);

    std::fill(& coeff[0], & coeff[0] + FluidSystem::numActivePhases(), 0.0);

    // Determinant of 'R' matrix
    const Scalar detRw = 1.0 - (Rsw * Rvw);

    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        // q[w]_r = 1/(bw * (1 - rsw*rvw)) * (q[w]_s - rvw*q[g]_s)

        const Scalar bw = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rsw, saltConcentration);

        const Scalar den = bw * detRw;

        coeff[iw] += 1.0 / den;

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            coeff[ig] -= Rvw / den;
        }
    }

    // Determinant of 'R' matrix
    const Scalar detR = 1.0 - (Rs * Rv);

    // Currently we only support either gas in water or gas in oil
    // not both
    if (detR != 1 && detRw != 1) {
        std::string msg = "We currently support either gas in water or gas in oil, not both."
        "i.e. detR = " + std::to_string(detR) + " and detRw " + std::to_string(detRw) +
        "can not both be different from 1";
        throw std::range_error(msg);
    }

    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        // q[o]_r = 1/(bo * (1 - rs*rv)) * (q[o]_s - rv*q[g]_s)

        const Scalar bo = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rs);
        const Scalar den = bo * detR;

        coeff[io] += 1.0 / den;

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            coeff[ig] -= Rv / den;
        }
    }

    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        // q[g]_r = 1/(bg * (1 - rs*rv)) * (q[g]_s - rs*q[o]_s)
        const Scalar bg  = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rv, Rvw);
        if (FluidSystem::enableDissolvedGasInWater()) {
            const Scalar denw = bg * detRw;
            coeff[ig] += 1.0 / denw;
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                coeff[iw] -= Rsw / denw;
            }
        } else {
            const Scalar den = bg * detR;
            coeff[ig] += 1.0 / den;
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                coeff[io] -= Rs / den;
            }
        }
    }
}


template <class FluidSystem, class Region>
template <typename SurfaceRates, typename VoidageRates>
void SurfaceToReservoirVoidage<FluidSystem,  Region>::
calcReservoirVoidageRates(const int           pvtRegionIdx,
                          const Scalar        p,
                          const Scalar        rs,
                          const Scalar        rv,
                          const Scalar        rsw,
                          const Scalar        rvw,
                          const Scalar        T,
                          const Scalar        saltConcentration,
                          const SurfaceRates& surface_rates,
                          VoidageRates&       voidage_rates) const
{
    const int   iw = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
    const int   io = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
    const int   ig = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);

    const auto [Rs, Rv] =
        dissolvedVaporisedRatio(io, ig, rs, rv, surface_rates);

    const auto [Rsw, Rvw] =
        dissolvedVaporisedRatio(iw, ig, rsw, rvw, surface_rates);


    std::fill_n(&voidage_rates[0], FluidSystem::numActivePhases(), 0.0);


    // Determinant of 'R' matrix
    const auto detRw = 1.0 - (Rsw * Rvw);

    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
        // q[w]_r = 1/(bw * (1 - rsw*rvw)) * (q[w]_s - rvw*q[g]_s)
        voidage_rates[iw] = surface_rates[iw];

        const auto bw = FluidSystem::waterPvt()
            .inverseFormationVolumeFactor(pvtRegionIdx, T, p,
                                          Rsw,
                                          saltConcentration);

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            voidage_rates[iw] -= Rvw * surface_rates[ig];
        }
        voidage_rates[iw] /= bw * detRw;
    }

    // Determinant of 'R' matrix
    const auto detR = 1.0 - (Rs * Rv);

    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        // q[o]_r = 1/(bo * (1 - rs*rv)) * (q[o]_s - rv*q[g]_s)
        voidage_rates[io] = surface_rates[io];

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            voidage_rates[io] -= Rv * surface_rates[ig];
        }

        const auto bo = FluidSystem::oilPvt()
            .inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rs);

        voidage_rates[io] /= bo * detR;
    }

    // we only support either gas in water
    // or gas in oil
    if (detR != 1 && detRw != 1) {
        std::string msg = "only support " + std::to_string(detR) + " " + std::to_string(detR);
        throw std::range_error(msg);
    }
    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        // q[g]_r = 1/(bg * (1 - rs*rv)) * (q[g]_s - rs*q[o]_s)
        voidage_rates[ig] = surface_rates[ig];
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            voidage_rates[ig] -= Rs * surface_rates[io];
        }
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            voidage_rates[ig] -= Rsw * surface_rates[iw];
        }

        const auto bg = FluidSystem::gasPvt()
            .inverseFormationVolumeFactor(pvtRegionIdx, T, p,
                                          Rv, Rvw);

        // we only support either gas in water or gas in oil
        if (detRw == 1) {
            voidage_rates[ig] /= bg * detR;
        } else {
            voidage_rates[ig] /= bg * detRw;
        }
    }
}

template <class FluidSystem, class Region>
template <class Rates>
void SurfaceToReservoirVoidage<FluidSystem,  Region>::
calcReservoirVoidageRates(const RegionId r,
                          const int      pvtRegionIdx,
                          const Rates&   surface_rates,
                          Rates&         voidage_rates) const
{
    const auto& ra = this->attr_.attributes(r);

    this->calcReservoirVoidageRates(pvtRegionIdx,
                                    ra.pressure, ra.rs, ra.rv,
                                    ra.rsw, ra.rvw,
                                    ra.temperature,
                                    ra.saltConcentration,
                                    surface_rates,
                                    voidage_rates);
}

template <class FluidSystem, class Region>
template <class Rates>
std::pair<typename FluidSystem::Scalar, typename FluidSystem::Scalar>
SurfaceToReservoirVoidage<FluidSystem, Region>::
inferDissolvedVaporisedRatio(const Scalar rsMax,
                             const Scalar rvMax,
                             const Rates& surface_rates) const
{
    // we probably should add checking whether the phases are active
    const auto io = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
    const auto ig = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
    return dissolvedVaporisedRatio(io, ig, rsMax, rvMax, surface_rates);
}

template<class Scalar>
using FS = BlackOilFluidSystem<Scalar, BlackOilDefaultFluidSystemIndices>;

#define INSTANTIATE_TYPE(T)                                              \
    template void SurfaceToReservoirVoidage<FS<T>,std::vector<int>>::    \
        sumRates(std::unordered_map<int,Attributes>&,                    \
                 std::unordered_map<int,Attributes>&,                    \
                 Parallel::Communication);                               \
    template void SurfaceToReservoirVoidage<FS<T>,std::vector<int>>::    \
        calcInjCoeff(const int, const int,                               \
                     std::vector<T>&) const;                             \
    template void SurfaceToReservoirVoidage<FS<T>,std::vector<int>>::    \
        calcCoeff(const int, const int, std::vector<T>&) const;          \
    template void SurfaceToReservoirVoidage<FS<T>,std::vector<int>>::    \
        calcCoeff(const int,                                             \
                  const int,                                             \
                  const std::vector<T>&,                                 \
                  std::vector<T>&) const;                                \
    template void SurfaceToReservoirVoidage<FS<T>,std::vector<int>>::    \
        calcReservoirVoidageRates(const int,                             \
                                  const T,                               \
                                  const T,                               \
                                  const T,                               \
                                  const T,                               \
                                  const T,                               \
                                  const T,                               \
                                  const T,                               \
                                  const std::vector<T>&,                 \
                                  std::vector<T>&) const;                \
    template void SurfaceToReservoirVoidage<FS<T>,std::vector<int>>::    \
        calcReservoirVoidageRates(const int,                             \
                                  const T,                               \
                                  const T,                               \
                                  const T,                               \
                                  const T,                               \
                                  const T,                               \
                                  const T,                               \
                                  const T,                               \
                                  T const* const&,                       \
                                  T*&) const;                            \
    template void SurfaceToReservoirVoidage<FS<T>,std::vector<int>>::    \
        calcReservoirVoidageRates(const int,                             \
                                  const int,                             \
                                  const std::vector<T>&,                 \
                                  std::vector<T>&) const;                \
    template std::pair<T,T>                                              \
    SurfaceToReservoirVoidage<FS<T>,std::vector<int>>::                  \
        inferDissolvedVaporisedRatio(const T,                            \
                                     const T,                            \
                                     const std::vector<T>::iterator&) const;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace RateConverter
} // namespace Opm
