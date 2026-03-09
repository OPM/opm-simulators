/*
  Copyright 2020 Equinor ASA.

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

#include <config.h>
#include <opm/simulators/wells/TargetCalculator.hpp>

#include <opm/material/fluidsystems/PhaseUsageInfo.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/GroupState.hpp>
#include <opm/simulators/wells/GroupStateHelper.hpp>

#include <cassert>
#include <stdexcept>

namespace Opm::GroupStateHelpers
{
template<typename Scalar, typename IndexTraits>
TargetCalculator<Scalar, IndexTraits>::
TargetCalculator(const Opm::GroupStateHelper<Scalar, IndexTraits>& groupStateHelper,
                 const std::vector<Scalar>& resv_coeff,
                 const Group& group)
    : cmode_{groupStateHelper.groupState().production_control(group.name())}
    , groupStateHelper_{groupStateHelper}
    , resv_coeff_{resv_coeff}
{
}

template<typename Scalar, typename IndexTraits>
TargetCalculator<Scalar, IndexTraits>::
TargetCalculator(const Opm::GroupStateHelper<Scalar, IndexTraits>& groupStateHelper,
                 const std::vector<Scalar>& resv_coeff,
                 Group::ProductionCMode cmode)
    : cmode_{cmode}
    , groupStateHelper_{groupStateHelper}
    , resv_coeff_{resv_coeff}
{
}

template<typename Scalar, typename IndexTraits>
template <typename RateType>
RateType TargetCalculator<Scalar, IndexTraits>::calcModeRateFromRates(const RateType* rates) const
{
    const auto& pu = this->groupStateHelper_.phaseUsage();
    switch (this->cmode_) {
    case Group::ProductionCMode::ORAT: {
        assert(pu.phaseIsActive(IndexTraits::oilPhaseIdx));
        const int pos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        return rates[pos];
    }
    case Group::ProductionCMode::WRAT: {
        assert(pu.phaseIsActive(IndexTraits::waterPhaseIdx));
        const int pos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
        return rates[pos];
    }
    case Group::ProductionCMode::GRAT: {
        assert(pu.phaseIsActive(IndexTraits::gasPhaseIdx));
        const int pos = pu.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
        return rates[pos];
    }
    case Group::ProductionCMode::LRAT: {
        // TODO: do need both oil and water activate to use LRAT?
        assert(pu.phaseIsActive(IndexTraits::oilPhaseIdx));
        assert(pu.phaseIsActive(IndexTraits::waterPhaseIdx));
        const int opos = pu.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        const int wpos = pu.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx);
        return rates[opos] + rates[wpos];
    }
    case Group::ProductionCMode::RESV: {
        auto mode_rate = rates[0] * this->resv_coeff_[0];
        const int num_phases = pu.numActivePhases();
        for (int phase = 1; phase < num_phases; ++phase) {
            mode_rate += rates[phase] * this->resv_coeff_[phase];
        }
        return mode_rate;
    }
    default:
        // Should never be here.
        assert(false);
        return rates[0];
    }
}

template<typename Scalar, typename IndexTraits>
InjectionTargetCalculator<Scalar, IndexTraits>::
InjectionTargetCalculator(const GroupStateHelperType& groupStateHelper,
                          const Phase& injection_phase)
{
    pos_ = groupStateHelper.phaseToActivePhaseIdx(injection_phase);

    switch (injection_phase) {
    case Phase::WATER:
    case Phase::OIL:
    case Phase::GAS:
        break;
    default:
        OPM_DEFLOG_THROW(std::logic_error,
                         "Invalid injection phase in InjectionTargetCalculator",
                         groupStateHelper.deferredLogger());
    }
}

#define INSTANTIATE_TARGET_CALCULATOR(T,...) \
    template __VA_ARGS__                     \
    TargetCalculator<T, BlackOilDefaultFluidSystemIndices>::calcModeRateFromRates(const __VA_ARGS__* rates) const;

#define INSTANTIATE_TYPE(T)                                       \
    template class TargetCalculator<T, BlackOilDefaultFluidSystemIndices>;                           \
    template class InjectionTargetCalculator<T, BlackOilDefaultFluidSystemIndices>;                  \
    INSTANTIATE_TARGET_CALCULATOR(T,T)                            \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,3,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,4,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,5,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,6,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,7,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,8,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,9,0>)   \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,10,0>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,4>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,5>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,6>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,7>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,8>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,9>)  \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,10>) \
    INSTANTIATE_TARGET_CALCULATOR(T,DenseAd::Evaluation<T,-1,11>)

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm::GroupStateHelpers
