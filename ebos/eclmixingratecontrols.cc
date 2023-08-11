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
 * \copydoc Opm::EclProblem
 */

#include <config.h>
#include <ebos/eclmixingratecontrols.hh>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <algorithm>
#include <limits>

namespace Opm {

template<class FluidSystem, class Scalar>
EclMixingRateControls<FluidSystem,Scalar>::
EclMixingRateControls(const Schedule& schedule)
    : schedule_(schedule)
{}

template<class FluidSystem, class Scalar>
EclMixingRateControls<FluidSystem,Scalar>::
EclMixingRateControls(const EclMixingRateControls& rhs)
    : schedule_(rhs.schedule_)
{
    *this = rhs;
}

template<class FluidSystem, class Scalar>
EclMixingRateControls<FluidSystem,Scalar>
EclMixingRateControls<FluidSystem,Scalar>::
serializationTestObject(const Schedule& schedule)
{
    EclMixingRateControls<FluidSystem,Scalar> result(schedule);
    result.lastRv_ = {21.0};
    result.maxDRv_ = {22.0, 23.0};
    result.convectiveDrs_ = {24.0, 25.0, 26.0};
    result.lastRs_ = {27.0};
    result.maxDRs_ = {28.0};
    result.dRsDtOnlyFreeGas_ = {false, true};

    return result;
}

template<class FluidSystem, class Scalar>
bool EclMixingRateControls<FluidSystem,Scalar>::
operator==(const EclMixingRateControls& rhs) const
{
    return this->lastRv_ == rhs.lastRv_ &&
           this->maxDRv_ == rhs.maxDRv_ &&
           this->convectiveDrs_ == rhs.convectiveDrs_ &&
           this->lastRs_ == rhs.lastRs_ &&
           this->maxDRs_ == rhs.maxDRs_ &&
           this->dRsDtOnlyFreeGas_ == rhs.dRsDtOnlyFreeGas_;
}

template<class FluidSystem, class Scalar>
EclMixingRateControls<FluidSystem,Scalar>&
EclMixingRateControls<FluidSystem,Scalar>::
operator=(const EclMixingRateControls& rhs)
{
    this->lastRv_ = rhs.lastRv_;
    this->maxDRv_ = rhs.maxDRv_;
    this->convectiveDrs_ = rhs.convectiveDrs_;
    this->lastRs_ = rhs.lastRs_;
    this->maxDRs_ = rhs.maxDRs_;
    this->dRsDtOnlyFreeGas_ = rhs.dRsDtOnlyFreeGas_;

    return *this;
}

template<class FluidSystem, class Scalar>
void EclMixingRateControls<FluidSystem,Scalar>::
init(std::size_t numDof, int episodeIdx, const unsigned ntpvt)
{
    // deal with DRSDT
    //TODO We may want to only allocate these properties only if active.
    //But since they may be activated at later time we need some more
    //intrastructure to handle it
    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
        FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
        maxDRv_.resize(ntpvt, 1e30);
        lastRv_.resize(numDof, 0.0);
        maxDRs_.resize(ntpvt, 1e30);
        dRsDtOnlyFreeGas_.resize(ntpvt, false);
        lastRs_.resize(numDof, 0.0);
        maxDRv_.resize(ntpvt, 1e30);
        lastRv_.resize(numDof, 0.0);
        if (this->drsdtConvective(episodeIdx)) {
            convectiveDrs_.resize(numDof, 1.0);
        }
    }
}

template<class FluidSystem, class Scalar>
bool EclMixingRateControls<FluidSystem,Scalar>::
drsdtActive(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    const bool bothOilGasActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                                  FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    return (oilVaporizationControl.drsdtActive() && bothOilGasActive);
}

template<class FluidSystem, class Scalar>
bool EclMixingRateControls<FluidSystem,Scalar>::
drvdtActive(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    const bool bothOilGasActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                                  FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    return (oilVaporizationControl.drvdtActive() && bothOilGasActive);
}

template<class FluidSystem, class Scalar>
bool EclMixingRateControls<FluidSystem,Scalar>::
drsdtConvective(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    const bool bothOilGasActive = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                                  FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
    return (oilVaporizationControl.drsdtConvective() && bothOilGasActive);
}

template<class FluidSystem, class Scalar>
void EclMixingRateControls<FluidSystem,Scalar>::
updateExplicitQuantities(const int episodeIdx,
                         const Scalar timeStepSize)
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    if (this->drsdtActive(episodeIdx)) {
        // DRSDT is enabled
        for (std::size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRs_.size(); ++pvtRegionIdx)
            maxDRs_[pvtRegionIdx] = oilVaporizationControl.getMaxDRSDT(pvtRegionIdx) * timeStepSize;
    }

    if (this->drvdtActive(episodeIdx)) {
        // DRVDT is enabled
        for (std::size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRv_.size(); ++pvtRegionIdx)
            maxDRv_[pvtRegionIdx] = oilVaporizationControl.getMaxDRVDT(pvtRegionIdx) * timeStepSize;
    }
}

template<class FluidSystem, class Scalar>
void EclMixingRateControls<FluidSystem,Scalar>::
updateLastValues(const unsigned elemIdx,
                 const Scalar Rs,
                 const Scalar Rv)
{
    if (!lastRs_.empty()) {
        lastRs_[elemIdx] = Rs;
    }

    if (!lastRv_.empty()) {
        lastRv_[elemIdx] = Rv;
    }
}

template<class FluidSystem, class Scalar>
void EclMixingRateControls<FluidSystem,Scalar>::
updateMaxValues(const int episodeIdx,
                const Scalar timeStepSize)
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    if (this->drsdtActive(episodeIdx)) {
        // DRSDT is enabled
        for (std::size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRs_.size(); ++pvtRegionIdx) {
            maxDRs_[pvtRegionIdx] = oilVaporizationControl.getMaxDRSDT(pvtRegionIdx) * timeStepSize;
        }
    }

    if (this->drvdtActive(episodeIdx)) {
        // DRVDT is enabled
        for (std::size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRv_.size(); ++pvtRegionIdx) {
            maxDRv_[pvtRegionIdx] = oilVaporizationControl.getMaxDRVDT(pvtRegionIdx) * timeStepSize;
        }
    }
}

template<class FluidSystem, class Scalar>
Scalar EclMixingRateControls<FluidSystem,Scalar>::
drsdtcon(const unsigned elemIdx,
         int episodeIdx,
         const int pvtRegionIdx) const
{
    if (convectiveDrs_.empty()) {
        return 0;
    }

    // The episode index is set to -1 in the initialization phase.
    // Output drsdt value for index 0
    episodeIdx = std::max(episodeIdx, 0);
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    return oilVaporizationControl.getMaxDRSDT(pvtRegionIdx) * convectiveDrs_[elemIdx];
}

template<class FluidSystem, class Scalar>
Scalar EclMixingRateControls<FluidSystem,Scalar>::
maxGasDissolutionFactor(const unsigned timeIdx,
                        const unsigned globalDofIdx,
                        const int episodeIdx,
                        const int pvtRegionIdx) const
{
    if (!this->drsdtActive(episodeIdx) || maxDRs_[pvtRegionIdx] < 0.0) {
        return std::numeric_limits<Scalar>::max() / 2.0;
    }

    Scalar scaling = 1.0;
    if (this->drsdtConvective(episodeIdx)) {
       scaling = convectiveDrs_[globalDofIdx];
    }

    // this is a bit hacky because it assumes that a time discretization with only
    // two time indices is used.
    if (timeIdx == 0) {
        return lastRs_[globalDofIdx] + maxDRs_[pvtRegionIdx] * scaling;
    } else {
        return lastRs_[globalDofIdx];
    }
}

template<class FluidSystem, class Scalar>
Scalar EclMixingRateControls<FluidSystem,Scalar>::
maxOilVaporizationFactor(const unsigned timeIdx,
                         const unsigned globalDofIdx,
                         const int episodeIdx,
                         const int pvtRegionIdx) const
{
    if (!this->drvdtActive(episodeIdx) || maxDRv_[pvtRegionIdx] < 0.0) {
        return std::numeric_limits<Scalar>::max() / 2.0;
    }

    // this is a bit hacky because it assumes that a time discretization with only
    // two time indices is used.
    if (timeIdx == 0) {
        return lastRv_[globalDofIdx] + maxDRv_[pvtRegionIdx];
    } else {
        return lastRv_[globalDofIdx];
    }
}
template<class FluidSystem, class Scalar>
void EclMixingRateControls<FluidSystem,Scalar>::
updateConvectiveDRsDt_(const unsigned compressedDofIdx,
                       const Scalar t,
                       const Scalar p,
                       const Scalar rs,
                       const Scalar so,
                       const Scalar poro,
                       const Scalar permz,
                       const Scalar distZ,
                       const Scalar gravity,
                       const int pvtRegionIndex)
{
    const Scalar rssat = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIndex, t, p);
    const Scalar saturatedInvB
        = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvtRegionIndex, t, p);
    const Scalar rsZero = 0.0;
    const Scalar pureDensity
        = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIndex, t, p, rsZero)
        * FluidSystem::oilPvt().oilReferenceDensity(pvtRegionIndex);
    const Scalar saturatedDensity = saturatedInvB
        * (FluidSystem::oilPvt().oilReferenceDensity(pvtRegionIndex)
           + rssat * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIndex));
    const Scalar deltaDensity = saturatedDensity - pureDensity;
    const Scalar visc = FluidSystem::oilPvt().viscosity(pvtRegionIndex, t, p, rs);
    // Note that for so = 0 this gives no limits (inf) for the dissolution rate
    // Also we restrict the effect of convective mixing to positive density differences
    // i.e. we only allow for fingers moving downward
    convectiveDrs_[compressedDofIdx]
        = permz * rssat * max(0.0, deltaDensity) * gravity / (so * visc * distZ * poro);
}

template class EclMixingRateControls<BlackOilFluidSystem<double,BlackOilDefaultIndexTraits>, double>;

} // namespace Opm
