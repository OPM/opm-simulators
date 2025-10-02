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
 * \copydoc Opm::FlowProblemBlackoil
 */

#include <config.h>
#include <opm/simulators/flow/MixingRateControls.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <algorithm>
#include <limits>

namespace Opm {

template<class FluidSystem>
MixingRateControls<FluidSystem>::
MixingRateControls(const Schedule& schedule)
    : schedule_(schedule)
{}

template<class FluidSystem>
MixingRateControls<FluidSystem>::
MixingRateControls(const MixingRateControls& rhs)
    : schedule_(rhs.schedule_)
{
    *this = rhs;
}

template<class FluidSystem>
MixingRateControls<FluidSystem>
MixingRateControls<FluidSystem>::
serializationTestObject(const Schedule& schedule)
{
    MixingRateControls<FluidSystem> result(schedule);
    result.lastRv_ = {21.0};
    result.maxDRv_ = {22.0, 23.0};
    result.convectiveDrs_ = {24.0, 25.0, 26.0};
    result.lastRs_ = {27.0};
    result.maxDRs_ = {28.0};
    result.dRsDtOnlyFreeGas_ = {false, true};

    return result;
}

template<class FluidSystem>
bool MixingRateControls<FluidSystem>::
operator==(const MixingRateControls& rhs) const
{
    return this->lastRv_ == rhs.lastRv_ &&
           this->maxDRv_ == rhs.maxDRv_ &&
           this->convectiveDrs_ == rhs.convectiveDrs_ &&
           this->lastRs_ == rhs.lastRs_ &&
           this->maxDRs_ == rhs.maxDRs_ &&
           this->dRsDtOnlyFreeGas_ == rhs.dRsDtOnlyFreeGas_;
}

template<class FluidSystem>
MixingRateControls<FluidSystem>&
MixingRateControls<FluidSystem>::
operator=(const MixingRateControls& rhs)
{
    this->lastRv_ = rhs.lastRv_;
    this->maxDRv_ = rhs.maxDRv_;
    this->convectiveDrs_ = rhs.convectiveDrs_;
    this->lastRs_ = rhs.lastRs_;
    this->maxDRs_ = rhs.maxDRs_;
    this->dRsDtOnlyFreeGas_ = rhs.dRsDtOnlyFreeGas_;

    return *this;
}

template<class FluidSystem>
void MixingRateControls<FluidSystem>::
init(std::size_t numDof, int episodeIdx, const unsigned ntpvt)
{
    // allocate DRSDT related vectors
    if (this->drsdtActive(episodeIdx) && maxDRs_.empty()) {
        maxDRs_.resize(ntpvt, 1e30);
        dRsDtOnlyFreeGas_.resize(ntpvt, false);
        lastRs_.resize(numDof, 0.0);
    }
    if (this->drvdtActive(episodeIdx) && maxDRv_.empty()) {
        lastRv_.resize(numDof, 0.0);
        maxDRv_.resize(ntpvt, 1e30);
    }
    if (this->drsdtConvective(episodeIdx) && convectiveDrs_.empty()) {
        convectiveDrs_.resize(numDof, 1.0);
    }
}

template<class FluidSystem>
bool MixingRateControls<FluidSystem>::
drsdtActive(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    return (oilVaporizationControl.drsdtActive());
}

template<class FluidSystem>
bool MixingRateControls<FluidSystem>::
drvdtActive(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    return (oilVaporizationControl.drvdtActive());
}

template<class FluidSystem>
bool MixingRateControls<FluidSystem>::
drsdtConvective(int episodeIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    return (oilVaporizationControl.drsdtConvective());
}

template<class FluidSystem>
bool MixingRateControls<FluidSystem>::
drsdtActive(int episodeIdx, std::size_t pvtRegionIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    return (oilVaporizationControl.drsdtActive(pvtRegionIdx));
}

template<class FluidSystem>
bool MixingRateControls<FluidSystem>::
drvdtActive(int episodeIdx, std::size_t pvtRegionIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    return (oilVaporizationControl.drvdtActive(pvtRegionIdx));
}

template<class FluidSystem>
bool MixingRateControls<FluidSystem>::
drsdtConvective(int episodeIdx, std::size_t pvtRegionIdx) const
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    return (oilVaporizationControl.drsdtConvective(pvtRegionIdx));
}


template<class FluidSystem>
void MixingRateControls<FluidSystem>::
updateExplicitQuantities(const int episodeIdx,
                         const Scalar timeStepSize)
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    // DRSDT is enabled
    for (std::size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRs_.size(); ++pvtRegionIdx) {
        if (oilVaporizationControl.drsdtActive(pvtRegionIdx)) {
            maxDRs_[pvtRegionIdx] = oilVaporizationControl.getMaxDRSDT(pvtRegionIdx) * timeStepSize;
        }
    }

    // DRVDT is enabled
    for (std::size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRv_.size(); ++pvtRegionIdx) {
        if (oilVaporizationControl.drvdtActive(pvtRegionIdx)) {
            maxDRv_[pvtRegionIdx] = oilVaporizationControl.getMaxDRVDT(pvtRegionIdx) * timeStepSize;
        }
    }
}

template<class FluidSystem>
void MixingRateControls<FluidSystem>::
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

template<class FluidSystem>
void MixingRateControls<FluidSystem>::
updateMaxValues(const int episodeIdx,
                const Scalar timeStepSize)
{
    const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();
    // DRSDT is enabled
    for (std::size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRs_.size(); ++pvtRegionIdx) {
        if (this->drsdtActive(episodeIdx, pvtRegionIdx)) {
            maxDRs_[pvtRegionIdx] = oilVaporizationControl.getMaxDRSDT(pvtRegionIdx) * timeStepSize;
        }
    }

       // DRVDT is enabled
    for (std::size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRv_.size(); ++pvtRegionIdx) {
        if (this->drvdtActive(episodeIdx, pvtRegionIdx)) {
            maxDRv_[pvtRegionIdx] = oilVaporizationControl.getMaxDRVDT(pvtRegionIdx) * timeStepSize;
        }
    }
}

template<class FluidSystem>
typename MixingRateControls<FluidSystem>::Scalar
MixingRateControls<FluidSystem>::
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

    if (!oilVaporizationControl.drsdtConvective(pvtRegionIdx))
        return 0;

    return oilVaporizationControl.getMaxDRSDT(pvtRegionIdx) * convectiveDrs_[elemIdx];
}

template<class FluidSystem>
typename MixingRateControls<FluidSystem>::Scalar
MixingRateControls<FluidSystem>::
maxGasDissolutionFactor(const unsigned timeIdx,
                        const unsigned globalDofIdx,
                        const int episodeIdx,
                        const int pvtRegionIdx) const
{
    if (!this->drsdtActive(episodeIdx, pvtRegionIdx)) {
        return std::numeric_limits<Scalar>::max() / 2.0;
    }

    Scalar scaling = 1.0;
    if (this->drsdtConvective(episodeIdx, pvtRegionIdx)) {
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

template<class FluidSystem>
typename MixingRateControls<FluidSystem>::Scalar
MixingRateControls<FluidSystem>::
maxOilVaporizationFactor(const unsigned timeIdx,
                         const unsigned globalDofIdx,
                         const int episodeIdx,
                         const int pvtRegionIdx) const
{
    if (!this->drvdtActive(episodeIdx, pvtRegionIdx)) {
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

template<class FluidSystem>
void MixingRateControls<FluidSystem>::
updateConvectiveDRsDt_(const unsigned compressedDofIdx,
                       const Scalar t,
                       const Scalar p,
                       const Scalar pg,
                       const Scalar rs,
                       const Scalar sg,
                       const Scalar poro,
                       const Scalar permz,
                       const Scalar distZ,
                       const Scalar gravity,
                       const Scalar salt,
                       const Scalar Xhi,
                       const Scalar Psi,
                       const Scalar omegainn,
                       const int pvtRegionIndex)
{
    const Scalar rssat = (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) ?
        FluidSystem::waterPvt().saturatedGasDissolutionFactor(pvtRegionIndex, t, p, salt) :
        FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIndex, t, p);
    const Scalar saturatedInvB = (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) ?
        FluidSystem::waterPvt().saturatedInverseFormationVolumeFactor(pvtRegionIndex, t, p, salt) :
        FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvtRegionIndex, t, p);
    const Scalar rsZero = 0.0;
    const Scalar sg_max = 1.0;
    const Scalar pureDensity = (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) ?
        (FluidSystem::waterPvt().inverseFormationVolumeFactor(pvtRegionIndex, t, p, rsZero, salt)
        * FluidSystem::waterPvt().waterReferenceDensity(pvtRegionIndex)) :
         (FluidSystem::oilPvt().inverseFormationVolumeFactor(pvtRegionIndex, t, p, rsZero)
        * FluidSystem::oilPvt().oilReferenceDensity(pvtRegionIndex));
    const Scalar saturatedDensity = (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) ?
        (saturatedInvB * (FluidSystem::waterPvt().waterReferenceDensity(pvtRegionIndex)
        + rssat * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIndex))) :
        (saturatedInvB * (FluidSystem::oilPvt().oilReferenceDensity(pvtRegionIndex)
        + rssat * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIndex)));
    Scalar deltaDensity = saturatedDensity - pureDensity;
    const Scalar visc = (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) ?
        FluidSystem::waterPvt().viscosity(pvtRegionIndex, t, p, rs, salt) :
        FluidSystem::oilPvt().viscosity(pvtRegionIndex, t, p, rs);

    // Note that for sLiquid = 0 this gives no limits (inf) for the dissolution rate
    // Also we restrict the effect of convective mixing to positive density differences
    // i.e. we only allow for fingers moving downward

    Scalar co2Density =
        FluidSystem::gasPvt().inverseFormationVolumeFactor(pvtRegionIndex,t,p,Scalar{0.0} /*=Rv*/, Scalar{0.0} /*=Rvw*/) *
        FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIndex);
    Scalar factor = 1.0;
    Scalar X = (rs - rssat * sg) / (rssat * ( 1.0 - sg));
    Scalar omega = 0.0;
    const Scalar pCap = Opm::abs(pg - p);
    if ((rs >= (rssat * sg)) || (pCap < 1e-12)) {
        if (X > Psi) {
            factor = 0.0;
            omega = omegainn;
        }
    } else {
        factor /= Xhi;
        deltaDensity = (saturatedDensity - co2Density);
    }

    convectiveDrs_[compressedDofIdx] =
        factor * permz * rssat * max(Scalar{0.0}, deltaDensity) *
        gravity / ( std::max(sg_max - sg, Scalar{0.0}) *
        visc * distZ * poro) + (omega/Xhi);
}

#define INSTANTIATE_TYPE(T) \
    template class MixingRateControls<BlackOilFluidSystem<T, BlackOilDefaultFluidSystemIndices>>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
