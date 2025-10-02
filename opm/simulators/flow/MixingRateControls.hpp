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
#ifndef OPM_MIXING_RATE_CONTROLS_HPP
#define OPM_MIXING_RATE_CONTROLS_HPP

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <limits>
#include <vector>

namespace Opm {

class EclipseState;

//! \brief Class handling mixing rate controls for a FlowProblemBlackoil.
template<class FluidSystem>
class MixingRateControls
{
public:
    using Scalar = typename FluidSystem::Scalar;

    explicit MixingRateControls(const Schedule& schedule);
    MixingRateControls(const MixingRateControls& rhs);

    static MixingRateControls serializationTestObject(const Schedule& schedule);

    bool operator==(const MixingRateControls& rhs) const;
    MixingRateControls& operator=(const MixingRateControls& rhs);

    void init(std::size_t numDof, int episodeIdx, const unsigned ntpvt);

    bool drsdtActive(int episodeIdx) const;
    bool drvdtActive(int episodeIdx) const;
    bool drsdtConvective(int episodeIdx) const;

    bool drsdtActive(int episodeIdx, std::size_t pvtRegionIdx) const;
    bool drvdtActive(int episodeIdx, std::size_t pvtRegionIdx) const;
    bool drsdtConvective(int episodeIdx, std::size_t pvtRegionIdx) const;
    /*!
     * \brief Returns the dynamic drsdt convective mixing value
     */
    Scalar drsdtcon(const unsigned elemIdx,
                    int episodeIdx,
                    const int pvtRegionIdx) const;

    /*!
     * \brief Returns the maximum value of the gas dissolution factor at the current time
     *        for a given degree of freedom.
     */
    Scalar maxGasDissolutionFactor(unsigned timeIdx,
                                   unsigned globalDofIdx,
                                   const int episodeIdx,
                                   const int pvtRegionIdx) const;

    /*!
     * \brief Returns the maximum value of the oil vaporization factor at the current
     *        time for a given degree of freedom.
     */
    Scalar maxOilVaporizationFactor(const unsigned timeIdx,
                                    const unsigned globalDofIdx,
                                    const int episodeIdx,
                                    const int pvtRegionIdx) const;

    void updateExplicitQuantities(const int episodeIdx,
                                  const Scalar timeStepSize);

    void updateLastValues(const unsigned elemIdx,
                          const Scalar Rs,
                          const Scalar Rv);

    void updateMaxValues(const int episodeIdx,
                         const Scalar timeStepSize);

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(lastRv_);
        serializer(maxDRv_);
        serializer(convectiveDrs_);
        serializer(lastRs_);
        serializer(maxDRs_);
        serializer(dRsDtOnlyFreeGas_);
    }

    template<class IntensiveQuantities>
    void update(unsigned compressedDofIdx,
                const IntensiveQuantities& iq,
                const int episodeIdx,
                const Scalar gravity,
                const Scalar permZ,
                const Scalar distZ,
                const int pvtRegionIdx)
    {
        const auto& oilVaporizationControl = schedule_[episodeIdx].oilvap();

        if (oilVaporizationControl.drsdtConvective(pvtRegionIdx)) {
            // This implements the convective DRSDT as described in
            // Sandve et al. "Convective dissolution in field scale CO2 storage simulations using the OPM Flow
            // simulator" Submitted to TCCS 11, 2021
            // modification and introduction of regimes following Mykkeltvedt et al. Submitted to TIMP 2024
            const auto& fs = iq.fluidState();

            const auto& temperature = (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) ?
                getValue(fs.temperature(FluidSystem::waterPhaseIdx)) :
                getValue(fs.temperature(FluidSystem::oilPhaseIdx));
            const auto& pressure = (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) ?
                getValue(fs.pressure(FluidSystem::waterPhaseIdx)) :
                getValue(fs.pressure(FluidSystem::oilPhaseIdx));
            const auto& pressuregas = getValue(fs.pressure(FluidSystem::gasPhaseIdx));
            const auto& rs = (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) ?
                getValue(fs.Rsw()) :
                getValue(fs.Rs());

            const auto& salt = getValue(fs.saltSaturation());

            this->updateConvectiveDRsDt_(compressedDofIdx,
                                         temperature,
                                         pressure,
                                         pressuregas,
                                         rs,
                                         getValue(fs.saturation(FluidSystem::gasPhaseIdx)),
                                         getValue(iq.porosity()),
                                         permZ,
                                         distZ,
                                         gravity,
                                         salt,
                                         oilVaporizationControl.getMaxDRSDT(fs.pvtRegionIndex()),
                                         oilVaporizationControl.getPsi(fs.pvtRegionIndex()),
                                         oilVaporizationControl.getOmega(fs.pvtRegionIndex()),
                                         fs.pvtRegionIndex());
        }

        if (oilVaporizationControl.drsdtActive(pvtRegionIdx)) {
            const auto& fs = iq.fluidState();

            using FluidState = typename std::decay<decltype(fs)>::type;
            constexpr Scalar freeGasMinSaturation_ = 1e-7;
            if (oilVaporizationControl.getOption(pvtRegionIdx) ||
                fs.saturation(FluidSystem::gasPhaseIdx) > freeGasMinSaturation_) {
                lastRs_[compressedDofIdx]
                    = ((FluidSystem::enableDissolvedGasInWater())) ?
                    BlackOil::template getRsw_<FluidSystem, FluidState, Scalar>(fs, iq.pvtRegionIndex()) :
                    BlackOil::template getRs_<FluidSystem, FluidState, Scalar>(fs, iq.pvtRegionIndex());
            }
            else
                lastRs_[compressedDofIdx] = std::numeric_limits<Scalar>::infinity();
        }

        if (oilVaporizationControl.drvdtActive(pvtRegionIdx)) {
            const auto& fs = iq.fluidState();
            using FluidState = typename std::decay<decltype(fs)>::type;
            lastRv_[compressedDofIdx]
                = BlackOil::template getRv_<FluidSystem, FluidState, Scalar>(fs, iq.pvtRegionIndex());
        }
    }

private:
    void updateConvectiveDRsDt_(const unsigned compressedDofIdx,
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
                                const int pvtRegionIndex);

    std::vector<Scalar> lastRv_;
    std::vector<Scalar> maxDRv_;

    std::vector<Scalar> convectiveDrs_;
    std::vector<Scalar> lastRs_;
    std::vector<Scalar> maxDRs_;
    std::vector<bool> dRsDtOnlyFreeGas_; // apply the DRSDT rate limit only to cells that exhibit free gas

    const Schedule& schedule_;
};

} // namespace Opm

#endif // OPM_MIXING_RATE_CONTROLS_HPP
