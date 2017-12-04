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
 * \copydoc Opm::BlackOilFluidState
 */
#ifndef OPM_BLACK_OIL_FLUID_STATE_HH
#define OPM_BLACK_OIL_FLUID_STATE_HH

#include <opm/common/Valgrind.hpp>
#include <opm/common/Unused.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

namespace Opm {
/*!
 * \brief Implements a "taylor-made" fluid state class for the black-oil model.
 *
 * I.e., it only uses quantities which are available in the ECL blackoil model. Further
 * quantities are computed "on the fly" and are accessing them is thus relatively slow.
 */
template <class ScalarT, class FluidSystem>
class BlackOilFluidState
{
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };

    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };

public:
    typedef ScalarT Scalar;
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    /*!
     * \brief Make sure that all attributes are defined.
     *
     * This method does not do anything if the program is not run
     * under valgrind. If it is, then valgrind will print an error
     * message if some attributes of the object have not been properly
     * defined.
     */
    void checkDefined() const
    {
#ifndef NDEBUG
        Opm::Valgrind::CheckDefined(pvtRegionIdx_);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            Opm::Valgrind::CheckDefined(saturation_[phaseIdx]);
            Opm::Valgrind::CheckDefined(pressure_[phaseIdx]);
            Opm::Valgrind::CheckDefined(invB_[phaseIdx]);
        }

        Opm::Valgrind::CheckDefined(Rs_);
        Opm::Valgrind::CheckDefined(Rv_);
#endif // NDEBUG
    }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    void assign(const FluidState& fs OPM_UNUSED)
    {
        assert(false); // not yet implemented
    }

    void setPvtRegionIndex(unsigned newPvtRegionIdx)
    { pvtRegionIdx_ = static_cast<unsigned short>(newPvtRegionIdx); }

    void setPressure(unsigned phaseIdx, const Scalar& p)
    { pressure_[phaseIdx] = p; }

    void setSaturation(unsigned phaseIdx, const Scalar& S)
    { saturation_[phaseIdx] = S; }

    void setInvB(unsigned phaseIdx, const Scalar& b)
    { invB_[phaseIdx] = b; }

    void setDensity(unsigned phaseIdx, const Scalar& rho)
    { density_[phaseIdx] = rho; }

    void setRs(const Scalar& newRs)
    { Rs_ = newRs; }

    void setRv(const Scalar& newRv)
    { Rv_ = newRv; }

    const Scalar& pressure(unsigned phaseIdx) const
    { return pressure_[phaseIdx]; }

    const Scalar& saturation(unsigned phaseIdx) const
    { return saturation_[phaseIdx]; }

    const Scalar& temperature(unsigned phaseIdx OPM_UNUSED) const
    { return temperature_; }

    const Scalar& invB(unsigned phaseIdx) const
    { return invB_[phaseIdx]; }

    const Scalar& Rs() const
    { return Rs_; }

    const Scalar& Rv() const
    { return Rv_; }

    unsigned short pvtRegionIndex() const
    { return pvtRegionIdx_; }

    //////
    // slow methods
    //////
    Scalar density(unsigned phaseIdx) const
    { return density_[phaseIdx]; }

    Scalar molarDensity(unsigned phaseIdx) const
    {
        const auto& rho = density(phaseIdx);

        if (phaseIdx == waterPhaseIdx)
            return rho/FluidSystem::molarMass(waterCompIdx, pvtRegionIdx_);

        return
            rho*(moleFraction(phaseIdx, gasCompIdx)/FluidSystem::molarMass(gasCompIdx, pvtRegionIdx_)
                 + moleFraction(phaseIdx, oilCompIdx)/FluidSystem::molarMass(oilCompIdx, pvtRegionIdx_));

    }

    Scalar molarVolume(unsigned phaseIdx) const
    { return 1.0/molarDensity(phaseIdx); }

    Scalar viscosity(unsigned phaseIdx) const
    { return FluidSystem::viscosity(*this, phaseIdx, pvtRegionIdx_); }

    Scalar enthalpy(unsigned phaseIdx OPM_UNUSED) const
    {
        OPM_THROW(Opm::NotImplemented,
                  "The black-oil model does not support energy conservation yet.");
    }

    Scalar internalEnergy(unsigned phaseIdx OPM_UNUSED) const
    {
        OPM_THROW(Opm::NotImplemented,
                  "The black-oil model does not support energy conservation yet.");
    }

    Scalar massFraction(unsigned phaseIdx, unsigned compIdx) const
    {
        switch (phaseIdx) {
        case waterPhaseIdx:
            if (compIdx == waterCompIdx)
                return 1.0;
            return 0.0;

        case oilPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return 1.0 - FluidSystem::convertRsToXoG(Rs_, pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return FluidSystem::convertRsToXoG(Rs_, pvtRegionIdx_);
            }
            break;

        case gasPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return FluidSystem::convertRvToXgO(Rv_, pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return 1.0 - FluidSystem::convertRvToXgO(Rv_, pvtRegionIdx_);
            }
            break;
        }

        OPM_THROW(std::logic_error,
                  "Invalid phase or component index!");
    }

    Scalar moleFraction(unsigned phaseIdx, unsigned compIdx) const
    {
        switch (phaseIdx) {
        case waterPhaseIdx:
            if (compIdx == waterCompIdx)
                return 1.0;
            return 0.0;

        case oilPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return 1.0 - FluidSystem::convertXoGToxoG(FluidSystem::convertRsToXoG(Rs_, pvtRegionIdx_),
                                                          pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return FluidSystem::convertXoGToxoG(FluidSystem::convertRsToXoG(Rs_, pvtRegionIdx_),
                                                    pvtRegionIdx_);
            }
            break;

        case gasPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return FluidSystem::convertXgOToxgO(FluidSystem::convertRvToXgO(Rv_, pvtRegionIdx_),
                                                    pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return 1.0 - FluidSystem::convertXgOToxgO(FluidSystem::convertRvToXgO(Rv_, pvtRegionIdx_),
                                                          pvtRegionIdx_);
            }
            break;
        }

        OPM_THROW(std::logic_error,
                  "Invalid phase or component index!");
    }

    Scalar molarity(unsigned phaseIdx, unsigned compIdx) const
    { return moleFraction(phaseIdx, compIdx)*molarDensity(phaseIdx); }

    Scalar averageMolarMass(unsigned phaseIdx) const
    {
        Scalar result(0.0);
        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
            result += FluidSystem::molarMass(compIdx, pvtRegionIdx_)*moleFraction(phaseIdx, compIdx);
        return result;
    }

    Scalar fugacityCoefficient(unsigned phaseIdx, unsigned compIdx) const
    { return FluidSystem::fugacityCoefficient(*this, phaseIdx, compIdx, pvtRegionIdx_); }

    Scalar fugacity(unsigned phaseIdx, unsigned compIdx) const
    {
        return
            fugacityCoefficient(phaseIdx, compIdx)
            *moleFraction(phaseIdx, compIdx)
            *pressure(phaseIdx);
    }

private:
    static const Scalar temperature_;
    std::array<Scalar, numPhases> pressure_;
    std::array<Scalar, numPhases> saturation_;
    std::array<Scalar, numPhases> invB_;
    std::array<Scalar, numPhases> density_;
    Scalar Rs_;
    Scalar Rv_;
    unsigned short pvtRegionIdx_;
};

template <class Scalar, class FluidSystem>
const Scalar BlackOilFluidState<Scalar, FluidSystem>::temperature_ =
    FluidSystem::surfaceTemperature;

} // namespace Opm

#endif
