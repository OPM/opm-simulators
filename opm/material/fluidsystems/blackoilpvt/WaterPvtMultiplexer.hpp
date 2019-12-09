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
 * \copydoc Opm::WaterPvtMultiplexer
 */
#ifndef OPM_WATER_PVT_MULTIPLEXER_HPP
#define OPM_WATER_PVT_MULTIPLEXER_HPP

#include "ConstantCompressibilityWaterPvt.hpp"
#include "WaterPvtThermal.hpp"

#define OPM_WATER_PVT_MULTIPLEXER_CALL(codeToCall)                      \
    switch (approach_) {                                                \
    case ConstantCompressibilityWaterPvt: {                             \
        auto& pvtImpl = getRealPvt<ConstantCompressibilityWaterPvt>();  \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case ThermalWaterPvt: {                                             \
        auto& pvtImpl = getRealPvt<ThermalWaterPvt>();                  \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case NoWaterPvt:                                                    \
        throw std::logic_error("Not implemented: Water PVT of this deck!"); \
    }

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the water
 *        phase in the black-oil model.
 */
template <class Scalar, bool enableThermal = true>
class WaterPvtMultiplexer
{
public:
    typedef Opm::WaterPvtThermal<Scalar> WaterPvtThermal;

    enum WaterPvtApproach {
        NoWaterPvt,
        ConstantCompressibilityWaterPvt,
        ThermalWaterPvt
    };

    WaterPvtMultiplexer()
    {
        approach_ = NoWaterPvt;
        realWaterPvt_ = nullptr;
    }

    WaterPvtMultiplexer(WaterPvtApproach approach, void* realWaterPvt)
        : approach_(approach)
        , realWaterPvt_(realWaterPvt)
    { }

    WaterPvtMultiplexer(const WaterPvtMultiplexer<Scalar,enableThermal>& data)
    {
        *this = data;
    }

    ~WaterPvtMultiplexer()
    {
        switch (approach_) {
        case ConstantCompressibilityWaterPvt: {
            delete &getRealPvt<ConstantCompressibilityWaterPvt>();
            break;
        }
        case ThermalWaterPvt: {
            delete &getRealPvt<ThermalWaterPvt>();
            break;
        }
        case NoWaterPvt:
            break;
        }
    }

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for water using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVDG keywords.
     */
    void initFromDeck(const Deck& deck, const EclipseState& eclState)
    {
        bool enableWater = deck.hasKeyword("WATER");
        if (!enableWater)
            return;

        if (enableThermal && (deck.hasKeyword("THERMAL") || deck.hasKeyword("TEMP")))
            setApproach(ThermalWaterPvt);
        else if (deck.hasKeyword("PVTW"))
            setApproach(ConstantCompressibilityWaterPvt);

        OPM_WATER_PVT_MULTIPLEXER_CALL(pvtImpl.initFromDeck(deck, eclState));
    }
#endif // HAVE_ECL_INPUT

    void initEnd()
    { OPM_WATER_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd()); }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions()); return 1; }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of gas given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                        const Evaluation& temperature,
                        const Evaluation& pressure) const
    { OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.internalEnergy(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure) const
    { OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.viscosity(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure) const
    { OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.inverseFormationVolumeFactor(regionIdx, temperature, pressure)); return 0; }

    void setApproach(WaterPvtApproach appr)
    {
        switch (appr) {
        case ConstantCompressibilityWaterPvt:
            realWaterPvt_ = new Opm::ConstantCompressibilityWaterPvt<Scalar>;
            break;

        case ThermalWaterPvt:
            realWaterPvt_ = new Opm::WaterPvtThermal<Scalar>;
            break;

        case NoWaterPvt:
            throw std::logic_error("Not implemented: Water PVT of this deck!");
        }

        approach_ = appr;
    }

    /*!
     * \brief Returns the concrete approach for calculating the PVT relations.
     *
     * (This is only determined at runtime.)
     */
    WaterPvtApproach approach() const
    { return approach_; }

    // get the concrete parameter object for the water phase
    template <WaterPvtApproach approachV>
    typename std::enable_if<approachV == ConstantCompressibilityWaterPvt, Opm::ConstantCompressibilityWaterPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Opm::ConstantCompressibilityWaterPvt<Scalar>* >(realWaterPvt_);
    }

    template <WaterPvtApproach approachV>
    typename std::enable_if<approachV == ConstantCompressibilityWaterPvt, const Opm::ConstantCompressibilityWaterPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<Opm::ConstantCompressibilityWaterPvt<Scalar>* >(realWaterPvt_);
    }

    template <WaterPvtApproach approachV>
    typename std::enable_if<approachV == ThermalWaterPvt, Opm::WaterPvtThermal<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Opm::WaterPvtThermal<Scalar>* >(realWaterPvt_);
    }

    template <WaterPvtApproach approachV>
    typename std::enable_if<approachV == ThermalWaterPvt, const Opm::WaterPvtThermal<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<Opm::WaterPvtThermal<Scalar>* >(realWaterPvt_);
    }

    const void* realWaterPvt() const { return realWaterPvt_; }

    bool operator==(const WaterPvtMultiplexer<Scalar,enableThermal>& data) const
    {
        if (this->approach() != data.approach())
            return false;

        switch (approach_) {
        case ConstantCompressibilityWaterPvt:
            return *static_cast<const Opm::ConstantCompressibilityWaterPvt<Scalar>*>(realWaterPvt_) ==
                   *static_cast<const Opm::ConstantCompressibilityWaterPvt<Scalar>*>(data.realWaterPvt_);
        case ThermalWaterPvt:
            return *static_cast<const Opm::WaterPvtThermal<Scalar>*>(realWaterPvt_) ==
                   *static_cast<const Opm::WaterPvtThermal<Scalar>*>(data.realWaterPvt_);
        default:
            return true;
        }
    }

    WaterPvtMultiplexer<Scalar,enableThermal>& operator=(const WaterPvtMultiplexer<Scalar,enableThermal>& data)
    {
        approach_ = data.approach_;
        switch (approach_) {
        case ConstantCompressibilityWaterPvt:
            realWaterPvt_ = new Opm::ConstantCompressibilityWaterPvt<Scalar>(*static_cast<const Opm::ConstantCompressibilityWaterPvt<Scalar>*>(data.realWaterPvt_));
            break;
        case ThermalWaterPvt:
            realWaterPvt_ = new Opm::WaterPvtThermal<Scalar>(*static_cast<const Opm::WaterPvtThermal<Scalar>*>(data.realWaterPvt_));
            break;
        default:
            break;
        }

        return *this;
    }

private:
    WaterPvtApproach approach_;
    void* realWaterPvt_;
};

#undef OPM_WATER_PVT_MULTIPLEXER_CALL

} // namespace Opm

#endif
