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
 * \copydoc Opm::OilPvtMultiplexer
 */
#ifndef OPM_OIL_PVT_MULTIPLEXER_HPP
#define OPM_OIL_PVT_MULTIPLEXER_HPP

#include "ConstantCompressibilityOilPvt.hpp"
#include "DeadOilPvt.hpp"
#include "LiveOilPvt.hpp"
#include "OilPvtThermal.hpp"

namespace Opm {
#define OPM_OIL_PVT_MULTIPLEXER_CALL(codeToCall)                        \
    switch (approach_) {                                                \
    case ConstantCompressibilityOilPvt: {                               \
        auto& pvtImpl = getRealPvt<ConstantCompressibilityOilPvt>();    \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case DeadOilPvt: {                                                  \
        auto& pvtImpl = getRealPvt<DeadOilPvt>();                       \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case LiveOilPvt: {                                                  \
        auto& pvtImpl = getRealPvt<LiveOilPvt>();                       \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case ThermalOilPvt: {                                               \
        auto& pvtImpl = getRealPvt<ThermalOilPvt>();                    \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case NoOilPvt:                                                      \
        throw std::logic_error("Not implemented: Oil PVT of this deck!"); \
    }

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil
 *        phase in the black-oil model.
 *
 * This is the base class which which provides an API for the actual PVT implementation
 * classes which based on dynamic polymorphism. The rationale to use dynamic polymorphism
 * here is that this enables the fluid system to easily switch the used PVT relations for
 * the individual fluid phases.
 *
 * Note that, since the application for this class is the black-oil fluid system, the API
 * exposed by this class is pretty specific to the black-oil model.
 */
template <class Scalar, bool enableThermal = true>
class OilPvtMultiplexer
{
public:
    typedef Opm::OilPvtThermal<Scalar> OilPvtThermal;

    enum OilPvtApproach {
        NoOilPvt,
        LiveOilPvt,
        DeadOilPvt,
        ConstantCompressibilityOilPvt,
        ThermalOilPvt
    };

    OilPvtMultiplexer()
    {
        approach_ = NoOilPvt;
        realOilPvt_ = nullptr;
    }

    OilPvtMultiplexer(OilPvtApproach approach, void* realOilPvt)
        : approach_(approach)
        , realOilPvt_(realOilPvt)
    { }

    OilPvtMultiplexer(const OilPvtMultiplexer<Scalar,enableThermal>& data)
    {
        *this = data;
    }

    ~OilPvtMultiplexer()
    {
        switch (approach_) {
        case LiveOilPvt: {
            delete &getRealPvt<LiveOilPvt>();
            break;
        }
        case DeadOilPvt: {
            delete &getRealPvt<DeadOilPvt>();
            break;
        }
        case ConstantCompressibilityOilPvt: {
            delete &getRealPvt<ConstantCompressibilityOilPvt>();
            break;
        }
        case ThermalOilPvt: {
            delete &getRealPvt<ThermalOilPvt>();
            break;
        }

        case NoOilPvt:
            break;
        }
    }

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the parameters for water using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVTO/PVDO/PVCDO keywords.
     */
    void initFromDeck(const Deck& deck, const EclipseState& eclState)
    {
        bool enableOil = deck.hasKeyword("OIL");
        if (!enableOil)
            return;

        if (enableThermal && (deck.hasKeyword("THERMAL") || deck.hasKeyword("TEMP")))
            setApproach(ThermalOilPvt);
        else if (deck.hasKeyword("PVCDO"))
            setApproach(ConstantCompressibilityOilPvt);
        else if (deck.hasKeyword("PVDO"))
            setApproach(DeadOilPvt);
        else if (deck.hasKeyword("PVTO"))
            setApproach(LiveOilPvt);

        OPM_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.initFromDeck(deck, eclState));
    }
#endif // HAVE_ECL_INPUT


    void initEnd()
    { OPM_OIL_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd()); }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.numRegions()); return 1; }

    /*!
     * \brief Returns the specific enthalpy [J/kg] oil given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx,
                        const Evaluation& temperature,
                        const Evaluation& pressure,
                        const Evaluation& Rs) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.internalEnergy(regionIdx, temperature, pressure, Rs)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& Rs) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.viscosity(regionIdx, temperature, pressure, Rs)); return 0; }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedViscosity(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& temperature,
                                            const Evaluation& pressure,
                                            const Evaluation& Rs) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.inverseFormationVolumeFactor(regionIdx, temperature, pressure, Rs)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& temperature,
                                                     const Evaluation& pressure) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedInverseFormationVolumeFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of saturated oil.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedGasDissolutionFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of saturated oil.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& temperature,
                                             const Evaluation& pressure,
                                             const Evaluation& oilSaturation,
                                             const Evaluation& maxOilSaturation) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedGasDissolutionFactor(regionIdx, temperature, pressure, oilSaturation, maxOilSaturation)); return 0; }

    /*!
     * \brief Returns the saturation pressure [Pa] of oil given the mass fraction of the
     *        gas component in the oil phase.
     *
     * Calling this method only makes sense for live oil. All other implementations of
     * the black-oil PVT interface will just throw an exception...
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& Rs) const
    { OPM_OIL_PVT_MULTIPLEXER_CALL(return pvtImpl.saturationPressure(regionIdx, temperature, Rs)); return 0; }

    void setApproach(OilPvtApproach appr)
    {
        switch (appr) {
        case LiveOilPvt:
            realOilPvt_ = new Opm::LiveOilPvt<Scalar>;
            break;

        case DeadOilPvt:
            realOilPvt_ = new Opm::DeadOilPvt<Scalar>;
            break;

        case ConstantCompressibilityOilPvt:
            realOilPvt_ = new Opm::ConstantCompressibilityOilPvt<Scalar>;
            break;

        case ThermalOilPvt:
            realOilPvt_ = new Opm::OilPvtThermal<Scalar>;
            break;

        case NoOilPvt:
            throw std::logic_error("Not implemented: Oil PVT of this deck!");
        }

        approach_ = appr;
    }

    /*!
     * \brief Returns the concrete approach for calculating the PVT relations.
     *
     * (This is only determined at runtime.)
     */
    OilPvtApproach approach() const
    { return approach_; }

    // get the concrete parameter object for the oil phase
    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == LiveOilPvt, Opm::LiveOilPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Opm::LiveOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == LiveOilPvt, const Opm::LiveOilPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<Opm::LiveOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == DeadOilPvt, Opm::DeadOilPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Opm::DeadOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == DeadOilPvt, const Opm::DeadOilPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<Opm::DeadOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == ConstantCompressibilityOilPvt, Opm::ConstantCompressibilityOilPvt<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Opm::ConstantCompressibilityOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == ConstantCompressibilityOilPvt, const Opm::ConstantCompressibilityOilPvt<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<Opm::ConstantCompressibilityOilPvt<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == ThermalOilPvt, Opm::OilPvtThermal<Scalar> >::type& getRealPvt()
    {
        assert(approach() == approachV);
        return *static_cast<Opm::OilPvtThermal<Scalar>* >(realOilPvt_);
    }

    template <OilPvtApproach approachV>
    typename std::enable_if<approachV == ThermalOilPvt, const Opm::OilPvtThermal<Scalar> >::type& getRealPvt() const
    {
        assert(approach() == approachV);
        return *static_cast<const Opm::OilPvtThermal<Scalar>* >(realOilPvt_);
    }

    const void* realOilPvt() const { return realOilPvt_; }

    bool operator==(const OilPvtMultiplexer<Scalar,enableThermal>& data) const
    {
        if (this->approach() != data.approach())
            return false;

        switch (approach_) {
        case ConstantCompressibilityOilPvt:
            return *static_cast<const Opm::ConstantCompressibilityOilPvt<Scalar>*>(realOilPvt_) ==
                   *static_cast<const Opm::ConstantCompressibilityOilPvt<Scalar>*>(data.realOilPvt_);
        case DeadOilPvt:
            return *static_cast<const Opm::DeadOilPvt<Scalar>*>(realOilPvt_) ==
                   *static_cast<const Opm::DeadOilPvt<Scalar>*>(data.realOilPvt_);
        case LiveOilPvt:
            return *static_cast<const Opm::LiveOilPvt<Scalar>*>(realOilPvt_) ==
                   *static_cast<const Opm::LiveOilPvt<Scalar>*>(data.realOilPvt_);
        case ThermalOilPvt:
            return *static_cast<const Opm::OilPvtThermal<Scalar>*>(realOilPvt_) ==
                   *static_cast<const Opm::OilPvtThermal<Scalar>*>(data.realOilPvt_);
        default:
            return true;
        }
    }

    OilPvtMultiplexer<Scalar,enableThermal>& operator=(const OilPvtMultiplexer<Scalar,enableThermal>& data)
    {
        approach_ = data.approach_;
        switch (approach_) {
        case ConstantCompressibilityOilPvt:
            realOilPvt_ = new Opm::ConstantCompressibilityOilPvt<Scalar>(*static_cast<const Opm::ConstantCompressibilityOilPvt<Scalar>*>(data.realOilPvt_));
            break;
        case DeadOilPvt:
            realOilPvt_ = new Opm::DeadOilPvt<Scalar>(*static_cast<const Opm::DeadOilPvt<Scalar>*>(data.realOilPvt_));
            break;
        case LiveOilPvt:
            realOilPvt_ = new Opm::LiveOilPvt<Scalar>(*static_cast<const Opm::LiveOilPvt<Scalar>*>(data.realOilPvt_));
            break;
        case ThermalOilPvt:
            realOilPvt_ = new Opm::OilPvtThermal<Scalar>(*static_cast<const Opm::OilPvtThermal<Scalar>*>(data.realOilPvt_));
            break;
        default:
            break;
        }

        return *this;
    }

private:
    OilPvtApproach approach_;
    void* realOilPvt_;
};

#undef OPM_OIL_PVT_MULTIPLEXER_CALL

} // namespace Opm

#endif
