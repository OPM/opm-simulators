// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser

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
*/
/*!
 * \file
 * \copydoc Opm::OilPvtMultiplexer
 */
#ifndef OPM_GAS_PVT_MULTIPLEXER_HPP
#define OPM_GAS_PVT_MULTIPLEXER_HPP

#include "DryGasPvt.hpp"
#include "WetGasPvt.hpp"

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#endif

namespace Opm {
template <class Scalar>
class OilPvtMultiplexer;

#define OPM_GAS_PVT_MULTIPLEXER_CALL(codeToCall)                        \
    switch (gasPvtApproach_) {                                          \
    case DryGasPvt: {                                                   \
        auto &pvtImpl = getRealGasPvt<DryGasPvt>();                     \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case WetGasPvt: {                                                   \
        auto &pvtImpl = getRealGasPvt<WetGasPvt>();                     \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case NoGasPvt:                                                      \
        OPM_THROW(std::logic_error, "Not implemented: Gas PVT of this deck!"); \
    }


/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas
 *        phase in the black-oil model.
 *
 * This is a multiplexer class which forwards all calls to the real implementation.
 *
 * Note that, since the main application for this class is the black oil fluid system,
 * the API exposed by this class is pretty specific to the assumptions made by the black
 * oil model.
 */
template <class Scalar>
class GasPvtMultiplexer
{
    typedef Opm::OilPvtMultiplexer<Scalar> OilPvtMultiplexer;

public:
    enum GasPvtApproach {
        NoGasPvt,
        DryGasPvt,
        WetGasPvt
    };

    GasPvtMultiplexer()
    {
        gasPvtApproach_ = NoGasPvt;
    }

    ~GasPvtMultiplexer()
    {
        switch (gasPvtApproach_) {
        case DryGasPvt: {
            delete &getRealGasPvt<DryGasPvt>();
            break;
        }
        case WetGasPvt: {
            delete &getRealGasPvt<WetGasPvt>();
            break;
        }
        case NoGasPvt:
            break;
        }
    }

#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize the parameters for gas using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVDG keywords.
     */
    void initFromDeck(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        if (deck->hasKeyword("PVTG"))
            setApproach(WetGasPvt);
        else if (deck->hasKeyword("PVDG"))
            setApproach(DryGasPvt);

        OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initFromDeck(deck, eclState));
    }
#endif // HAVE_OPM_PARSER

    void setApproach(GasPvtApproach gasPvtApproach)
    {
        switch (gasPvtApproach) {
        case DryGasPvt:
            realGasPvt_ = new Opm::DryGasPvt<Scalar>;
            break;

        case WetGasPvt:
            realGasPvt_ = new Opm::WetGasPvt<Scalar>;
            break;

        case NoGasPvt:
            OPM_THROW(std::logic_error, "Not implemented: Gas PVT of this deck!");
        }

        gasPvtApproach_ = gasPvtApproach;
    }

    void initEnd(const OilPvtMultiplexer *oilPvt)
    { OPM_GAS_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd(oilPvt)); }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation = Scalar>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& XgO) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.viscosity(regionIdx, temperature, pressure, XgO)); return 0; }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation = Scalar>
    Evaluation formationVolumeFactor(unsigned regionIdx,
                                     const Evaluation& temperature,
                                     const Evaluation& pressure,
                                     const Evaluation& XgO) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.formationVolumeFactor(regionIdx, temperature, pressure, XgO)); return 0; }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    template <class Evaluation = Scalar>
    Evaluation density(unsigned regionIdx,
                       const Evaluation& temperature,
                       const Evaluation& pressure,
                       const Evaluation& XgO) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.density(regionIdx, temperature, pressure, XgO)); return 0; }

    /*!
     * \brief Returns the fugacity coefficient [Pa] of the gas component in the gas phase
     *        given a set of parameters.
     */
    template <class Evaluation = Scalar>
    Evaluation fugacityCoefficientGas(unsigned regionIdx,
                                      const Evaluation& temperature,
                                      const Evaluation& pressure) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.fugacityCoefficientGas(regionIdx, temperature, pressure)); return 0; }

    template <class Evaluation = Scalar>
    Evaluation fugacityCoefficientOil(unsigned regionIdx,
                                      const Evaluation& temperature,
                                      const Evaluation& pressure) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.fugacityCoefficientOil(regionIdx, temperature, pressure)); return 0; }

    template <class Evaluation = Scalar>
    Evaluation fugacityCoefficientWater(unsigned regionIdx,
                                        const Evaluation& temperature,
                                        const Evaluation& pressure) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.fugacityCoefficientWater(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ [m^3/m^3] of oil saturated gas.
     */
    template <class Evaluation = Scalar>
    Evaluation oilVaporizationFactor(unsigned regionIdx,
                                     const Evaluation& temperature,
                                     const Evaluation& pressure) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.oilVaporizationFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the oil component
     *
     * \param XgO The mass fraction of the oil component in the gas phase [-]
     */
    template <class Evaluation = Scalar>
    Evaluation gasSaturationPressure(unsigned regionIdx,
                                     const Evaluation& temperature,
                                     const Evaluation& XgO) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.gasSaturationPressure(regionIdx, temperature, XgO)); return 0; }

    /*!
     * \brief Returns the gas mass fraction of oil-saturated gas at a given temperatire
     *        and pressure [-].
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedGasOilMassFraction(unsigned regionIdx,
                                           const Evaluation& temperature,
                                           const Evaluation& pressure) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedGasOilMassFraction(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the gas mole fraction of oil-saturated gas at a given temperatire
     *        and pressure [-].
     */
    template <class Evaluation = Scalar>
    Evaluation saturatedGasOilMoleFraction(unsigned regionIdx,
                                           const Evaluation& temperature,
                                           const Evaluation& pressure) const
    { OPM_GAS_PVT_MULTIPLEXER_CALL(return pvtImpl.saturatedGasOilMoleFraction(regionIdx, temperature, pressure)); return 0; }


    GasPvtApproach gasPvtApproach() const
    { return gasPvtApproach_; }

    // get the parameter object for the dry gas case
    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == DryGasPvt, Opm::DryGasPvt<Scalar> >::type& getRealGasPvt()
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<Opm::DryGasPvt<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == DryGasPvt, const Opm::DryGasPvt<Scalar> >::type& getRealGasPvt() const
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<const Opm::DryGasPvt<Scalar>* >(realGasPvt_);
    }

    // get the parameter object for the wet gas case
    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == WetGasPvt, Opm::WetGasPvt<Scalar> >::type& getRealGasPvt()
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<Opm::WetGasPvt<Scalar>* >(realGasPvt_);
    }

    template <GasPvtApproach approachV>
    typename std::enable_if<approachV == WetGasPvt, const Opm::WetGasPvt<Scalar> >::type& getRealGasPvt() const
    {
        assert(gasPvtApproach() == approachV);
        return *static_cast<const Opm::WetGasPvt<Scalar>* >(realGasPvt_);
    }

private:
    GasPvtApproach gasPvtApproach_;
    void *realGasPvt_;
};

#undef OPM_GAS_MULTIPLEXER_CALL

} // namespace Opm

#endif
