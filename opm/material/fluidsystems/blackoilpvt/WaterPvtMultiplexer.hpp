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
 * \copydoc Opm::WaterPvtMultiplexer
 */
#ifndef OPM_WATER_PVT_MULTIPLEXER_HPP
#define OPM_WATER_PVT_MULTIPLEXER_HPP

#include "ConstantCompressibilityWaterPvt.hpp"

#define OPM_WATER_PVT_MULTIPLEXER_CALL(codeToCall)                      \
    switch (approach_) {                                                \
    case ConstantCompressibilityWaterPvt: {                             \
        auto &pvtImpl = getRealPvt<ConstantCompressibilityWaterPvt>();  \
        codeToCall;                                                     \
        break;                                                          \
    }                                                                   \
    case NoWaterPvt:                                                    \
        OPM_THROW(std::logic_error, "Not implemented: Water PVT of this deck!"); \
    }

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the water
 *        phase in the black-oil model.
 */
template <class Scalar>
class WaterPvtMultiplexer
{
public:
    enum WaterPvtApproach {
        NoWaterPvt,
        ConstantCompressibilityWaterPvt
    };

    WaterPvtMultiplexer()
    {
        approach_ = NoWaterPvt;
    }

    ~WaterPvtMultiplexer()
    {
        switch (approach_) {
        case ConstantCompressibilityWaterPvt: {
            delete &getRealPvt<ConstantCompressibilityWaterPvt>();
            break;
        }
        case NoWaterPvt:
            break;
        }
    }

#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize the parameters for water using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVDG keywords.
     */
    void initFromDeck(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        if (deck->hasKeyword("PVTW"))
            setApproach(ConstantCompressibilityWaterPvt);

        OPM_WATER_PVT_MULTIPLEXER_CALL(pvtImpl.initFromDeck(deck, eclState));
    }
#endif // HAVE_OPM_PARSER

    void initEnd()
    { OPM_WATER_PVT_MULTIPLEXER_CALL(pvtImpl.initEnd()); }

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
    Evaluation formationVolumeFactor(unsigned regionIdx,
                                     const Evaluation& temperature,
                                     const Evaluation& pressure) const
    { OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.formationVolumeFactor(regionIdx, temperature, pressure)); return 0; }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation density(unsigned regionIdx,
                       const Evaluation& temperature,
                       const Evaluation& pressure) const
    { OPM_WATER_PVT_MULTIPLEXER_CALL(return pvtImpl.density(regionIdx, temperature, pressure)); return 0; }

    void setApproach(WaterPvtApproach approach)
    {
        switch (approach) {
        case ConstantCompressibilityWaterPvt:
            realWaterPvt_ = new Opm::ConstantCompressibilityWaterPvt<Scalar>;
            break;

        case NoWaterPvt:
            OPM_THROW(std::logic_error, "Not implemented: Water PVT of this deck!");
        }

        approach_ = approach;
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

private:
    WaterPvtApproach approach_;
    void *realWaterPvt_;
};

#undef OPM_WATER_PVT_MULTIPLEXER_CALL

} // namespace Opm

#endif
