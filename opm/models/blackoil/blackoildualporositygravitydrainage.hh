// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 TNO

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
/*!
 * \file
 *
 * \brief Gravity drainage between a matrix cell and its fracture twin in the
 *        dual-porosity models.
 *
 * Each continuum holds, per cell, a fraction of its height filled with mobile
 * water and a fraction filled with mobile gas, estimated from the saturation
 * with a vertical equilibrium picture.  The difference of these fractions
 * between a matrix cell and its fracture twin, times the matrix block height
 * and the density difference of the two phases, is a gravity head that joins
 * the matrix-fracture flow of each phase.
 */
#ifndef OPM_BLACKOIL_DUAL_POROSITY_GRAVITY_DRAINAGE_HH
#define OPM_BLACKOIL_DUAL_POROSITY_GRAVITY_DRAINAGE_HH

#include <vector>

namespace Opm {

namespace DualPorosityFractions {

//! Saturation end points that bound the mobile water fraction.
struct WaterFractionEndPoints
{
    double swco;    //!< connate water saturation
    double swcr;    //!< critical water saturation
    double scohy;   //!< connate hydrocarbon saturation
    double scrhy;   //!< critical hydrocarbon saturation
};

//! Saturation end points that bound the mobile gas fraction.
struct GasFractionEndPoints
{
    double sgco;    //!< connate gas saturation
    double sgcr;    //!< critical gas saturation
    double slco;    //!< connate liquid saturation
    double slcr;    //!< critical liquid saturation
};

//! Clamp to [0, 1].  The current-fraction functions are templated on the
//! saturation type, so automatic-differentiation types flow through them.
template <class Value>
constexpr Value clampFraction(const Value& x) noexcept
{
    if (x < 0.0) {
        return Value(0.0);
    }
    if (x > 1.0) {
        return Value(1.0);
    }
    return x;
}

//! Fraction of the cell initially below the water contact.
constexpr double initialWaterFraction(const double swi,
                                      const WaterFractionEndPoints& ep) noexcept
{
    const double denom = 1.0 - ep.scohy - ep.swco;
    if (denom <= 0.0) {
        return 0.0;
    }
    return clampFraction((swi - ep.swco) / denom);
}

//! Current fraction of the cell filled with mobile water.  The water level
//! rises from the initial contact between the connate water and the critical
//! hydrocarbon saturations, and falls between the critical water and the
//! connate hydrocarbon saturations.
template <class Value>
constexpr Value waterFraction(const Value& sw,
                              const double swi,
                              const double xwi,
                              const WaterFractionEndPoints& ep) noexcept
{
    const double denom = (sw >= swi)
        ? 1.0 - ep.scrhy - ep.swco
        : 1.0 - ep.scohy - ep.swcr;
    if (denom <= 0.0) {
        return Value(xwi);
    }
    const double offset = (sw >= swi)
        ? xwi * (ep.scrhy - ep.scohy)
        : xwi * (ep.swcr - ep.swco);
    return clampFraction(Value((sw - offset - ep.swco) / denom));
}

//! Fraction of the cell initially above the gas contact.
constexpr double initialGasFraction(const double sgi,
                                    const GasFractionEndPoints& ep) noexcept
{
    const double denom = 1.0 - ep.slco - ep.sgco;
    if (denom <= 0.0) {
        return 0.0;
    }
    return clampFraction((sgi - ep.sgco) / denom);
}

//! Current fraction of the cell filled with mobile gas, as waterFraction()
//! for the gas level.
template <class Value>
constexpr Value gasFraction(const Value& sg,
                            const double sgi,
                            const double xgi,
                            const GasFractionEndPoints& ep) noexcept
{
    const double denom = (sg >= sgi)
        ? 1.0 - ep.slcr - ep.sgco
        : 1.0 - ep.slco - ep.sgcr;
    if (denom <= 0.0) {
        return Value(xgi);
    }
    const double offset = (sg >= sgi)
        ? xgi * (ep.slcr - ep.slco)
        : xgi * (ep.sgcr - ep.sgco);
    return clampFraction(Value((sg - offset - ep.sgco) / denom));
}

//! Gas-oil gravity head between two cells: positive when cell a holds the
//! larger mobile gas fraction.
template <class Value>
constexpr Value gasOilGravityHead(const double gravity,
                                  const double dzMatrix,
                                  const Value& rhoOil,
                                  const Value& rhoGas,
                                  const Value& xgA,
                                  const Value& xgB) noexcept
{
    return gravity * dzMatrix * (rhoOil - rhoGas) * (xgA - xgB);
}

//! Water-oil gravity head between two cells: positive when cell a holds the
//! larger mobile water fraction.
template <class Value>
constexpr Value waterOilGravityHead(const double gravity,
                                    const double dzMatrix,
                                    const Value& rhoWater,
                                    const Value& rhoOil,
                                    const Value& xwA,
                                    const Value& xwB) noexcept
{
    return gravity * dzMatrix * (rhoWater - rhoOil) * (xwA - xwB);
}

//! Gravity head of each phase.
template <class Value>
struct PhaseGravityHeads
{
    Value oil;
    Value gas;
    Value water;
};

//! Share the two pair heads among the phases: each pair head is shared half
//! and half between its two phases, and the oil phase takes part in both.
template <class Value>
constexpr PhaseGravityHeads<Value> phaseGravityHeads(const Value& gasOilHead,
                                                     const Value& waterOilHead) noexcept
{
    return { -0.5 * (gasOilHead + waterOilHead),
             0.5 * gasOilHead,
             0.5 * waterOilHead };
}

} // namespace DualPorosityFractions

/*!
 * \brief Per-cell data of the dual-porosity gravity drainage: the end points
 *        of the mobile fractions, the initial saturations and fractions, and
 *        the matrix block height.
 *
 * Filled on the host when the run selects a gravity drainage model, over the
 * cells of the process (a matrix cell and its fracture twin are on the same
 * process); empty otherwise.
 */
template<class Scalar>
struct DualPorosityGravityDrainageParam
{
    struct CellData
    {
        DualPorosityFractions::WaterFractionEndPoints water{};
        DualPorosityFractions::GasFractionEndPoints gas{};
        Scalar swi = 0.0;       //!< initial water saturation
        Scalar sgi = 0.0;       //!< initial gas saturation
        Scalar xwi = 0.0;       //!< initial mobile water fraction
        Scalar xgi = 0.0;       //!< initial mobile gas fraction
        Scalar dzMatrix = 0.0;  //!< matrix block height (DZMTRXV)
    };

    Scalar gravity = 0.0;
    std::vector<CellData> cell;
};

//! The module parameters of problems without dual-porosity gravity drainage
//! (and of the device copies of the problems) carry no data for it.
struct NoDualPorosityGravityDrainageParam
{};

} // namespace Opm

#endif // OPM_BLACKOIL_DUAL_POROSITY_GRAVITY_DRAINAGE_HH
