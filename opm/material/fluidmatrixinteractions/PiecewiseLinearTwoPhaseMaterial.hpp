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
 * \copydoc Opm::PiecewiseLinearTwoPhaseMaterial
 */
#ifndef OPM_PIECEWISE_LINEAR_TWO_PHASE_MATERIAL_HPP
#define OPM_PIECEWISE_LINEAR_TWO_PHASE_MATERIAL_HPP

#include "PiecewiseLinearTwoPhaseMaterialParams.hpp"

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <algorithm>
#include <cmath>
#include <cassert>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Implementation of a tabulated, piecewise linear capillary
 *        pressure law.
 *
 * It would be equally possible to use cubic splines, but since the
 * ECLIPSE reservoir simulator uses linear interpolation for capillary
 * pressure and relperm curves, we do the same.
 */
template <class TraitsT, class ParamsT = PiecewiseLinearTwoPhaseMaterialParams<TraitsT> >
class PiecewiseLinearTwoPhaseMaterial : public TraitsT
{
    typedef typename ParamsT::ValueVector ValueVector;

public:
    //! The traits class for this material law
    typedef TraitsT Traits;

    //! The type of the parameter objects for this law
    typedef ParamsT Params;

    //! The type of the scalar values for this law
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The piecewise linear two-phase capillary pressure law only"
                  "applies to the case of two fluid phases");

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = true;

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = true;

    //! Specify whether the quantities defined by this material law
    //! are saturation dependent
    static const bool isSaturationDependent = true;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the absolute pressure
    static const bool isPressureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are temperature dependent
    static const bool isTemperatureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the phase composition
    static const bool isCompositionDependent = false;

    /*!
     * \brief The capillary pressure-saturation curve.
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = pcnw<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief The saturations of the fluid phases starting from their
     *        pressure differences.
     */
    template <class Container, class FluidState>
    static void saturations(Container& /* values */, const Params& /* params */, const FluidState& /* fs */)
    { throw std::logic_error("Not implemented: saturations()"); }

    /*!
     * \brief The relative permeabilities
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = krw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = krn<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params& params, const FluidState& fs)
    {
        const auto& Sw =
            Opm::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));

        return twoPhaseSatPcnw(params, Sw);
    }

    /*!
     * \brief The saturation-capillary pressure curve
     */
    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params& params, const Evaluation& Sw)
    { return eval_(params.SwPcwnSamples(), params.pcnwSamples(), Sw); }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnwInv(const Params& params, const Evaluation& pcnw)
    { return eval_(params.pcnwSamples(), params.SwPcwnSamples(), pcnw); }

    /*!
     * \brief The saturation-capillary pressure curve
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params& /* params */, const FluidState& /* fs */)
    { throw std::logic_error("Not implemented: Sw()"); }

    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params& /* params */, const Evaluation& /* pC */)
    { throw std::logic_error("Not implemented: twoPhaseSatSw()"); }

    /*!
     * \brief Calculate the non-wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params& params, const FluidState& fs)
    { return 1 - Sw<FluidState, Scalar>(params, fs); }

    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params& params, const Evaluation& pC)
    { return 1 - twoPhaseSatSw(params, pC); }

    /*!
     * \brief The relative permeability for the wetting phase of the
     *        porous medium
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params& params, const FluidState& fs)
    {
        const auto& Sw =
            Opm::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));

        return twoPhaseSatKrw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params& params, const Evaluation& Sw)
    { return eval_(params.SwKrwSamples(), params.krwSamples(), Sw); }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrwInv(const Params& params, const Evaluation& krw)
    { return eval_(params.krwSamples(), params.SwKrwSamples(), krw); }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the porous medium
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params& params, const FluidState& fs)
    {
        const auto& Sw =
            Opm::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));

        return twoPhaseSatKrn(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params& params, const Evaluation& Sw)
    { return eval_(params.SwKrnSamples(), params.krnSamples(), Sw); }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrnInv(const Params& params, const Evaluation& krn)
    { return eval_(params.krnSamples(), params.SwKrnSamples(), krn); }

private:
    template <class Evaluation>
    static Evaluation eval_(const ValueVector& xValues,
                            const ValueVector& yValues,
                            const Evaluation& x)
    {
        if (xValues.front() < xValues.back())
            return evalAscending_(xValues, yValues, x);
        return evalDescending_(xValues, yValues, x);
    }

    template <class Evaluation>
    static Evaluation evalAscending_(const ValueVector& xValues,
                                     const ValueVector& yValues,
                                     const Evaluation& x)
    {
        if (x <= xValues.front())
            return yValues.front();
        if (x >= xValues.back())
            return yValues.back();

        size_t segIdx = findSegmentIndex_(xValues, Opm::scalarValue(x));

        Scalar x0 = xValues[segIdx];
        Scalar x1 = xValues[segIdx + 1];

        Scalar y0 = yValues[segIdx];
        Scalar y1 = yValues[segIdx + 1];

        Scalar m = (y1 - y0)/(x1 - x0);

        return y0 + (x - x0)*m;
    }

    template <class Evaluation>
    static Evaluation evalDescending_(const ValueVector& xValues,
                                      const ValueVector& yValues,
                                      const Evaluation& x)
    {
        if (x >= xValues.front())
            return yValues.front();
        if (x <= xValues.back())
            return yValues.back();

        size_t segIdx = findSegmentIndexDescending_(xValues, Opm::scalarValue(x));

        Scalar x0 = xValues[segIdx];
        Scalar x1 = xValues[segIdx + 1];

        Scalar y0 = yValues[segIdx];
        Scalar y1 = yValues[segIdx + 1];

        Scalar m = (y1 - y0)/(x1 - x0);

        return y0 + (x - x0)*m;
    }

    template <class Evaluation>
    static Evaluation evalDeriv_(const ValueVector& xValues,
                                 const ValueVector& yValues,
                                 const Evaluation& x)
    {
        if (x <= xValues.front())
            return 0.0;
        if (x >= xValues.back())
            return 0.0;

        size_t segIdx = findSegmentIndex_(xValues, Opm::scalarValue(x));

        Scalar x0 = xValues[segIdx];
        Scalar x1 = xValues[segIdx + 1];

        Scalar y0 = yValues[segIdx];
        Scalar y1 = yValues[segIdx + 1];

        return (y1 - y0)/(x1 - x0);
    }

    static size_t findSegmentIndex_(const ValueVector& xValues, Scalar x)
    {
        assert(xValues.size() > 1); // we need at least two sampling points!
        size_t n = xValues.size() - 1;
        if (xValues.back() <= x)
            return n - 1;
        else if (x <= xValues.front())
            return 0;

        // bisection
        size_t lowIdx = 0, highIdx = n;
        while (lowIdx + 1 < highIdx) {
            size_t curIdx = (lowIdx + highIdx)/2;
            if (xValues[curIdx] < x)
                lowIdx = curIdx;
            else
                highIdx = curIdx;
        }

        return lowIdx;
    }

    static size_t findSegmentIndexDescending_(const ValueVector& xValues, Scalar x)
    {
        assert(xValues.size() > 1); // we need at least two sampling points!
        size_t n = xValues.size() - 1;
        if (x <= xValues.back())
            return n;
        else if (xValues.front() <= x)
            return 0;

        // bisection
        size_t lowIdx = 0, highIdx = n;
        while (lowIdx + 1 < highIdx) {
            size_t curIdx = (lowIdx + highIdx)/2;
            if (xValues[curIdx] >= x)
                lowIdx = curIdx;
            else
                highIdx = curIdx;
        }

        return lowIdx;
    }
};
} // namespace Opm

#endif
