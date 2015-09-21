// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2009-2013 by Andreas Lauser

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
 * \copydoc Opm::SplineTwoPhaseMaterial
 */
#ifndef OPM_SPLINE_TWO_PHASE_MATERIAL_HPP
#define OPM_SPLINE_TWO_PHASE_MATERIAL_HPP

#include "SplineTwoPhaseMaterialParams.hpp"

#include <opm/material/common/ErrorMacros.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <algorithm>
#include <cmath>
#include <cassert>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Implementation of a tabulated capillary pressure and relperm law which uses
 *        spline curves as interpolation functions.
 */
template <class TraitsT, class ParamsT = SplineTwoPhaseMaterialParams<TraitsT> >
class SplineTwoPhaseMaterial : public TraitsT
{
    typedef typename ParamsT::SamplePoints SamplePoints;

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
    static void capillaryPressures(Container &values, const Params &params, const FluidState &fluidState)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = pcnw<FluidState, Evaluation>(params, fluidState);
    }

    /*!
     * \brief The saturations of the fluid phases starting from their
     *        pressure differences.
     */
    template <class Container, class FluidState>
    static void saturations(Container &values, const Params &params, const FluidState &fluidState)
    { OPM_THROW(std::logic_error, "Not implemented: saturations()"); }

    /*!
     * \brief The relative permeabilities
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container &values, const Params &params, const FluidState &fluidState)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = krw<FluidState, Evaluation>(params, fluidState);
        values[Traits::nonWettingPhaseIdx] = krn<FluidState, Evaluation>(params, fluidState);
    }

    /*!
     * \brief The capillary pressure-saturation curve
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params &params, const FluidState &fluidState)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        const Evaluation& Sw =
            FsToolbox::template toLhs<Evaluation>(fluidState.saturation(Traits::wettingPhaseIdx));

        return twoPhaseSatPcnw(params, Sw);
    }

    /*!
     * \brief The saturation-capillary pressure curve
     */
    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params &params, const Evaluation& Sw)
    {
        // this assumes that the capillary pressure is monotonically decreasing
        if (Sw <= params.pcnwSpline().xMin())
            return Evaluation(params.pcnwSpline().yFirst());
        if (Sw >= params.pcnwSpline().xMax())
            return Evaluation(params.pcnwSpline().yLast());

        return params.pcnwSpline().eval(Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnwInv(const Params &params, const Evaluation& pcnw)
    {
        static const Evaluation nil(0.0);

        // this assumes that the capillary pressure is monotonically decreasing
        if (pcnw >= params.pcnwSpline().yFirst())
            return Evaluation(params.pcnwSpline().xMin());
        if (pcnw <= params.pcnwSpline().yLast())
            return Evaluation(params.pcnwSpline().xMax());

        // the intersect() method of splines is a bit slow, but this code path is not too
        // time critical...
        return params.pcnwSpline().intersect(/*a=*/nil, /*b=*/nil, /*c=*/nil, /*d=*/pcnw);
    }

    /*!
     * \brief The saturation-capillary pressure curve
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params &params, const FluidState &fluidState)
    { OPM_THROW(std::logic_error, "Not implemented: Sw()"); }

    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params &params, const Evaluation& pC)
    { OPM_THROW(std::logic_error, "Not implemented: twoPhaseSatSw()"); }

    /*!
     * \brief Calculate the non-wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params &params, const FluidState &fluidState)
    { return 1 - Sw<FluidState, Evaluation>(params, fluidState); }

    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params &params, const Evaluation& pC)
    { return 1 - twoPhaseSatSw(params, pC); }

    /*!
     * \brief The relative permeability for the wetting phase of the
     *        porous medium
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params &params, const FluidState &fluidState)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        const Evaluation& Sw =
            FsToolbox::template toLhs<Evaluation>(fluidState.saturation(Traits::wettingPhaseIdx));

        return twoPhaseSatKrw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params &params, const Evaluation& Sw)
    {
        if (Sw <= params.krnSpline().xMin())
            return Evaluation(params.krwSpline().yFirst());
        if (Sw >= params.krnSpline().xMax())
            return Evaluation(params.krwSpline().yLast());

        return params.krwSpline().eval(Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrwInv(const Params &params, const Evaluation& krw)
    {
        static const Evaluation nil(0.0);

        if (krw <= params.krwSpline().yFirst())
            return Evaluation(params.krwSpline().xMin());
        if (krw >= params.krwSpline().yLast())
            return Evaluation(params.krwSpline().xMax());

        return params.krwSpline().intersect(/*a=*/nil, /*b=*/nil, /*c=*/nil, /*d=*/krw);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the porous medium
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params &params, const FluidState &fluidState)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        const Evaluation& Sn =
            FsToolbox::template toLhs<Evaluation>(fluidState.saturation(Traits::nonWettingPhaseIdx));

        return twoPhaseSatKrn(params, 1.0 - Sn);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params &params, const Evaluation& Sw)
    {
        if (Sw <= params.krnSpline().xMin())
            return Evaluation(params.krnSpline().yFirst());
        if (Sw >= params.krnSpline().xMax())
            return Evaluation(params.krnSpline().yLast());

        return params.krnSpline().eval(Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrnInv(const Params &params, const Evaluation& krn)
    {
        static const Evaluation nil(0.0);

        if (krn >= params.krnSpline().yFirst())
            return Evaluation(params.krnSpline().xMin());
        if (krn <= params.krnSpline().yLast())
            return Evaluation(params.krnSpline().xMax());

        return params.krnSpline().intersect(/*a=*/nil, /*b=*/nil, /*c=*/nil, /*d=*/krn);
    }
};
} // namespace Opm

#endif
