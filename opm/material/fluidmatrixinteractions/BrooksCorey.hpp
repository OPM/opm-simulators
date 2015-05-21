/*
  Copyright (C) 2009-2013 by Andreas Lauser
  Copyright (C) 2011 by Holger Class
  Copyright (C) 2010 by Philipp Nuske

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
 * \copydoc Opm::BrooksCorey
 */
#ifndef OPM_BROOKS_COREY_HPP
#define OPM_BROOKS_COREY_HPP

#include "BrooksCoreyParams.hpp"

#include <opm/material/common/MathToolbox.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Implementation of the Brooks-Corey capillary pressure <->
 *        saturation relation.
 *
 * This class provides the "raw" curves as static members and doesn't
 * concern itself converting absolute to effective saturations and
 * vice versa.
 *
 *\see BrooksCoreyParams
 */
template <class TraitsT, class ParamsT = BrooksCoreyParams<TraitsT> >
class BrooksCorey : public TraitsT
{
public:
    typedef TraitsT Traits;
    typedef ParamsT Params;
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases to which this material law applies.
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The Brooks-Corey capillary pressure law only applies "
                  "to the case of two fluid phases");

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

    static_assert(Traits::numPhases == 2,
                  "The number of fluid phases must be two if you want to use "
                  "this material law!");

    /*!
     * \brief The capillary pressure-saturation curves.
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container &values, const Params &params, const FluidState &fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = pcnw<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief Calculate the saturations of the phases starting from
     *        their pressure differences.
     */
    template <class Container, class FluidState>
    static void saturations(Container &values, const Params &params, const FluidState &fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = Sw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = 1 - values[Traits::wettingPhaseIdx];
    }

    /*!
     * \brief The relative permeability-saturation curves.
     *
     * \param values A random access container which stores the
     *               relative permeability of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the material law.
     * \param fs The fluid state for which the relative permeabilities
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container &values, const Params &params, const FluidState &fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = krw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = krn<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief The capillary pressure-saturation curve according to
     *        Brooks and Corey.
     *
     * The empirical Brooks-Corey capillary pressure-saturation
     * function is defined as
     * \f[
     * p_C = p_e\overline{S}_w^{-1/\lambda}
     * \f]
     *
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params &params, const FluidState &fs)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        const Evaluation& Sw =
            FsToolbox::template toLhs<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));

        assert(0 <= Sw && Sw <= 1);

        return twoPhaseSatPcnw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params &params, const Evaluation& Sw)
    {
        typedef MathToolbox<Evaluation> Toolbox;

        assert(0 <= Sw && Sw <= 1);

        return params.entryPressure()*Toolbox::pow(Sw, -1.0/params.lambda());
    }

    /*!
     * \brief The saturation-capillary pressure curve according to
     *        Brooks & Corey.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     \overline{S}_w = (\frac{p_C}{p_e})^{-\lambda}
     \f]
     *
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params &params, const FluidState &fs)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        Evaluation pC =
            FsToolbox::template toLhs<Evaluation>(fs.pressure(Traits::nonWettingPhaseIdx))
            - FsToolbox::template toLhs<Evaluation>(fs.pressure(Traits::wettingPhaseIdx));
        return twoPhaseSatSw(params, pC);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params &params, const Evaluation& pc)
    {
        typedef MathToolbox<Evaluation> Toolbox;

        assert(pc > 0); // if we don't assume that, std::pow will screw up!

        return Toolbox::pow(pc/params.entryPressure(), -params.lambda());
    }

    /*!
     * \brief Calculate the non-wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params &params, const FluidState &fs)
    { return 1 - Sw<FluidState, Evaluation>(params, fs); }

    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params &params, const Evaluation& pc)
    { return 1 - twoPhaseSatSw(params, pc); }
    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by the Brooks-Corey
     *        parameterization.
     *
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params &params, const FluidState &fs)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& Sw =
            FsToolbox::template toLhs<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));

        return twoPhaseSatKrw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params &params, const Evaluation& Sw)
    {
        typedef MathToolbox<Evaluation> Toolbox;

        assert(0 <= Sw && Sw <= 1);

        return Toolbox::pow(Sw, 2.0/params.lambda() + 3);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params &params, const FluidState &fs)
    {
        typedef MathToolbox<typename FluidState::Scalar> FsToolbox;

        const Evaluation& Sw =
            1.0 - FsToolbox::template toLhs<Evaluation>(fs.saturation(Traits::nonWettingPhaseIdx));

        return twoPhaseSatKrn(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params &params, const Evaluation& Sw)
    {
        typedef MathToolbox<Evaluation> Toolbox;

        assert(0 <= Sw && Sw <= 1);

        Scalar exponent = 2.0/params.lambda() + 1;
        const Evaluation Sn = 1. - Sw;
        return Sn*Sn*(1. - Toolbox::pow(Sw, exponent));
    }
};
} // namespace Opm

#endif // BROOKS_COREY_HPP
