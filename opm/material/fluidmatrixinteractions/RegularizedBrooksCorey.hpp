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
 * \copydoc Opm::RegularizedBrooksCorey
 */
#ifndef REGULARIZED_BROOKS_COREY_HPP
#define REGULARIZED_BROOKS_COREY_HPP

#include "BrooksCorey.hpp"
#include "RegularizedBrooksCoreyParams.hpp"

#include <opm/material/common/Spline.hpp>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 * \brief Implementation of the regularized Brooks-Corey capillary
 *        pressure / relative permeability <-> saturation relation.
 *
 * This class bundles the "raw" curves as static members and doesn't
 * concern itself converting absolute to effective saturations and
 * vice versa.
 *
 * In order to avoid very steep gradients the marginal values are
 * "regularized".  This means that in stead of following the curve of
 * the material law in these regions, some linear approximation is
 * used.  Doing this is not worse than following the material
 * law. E.g. for very low wetting phase values the material laws
 * predict infinite values for \f$p_c\f$ which is completely
 * unphysical. In case of very high wetting phase saturations the
 * difference between regularized and "pure" material law is not big.
 *
 * Regularizing has the additional benefit of being numerically
 * friendly: Newton's method does not like infinite gradients.
 *
 * The implementation is accomplished as follows:
 * - check whether we are in the range of regularization
 *   - yes: use the regularization
 *   - no: forward to the standard material law.
 *
 * \see BrooksCorey
 */
template <class TraitsT, class ParamsT = RegularizedBrooksCoreyParams<TraitsT> >
class RegularizedBrooksCorey : public TraitsT
{
    typedef ::Opm::BrooksCorey<TraitsT, ParamsT> BrooksCorey;

public:
    typedef TraitsT Traits;
    typedef ParamsT Params;
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The regularized Brooks-Corey capillary pressure law only "
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

    static_assert(Traits::numPhases == 2,
                  "The number of fluid phases must be two if you want to use "
                  "this material law!");

    /*!
     * \brief The capillary pressure-saturation curves depending on absolute saturations.
     *
     * \param values A random access container which stores the
     *               relative pressure of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the material law.
     * \param fs The fluid state for which the capillary pressure
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container& values, const Params& params, const FluidState& fs)
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
    static void saturations(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = Sw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = 1 - values[Traits::wettingPhaseIdx];
    }

    /*!
     * \brief The relative permeability-saturation curves depending on absolute saturations.
     *
     * \param values A random access container which stores the
     *               relative permeability of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the material law.
     * \param fs The fluid state for which the relative permeabilities
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = krw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = krn<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief A regularized Brooks-Corey capillary pressure-saturation curve.
     *
     * This is a regularized variant of the Brooks-Corey curve:
     *
     * - For wetting phase saturations lower than the threshold saturation, the
     *   \f$p_c(S_w)\f$ curve is extrapolated using a straight line exhibiting the slope
     *   unregularized capillary pressure curve at the threshold saturation.
     * - For wetting phase saturations larger than 1, the curve is extrapolated using a
     *   straight line that exhibits the slope of the unregularized curve at \f$S_w =
     *   1\f$
     *
     * \sa BrooksCorey::pcnw
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params& params, const FluidState& fs)
    {
        const auto& Sw = decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));
        return twoPhaseSatPcnw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params& params, const Evaluation& Sw)
    {
        const Scalar Sthres = params.pcnwLowSw();

        if (Sw <= Sthres) {
            Scalar m = params.pcnwSlopeLow();
            Scalar pcnw_SwLow = params.pcnwLow();
            return pcnw_SwLow + m*(Sw - Sthres);
        }
        else if (Sw >= 1.0) {
            Scalar m = params.pcnwSlopeHigh();
            Scalar pcnw_SwHigh = params.pcnwHigh();
            return pcnw_SwHigh + m*(Sw - 1.0);
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real Brooks-Corey law...
        return BrooksCorey::twoPhaseSatPcnw(params, Sw);
    }

    /*!
     * \brief A regularized Brooks-Corey saturation-capillary pressure
     *        curve.
     *
     * This is the inverse of the pcnw() method.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params& params, const FluidState& fs)
    {
        const Evaluation& pC =
            decay<Evaluation>(fs.pressure(Traits::nonWettingPhaseIdx))
            - decay<Evaluation>(fs.pressure(Traits::wettingPhaseIdx));
        return twoPhaseSatSw(params, pC);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params& params, const Evaluation& pcnw)
    {
        const Scalar Sthres = params.pcnwLowSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized version of the
        // Brooks-Corey law. If the input capillary pressure is
        // smaller than the entry pressure, make sure that we will
        // regularize.
        Evaluation Sw = 1.5;
        if (pcnw >= params.entryPressure())
            Sw = BrooksCorey::twoPhaseSatSw(params, pcnw);

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Sw <= Sthres) {
            // invert the low saturation regularization of pcnw()
            Scalar m = params.pcnwSlopeLow();
            Scalar pcnw_SwLow = params.pcnwLow();
            return Sthres + (pcnw - pcnw_SwLow)/m;
        }
        else if (Sw > 1.0) {
            Scalar m = params.pcnwSlopeHigh();
            Scalar pcnw_SwHigh = params.pcnwHigh();
            return 1.0 + (pcnw - pcnw_SwHigh)/m;;
        }

        return BrooksCorey::twoPhaseSatSw(params, pcnw);
    }

    /*!
     * \brief Calculate the non-wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params& params, const FluidState& fs)
    { return 1 - Sw<FluidState, Evaluation>(params, fs); }

    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params& params, const Evaluation& pcnw)
    { return 1 - twoPhaseSatSw(params, pcnw); }

    /*!
     * \brief Regularized version of the relative permeability of the
     *        wetting phase of the Brooks-Corey curves.
     *
     * The approach for regularization is very similar to the one of
     * the capillary pressure, but it does not avoid kinks:
     * - For wetting phase saturations between 0 and 1, use the
     *   unregularized Brooks-Corey wetting phase relative
     *   permeability
     * - For wetting phase saturations smaller than 0, return 0
     * - For wetting phase saturations larger than 1, return 1
     *
     * \sa BrooksCorey::krw
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params& params, const FluidState& fs)
    {
        const auto& Sw = decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));
        return twoPhaseSatKrw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params& params, const Evaluation& Sw)
    {
        if (Sw <= 0.0)
            return 0.0;
        else if (Sw >= 1.0)
            return 1.0;

        return BrooksCorey::twoPhaseSatKrw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrwInv(const Params& params, const Evaluation& krw)
    {
        if (krw <= 0.0)
            return 0.0;
        else if (krw >= 1.0)
            return 1.0;

        return BrooksCorey::twoPhaseSatKrwInv(params, krw);
    }

    /*!
     * \brief Regularized version of the relative permeability of the
     *        non-wetting phase of the Brooks-Corey curves.
     *
     * The approach for regularization is very similar to the one of
     * the capillary pressure, but it does not avoid kinks:
     * - For wetting phase saturations between 0 and 1, use the
     *   unregularized Brooks-Corey non-wetting phase relative
     *   permeability
     * - For wetting phase saturations smaller than 0, return 1
     * - For wetting phase saturations larger than 1, return 0
     *
     * \sa BrooksCorey::krn
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params& params, const FluidState& fs)
    {
        const Evaluation& Sw =
            1.0 - decay<Evaluation>(fs.saturation(Traits::nonWettingPhaseIdx));
        return twoPhaseSatKrn(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params& params, const Evaluation& Sw)
    {
        if (Sw >= 1.0)
            return 0.0;
        else if (Sw <= 0.0)
            return 1.0;

        return BrooksCorey::twoPhaseSatKrn(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrnInv(const Params& params, const Evaluation& krn)
    {
        if (krn <= 0.0)
            return 1.0;
        else if (krn >= 1.0)
            return 0.0;

        return BrooksCorey::twoPhaseSatKrnInv(params, krn);
    }
};
} // namespace Opm

#endif
