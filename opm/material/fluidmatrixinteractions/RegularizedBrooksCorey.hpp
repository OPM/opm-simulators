/*
  Copyright (C) 2008-2013 by Andreas Lauser
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
 * \copydoc Opm::RegularizedBrooksCorey
 */
#ifndef REGULARIZED_BROOKS_COREY_HPP
#define REGULARIZED_BROOKS_COREY_HPP

#include "BrooksCorey.hpp"
#include "RegularizedBrooksCoreyParams.hpp"

#include <opm/core/utility/Spline.hpp>

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
    typedef Opm::BrooksCorey<TraitsT, ParamsT> BrooksCorey;

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
    static void capillaryPressures(Container &values, const Params &params, const FluidState &fs)
    {
        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = pcnw(params, fs);
    }

    /*!
     * \brief Calculate the saturations of the phases starting from
     *        their pressure differences.
     */
    template <class Container, class FluidState>
    static void saturations(Container &values, const Params &params, const FluidState &fs)
    {
        values[Traits::wettingPhaseIdx] = Sw(params, fs);
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
    static void relativePermeabilities(Container &values, const Params &params, const FluidState &fs)
    {
        values[Traits::wettingPhaseIdx] = krw(params, fs);
        values[Traits::nonWettingPhaseIdx] = krn(params, fs);
    }

    /*!
     * \brief The derivative of all capillary pressures in regard to
     *        a given phase saturation.
     */
    template <class ContainerT, class FluidState>
    static void dCapillaryPressures_dSaturation(ContainerT &values,
                                                const Params &params,
                                                const FluidState &state,
                                                int satPhaseIdx)
    {
        values[Traits::wettingPhaseIdx] = 0;
        values[Traits::nonWettingPhaseIdx] = 0;
        if (satPhaseIdx == Traits::wettingPhaseIdx)
            values[Traits::nonWettingPhaseIdx] = twoPhaseSatDPcnw_dSw(params, state.saturation(Traits::wettingPhaseIdx));
    }

    /*!
     * \brief The derivative of all capillary pressures in regard to
     *        a given phase pressure.
     */
    template <class ContainerT, class FluidState>
    static void dCapillaryPressures_dPressure(ContainerT &values,
                                              const Params &params,
                                              const FluidState &state,
                                              int pPhaseIdx)
    {
        // -> not pressure dependent
        for (int pcPhaseIdx = 0; pcPhaseIdx < numPhases; ++pcPhaseIdx)
            values[pcPhaseIdx] = 0.0;
    }

    /*!
     * \brief The derivative of all capillary pressures in regard to
     *        temperature.
     */
    template <class ContainerT, class FluidState>
    static void dCapillaryPressures_dTemperature(ContainerT &values,
                                                 const Params &params,
                                                 const FluidState &state)
    {
        // -> not temperature dependent
        for (int pcPhaseIdx = 0; pcPhaseIdx < numPhases; ++pcPhaseIdx)
            values[pcPhaseIdx] = 0.0;
    }

    /*!
     * \brief The derivative of all capillary pressures in regard to
     *        a given mole fraction of a component in a phase.
     */
    template <class ContainerT, class FluidState>
    static void dCapillaryPressures_dMoleFraction(ContainerT &values,
                                                  const Params &params,
                                                  const FluidState &state,
                                                  int phaseIdx,
                                                  int compIdx)
    {
        // -> not composition dependent
        for (int pcPhaseIdx = 0; pcPhaseIdx < numPhases; ++pcPhaseIdx)
            values[pcPhaseIdx] = 0.0;
    }

    /*!
     * \brief The derivative of all relative permeabilities in regard to
     *        a given phase saturation.
     */
    template <class ContainerT, class FluidState>
    static void dRelativePermeabilities_dSaturation(ContainerT &values,
                                                    const Params &params,
                                                    const FluidState &state,
                                                    int satPhaseIdx)
    {
        if (satPhaseIdx == Traits::wettingPhaseIdx) {
            values[Traits::wettingPhaseIdx] = twoPhaseSatDKrw_dSw(params, state.saturation(Traits::wettingPhaseIdx));
            values[Traits::nonWettingPhaseIdx] = 0;
        }
        else {
            values[Traits::wettingPhaseIdx] = 0;
            values[Traits::nonWettingPhaseIdx] = - twoPhaseSatDKrn_dSw(params, 1 - state.saturation(Traits::nonWettingPhaseIdx));
        }
    }

    /*!
     * \brief The derivative of all relative permeabilities in regard to
     *        a given phase pressure.
     */
    template <class ContainerT, class FluidState>
    static void dRelativePermeabilities_dPressure(ContainerT &values,
                                                  const Params &params,
                                                  const FluidState &state,
                                                  int pPhaseIdx)
    {
        // -> not pressure dependent
        for (int krPhaseIdx = 0; krPhaseIdx < numPhases; ++krPhaseIdx)
            values[krPhaseIdx] = 0.0;
    }

    /*!
     * \brief The derivative of all relative permeabilities in regard to
     *        temperature.
     */
    template <class ContainerT, class FluidState>
    static void dRelativePermeabilities_dTemperature(ContainerT &values,
                                                     const Params &params,
                                                     const FluidState &state)
    {
        // -> not temperature dependent
        for (int krPhaseIdx = 0; krPhaseIdx < numPhases; ++krPhaseIdx)
            values[krPhaseIdx] = 0.0;
    }

    /*!
     * \brief The derivative of all relative permeabilities in regard to
     *        a given mole fraction of a component in a phase.
     */
    template <class ContainerT, class FluidState>
    static void dRelativePermeabilities_dMoleFraction(ContainerT &values,
                                                      const Params &params,
                                                      const FluidState &state,
                                                      int phaseIdx,
                                                      int compIdx)
    {
        // -> not composition dependent
        for (int krPhaseIdx = 0; krPhaseIdx < numPhases; ++krPhaseIdx)
            values[krPhaseIdx] = 0.0;
    }

    /*!
     * \brief A regularized Brooks-Corey capillary pressure-saturation
     *        curve.
     *
     * This is a regularized variant of the Brooks-Corey curve. For
     * wetting phase saturations between a lower threshold saturation
     * and \f$S_w=1\f$, for other wetting phase saturations it is
     * regularized in a way which removes the singularity at
     * \f$S_w=0\f$, avoids kinks and allows the capillary pressure to
     * reach arbitrary values. (Albeit, to reach a given capillary
     * pressure, the saturations can become unphysical). The
     * regularization is done in the following way:
     *
     * - For wetting phase saturations lower than the threshold
     *   saturation, the \f$p_c(S_w)\f$ curve is extrapolated using a
     *   straight line exhibiting the slope unregularized capillary
     *   pressure curve at the threshold saturation.
     * - For wetting phase saturations larger than 1, the Brooks-Corey
     *   curve is extrapolated using a straight line that exhibits the
     *   slope of the unregularized Brooks-Corey curve at \f$S_w =
     *   1\f$
     *
     * \sa BrooksCorey::pcnw
     */
    template <class FluidState>
    static Scalar pcnw(const Params &params, const FluidState &fs)
    { return twoPhaseSatPcnw(params, fs.saturation(Traits::wettingPhaseIdx)); }

    static Scalar twoPhaseSatPcnw(const Params &params, Scalar Sw)
    {
        const Scalar Sthres = params.thresholdSw();

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative is calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Sw <= Sthres) {
            Scalar m = BrooksCorey::twoPhaseSatDPcnw_dSw(params, Sthres);
            Scalar pcnw_SwLow = BrooksCorey::twoPhaseSatPcnw(params, Sthres);
            return pcnw_SwLow + m*(Sw - Sthres);
        }
        else if (Sw > 1.0) {
            Scalar m = BrooksCorey::twoPhaseSatDPcnw_dSw(params, 1.0);
            Scalar pcnw_SwHigh = params.entryPressure();
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
    template <class FluidState>
    static Scalar Sw(const Params &params, const FluidState &fs)
    {
        Scalar pcnw = fs.pressure(Traits::nonWettingPhaseIdx) - fs.pressure(Traits::wettingPhaseIdx);
        return twoPhaseSatSw(params, pcnw);
    }

    static Scalar twoPhaseSatSw(const Params &params, Scalar pcnw)
    {
        const Scalar Sthres = params.thresholdSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized version of the
        // Brooks-Corey law. If the input capillary pressure is
        // smaller than the entry pressure, make sure that we will
        // regularize.
        Scalar Sw = 1.5;
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
            Scalar m = BrooksCorey::twoPhaseSatDPcnw_dSw(params, Sthres);
            Scalar pcnw_SwLow = BrooksCorey::twoPhaseSatPcnw(params, Sthres);
            return Sthres + (pcnw - pcnw_SwLow)/m;
        }
        else if (Sw > 1.0) {
            Scalar m = BrooksCorey::twoPhaseSatDPcnw_dSw(params, 1.0);
            Scalar pcnw_SwHigh = BrooksCorey::twoPhaseSatPcnw(params, 1.0);
            return 1.0 + (pcnw - pcnw_SwHigh)/m;;
        }

        return BrooksCorey::twoPhaseSatSw(params, pcnw);
    }

    /*!
     * \brief Calculate the non-wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState>
    static Scalar Sn(const Params &params, const FluidState &fs)
    { return 1 - Sw(params, fs); }

    static Scalar twoPhaseSatSn(const Params &params, Scalar pcnw)
    { return 1 - twoPhaseSatSw(params, pcnw); }

    /*!
     * \brief The derivative of the regularized Brooks-Corey capillary
     *        pressure-saturation curve.
     */
    static Scalar twoPhaseSatDPcnw_dSw(const Params &params, Scalar Sw)
    {
        const Scalar Sthres = params.thresholdSw();

        // derivative of the regualarization
        if (Sw <= Sthres) {
            // calculate the slope of the straight line used in pcnw()
            Scalar m = BrooksCorey::twoPhaseSatDPcnw_dSw(params, Sthres);
            return m;
        }
        else if (Sw > 1.0) {
            // calculate the slope of the straight line used in pcnw()
            Scalar m = BrooksCorey::twoPhaseSatDPcnw_dSw(params, 1.0);
            return m;
        }

        return BrooksCorey::twoPhaseSatDPcnw_dSw(params, Sw);
    }

    /*!
     * \brief The derivative of the regularized Brooks-Corey
     *        saturation-capillary pressure curve.
     */
    static Scalar twoPhaseSatDSw_dpcnw(const Params &params, Scalar pcnw)
    {
        const Scalar Sthres = params.thresholdSw();

        // calculate the saturation which corresponds to the
        // saturation in the non-regularized version of the
        // Brooks-Corey law
        Scalar Sw;
        if (pcnw < params.entryPressure())
            Sw = 1.5; // make sure we regularize (see below)
        else
            Sw = BrooksCorey::Sw(params, pcnw);

        // derivative of the regularization
        if (Sw <= Sthres) {
            // calculate the slope of the straight line used in pcnw()
            Scalar m = BrooksCorey::dPcnw_dSw(params, Sthres);
            return 1/m;
        }
        else if (Sw > 1.0) {
            // calculate the slope of the straight line used in pcnw()
            Scalar m = BrooksCorey::dPcnw_dSw(params, 1.0);
            return 1/m;
        }
        return 1.0/BrooksCorey::dPcnw_dSw(params, Sw);
    }

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
    template <class FluidState>
    static Scalar krw(const Params &params, const FluidState &fs)
    { return twoPhaseSatKrw(params, fs.saturation(Traits::wettingPhaseIdx)); }

    static Scalar twoPhaseSatKrw(const Params &params, Scalar Sw)
    {
        if (Sw <= 0.0)
            return 0.0;
        else if (Sw >= 1.0)
            return 1.0;

        return BrooksCorey::twoPhaseSatKrw(params, Sw);
    }

    static Scalar twoPhaseSatDKrw_dSw(const Params &params, Scalar Sw)
    {
        if (Sw <= 0.0 || Sw >= 1.0)
            return 0.0;

        return BrooksCorey::twoPhaseSatDKrw_dSw(params, Sw);
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
    template <class FluidState>
    static Scalar krn(const Params &params, const FluidState &fs)
    { return twoPhaseSatKrn(params, 1.0 - fs.saturation(Traits::nonWettingPhaseIdx)); }

    static Scalar twoPhaseSatKrn(const Params &params, Scalar Sw)
    {
        if (Sw >= 1.0)
            return 0.0;
        else if (Sw <= 0.0)
            return 1.0;

        return BrooksCorey::twoPhaseSatKrn(params, Sw);
    }

    static Scalar twoPhaseSatDKrn_dSw(const Params &params, Scalar Sw)
    {
        if (Sw <= 0.0 || Sw >= 1.0)
            return 0.0;

        return BrooksCorey::twoPhaseSatDKrn_dSw(params, Sw);
    }

};
} // namespace Opm

#endif
