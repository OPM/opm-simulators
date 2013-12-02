// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2008-2013 by Andreas Lauser
  Copyright (C) 2010 by Philipp Nuske
  Copyright (C) 2010 by Bernd Flemisch

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
 * \copydoc Opm::RegularizedVanGenuchten
 */
#ifndef OPM_REGULARIZED_VAN_GENUCHTEN_HPP
#define OPM_REGULARIZED_VAN_GENUCHTEN_HPP

#include "VanGenuchten.hpp"
#include "RegularizedVanGenuchtenParams.hpp"

#include <opm/core/utility/Spline.hpp>

#include <algorithm>

namespace Opm {

/*!
 * \ingroup FluidMatrixInteractions
 * \brief Implementation of the regularized  van Genuchten's
       capillary pressure / relative permeability  <-> saturation relation.
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
 *  - yes: use the regularization
 *  - no: forward to the standard material law.
 *
 * An example of the regularization of the capillary pressure curve is
 * shown below: \image html regularizedVanGenuchten.png
 *
 * \see VanGenuchten
 */
template <class TraitsT, class ParamsT = RegularizedVanGenuchtenParams<TraitsT> >
class RegularizedVanGenuchten : public TraitsT
{
    typedef Opm::VanGenuchten<TraitsT, ParamsT> VanGenuchten;

public:
    typedef TraitsT Traits;
    typedef ParamsT Params;

    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The regularized van Genuchten capillary pressure law only "
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
     * \brief Calculate the pressure difference of the phases in the
     *        most generic way.
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container &values, const Params &params, const FluidState &fs)
    {
        values[Traits::wPhaseIdx] = 0.0; // reference phase
        values[Traits::nPhaseIdx] = pcnw(params, fs);
    }

    /*!
     * \brief Calculate the saturations of the phases starting from
     *        their pressure differences.
     */
    template <class Container, class FluidState>
    static void saturations(Container &values, const Params &params, const FluidState &fs)
    {
        values[Traits::wPhaseIdx] = Sw(params, fs);
        values[Traits::nPhaseIdx] = 1 - values[Traits::wPhaseIdx];
    }

    /*!
     * \brief Returns the relative permeabilities of the phases
     *        dependening on the phase saturations.
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container &values, const Params &params, const FluidState &fs)
    {
        values[Traits::wPhaseIdx] = krw(params, fs);
        values[Traits::nPhaseIdx] = krn(params, fs);
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
        values[Traits::wPhaseIdx] = 0;
        values[Traits::nPhaseIdx] = 0;
        if (satPhaseIdx == Traits::wPhaseIdx)
            values[Traits::nPhaseIdx] = twoPhaseSatDPcnw_dSw(params, state.saturation(Traits::wPhaseIdx));
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
        if (satPhaseIdx == Traits::wPhaseIdx) {
            values[Traits::wPhaseIdx] = twoPhaseSatDKrw_dSw(params, state.saturation(Traits::wPhaseIdx));
            values[Traits::nPhaseIdx] = 0;
        }
        else {
            values[Traits::wPhaseIdx] = 0;
            values[Traits::nPhaseIdx] = - twoPhaseSatDKrn_dSw(params, 1 - state.saturation(Traits::nPhaseIdx));
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
     * \brief A regularized van Genuchten capillary pressure-saturation
     *          curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$p_c(S_w)\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$ \overline S_w =1\f$ by a straight line (yes, there is a kink :-( ).
     *
     *  For not-regularized part:
     *
     * \copydetails VanGenuchten::pC()
     */
    template <class FluidState>
    static Scalar pcnw(const Params &params, const FluidState &fs)
    { return twoPhaseSatPcnw(params, fs.saturation(Traits::wPhaseIdx)); }

    static Scalar twoPhaseSatPcnw(const Params &params, Scalar Sw)
    {
        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar SwThLow = params.pcnwLowSw();
        const Scalar SwThHigh = params.pcnwHighSw();

        // make sure that the capillary pressure observes a derivative
        // != 0 for 'illegal' saturations. This is favourable for the
        // newton solver (if the derivative is calculated numerically)
        // in order to get the saturation moving to the right
        // direction if it temporarily is in an 'illegal' range.
        if (Sw < SwThLow) {
            return params.pcnwLow() + params.pcnwSlopeLow()*(Sw - SwThLow);
        }
        else if (Sw > SwThHigh)
        {
            Scalar yTh = params.pcnwHigh();
            Scalar m1 = (0.0 - yTh)/(1.0 - SwThHigh)*2;

            if (Sw < 1.0) {
                // use spline between threshold Sw and 1.0
                const Spline<Scalar> &sp = params.pcnwHighSpline();

                return sp.eval(Sw);
            }
            else {
                // straight line for Sw > 1.0
                return m1*(Sw - 1.0) + 0.0;
            }
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real van genuchten law...
        return VanGenuchten::twoPhaseSatPcnw(params, Sw);
    }

    /*!
     * \brief   A regularized van Genuchten saturation-capillary pressure curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$p_c(S_w)\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$ \overline S_w =1\f$ by a straight line (yes, there is a kink :-( ).
     *
     *  The according quantities are obtained by exploiting theorem of intersecting lines.
     *
     *  For not-regularized part:
     *
         \copydetails VanGenuchten::Sw()
     *
     */
    template <class FluidState>
    static Scalar Sw(const Params &params, const FluidState &fs)
    {
        Scalar pC = fs.pressure(Traits::nPhaseIdx) - fs.pressure(Traits::wPhaseIdx);
        return twoPhaseSatSw(params, pC);
    }

    static Scalar twoPhaseSatSw(const Params &params, Scalar pC)
    {
        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar SwThLow = params.pcnwLowSw();
        const Scalar SwThHigh = params.pcnwHighSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized verision of van
        // Genuchten's law
        Scalar Sw;
        if (pC <= 0) {
            // invert straight line for Sw > 1.0
            Scalar m1 = params.pcnwSlopeHigh();
            return pC/m1 + 1.0;
        }
        else
            Sw = VanGenuchten::twoPhaseSatSw(params, pC);

        // invert the regularization if necessary
        if (Sw <= SwThLow) {
            // invert the low saturation regularization of pC()
            Scalar pC_SwLow = VanGenuchten::twoPhaseSatPcnw(params, SwThLow);
            return (pC - pC_SwLow)/params.pcnwSlopeLow() + SwThLow;
        }
        else if (SwThHigh < Sw /* && Sw < 1.0*/)
        {
            // invert spline between threshold saturation and 1.0
            const Spline<Scalar>& spline = params.pcnwHighSpline();

            return spline.intersectInterval(/*x0=*/SwThHigh, /*x1=*/1.0,
                                            /*a=*/0, /*b=*/0, /*c=*/0, /*d=*/pC);
        }

        return Sw;
    }

    /*!
     * \brief Calculate the non-wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState>
    static Scalar Sn(const Params &params, const FluidState &fs)
    { return 1 - Sw(params, fs); }

    static Scalar twoPhaseSatSn(const Params &params, Scalar pc)
    { return 1 - twoPhaseSatSw(params, pc); }

    /*!
     * \brief A regularized version of the partial derivative
     *        of the \f$p_c(\overline S_w)\f$ w.r.t. effective saturation
     *        according to van Genuchten.
     *
     * regularized part:
     *    - low saturation:  use the slope of the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$ \overline S_w =1\f$ by a straight line and use that slope (yes, there is a kink :-( ).
     *
     *        For not-regularized part:
     *
     * \copydetails VanGenuchten::dpC_dSw()
     *
     */
    template <class FluidState>
    static Scalar dPcnw_dSw(const Params &params, const FluidState &fs)
    { return twoPhaseSatDPcnw_dSw(params, fs.saturation(Traits::wPhaseIdx)); }

    static Scalar twoPhaseSatDPcnw_dSw(const Params &params, Scalar Sw)
    {
        // derivative of the regualarization
        if (Sw < params.pcnwLowSw()) {
            // the slope of the straight line used in pC()
            return params.pcnwSlopeLow();
        }
        else if (params.pcnwHighSw() <= Sw) {
            if (Sw < 1)
                return params.pcnwHighSpline().evalDerivative(Sw);
            else
                // the slope of the straight line used for the right
                // side of the capillary pressure function
                return params.pcnwSlopeHigh();
        }

        return VanGenuchten::twoPhaseSatDPcnw_dSw(params, Sw);
    }

    /*!
     * \brief A regularized version of the partial derivative
     *        of the \f$\overline S_w(p_c)\f$ w.r.t. cap.pressure
     *        according to van Genuchten.
     *
     *  regularized part:
     *    - low saturation:  use the slope of the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$ \overline S_w =1\f$ by a straight line and use that slope (yes, there is a kink :-( ).
     *
     *        For not-regularized part:
     * \copydetails VanGenuchten::dSw_dpC()
     */
    template <class FluidState>
    static Scalar dSw_dpC(const Params &params, const FluidState &fs)
    {
        Scalar pC = fs.pressure(Traits::nPhaseIdx) - fs.pressure(Traits::wPhaseIdx);

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized verision of van
        // Genuchten's law
        Scalar Sw;
        if (pC < 0)
            Sw = 1.5; // make sure we regularize below
        else
            Sw = VanGenuchten::Sw_raw(params, pC);

        // derivative of the regularization
        if (Sw < params.pcnwLowSw()) {
            // same as in dpC_dSw() but inverted
            return 1/params.pcnwSlopeLow();
        }
        if (Sw > params.pcnwHighSw()) {
            if (Sw < 1)
                return 1/params.pcnwHighSpline().evalDerivative(Sw);

            // same as in dpC_dSw() but inverted
            return 1/params.pcnwSlopHigh();
        }

        return VanGenuchten::dSw_dpnw(params, fs);
    }

    /*!
     * \brief   Regularized version of the  relative permeability
     *          for the wetting phase of
     *          the medium implied by the van Genuchten
     *          parameterization.
     *
     *  regularized part:
     *    - below \f$ \overline S_w =0\f$:                  set relative permeability to zero
     *    - above \f$ \overline S_w =1\f$:                  set relative permeability to one
     *    - between \f$ 0.95 \leq \overline S_w \leq 1\f$:  use a spline as interpolation
     *
     *  For not-regularized part:
        \copydetails VanGenuchten::krw()
     */
    template <class FluidState>
    static Scalar krw(const Params &params, const FluidState &fs)
    { return twoPhaseSatKrw(params, fs.saturation(Traits::wPhaseIdx)); }

    static Scalar twoPhaseSatKrw(const Params &params, Scalar Sw)
    {
        // regularize
        if (Sw < 0)
            return 0;
        else if (Sw > 1)
            return 1;

        return VanGenuchten::twoPhaseSatKrw(params, Sw);
    }

    static Scalar twoPhaseSatDKrw_dSw(const Params &params, Scalar Sw)
    {
        if (Sw <= 0.0 || Sw >= 1.0)
            return 0.0;

        return VanGenuchten::twoPhaseSatDKrw_dSw(params, Sw);
    }

    /*!
     * \brief   Regularized version of the  relative permeability
     *          for the non-wetting phase of
     *          the medium implied by the van Genuchten
     *          parameterization.
     *
     * regularized part:
     *    - below \f$ \overline S_w =0\f$:                  set relative permeability to zero
     *    - above \f$ \overline S_w =1\f$:                  set relative permeability to one
     *    - for \f$ 0 \leq \overline S_w \leq 0.05 \f$:     use a spline as interpolation
     *
         \copydetails VanGenuchten::krn()
     *
     */
    template <class FluidState>
    static Scalar krn(const Params &params, const FluidState &fs)
    { return twoPhaseSatKrn(params, 1.0 - fs.saturation(Traits::nPhaseIdx)); }

    static Scalar twoPhaseSatKrn(const Params &params, Scalar Sw)
    {
        // regularize
        if (Sw <= 0)
            return 1;
        else if (Sw >= 1)
            return 0;

        return VanGenuchten::twoPhaseSatKrn(params, Sw);
    }

    static Scalar twoPhaseSatDKrn_dSw(const Params &params, Scalar Sw)
    {
        if (Sw <= 0.0 || Sw >= 1.0)
            return 0.0;

        return VanGenuchten::twoPhaseSatDKrn_dSw(params, Sw);
    }
};

} // namespace Opm

#endif
