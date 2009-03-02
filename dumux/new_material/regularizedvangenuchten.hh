/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Implementation of a regularized version of van Genuchten's
 *        capillary pressure <-> saturation relation.
 */
#ifndef REGULARIZED_VAN_GENUCHTEN_HH
#define REGULARIZED_VAN_GENUCHTEN_HH

#include <dumux/auxiliary/apis.hh>

#include <dumux/new_material/vangenuchten.hh>
#include <dumux/new_material/regularizedvangenuchtenstate.hh>

#include <algorithm>

#include <math.h>
#include <assert.h>

namespace Dune
{
/*!\ingroup material
 *
 * \brief Implementation of van Genuchten's capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vince versa.
 *
 * \sa VanGenuchten, VanGenuchtenTwophase
 */
template <class StateT>
class RegularizedVanGenuchten
{
public:
    typedef StateT State;
    typedef typename State::Scalar Scalar;

    typedef Dune::VanGenuchten<StateT> VanGenuchten;

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     * Van Genuchten's empirical capillary pressure <-> saturation
     * function is given by
     * \f[
     p_C = (\overline{S}_w^{-1/m} - 1)^{1/n}/\alpha
     \f]
     * \param Swe Effective saturation of of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar pC(const State &state, Scalar Swe)
    {
        Api::require<Api::RegularizedVanGenuchtenParams>(state);

        // make sure that the capilarry pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Swe <= SweLow_) {
            Scalar pC_SweLow = VanGenuchten::pC(state, SweLow_);
            Scalar m = (state.vgMaxPC() - pC_SweLow)
                / (state.vgMinSw() - SweLow_);
            return pC_SweLow + m*(Swe - SweLow_);
        }
        else if (Swe >= SweHigh_) {
            Scalar pC_SweHigh = VanGenuchten::pC(state, SweHigh_);
            Scalar m = (pC_SweHigh - 0)/(SweHigh_ - 1.0);
            return pC_SweHigh + m*(Swe - SweHigh_);
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real van genuchten law...
        return VanGenuchten::pC(state, Swe);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     \overline{S}_w = {p_C}^{-1} = ((\alpha p_C)^n + 1)^{-m}
     \f]
     *
     * \param pC Capillary pressure \f$\p_C\f$
     * \return The effective saturaion of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar Sw(const State &state, Scalar pC)
    {
        Api::require<Api::RegularizedVanGenuchtenParams>(state);

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized verision of van
        // Genuchten's law
        Scalar Swe;
        if (pC < 0)
            Swe = 1.5; // make sure we regularize below
        else
            Swe = VanGenuchten::Sw(state, pC);

        // make sure that the capilarry pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Swe <= SweLow_) {
            // invert the low saturation regularization of pC()
            Scalar pC_SweLow = VanGenuchten::pC(state, SweLow_);
            Scalar m = (state.vgMaxPC() - pC_SweLow)
                / (state.vgMinSw() - SweLow_);
            return SweLow_ + (pC - pC_SweLow)/m;
        }
        else if (Swe >= SweHigh_) {
            // invert the high saturation regularization of pC()
            Scalar pC_SweHigh = VanGenuchten::pC(state, SweHigh_);
            Scalar m = (pC_SweHigh - 0)/(SweHigh_ - 1.0);
            return SweHigh_ + (pC - pC_SweHigh)/m;
        }

        return VanGenuchten::Sw(state, pC);
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{1}{\alpha} (\overline{S}_w^{-1/m} - 1)^{1/n - }
     \overline{S}_w^{-1/m} / \overline{S}_w / m
     \f]
    */
    static Scalar dpC_dSw(const State &state, Scalar Swe)
    {
        Api::require<Api::RegularizedVanGenuchtenParams>(state);

        // derivative of the regualarization
        if (Swe <= SweLow_) {
            // calculate the slope of the straight line used in pC()
            Scalar pC_SweLow = VanGenuchten::pC(state, SweLow_);
            Scalar m = (state.vgMaxPC() - pC_SweLow)
                / (state.vgMinSw() - SweLow_);
            return m;
        }
        else if (Swe >= SweHigh_) {
            // calculate the slope of the straight line used in pC()
            Scalar pC_SweHigh = VanGenuchten::pC(state, SweHigh_);
            Scalar m = (pC_SweHigh - 0)/(SweHigh_ - 1.0);
            return m;
        }

        return VanGenuchten::dpC_dSw(state, Swe);
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const State &state, Scalar pC)
    {
        Api::require<Api::RegularizedVanGenuchtenParams>(state);

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized verision of van
        // Genuchten's law
        Scalar Swe;
        if (pC < 0)
            Swe = 1.5; // make sure we regularize below
        else
            Swe = VanGenuchten::Sw(state, pC);

        // derivative of the regularization
        if (Swe <= SweLow_) {
            // same as in dpC_dSw() but inverted
            Scalar pC_SweLow = VanGenuchten::pC(state, SweLow_);
            Scalar m = (state.vgMaxPC() - pC_SweLow)
                / (state.vgMinSw() - SweLow_);
            return 1/m;
        }
        else if (Swe >= SweHigh_) {
            // same as in dpC_dSw() but inverted
            Scalar pC_SweHigh = VanGenuchten::pC(state, SweHigh_);
            Scalar m = (pC_SweHigh - 0)/(SweHigh_ - 1.0);
            return 1/m;
        }

        return VanGenuchten::dpC_dSw(state, pC);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     */
    static Scalar krw(const State &state, Scalar Sw_mob)
    {
        Api::require<Api::RegularizedVanGenuchtenParams>(state);
        if (Sw_mob < 0)
            return 0;
        else if (Sw_mob > 1)
            return 1;

        return VanGenuchten::krw(state, Sw_mob);
    };

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     */
    static Scalar krn(const State &state, Scalar Sw_mob)
    {
        Api::require<Api::RegularizedVanGenuchtenParams>(state);

        if (Sw_mob <= 0)
            return 1;
        else if (Sw_mob >= 1)
            return 0;

        return VanGenuchten::krn(state, Sw_mob);
    }

    //! Effective saturation below which we regularize
    static const Scalar SweLow_;
    //! Effective saturation above which we regularize
    static const Scalar SweHigh_;
};

template <class StateT>
const typename StateT::Scalar RegularizedVanGenuchten<StateT>::SweLow_(0.03);
template <class StateT>
const typename StateT::Scalar RegularizedVanGenuchten<StateT>::SweHigh_(0.97);
}

#endif
