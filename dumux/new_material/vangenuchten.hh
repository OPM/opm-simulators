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
 * \file VanGenuchten.hh Implementation of van Genuchten's capillary
 *                       pressure <-> saturation relation
 */
#ifndef VAN_GENUCHTEN_HH
#define VAN_GENUCHTEN_HH

#include <dumux/auxiliary/apis.hh>
#include <dumux/new_material/vangenuchtenstate.hh>

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
    class VanGenuchten
    {
    public:
        typedef StateT State;
        typedef typename State::Scalar Scalar;

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
            Api::require<Api::VanGenuchtenParams>(state);
            assert(0 <= Swe && Swe <= 1);
            return pow(pow(Swe, -1.0/state.vgM()) - 1, 1.0/state.vgN())/state.vgAlpha();
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
            Api::require<Api::VanGenuchtenParams>(state);
            assert(pC >= 0);

            return pow(pow(state.vgAlpha()*pC, state.vgN()) + 1, -state.vgM());
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
            Api::require<Api::VanGenuchtenParams>(state);
            assert(0 <= Swe && Swe <= 1);

            Scalar powSwe = pow(Swe, -1/state.vgM());
            return - 1/state.vgAlpha() * pow(powSwe - 1, 1/state.vgN() - 1)/state.vgN()
                                    * powSwe/Swe/state.vgM();
        }

        /*!
         * \brief Returns the partial derivative of the effective
         *        saturation to the capillary pressure.
         */
        static Scalar dSw_dpC(const State &state, Scalar pC)
        {
            Api::require<Api::VanGenuchtenParams>(state);
            assert(pC >= 0);

            Scalar powAlphaPc = pow(state.vgAlpha()*pC, state.vgN());
            return -pow(powAlphaPc + 1, -state.vgM()-1)*
                        state.vgM()*powAlphaPc/pC*state.vgN();
        }

        /*!
         * \brief The relative permeability for the wetting phase of
         *        the medium implied by van Genuchten's
         *        parameterization.
         *
         * \param Sw_mob The mobile saturation of the wetting phase.
         * \param m      The "m" shape parameter according to van Genuchten
         */
        static Scalar krw(const State &state, Scalar Sw_mob)
        {
            Api::require<Api::VanGenuchtenParams>(state);
            assert(0 <= Sw_mob && Sw_mob <= 1);

            Scalar r = 1. - pow(1 - pow(Sw_mob, 1/state.vgM()), state.vgM());
            return sqrt(Sw_mob)*r*r;
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
            Api::require<Api::VanGenuchtenParams>(state);
            assert(0 <= Sw_mob && Sw_mob <= 1);

            Scalar r = pow(1 - pow(Sw_mob, 1/state.vgM()), state.vgM());
            return sqrt(1 - Sw_mob)*r*r;
        }
    };
}

#endif
