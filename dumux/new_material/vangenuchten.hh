/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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

#include <dumux/new_material/vangenuchtencontext.hh>

#include <algorithm>

#include <math.h>
#include <assert.h>

namespace Dune
{
/*!
 * \ingroup material
 *
 * \brief Implementation of van Genuchten's capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vince versa.
 *
 * \sa VanGenuchten, VanGenuchtenTwophase
 */
template <class ContextT>
class VanGenuchten
{
public:
    typedef ContextT Context;
    typedef typename Context::Scalar Scalar;

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
    static Scalar pC(const Context &context, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);
        return pow(pow(Swe, -1.0/context.vgM()) - 1, 1.0/context.vgN())/context.vgAlpha();
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
    static Scalar Sw(const Context &context, Scalar pC)
    {
        assert(pC >= 0);

        return pow(pow(context.vgAlpha()*pC, context.vgN()) + 1, -context.vgM());
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
    static Scalar dpC_dSw(const Context &context, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        Scalar powSwe = pow(Swe, -1/context.vgM());
        return - 1/context.vgAlpha() * pow(powSwe - 1, 1/context.vgN() - 1)/context.vgN()
            * powSwe/Swe/context.vgM();
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Context &context, Scalar pC)
    {
        assert(pC >= 0);

        Scalar powAlphaPc = pow(context.vgAlpha()*pC, context.vgN());
        return -pow(powAlphaPc + 1, -context.vgM()-1)*
            context.vgM()*powAlphaPc/pC*context.vgN();
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     * \param m      The "m" shape parameter according to van Genuchten
     */
    static Scalar krw(const Context &context, Scalar Sw_mob)
    {
        assert(0 <= Sw_mob && Sw_mob <= 1);

        Scalar r = 1. - pow(1 - pow(Sw_mob, 1/context.vgM()), context.vgM());
        return sqrt(Sw_mob)*r*r;
    };

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     */
    static Scalar krn(const Context &context, Scalar Sw_mob)
    {
        assert(0 <= Sw_mob && Sw_mob <= 1);

        Scalar r = pow(1 - pow(Sw_mob, 1/context.vgM()), context.vgM());
        return sqrt(1 - Sw_mob)*r*r;
    }
};
}

#endif
