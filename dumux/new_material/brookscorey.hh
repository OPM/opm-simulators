/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
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
 * \file
 *
 * \brief Implementation of the capillary pressure <-> saturation
 *         relation due to Brooks and Corey.
 */
#ifndef DUMUX_BROOKS_COREY_HH
#define DUMUX_BROOKS_COREY_HH

#include <dumux/new_material/brookscoreycontext.hh>

#include <algorithm>

#include <math.h>
#include <assert.h>


namespace Dune
{
/*!\ingroup material
 *
 * \brief Implementation of the Brooks-Corey capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vince versa.
 *
 * \sa BrooksCorey, BrooksCoreyTwophase
 */
template <class ContextT>
class BrooksCorey
{
public:
    typedef ContextT Context;
    typedef typename Context::Scalar Scalar;

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     * The Brooks-Corey empirical capillary pressure <-> saturation
     * function is given by
     * \f[
     p_C = p_e\overline{S}_w^{-1/\alpha}
     \f]
     * \param Swe   Effective saturation of of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar pC(const Context &context, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return context.pe()*pow(Swe, -1.0/context.alpha());
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     \overline{S}_w = (\frac{p_C}{p_e})^{-\alpha}
     \f]
     *
     * \param pC Capillary pressure \f$\p_C\f$
     * \return The effective saturaion of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar Sw(const Context &context, Scalar pC)
    {
        assert(pC >= 0);

        Scalar tmp = pow(pC/context.pe(), -context.alpha());
        return std::min(std::max(tmp, Scalar(0.0)), Scalar(1.0));
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{p_e}{\alpha} \overline{S}_w^{-1/\alpha - 1}
     \f]
    */
    static Scalar dpC_dSw(const Context &context, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return - context.pe()/context.alpha() * pow(Swe, -1/context.alpha() - 1);
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Context &context, Scalar pC)
    {
        assert(pC >= 0);

        return -context.alpha()/context.pe() * pow(pC/context.pe(), - context.alpha() - 1);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     */
    static Scalar krw(const Context &context, Scalar Sw_mob)
    {
        assert(0 <= Sw_mob && Sw_mob <= 1);

        return pow(Sw_mob, (2. + 3*context.alpha()) / context.alpha());
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

        Scalar exponent = (2. + context.alpha())/context.alpha();
        Scalar tmp = 1. - Sw_mob;
        return tmp*tmp*(1. - pow(Sw_mob, exponent));
    }
};
}

#endif // BROOKS_COREY_HH
