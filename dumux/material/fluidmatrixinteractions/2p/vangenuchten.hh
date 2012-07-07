// $Id: vangenuchten.hh 3777 2010-06-24 06:46:46Z bernd $
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

#include "vangenuchtenparams.hh"

#include <algorithm>

#include <math.h>
#include <assert.h>

namespace Dumux
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
template <class ScalarT, class ParamsT = VanGenuchtenParams<ScalarT> >
class VanGenuchten
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     * Van Genuchten's empirical capillary pressure <-> saturation
     * function is given by
     * \f[
     p_C = (\overline{S}_w^{-1/m} - 1)^{1/n}/\alpha
     \f]
     * \param Sw Effective saturation of of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar pC(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);
        return pow(pow(Sw, -1.0/params.vgM()) - 1, 1.0/params.vgN())/params.vgAlpha();
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
    static Scalar Sw(const Params &params, Scalar pC)
    {
        assert(pC >= 0);

        return pow(pow(params.vgAlpha()*pC, params.vgN()) + 1, -params.vgM());
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
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        Scalar powSw = pow(Sw, -1/params.vgM());
        return - 1/params.vgAlpha() * pow(powSw - 1, 1/params.vgN() - 1)/params.vgN()
            * powSw/Sw/params.vgM();
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        assert(pC >= 0);

        Scalar powAlphaPc = pow(params.vgAlpha()*pC, params.vgN());
        return -pow(powAlphaPc + 1, -params.vgM()-1)*
            params.vgM()*powAlphaPc/pC*params.vgN();
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar krw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        Scalar r = 1. - pow(1 - pow(Sw, 1/params.vgM()), params.vgM());
        return sqrt(Sw)*r*r;
    };

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase in regard to the wetting saturation of the
     *        medium implied by the van Genuchten parameterization.
     *
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar dkrw_dSw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        const Scalar x = 1 - std::pow(Sw, 1.0/params.vgM());
        const Scalar xToM = std::pow(x, params.vgM());
        return (1 - xToM)/std::sqrt(Sw) * ( (1 - xToM)/2 + 2*xToM*(1-x)/x );
    };


    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar krn(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        return
            pow(1 - Sw, 1.0/3) *
            pow(1 - pow(Sw, 1/params.vgM()), 2*params.vgM());
    };

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the van Genuchten
     *        parameterization.
     *
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar dkrn_dSw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        const Scalar x = std::pow(Sw, 1.0/params.vgM());
        return
            -std::pow(1 - x, 2*params.vgM())
            *std::pow(1 - Sw, -2/3)
            *(1.0/3 + 2*x/Sw);
    }

};
}

#endif
