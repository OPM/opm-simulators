// $Id$
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

#include "brookscoreyparams.hh"

#include <algorithm>



namespace Dumux
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
template <class ScalarT, class ParamsT = BrooksCoreyParams<ScalarT> >
class BrooksCorey
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

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
    static Scalar pC(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return params.pe()*pow(Swe, -1.0/params.alpha());
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     \overline{S}_w = (\frac{p_C}{p_e})^{-\alpha}
     \f]
     *
     * \param pC Capillary pressure \f$p_C\f$
     * \return The effective saturaion of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        assert(pC >= 0);

        Scalar tmp = pow(pC/params.pe(), -params.alpha());
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
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return - params.pe()/params.alpha() * pow(Swe, -1/params.alpha() - 1);
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        assert(pC >= 0);

        return -params.alpha()/params.pe() * pow(pC/params.pe(), - params.alpha() - 1);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar krw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        return pow(Sw, 2.0/params.alpha() + 3);
    };

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase in regard to the wetting saturation of the
     *        medium implied by the Brooks-Corey parameterization.
     *
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar dkrw_dSw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        return (2.0/params.alpha() + 3)*pow(Sw, 2.0/params.alpha() + 2);
    };

    /*!
     * \brief The relative permeability for the non-wetting phase of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar krn(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        Scalar exponent = 2.0/params.alpha() + 1;
        Scalar tmp = 1. - Sw;
        return tmp*tmp*(1. - pow(Sw, exponent));
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar dkrn_dSw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        return
            2.0*(Sw - 1)*(
                1 +
                pow(Sw, 2.0/params.alpha())*(
                    1.0/params.alpha() + 1.0/2 -
                    Sw*(1.0/params.alpha() + 1.0/2)
                    )
                );
    };

};
}

#endif // BROOKS_COREY_HH
