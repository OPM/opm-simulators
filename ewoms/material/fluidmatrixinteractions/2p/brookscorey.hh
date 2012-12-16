// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Philipp Nuske                                     *
 *   Copyright (C) 2010 by Bernd Flemisch                                    *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Ewoms::BrooksCorey
 */
#ifndef EWOMS_BROOKS_COREY_HH
#define EWOMS_BROOKS_COREY_HH

#include "brookscoreyparams.hh"

#include <algorithm>

namespace Ewoms {
/*!
 * \ingroup fluidmatrixinteractionslaws
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
template <class ScalarT, class ParamsT = BrooksCoreyParams<ScalarT> >
class BrooksCorey
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

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
     * \param Swe Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    static Scalar pC(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return params.pe()*std::pow(Swe, -1.0/params.lambda());
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
     * \param pc Capillary pressure \f$[Pa]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    static Scalar Sw(const Params &params, Scalar pc)
    {
        assert(pc > 0); // if we don't assume that, std::pow will screw up!

        return std::pow(pc/params.pe(), -params.lambda());
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to Brooks & Corey.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{p_e}{\lambda} \overline{S}_w^{-1/\lambda - 1}
     \f]
     *
     * \param Swe Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
    */
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return - params.pe()/params.lambda() * std::pow(Swe, -1/params.lambda() - 1);
    }

    /*!
     * \brief The partial derivative of the effective saturation with
     *        regard to the capillary pressure according to Brooks and
     *        Corey.
     *
     * \param pc Capillary pressure \f$[Pa]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        assert(pC >= 0);

        return -params.lambda()/params.pe() * std::pow(pC/params.pe(), - params.lambda() - 1);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Swe Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    static Scalar krw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return std::pow(Swe, 2.0/params.lambda() + 3);
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase with regard to the wetting saturation of the
     *        medium implied by the Brooks-Corey parameterization.
     *
     * \param Swe Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    static Scalar dkrw_dSw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return (2.0/params.lambda() + 3)*std::pow(Swe, 2.0/params.lambda() + 2);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Swe Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    static Scalar krn(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        Scalar exponent = 2.0/params.lambda() + 1;
        Scalar tmp = 1. - Swe;
        return tmp*tmp*(1. - std::pow(Swe, exponent));
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Swe Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    static Scalar dkrn_dSw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return
            2.0*(Swe - 1)*(
                1 +
                std::pow(Swe, 2.0/params.lambda())*(
                    1.0/params.lambda() + 1.0/2 -
                    Swe*(1.0/params.lambda() + 1.0/2)
                    )
                );
    }

};
}

#endif // BROOKS_COREY_HH
