// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 *
 * \brief Implementation of the capillary pressure and
 * relative permeability <-> saturation relations according to Brooks and Corey.
 *
 */
#ifndef DUMUX_BROOKS_COREY_HH
#define DUMUX_BROOKS_COREY_HH

#include "brookscoreyparams.hh"

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Implementation of the Brooks-Corey capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vice versa.
 *
 * For general info: EffToAbsLaw
 *
 *\see BrooksCoreyParams
 */
template <class ScalarT, class ParamsT = BrooksCoreyParams<ScalarT> >
class BrooksCorey
{
public:
    typedef ParamsT     Params;
    typedef typename    Params::Scalar Scalar;

    /*!
     * \brief The capillary pressure-saturation curve according to Brooks & Corey.
     *
     * The Brooks-Corey empirical  capillary pressure <-> saturation
     * function is given by
     *
     *  \f[
        p_C = p_e\overline{S}_w^{-1/\alpha}
    *  \f]
    *
     * \param Swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Capillary pressure calculated by Brooks & Corey constitutive relation.
     */
    static Scalar pC(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return params.pe()*pow(Swe, -1.0/params.alpha());
    }

    /*!
     * \brief The saturation-capillary pressure curve according to Brooks & Corey.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     \overline{S}_w = (\frac{p_C}{p_e})^{-\alpha}
     \f]
     *
     * \param pC        Capillary pressure \f$p_C\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Effective wetting phase saturation calculated as inverse of BrooksCorey constitutive relation.
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        assert(pC >= 0);

        Scalar tmp = pow(pC/params.pe(), -params.alpha());
        return std::min(std::max(tmp, Scalar(0.0)), Scalar(1.0));
    }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to Brooks & Corey.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{p_e}{\alpha} \overline{S}_w^{-1/\alpha - 1}
     \f]
     *
     * \param Swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of \f$p_c\f$ w.r.t. effective saturation according to Brooks & Corey.
    */
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return - params.pe()/params.alpha() * pow(Swe, -1/params.alpha() - 1);
    }

    /*!
     * \brief The partial derivative of the effective
     *        saturation w.r.t. the capillary pressure according to Brooks & Corey.
     *
     * \param pC        Capillary pressure \f$p_C\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of effective saturation w.r.t. \f$p_c\f$ according to Brooks & Corey.
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
     * \param Swe       The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the wetting phase calculated as implied by Brooks & Corey.
     */
    static Scalar krw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return pow(Swe, 2.0/params.alpha() + 3);
    };

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase with regard to the wetting saturation of the
     *        medium implied by the Brooks-Corey parameterization.
     *
     * \param Swe       The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Derivative of the relative permeability of the wetting phase w.r.t. effective wetting phase saturation calculated as implied by Brooks & Corey.
     */
    static Scalar dkrw_dSw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return (2.0/params.alpha() + 3)*pow(Swe, 2.0/params.alpha() + 2);
    };

    /*!
     * \brief The relative permeability for the non-wetting phase of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Swe       The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the non-wetting phase calculated as implied by Brooks & Corey.
     */
    static Scalar krn(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        Scalar exponent = 2.0/params.alpha() + 1;
        Scalar tmp = 1. - Swe;
        return tmp*tmp*(1. - pow(Swe, exponent));
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Swe       The mobile saturation of the wetting phase.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Derivative of the relative permeability of the non-wetting phase w.r.t. effective wetting phase saturation calculated as implied by Brooks & Corey.
     */
    static Scalar dkrn_dSw(const Params &params, Scalar Swe)
    {
        assert(0 <= Swe && Swe <= 1);

        return
            2.0*(Swe - 1)*(
                1 +
                pow(Swe, 2.0/params.alpha())*(
                    1.0/params.alpha() + 1.0/2 -
                    Swe*(1.0/params.alpha() + 1.0/2)
                    )
                );
    };

};
}

#endif // BROOKS_COREY_HH
