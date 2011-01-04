/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief   Regularized linear capillary pressure and
 *          relative permeability <-> saturation relations.
 */
#ifndef DUMUX_REGULARIZED_LINEAR_MATERIAL_HH
#define DUMUX_REGULARIZED_LINEAR_MATERIAL_HH

#include "linearmaterial.hh"
#include "regularizedlinearmaterialparams.hh"



#include <dumux/common/spline.hh>

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Implements a linear saturation-capillary pressure relation
 *
 * The entry pressure is reached at \f$\overline S_w = 1\f$, the maximum
 * capillary pressure is observed at \f$\overline S_w = 0\f$.
 *
 * It may seem strange to regularize a linear material, here comes the rationale:
 *
 * The relative permeabilities are 0 or 1 outside of the range of effective saturation.
 * However, the transition between the linearly changing and the constant part is not smooth but with a kink.
 * The Newton scheme does not like that. Therefore a smooth transition is accomplished by interpolating these
 * regions with a spline.
 *
 * An example of the regularization of the relative permeability is shown below:
 * \image html regularizedLinearKr.png
 *

 *
 * \see LinearMaterial
 */
template <class ScalarT, class ParamsT = RegularizedLinearMaterialParams<ScalarT> >
class RegularizedLinearMaterial
{
    typedef Dumux::LinearMaterial<ScalarT, ParamsT> LinearMaterial;

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f[
     p_C = (1 - \overline{S}_w) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     \f]
     *
     * \param Swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$ conversion from absolute saturation happened in EffToAbsLaw.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar pC(const Params &params, Scalar Swe)
    {
        return LinearMaterial::pC(params, Swe);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     S_w = 1 - \frac{p_C - p_{C,entry}}{p_{C,max} - p_{C,entry}}
     \f]
     *
     * \param pC        Capillary pressure \f$p_C\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return The effective saturation of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        return LinearMaterial::Sw(params, pC);
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     - (p_{C,max} - p_{C,min})
     \f]
     * \param Swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$ conversion from absolute saturation happened in EffToAbsLaw.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
    */
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        return LinearMaterial::dpC_dSw(params, Swe);
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     *
     * \param pC        Capillary pressure \f$p_C\f$
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        return LinearMaterial::dSw_dpC(params, pC);
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \param Swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$ conversion from absolute saturation happened in EffToAbsLaw.
     */
    static Scalar krw(const Params &params, Scalar Swe)
    {
        return relperm_(params, Swe);
    };

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \param Swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$ conversion from absolute saturation happened in EffToAbsLaw.
     */
    static Scalar krn(const Params &params, Scalar Swe)
    {
        Scalar Sne = 1 - Swe;
        return relperm_(params, Sne);
    }

private:
    static Scalar relperm_(const Params &params, Scalar S)
    {
        const Scalar lowS = params.krLowS();
        const Scalar highS = params.krHighS();


        const Scalar m = (1 - ((1 - highS) + lowS)/2 ) / (1 - (1 - highS) - lowS);

        // check whether the saturation is unpyhsical
        if (S >= 1.0)
            return 1.0;
        else if (S <= 0.0)
            return 0;
        // check wether the permeability needs to be regularized
        else if (S < lowS) {
            typedef Dumux::Spline<Scalar> Spline;
            Spline sp(0, lowS,
                      0, lowS/2,
                      0, m);
            return sp.eval(S);
        }
        else if (S > highS) {
            typedef Dumux::Spline<Scalar> Spline;
            Spline sp(highS, 1,
                      1 - (1 - highS)/2, 1,
                      m, 0);
            return sp.eval(S);
        }

        // straight line for S \in [lowS, highS]
        return lowS/2 + m*(S - lowS);
    }
};
}

#endif
