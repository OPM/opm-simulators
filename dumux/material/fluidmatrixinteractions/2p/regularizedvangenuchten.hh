// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 * \brief Implementation of the regularized version of the van Genuchten's
 *        capillary pressure / relative permeability  <-> saturation relation.
 */
#ifndef REGULARIZED_VAN_GENUCHTEN_HH
#define REGULARIZED_VAN_GENUCHTEN_HH

#include "vangenuchten.hh"
#include "regularizedvangenuchtenparams.hh"

#include <algorithm>


#include <dumux/common/spline.hh>

namespace Dumux
{
/*!\ingroup fluidmatrixinteractionslaws


 * \brief Implementation of the regularized  van Genuchten's
 *        capillary pressure / relative permeability  <-> saturation relation.
 *
 *        This class bundles the "raw" curves as
 *        static members and doesn't concern itself converting
 *        absolute to effective saturations and vice versa.
 *
 *        In order to avoid very steep gradients the marginal values are "regularized".
 *        This means that in stead of following the curve of the material law in these regions, some linear approximation is used.
 *        Doing this is not worse than following the material law. E.g. for very low wetting phase values the material
 *        laws predict infinite values for \f$p_c\f$ which is completely unphysical. In case of very high wetting phase
 *        saturations the difference between regularized and "pure" material law is not big.
 *
 *        Regularizing has the additional benefit of being numerically friendly: Newton's method does not like infinite gradients.
 *
 *        The implementation is accomplished as follows:
 *        - check whether we are in the range of regularization
 *         - yes: use the regularization
 *         - no: forward to the standard material law.
 *
 *        An example of the regularization of the capillary pressure curve is shown below:
 *        \image html regularizedVanGenuchten.png
 *
 * \see VanGenuchten
 */
template <class ScalarT, class ParamsT = RegularizedVanGenuchtenParams<ScalarT> >
class RegularizedVanGenuchten
{
    typedef Dumux::VanGenuchten<ScalarT, ParamsT> VanGenuchten;

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

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
         \copydetails VanGenuchten::pC()
     */
    static Scalar pC(const Params &params, Scalar Swe)
    {
        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar SwThLow = params.pCLowSw();
        const Scalar SwThHigh = params.pCHighSw();

        // make sure that the capillary pressure observes a derivative
        // != 0 for 'illegal' saturations. This is favourable for the
        // newton solver (if the derivative is calculated numerically)
        // in order to get the saturation moving to the right
        // direction if it temporarily is in an 'illegal' range.
        if (Swe < SwThLow) {
            return VanGenuchten::pC(params, SwThLow) + mLow_(params)*(Swe - SwThLow);
        }
        else if (Swe > SwThHigh) {
            return VanGenuchten::pC(params, SwThHigh) + mHigh_(params)*(Swe - SwThHigh);
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real van genuchten law...
        return VanGenuchten::pC(params, Swe);
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
    static Scalar Sw(const Params &params, Scalar pC)
    {
        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized verision of van
        // Genuchten's law
        Scalar Sw;
        if (pC <= 0)
            // make sure we invert the regularization
            Sw = 1.5;
        else
            Sw = VanGenuchten::Sw(params, pC);

        // retrieve the low and the high threshold saturations for the
        // unregularized capillary pressure curve from the parameters
        const Scalar SwThLow = params.pCLowSw();
        const Scalar SwThHigh = params.pCHighSw();

        // invert the regularization if necessary
        if (Sw <= SwThLow) {
            // invert the low saturation regularization of pC()
            Scalar pC_SwLow = VanGenuchten::pC(params, SwThLow);
            return (pC - pC_SwLow)/mLow_(params) + SwThLow;
        }
        else if (Sw >= SwThHigh) {
            // invert the high saturation regularization of pC()
            Scalar pC_SwHigh = VanGenuchten::pC(params, SwThHigh);
            return (pC - pC_SwHigh)/mHigh_(params) + SwThHigh;
        }

        return Sw;
    }

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
      \copydetails VanGenuchten::dpC_dSw()
    *
    */
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        // derivative of the regualarization
        if (Swe < params.pCLowSw()) {
            // the slope of the straight line used in pC()
            return mLow_(params);
        }
        else if (Swe > params.pCHighSw()) {
            // the slope of the straight line used in pC()
            return mHigh_(params);
        }

        return VanGenuchten::dpC_dSw(params, Swe);
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
        \copydetails VanGenuchten::dSw_dpC()
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized verision of van
        // Genuchten's law
        Scalar Sw;
        if (pC < 0)
            Sw = 1.5; // make sure we regularize below
        else
            Sw = VanGenuchten::Sw(params, pC);

        // derivative of the regularization
        if (Sw < params.pCLowSw()) {
            // same as in dpC_dSw() but inverted
            return 1/mLow_(params);
        }
        if (Sw > params.pCHighSw()) {
            // same as in dpC_dSw() but inverted
            return 1/mHigh_(params);
        }

        return VanGenuchten::dSw_dpC(params, pC);
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
    static Scalar krw(const Params &params, Scalar Swe)
    {
        // retrieve the high threshold saturation for the
        // unregularized relative permeability curve of the wetting
        // phase from the parameters
        const Scalar SwThHigh = params.krwHighSw();

        if (Swe < 0)
            return 0;
        else if (Swe > 1)
            return 1;
        else if (Swe > SwThHigh) {
            typedef Dumux::Spline<Scalar> Spline;
            Spline sp(SwThHigh, 1.0, // x1, x2
                      VanGenuchten::krw(params, SwThHigh), 1.0, // y1, y2
                      VanGenuchten::dkrw_dSw(params, SwThHigh), 0); // m1, m2
            return sp.eval(Swe);
        }

        return VanGenuchten::krw(params, Swe);
    };

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
    static Scalar krn(const Params &params, Scalar Swe)
    {
        // retrieve the low threshold saturation for the unregularized
        // relative permeability curve of the non-wetting phase from
        // the parameters
        const Scalar SwThLow = params.krnLowSw();

        if (Swe <= 0)
            return 1;
        else if (Swe >= 1)
            return 0;
        else if (Swe < SwThLow) {
            typedef Dumux::Spline<Scalar> Spline;
            Spline sp(0.0, SwThLow, // x1, x2
                      1.0, VanGenuchten::krn(params, SwThLow), // y1, y2
                      0.0, VanGenuchten::dkrn_dSw(params, SwThLow)); // m1, m2
            return sp.eval(Swe);
        }

        return VanGenuchten::krn(params, Swe);
    }

private:
    // the slope of the straight line used to regularize saturations
    // below the minimum saturation

    /*!
     * \brief   The slope of the straight line used to regularize
     *          saturations below the minimum saturation.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar mLow_(const Params &params)
    {
        const Scalar SwThLow = params.pCLowSw();

        return VanGenuchten::dpC_dSw(params, SwThLow);
    }

    /*!
     * \brief   The slope of the straight line used to regularize
     *          saturations above the minimum saturation.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     */
    static Scalar mHigh_(const Params &params)
    {
        const Scalar SwThHigh = params.pCHighSw();

        Scalar pC_SwHigh = VanGenuchten::pC(params, SwThHigh);
        return (0 - pC_SwHigh)/(1.0 - SwThHigh);
    }
};

}

#endif
