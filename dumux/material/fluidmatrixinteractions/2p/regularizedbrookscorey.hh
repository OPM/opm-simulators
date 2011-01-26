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
 * \brief Implementation of a regularized version of the Brooks-Corey
 *        capillary pressure / relative permeability  <-> saturation relation.
 */
#ifndef REGULARIZED_BROOKS_COREY_HH
#define REGULARIZED_BROOKS_COREY_HH

#include "brookscorey.hh"
#include "regularizedbrookscoreyparams.hh"



#include <dumux/common/spline.hh>

namespace Dumux
{
/*!\ingroup fluidmatrixinteractionslaws
 *
 * \brief Implementation of the regularized  Brooks-Corey
 *        capillary pressure / relative permeability  <-> saturation relation.
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
 *         For an example figure of the regularization: RegularizedVanGenuchten
 *
 * \see BrooksCorey
 */
template <class ScalarT, class ParamsT = RegularizedBrooksCoreyParams<ScalarT> >
class RegularizedBrooksCorey
{
    typedef Dumux::BrooksCorey<ScalarT, ParamsT> BrooksCorey;

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief   A regularized Brooks-Corey capillary pressure-saturation
     *          curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$p_c(S_w)\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$ \overline S_w =1\f$ by a straight line (yes, there is a kink :-( ).
     *
     *  For not-regularized part:
     *
         \copydetails BrooksCorey::pC()
     */


    static Scalar pC(const Params &params, Scalar Swe)
    {
        const Scalar Sthres = params.thresholdSw();

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative is calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Swe <= Sthres) {
            Scalar m = BrooksCorey::dpC_dSw(params, Sthres);
            Scalar pC_SweLow = BrooksCorey::pC(params, Sthres);
            return pC_SweLow + m*(Swe - Sthres);
        }
        else if (Swe > 1) {
            Scalar m = BrooksCorey::dpC_dSw(params, 1.0);
            Scalar pC_SweHigh = BrooksCorey::pC(params, 1.0);
            return pC_SweHigh + m*(Swe - 1.0);
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real Brooks-Corey law...
        return BrooksCorey::pC(params, Swe);
    }

    /*!
     * \brief   A regularized Brooks-Corey saturation-capillary pressure curve.
     *
     * regularized part:
     *    - low saturation:  extend the \f$p_c(S_w)\f$ curve with the slope at the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$ \overline S_w =1\f$ by a straight line (yes, there is a kink :-( ).
     *
     *  The according quantities are obtained by exploiting theorem of intersecting lines.
     *
     *  For not-regularized part:
     *
         \copydetails BrooksCorey::Sw()
     *
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        const Scalar Sthres = params.thresholdSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized version of
        // the Brooks-Corey law
        Scalar Swe = BrooksCorey::Sw(params, pC);

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Swe <= Sthres) {
            // invert the low saturation regularization of pC()
            Scalar m = BrooksCorey::dpC_dSw(params, Sthres);
            Scalar pC_SweLow = BrooksCorey::pC(params, Sthres);
            return Sthres + (pC - pC_SweLow)/m;
        }
        else if (Swe > 1) {
            Scalar m = BrooksCorey::dpC_dSw(params, 1.0);
            Scalar pC_SweHigh = BrooksCorey::pC(params, 1.0);
            return 1.0 + (pC - pC_SweHigh)/m;;
        }

        return BrooksCorey::Sw(params, pC);
    }

    /*!
     * \brief A regularized version of the partial derivative
     *        of the \f$p_c(\overline S_w)\f$ w.r.t. effective saturation
     *        according to Brooks & Corey.
     *
     * regularized part:
     *    - low saturation:  use the slope of the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$ \overline S_w =1\f$ by a straight line and use that slope (yes, there is a kink :-( ).
     *
     *        For not-regularized part:
     *
       \copydetails BrooksCorey::dpC_dSw()
     *
    */
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        const Scalar Sthres = params.thresholdSw();

        // derivative of the regualarization
        if (Swe <= Sthres) {
            // calculate the slope of the straight line used in pC()
            Scalar m = BrooksCorey::dpC_dSw(params, Sthres);
            return m;
        }
        else if (Swe > 1.0) {
            // calculate the slope of the straight line used in pC()
            Scalar m = BrooksCorey::dpC_dSw(params, 1.0);
            return m;
        }

        return BrooksCorey::dpC_dSw(params, Swe);
    }

    /*!
     * \brief A regularized version of the partial derivative
     *        of the \f$\overline S_w(p_c)\f$ w.r.t. cap.pressure
     *        according to Brooks & Corey.
     *
     *  regularized part:
     *    - low saturation:  use the slope of the regularization point (i.e. no kink).
     *    - high saturation: connect the high regularization point with \f$ \overline S_w =1\f$ by a straight line and use that slope (yes, there is a kink :-( ).
     *
     *        For not-regularized part:
        \copydetails BrooksCorey::dSw_dpC()
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        const Scalar Sthres = params.thresholdSw();

        //instead of return value = inf, return a very large number
        if (params.pe() == 0.0)
        {
            return 1e100;
        }

        // calculate the saturation which corresponds to the
        // saturation in the non-regularized version of the
        // Brooks-Corey law
        Scalar Swe;
        if (pC < 0)
            Swe = 1.5; // make sure we regularize below
        else
            Swe = BrooksCorey::Sw(params, pC);

        // derivative of the regularization
        if (Swe <= Sthres) {
            // calculate the slope of the straight line used in pC()
            Scalar m = BrooksCorey::dpC_dSw(params, Sthres);
            return 1/m;
        }
        else if (Swe > 1.0) {
            // calculate the slope of the straight line used in pC()
            Scalar m = BrooksCorey::dpC_dSw(params, 1.0);
            return 1/m;
        }
        return 1.0/BrooksCorey::dpC_dSw(params, Swe);
    }

    /*!
     * \brief   Regularized version of the  relative permeability
     *          for the wetting phase of
     *          the medium implied by the Brooks-Corey
     *          parameterization.
     *
     *  regularized part:
     *    - below \f$ \overline S_w =0\f$:                  set relative permeability to zero
     *    - above \f$ \overline S_w =1\f$:                  set relative permeability to one
     *    - between \f$ 0.95 \leq \overline S_w \leq 1\f$:  use a spline as interpolation
     *
     *  For not-regularized part:
        \copydetails BrooksCorey::krw()
     */
    static Scalar krw(const Params &params, Scalar Swe)
    {
        if (Swe <= 0)
            return 0;
        else if (Swe >= 1)
            return 1.0;
        else if (Swe >= 1 - 0.05) {
            Scalar m1 = BrooksCorey::dkrw_dSw(params, 1.0 - 0.05);
            Scalar y1 = BrooksCorey::krw(params, 1.0 - 0.05);
            Dumux::Spline<Scalar> sp(1 - 0.05, 1.0,
                                     y1, 1.0,
                                     m1, 0);
            return sp.eval(Swe);
        }

        return BrooksCorey::krw(params, Swe);
    };

    /*!
     * \brief   Regularized version of the  relative permeability
     *          for the non-wetting phase of
     *          the medium implied by the Brooks-Corey
     *          parameterization.
     *
     * regularized part:
     *    - below \f$ \overline S_w =0\f$:                  set relative permeability to zero
     *    - above \f$ \overline S_w =1\f$:                  set relative permeability to one
     *    - for \f$ 0 \leq \overline S_w \leq 0.05 \f$:     use a spline as interpolation
     *
         \copydetails BrooksCorey::krn()
     *
     */
    static Scalar krn(const Params &params, Scalar Swe)
    {
        if (Swe >= 1)
            return 0;
        // check if we need to regularize the relative permeability
        else if (Swe <= 0)
            return 1.0;
        else if (Swe < 0.05) {
            Scalar m1 = BrooksCorey::dkrn_dSw(params, 0.05);
            Scalar y1 = BrooksCorey::krn(params, 0.05);
            Dumux::Spline<Scalar> sp(0.0, 0.05,
                                     1.0, y1,
                                     0, m1);
            return sp.eval(Swe);
        }
        return BrooksCorey::krn(params, Swe);
    }
};
}

#endif
