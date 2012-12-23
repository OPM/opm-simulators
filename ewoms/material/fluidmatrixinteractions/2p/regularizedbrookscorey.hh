// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Philipp Nuske                                     *
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
 * \copydoc Ewoms::RegularizedBrooksCorey
 */
#ifndef REGULARIZED_BROOKS_COREY_HH
#define REGULARIZED_BROOKS_COREY_HH

#include "brookscorey.hh"
#include "regularizedbrookscoreyparams.hh"

#include <ewoms/common/spline.hh>

namespace Ewoms {
/*!
 * \ingroup fluidmatrixinteractionslaws
 * \brief Implementation of the regularized Brooks-Corey capillary
 *        pressure / relative permeability <-> saturation relation.
 *
 * This class bundles the "raw" curves as static members and doesn't
 * concern itself converting absolute to effective saturations and
 * vice versa.
 *
 * In order to avoid very steep gradients the marginal values are
 * "regularized".  This means that in stead of following the curve of
 * the material law in these regions, some linear approximation is
 * used.  Doing this is not worse than following the material
 * law. E.g. for very low wetting phase values the material laws
 * predict infinite values for \f$p_c\f$ which is completely
 * unphysical. In case of very high wetting phase saturations the
 * difference between regularized and "pure" material law is not big.
 *
 * Regularizing has the additional benefit of being numerically
 * friendly: Newton's method does not like infinite gradients.
 *
 * The implementation is accomplished as follows:
 * - check whether we are in the range of regularization
 *   - yes: use the regularization
 *   - no: forward to the standard material law.
 *
 * \see BrooksCorey
 */
template <class ScalarT, class ParamsT = RegularizedBrooksCoreyParams<ScalarT> >
class RegularizedBrooksCorey
{
    typedef Ewoms::BrooksCorey<ScalarT, ParamsT> BrooksCorey;

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief A regularized Brooks-Corey capillary pressure-saturation
     *        curve.
     *
     * This is a regularized variant of the Brooks-Corey curve. For
     * wetting phase saturations between a lower threshold saturation
     * and \f$S_w=1\f$, for other wetting phase saturations it is
     * regularized in a way which removes the singularity at
     * \f$S_w=0\f$, avoids kinks and allows the capillary pressure to
     * reach arbitrary values. (Albeit, to reach a given capillary
     * pressure, the saturations can become unphysical). The
     * regularization is done in the following way:
     *
     * - For wetting phase saturations lower than the threshold
     *   saturation, the \f$p_c(S_w)\f$ curve is extrapolated using a
     *   straight line exhibiting the slope unregularized capillary
     *   pressure curve at the threshold saturation.
     * - For wetting phase saturations larger than 1, the Brooks-Corey
     *   curve is extrapolated using a straight line that exhibits the
     *   slope of the unregularized Brooks-Corey curve at \f$S_w =
     *   1\f$
     *
     * \sa BrooksCorey::pC
     */
    static Scalar pC(const Params &params, Scalar Sw)
    {
        const Scalar Sthres = params.thresholdSw();

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative is calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Sw <= Sthres) {
            Scalar m = BrooksCorey::dpC_dSw(params, Sthres);
            Scalar pC_SwLow = BrooksCorey::pC(params, Sthres);
            return pC_SwLow + m*(Sw - Sthres);
        }
        else if (Sw > 1.0) {
            Scalar m = BrooksCorey::dpC_dSw(params, 1.0);
            Scalar pC_SwHigh = BrooksCorey::pC(params, 1.0);
            return pC_SwHigh + m*(Sw - 1.0);
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real Brooks-Corey law...
        return BrooksCorey::pC(params, Sw);
    }

    /*!
     * \brief A regularized Brooks-Corey saturation-capillary pressure
     *        curve.
     *
     * This is the inverse of the pC() method.
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        const Scalar Sthres = params.thresholdSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized version of the
        // Brooks-Corey law. If the input capillary pressure is
        // smaller than the entry pressure, make sure that we will
        // regularize.
        Scalar Sw = 1.5;
        if (pC >= params.pe())
            Sw = BrooksCorey::Sw(params, pC);

        // make sure that the capilary pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Sw <= Sthres) {
            // invert the low saturation regularization of pC()
            Scalar m = BrooksCorey::dpC_dSw(params, Sthres);
            Scalar pC_SwLow = BrooksCorey::pC(params, Sthres);
            return Sthres + (pC - pC_SwLow)/m;
        }
        else if (Sw > 1.0) {
            Scalar m = BrooksCorey::dpC_dSw(params, 1.0);
            Scalar pC_SwHigh = BrooksCorey::pC(params, 1.0);
            return 1.0 + (pC - pC_SwHigh)/m;;
        }

        return BrooksCorey::Sw(params, pC);
    }

    /*!
     * \brief The derivative of the regularized Brooks-Corey capillary
     *        pressure-saturation curve.
     */
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        const Scalar Sthres = params.thresholdSw();

        // derivative of the regualarization
        if (Sw <= Sthres) {
            // calculate the slope of the straight line used in pC()
            Scalar m = BrooksCorey::dpC_dSw(params, Sthres);
            return m;
        }
        else if (Sw > 1.0) {
            // calculate the slope of the straight line used in pC()
            Scalar m = BrooksCorey::dpC_dSw(params, 1.0);
            return m;
        }

        return BrooksCorey::dpC_dSw(params, Sw);
    }

    /*!
     * \brief The derivative of the regularized Brooks-Corey
     *        saturation-capillary pressure curve.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        const Scalar Sthres = params.thresholdSw();

        // calculate the saturation which corresponds to the
        // saturation in the non-regularized version of the
        // Brooks-Corey law
        Scalar Sw;
        if (pC < params.pe())
            Sw = 1.5; // make sure we regularize (see below)
        else
            Sw = BrooksCorey::Sw(params, pC);

        // derivative of the regularization
        if (Sw <= Sthres) {
            // calculate the slope of the straight line used in pC()
            Scalar m = BrooksCorey::dpC_dSw(params, Sthres);
            return 1/m;
        }
        else if (Sw > 1.0) {
            // calculate the slope of the straight line used in pC()
            Scalar m = BrooksCorey::dpC_dSw(params, 1.0);
            return 1/m;
        }
        return 1.0/BrooksCorey::dpC_dSw(params, Sw);
    }

    /*!
     * \brief Regularized version of the relative permeability of the
     *        wetting phase of the Brooks-Corey curves.
     *
     * The approach for regularization is very similar to the one of
     * the capillary pressure, but it does not avoid kinks:
     * - For wetting phase saturations between 0 and 1, use the
     *   unregularized Brooks-Corey wetting phase relative
     *   permeability
     * - For wetting phase saturations smaller than 0, return 0
     * - For wetting phase saturations larger than 1, return 1
     *
     * \sa BrooksCorey::krw
     */
    static Scalar krw(const Params &params, Scalar Sw)
    {
        if (Sw <= 0.0)
            return 0.0;
        else if (Sw >= 1.0)
            return 1.0;

        return BrooksCorey::krw(params, Sw);
    }

    /*!
     * \brief Regularized version of the relative permeability of the
     *        non-wetting phase of the Brooks-Corey curves.
     *
     * The approach for regularization is very similar to the one of
     * the capillary pressure, but it does not avoid kinks:
     * - For wetting phase saturations between 0 and 1, use the
     *   unregularized Brooks-Corey non-wetting phase relative
     *   permeability
     * - For wetting phase saturations smaller than 0, return 1
     * - For wetting phase saturations larger than 1, return 0
     *
     * \sa BrooksCorey::krn
     */
    static Scalar krn(const Params &params, Scalar Sw)
    {
        if (Sw >= 1.0)
            return 0.0;
        else if (Sw <= 0.0)
            return 1.0;

        return BrooksCorey::krn(params, Sw);
    }
};
}

#endif
