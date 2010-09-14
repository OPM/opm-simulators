// $Id$
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
 * \file
 *
 * \brief Implementation of a regularized version of the Brooks-Corey
 *        capillary pressure <-> saturation relation.
 */
#ifndef REGULARIZED_BROOKS_COREY_HH
#define REGULARIZED_BROOKS_COREY_HH

#include "brookscorey.hh"
#include "regularizedbrookscoreyparams.hh"



#include <dumux/common/spline.hh>

namespace Dumux
{
/*!\ingroup material
 *
 * \brief Implementation of the Brooks-Corey capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves as
 *        static members and doesn't concern itself converting
 *        absolute to effective saturations and vince versa.
 *
 * \sa BrooksCorey, BrooksCoreyTwophase
 */
template <class ScalarT, class ParamsT = RegularizedBrooksCoreyParams<ScalarT> >
class RegularizedBrooksCorey
{
    typedef Dumux::BrooksCorey<ScalarT, ParamsT> BrooksCorey;

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief A regularized Brooks-Corey capillary pressure-saturation
     *        curve.
     */
    static Scalar pC(const Params &params, Scalar Swe)
    {
        const Scalar Sthres = params.thresholdSw();

        // make sure that the capilarry pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
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
     * \brief The saturation-capillary pressure curve.
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        const Scalar Sthres = params.thresholdSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized version of
        // the Brooks-Corey law
        Scalar Swe = BrooksCorey::Sw(params, pC);

        // make sure that the capilarry pressure observes a
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
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
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
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        const Scalar Sthres = params.thresholdSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized verision of the
        // Brooks-Corey law
        Scalar Swe;
        if (pC < 0)
            Swe = 1.5; // make sure we regularize below
        else
            Swe = BrooksCorey::Sw(params, pC);

        // derivative of the regualarization
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

        return BrooksCorey::dpC_dSw(params, pC);
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
        if (Sw <= 0)
            return 0;
        else if (Sw > 0)
            return 1.0;
        else if (Sw >= 1 - 0.01) {
            Scalar m1 = BrooksCorey::dkrw_dSw(params, 1.0 - 0.01);
            Scalar y1 = BrooksCorey::krw(params, 1.0 - 0.01);
            Dumux::Spline<Scalar> sp(1 - 0.01, 1.0,
                                     y1, 1.0,
                                     m1, 0);
            return sp.eval(Sw);
        }

        return BrooksCorey::krw(params, Sw);
    };

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Sw The mobile saturation of the wetting phase.
     */
    static Scalar krn(const Params &params, Scalar Sw)
    {
        if (Sw >= 1)
            return 0;
        // check if we need to regularize the relative permeability
        else if (Sw <= 0)
            return 1.0;
        else if (Sw < 0.01) {
            Scalar m1 = BrooksCorey::dkrn_dSw(params, 0.01);
            Scalar y1 = BrooksCorey::krn(params, 0.01);
            Dumux::Spline<Scalar> sp(0.0, 0.01,
                                     1.0, y1,
                                     0, m1);
            return sp.eval(Sw);
        }

        return BrooksCorey::krn(params, Sw);
    }
};
}

#endif
