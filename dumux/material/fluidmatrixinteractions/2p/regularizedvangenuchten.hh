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
 * \brief Implementation of a regularized version of van Genuchten's
 *        capillary pressure <-> saturation relation.
 */
#ifndef REGULARIZED_VAN_GENUCHTEN_HH
#define REGULARIZED_VAN_GENUCHTEN_HH

#include "vangenuchten.hh"
#include "regularizedvangenuchtenparams.hh"

#include <algorithm>


#include <dumux/common/spline.hh>

namespace Dumux
{
/*!\ingroup material
 *
 * \brief Implementation of a regularized version of van Genuchten's
 *        capillary pressure <-> saturation relation.
 *
 * This class bundles the "raw" curves as static members and doesn't
 * concern itself converting absolute to effective saturations and
 * vince versa.
 */
template <class ScalarT, class ParamsT = RegularizedVanGenuchtenParams<ScalarT> >
class RegularizedVanGenuchten
{
    typedef Dumux::VanGenuchten<ScalarT, ParamsT> VanGenuchten;

public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief \copybrief VanGenuchten::pC
     *
     * \copydetails VanGenuchten::pC
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
     * \brief \copybrief VanGenuchten::Sw
     *
     * \copydetails VanGenuchten::Sw
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
     * \brief \copybrief VanGenuchten::dpC_dSw
     *
     * \copydetails VanGenuchten::dpC_dSw
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
     * \brief \copybrief VanGenuchten::dSw_dpC
     *
     * \copydetails VanGenuchten::dSw_dpC
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
     * \brief \copybrief VanGenuchten::krw
     *
     * \copydetails VanGenuchten::krw
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
     * \brief \copybrief VanGenuchten::krn
     *
     * \copydetails VanGenuchten::krn
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
    static Scalar mLow_(const Params &params)
    {
        const Scalar SwThLow = params.pCLowSw();

        return VanGenuchten::dpC_dSw(params, SwThLow);
    }

    // the slope of the straight line used to regularize saturations
    // above the maximum saturation
    static Scalar mHigh_(const Params &params)
    {
        const Scalar SwThHigh = params.pCHighSw();

        Scalar pC_SwHigh = VanGenuchten::pC(params, SwThHigh);
        return (0 - pC_SwHigh)/(1.0 - SwThHigh);
    }
};

}

#endif
