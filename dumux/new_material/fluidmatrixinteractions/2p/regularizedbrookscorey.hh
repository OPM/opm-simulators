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

#include <algorithm>

#include <math.h>
#include <assert.h>

#include <dumux/auxiliary/spline.hh>

namespace Dune
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
template <class ParamsT>
class RegularizedBrooksCorey
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    typedef Dune::BrooksCorey<ParamsT> BrooksCorey;

    /*!
     * \brief A regularized Brooks-Corey capillary pressure-saturation
     *        curve.
     */
    static Scalar pC(const Params &params, Scalar Swe)
    {
        const Scalar Sthres = params.pCLowSw();

        // make sure that the capilarry pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Swe <= Sthres) {
            Scalar pC_SweLow  = BrooksCorey::pC(params, Sthres);
            Scalar pC_SweLow2 = BrooksCorey::pC(params, Sthres/2);
            Scalar m = (pC_SweLow2 - pC_SweLow)/(Sthres/2 - Sthres);
            return pC_SweLow + m*(Swe - Sthres);
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
        const Scalar Sthres = params.pCLowSw();

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
            Scalar pC_SweLow  = BrooksCorey::pC(params, Sthres);
            Scalar pC_SweLow2 = BrooksCorey::pC(params, Sthres/2);
            Scalar m = (pC_SweLow2 - pC_SweLow)/(Sthres/2 - Sthres);
            return Sthres + (pC - pC_SweLow)/m;
        }

        return BrooksCorey::Sw(params, pC);
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
    */
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        const Scalar Sthres = params.pCLowSw();

        // derivative of the regualarization
        if (Swe <= Sthres) {
            // calculate the slope of the straight line used in pC()
            Scalar pC_SweLow  = BrooksCorey::pC(params, Sthres);
            Scalar pC_SweLow2 = BrooksCorey::pC(params, Sthres/2);
            Scalar m = (pC_SweLow2 - pC_SweLow)/(Sthres/2 - Sthres);
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
        const Scalar Sthres = params.pCLowSw();

        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized verision of the
        // Brooks-Corey law
        Scalar Swe;
        if (pC < 0)
            Swe = 1.5; // make sure we regularize below
        else
            Swe = BrooksCorey::Sw(params, pC);

        // derivative of the regularization
        if (Swe <= Sthres) {
            // same as in dpC_dSw() but inverted
            Scalar pC_SweLow  = BrooksCorey::pC(params, Sthres);
            Scalar pC_SweLow2 = BrooksCorey::pC(params, Sthres/2);
            Scalar m = (pC_SweLow2 - pC_SweLow)/(Sthres/2 - Sthres);
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
        const Scalar Sthres = params.krwHighSw();
        if (Sw < 0)
            return 0;
        else if (Sw >= 1)
            return 1;       
        // check if we need to regularize the relative permeability
        else if (Sw > Sthres) {
            typedef Dune::Spline<Scalar> Spline;
            Spline sp(Sthres, 1.0, // x1, x2
                      BrooksCorey::krw(params, Sthres), 1.0, // y1, y2
                      BrooksCorey::dkrw_dSw(params, Sthres), 0); // m1, m2

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
        const Scalar Sthres = params.krnLowSw();
        if (Sw <= 0)
            return 1;
        else if (Sw >= 1)
            return 0;
        // check if we need to regularize the relative permeability
        else if (Sw < Sthres) {
            typedef Dune::Spline<Scalar> Spline;
            Spline sp(0.0, Sthres, // x1, x2
                      1.0, BrooksCorey::krn(params, Sthres), // y1, y2
                      0.0, BrooksCorey::dkrn_dSw(params, Sthres)); // m1, m2
            
            return sp.eval(Sw);
        }
        
        return BrooksCorey::krn(params, Sw);
    }
};
}

#endif
