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

#include <dumux/new_material/brookscorey.hh>
#include <dumux/new_material/regularizedbrookscoreycontext.hh>

#include <algorithm>

#include <math.h>
#include <assert.h>

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
template <class ContextT>
class RegularizedBrooksCorey
{
public:
    typedef ContextT Context;
    typedef typename Context::Scalar Scalar;

    typedef Dune::BrooksCorey<ContextT> BrooksCorey;

    /*!
     * \brief A regularized Brooks-Corey capillary pressure-saturation
     *        curve.
     */
    static Scalar pC(const Context &context, Scalar Swe)
    {
        // make sure that the capilarry pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Swe <= SweLow_) {
            Scalar pC_SweLow  = BrooksCorey::pC(context, SweLow_);
            Scalar pC_SweLow2 = BrooksCorey::pC(context, SweLow_/2);
            Scalar m = (pC_SweLow2 - pC_SweLow)/(SweLow_/2 - SweLow_);
            return pC_SweLow + m*(Swe - SweLow_);
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real Brooks-Corey law...
        return BrooksCorey::pC(context, Swe);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     */
    static Scalar Sw(const Context &context, Scalar pC)
    {
        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized version of 
        // the Brooks-Corey law
        Scalar Swe = BrooksCorey::Sw(context, pC);

        // make sure that the capilarry pressure observes a
        // derivative != 0 for 'illegal' saturations. This is
        // required for example by newton solvers (if the
        // derivative calculated numerically) in order to get the
        // saturation moving to the right direction if it
        // temporarily is in an 'illegal' range.
        if (Swe <= SweLow_) {
            // invert the low saturation regularization of pC()
            Scalar pC_SweLow  = BrooksCorey::pC(context, SweLow_);
            Scalar pC_SweLow2 = BrooksCorey::pC(context, SweLow_/2);
            Scalar m = (pC_SweLow2 - pC_SweLow)/(SweLow_/2 - SweLow_);
            return SweLow_ + (pC - pC_SweLow)/m;
        }

        return BrooksCorey::Sw(context, pC);
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
    */
    static Scalar dpC_dSw(const Context &context, Scalar Swe)
    {
        // derivative of the regualarization
        if (Swe <= SweLow_) {
            // calculate the slope of the straight line used in pC()
            Scalar pC_SweLow  = BrooksCorey::pC(context, SweLow_);
            Scalar pC_SweLow2 = BrooksCorey::pC(context, SweLow_/2);
            Scalar m = (pC_SweLow2 - pC_SweLow)/(SweLow_/2 - SweLow_);
            return m;
        }

        return BrooksCorey::dpC_dSw(context, Swe);
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Context &context, Scalar pC)
    {
        // calculate the saturation which corrosponds to the
        // saturation in the non-regularized verision of the
        // Brooks-Corey law
        Scalar Swe;
        if (pC < 0)
            Swe = 1.5; // make sure we regularize below
        else
            Swe = BrooksCorey::Sw(context, pC);

        // derivative of the regularization
        if (Swe <= SweLow_) {
            // same as in dpC_dSw() but inverted
            Scalar pC_SweLow  = BrooksCorey::pC(context, SweLow_);
            Scalar pC_SweLow2 = BrooksCorey::pC(context, SweLow_/2);
            Scalar m = (pC_SweLow2 - pC_SweLow)/(SweLow_/2 - SweLow_);
            return 1/m;
        }

        return BrooksCorey::dpC_dSw(context, pC);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     */
    static Scalar krw(const Context &context, Scalar Sw_mob)
    {
        if (Sw_mob < 0)
            return 0;
        else if (Sw_mob > 1)
            return 1;

        return BrooksCorey::krw(context, Sw_mob);
    };

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     */
    static Scalar krn(const Context &context, Scalar Sw_mob)
    {
        if (Sw_mob <= 0)
            return 1;
        else if (Sw_mob >= 1)
            return 0;

        return BrooksCorey::krn(context, Sw_mob);
    }

    //! Effective saturation below which we regularize
    static const Scalar SweLow_;
};

template <class ContextT>
const typename ContextT::Scalar RegularizedBrooksCorey<ContextT>::SweLow_(0.05);
}

#endif
