/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on absolute 
 *        saturations.
 */
#ifndef ABSOLUTE_SATURATIONS_LAW_HH
#define ABSOLUTE_SATURATIONS_LAW_HH

#include <dumux/new_material/absolutesaturationslawcontext.hh>

namespace Dune
{
/*!
 * \ingroup material
 *
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on absolute 
 *        saturations.
 */
template <class RawLawT>
class AbsoluteSaturationsLaw
{
    typedef RawLawT   RawLaw;

public:
    typedef typename RawLaw::Context  Context;
    typedef typename RawLaw::Scalar   Scalar;


    /*!
     * \brief The capillary pressure-saturation curve.
     */
    static Scalar pC(const Context &context, Scalar Sw)
    {
        return RawLaw::pC(context, SwToSwe_(context, Sw));
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * \param pC Capillary pressure \f$\p_C\f$
     * \return The absolute saturation of the wetting phase \f$S_w\f$
     */
    static Scalar Sw(const Context &context, Scalar pC)
    {
        return SweToSw_(context, RawLaw::Sw(context, pC));
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the absolute saturation.
    */
    static Scalar dpC_dSw(const Context &context, Scalar Sw)
    {
        return RawLaw::dpC_dSw(context, pC)/(1 - context.Swr() - context.Snr());
    }

    /*!
     * \brief Returns the partial derivative of the absolute
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Context &context, Scalar pC)
    {
        return RawLaw::dSw_dpC(context, pC)*(1 - context.Swr() - context.Snr());
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param Sw The absolute saturation of the wetting phase.
     */
    static Scalar krw(const Context &context, Scalar Sw)
    {
        return RawLaw::krw(context, SwToSwe_(context, Sw));
    };

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param Sw The absolute saturation of the wetting phase.
     */
    static Scalar krn(const Context &context, Scalar Sw)
    {
        return RawLaw::krn(context, SwToSwe_(context, Sw));
    }

private:
    // convert an absolute wetting saturation to an effective one
    static Scalar SwToSwe_(const Context &context, Scalar Sw)
    {
        return (Sw - context.Swr())/(1 - context.Swr() - context.Snr());
    }

    // convert an effective wetting saturation to an absolute one
    static Scalar SweToSw_(const Context &context, Scalar Swe)
    {
        return Swe*(1 - context.Swr() - context.Snr()) + context.Swr();
    }
};
}

#endif
