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

#include "absolutesaturationslawparams.hh"

namespace Dumux
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
    typedef typename RawLaw::Params   Params;
    typedef typename RawLaw::Scalar   Scalar;


    /*!
     * \brief The capillary pressure-saturation curve.
     */
    static Scalar pC(const Params &params, Scalar Sw)
    {
        return RawLaw::pC(params, SwToSwe_(params, Sw));
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * \param pC Capillary pressure \f$\p_C\f$
     * \return The absolute saturation of the wetting phase \f$S_w\f$
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        return SweToSw_(params, RawLaw::Sw(params, pC));
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the absolute saturation.
    */
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        return RawLaw::dpC_dSw(params, Sw)/(1 - params.Swr() - params.Snr());
    }

    /*!
     * \brief Returns the partial derivative of the absolute
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        return RawLaw::dSw_dpC(params, pC)*(1 - params.Swr() - params.Snr());
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param Sw The absolute saturation of the wetting phase.
     */
    static Scalar krw(const Params &params, Scalar Sw)
    {
        return RawLaw::krw(params, SwToSwe_(params, Sw));
    };

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param Sw The absolute saturation of the wetting phase.
     */
    static Scalar krn(const Params &params, Scalar Sw)
    {
        return RawLaw::krn(params, SwToSwe_(params, Sw));
    }

private:
    // convert an absolute wetting saturation to an effective one
    static Scalar SwToSwe_(const Params &params, Scalar Sw)
    {
        return (Sw - params.Swr())/(1 - params.Swr() - params.Snr());
    }

    // convert an effective wetting saturation to an absolute one
    static Scalar SweToSw_(const Params &params, Scalar Swe)
    {
        return Swe*(1 - params.Swr() - params.Snr()) + params.Swr();
    }
};
}

#endif
