// $Id$
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
 *        saturations and converts it to a material law defined on
 *        absolute saturations.
 */
#ifndef DUMUX_EFF_TO_ABS_LAW_HH
#define DUMUX_EFF_TO_ABS_LAW_HH

#include "efftoabslawparams.hh"

namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on absolute
 *        saturations.
 */
template <class EffLawT, class AbsParamsT = EffToAbsLawParams<typename EffLawT::Params> >
class EffToAbsLaw
{
    typedef EffLawT EffLaw;

public:
    typedef AbsParamsT Params;
    typedef typename EffLaw::Scalar Scalar;


    /*!
     * \brief The capillary pressure-saturation curve.
     */
    static Scalar pC(const Params &params, Scalar Sw)
    {
        return EffLaw::pC(params, SwToSwe(params, Sw));
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * \param pC Capillary pressure \f$p_C\f$
     * \return The absolute saturation of the wetting phase \f$S_w\f$
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        return SweToSw_(params, EffLaw::Sw(params, pC));
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the absolute saturation.
    */
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        return EffLaw::dpC_dSw(params, Sw)*dSwe_dSw_(params);
    }

    /*!
     * \brief Returns the partial derivative of the absolute
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        return EffLaw::dSw_dpC(params, pC)*dSw_dSwe_(params);
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param Sw The absolute saturation of the wetting phase.
     */
    static Scalar krw(const Params &params, Scalar Sw)
    {
        return EffLaw::krw(params, SwToSwe(params, Sw));
    };

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param Sw The absolute saturation of the wetting phase.
     */
    static Scalar krn(const Params &params, Scalar Sw)
    {
        return EffLaw::krn(params, SwToSwe(params, Sw));
    }

    // convert an absolute wetting saturation to an effective one
    static Scalar SwToSwe(const Params &params, Scalar Sw)
    {
        return (Sw - params.Swr())/(1 - params.Swr() - params.Snr());
    }

    // convert an absolute wetting saturation to an effective one
    static Scalar SnToSne(const Params &params, Scalar Sn)
    {
        return (Sn - params.Snr())/(1 - params.Swr() - params.Snr());
    }

private:
    // convert an effective wetting saturation to an absolute one
    static Scalar SweToSw_(const Params &params, Scalar Swe)
    {
        return Swe*(1 - params.Swr() - params.Snr()) + params.Swr();
    }

    // derivative of the effective saturation to the absolute
    // saturation.
    static Scalar dSwe_dSw_(const Params &params)
    { return 1.0/(1 - params.Swr() - params.Snr()); }

    // derivative of the absolute saturation to the effective
    // saturation.
    static Scalar dSw_dSwe_(const Params &params)
    { return 1 - params.Swr() - params.Snr(); }
};
}

#endif
