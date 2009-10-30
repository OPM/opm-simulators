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
 * \file RegularizedVanGenuchtenContext.hh

 * \brief Specification of the context API for the regualized van
 *        Genuchten capillary pressure model.
 */
#ifndef REGULARIZED_VAN_GENUCHTEN_CONTEXT_HH
#define REGULARIZED_VAN_GENUCHTEN_CONTEXT_HH

#include <dumux/auxiliary/apis.hh>
#include <dumux/new_material/statehelpermacros.hh>

#include <dumux/new_material/vangenuchten.hh>
#include <dumux/new_material/vangenuchtencontext.hh>

namespace Dune
{
/*!
 * \brief Reference implementation of a van Genuchten context
 */
template<class ScalarT>
class RegularizedVanGenuchtenContext
{
public:
    typedef ScalarT Scalar;
    typedef RegularizedVanGenuchtenContext<Scalar> Self;
    typedef Dune::VanGenuchten<Self> VanGenuchten;

    RegularizedVanGenuchtenContext()
    {}

    RegularizedVanGenuchtenContext(Scalar vgAlpha,
                                   Scalar vgN,
                                   Scalar minSw = 0.01)
    {
        setVgAlpha(vgAlpha);
        setVgN(vgN);
        setVgMinSw(minSw);
    };

    PROPERTY(Scalar, vgAlpha, setVgAlpha);

    // we also need to update vgM if vgN is changed (and vince versa),
    // so that we have to define the setter for Swr ourselfs
    PARAMETER(Scalar, vgM);
    void setVgM(Scalar vgM) { vgM_ = vgM; vgN_ = 1/(1 - vgM_); }
    PARAMETER(Scalar, vgN);
    void setVgN(Scalar vgN) { vgN_ = vgN; vgM_ = 1 - 1/vgN_; }

    PARAMETER(Scalar, vgMinSw);
    void setVgMinSw(Scalar vgMinSw)
    {
        vgMinSw_ = vgMinSw;
        vgMaxPC_ = VanGenuchten::pC(*this, vgMinSw_);
    }

    PARAMETER(Scalar, vgMaxPC);
    void setVgMaxPC(Scalar vgMaxPC)
    {
        vgMaxPC_ = vgMaxPC;
        vgMinSw_ = VanGenuchten::Sw(*this, vgMaxPC_);
    }
};
}; // namespace Dune

#endif
