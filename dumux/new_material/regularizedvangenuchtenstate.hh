/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
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
 * \file RegularizedVanGenuchtenState.hh

 * \brief Specification of the state API for the regualized van
 *        Genuchten capillary pressure model.
 */
#ifndef REGULARIZED_VAN_GENUCHTEN_STATE_HH
#define REGULARIZED_VAN_GENUCHTEN_STATE_HH

#include <dumux/auxiliary/apis.hh>
#include <dumux/new_material/statehelpermacros.hh>

#include <dumux/new_material/vangenuchten.hh>
#include <dumux/new_material/vangenuchtenstate.hh>

namespace Dune
{
namespace Api
{
    BEGIN_API_DEF(RegularizedVanGenuchtenParams)
    {
        require<VanGenuchtenParams>(impl);

        typedef typename Implementation::Scalar Scalar;
        Scalar tmp = 0.5;
        tmp = const_impl.vgMaxPC();
        tmp = const_impl.vgMinSw();
    }
    END_API_DEF;

    BEGIN_API_DEF(RegularizedVanGenuchtenState)
    {
        require<RegularizedVanGenuchtenParams>(impl);
        require<VanGenuchtenState>(impl);

        typedef typename Implementation::Scalar Scalar;
        Scalar tmp = 0.5f;
        impl.setVgMaxPC(tmp);
        impl.setVgMinSw(tmp);
    }
    END_API_DEF;
}; // namespace Api

    /*!
     * \brief Reference implementation of a van Genuchten state
     */
    template<class ScalarT>
    class RegularizedVanGenuchtenState
    {
    public:
        typedef ScalarT Scalar;
        typedef RegularizedVanGenuchtenState<Scalar> Self;
        typedef Dune::VanGenuchten<Self> VanGenuchten;

        RegularizedVanGenuchtenState()
            {}

        RegularizedVanGenuchtenState(Scalar vgAlpha,
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
        void setVgM(Scalar vgM) { _vgM = vgM; _vgN = 1/(1 - _vgM); }
        PARAMETER(Scalar, vgN);
        void setVgN(Scalar vgN) { _vgN = vgN; _vgM = 1 - 1/_vgN; }

        PARAMETER(Scalar, vgMinSw);
        void setVgMinSw(Scalar vgMinSw)
        {
            _vgMinSw = vgMinSw;
            _vgMaxPC = VanGenuchten::pC(*this, _vgMinSw);
        }

        PARAMETER(Scalar, vgMaxPC);
        void setVgMaxPC(Scalar vgMaxPC)
        {
            _vgMaxPC = vgMaxPC;
            _vgMinSw = VanGenuchten::Sw(*this, _vgMaxPC);
        }
    };
}; // namespace Dune

#endif
