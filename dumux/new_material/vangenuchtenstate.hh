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
 * \file VanGenuchtenState.hh Specification of the state API for the
 *                            van Genuchten capillary pressure model.
 */
#ifndef VAN_GENUCHTEN_STATE_HH
#define VAN_GENUCHTEN_STATE_HH

#include <dumux/auxiliary/apis.hh>
#include <dumux/new_material/statehelpermacros.hh>

namespace Dune
{
namespace Api
{
    BEGIN_API_DEF(VanGenuchtenParams)
    {
        typedef typename Implementation::Scalar Scalar;
        Scalar tmp;
        tmp = const_impl.vgAlpha();
        tmp = const_impl.vgN();
        tmp = const_impl.vgM();
    }
    END_API_DEF;

    BEGIN_API_DEF(VanGenuchtenState)
    {
        typedef typename Implementation::Scalar Scalar;
        Scalar tmp = 0.5f;
        impl.setVgAlpha(tmp);
        impl.setVgN(tmp);
        impl.setVgM(tmp);
    }
    END_API_DEF;
}; // namespace Api

    /*!
     * \brief Reference implementation of a van Genuchten state
     */
    template<class ScalarT>
    class VanGenuchtenState
    {
    public:
        typedef ScalarT Scalar;

        VanGenuchtenState()
            {}

        VanGenuchtenState(Scalar vgAlpha, Scalar vgN)
            {
                setVgAlpha(vgAlpha);
                setVgN(vgN);
            };

        PROPERTY(Scalar, vgAlpha, setVgAlpha);

        // we also need to update vgM if vgN is changed (and vince versa),
        // so that we have to define the setter for Swr ourselfs
        PARAMETER(Scalar, vgM);
        void setVgM(Scalar vgM) { _vgM = vgM; _vgN = 1/(1 - _vgM); }
        PARAMETER(Scalar, vgN);
        void setVgN(Scalar vgN) { _vgN = vgN; _vgM = 1 - 1/_vgN; }
    };
}; // namespace Dune

#endif
