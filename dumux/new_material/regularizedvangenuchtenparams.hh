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
 * \file RegularizedVanGenuchtenParams.hh

 * \brief Reference implementation of a regularized van Genuchten params
 */
#ifndef REGULARIZED_VAN_GENUCHTEN_PARAMS_HH
#define REGULARIZED_VAN_GENUCHTEN_PARAMS_HH

#include <dumux/new_material/vangenuchten.hh>
#include <dumux/new_material/vangenuchtenparams.hh>

namespace Dune
{
/*!
 * \brief Reference implementation of a regularized van Genuchten params
 */
template<class ScalarT>
class RegularizedVanGenuchtenParams : public VanGenuchtenParams<ScalarT>
{
public:
    typedef ScalarT Scalar;
    typedef VanGenuchtenParams<Scalar> Parent;
    typedef RegularizedVanGenuchtenParams<Scalar> Self;

    RegularizedVanGenuchtenParams()
    {}

    RegularizedVanGenuchtenParams(Scalar vgAlpha,
                                   Scalar vgN)
        : Parent(vgAlpha, vgN)
    {};
};
}; // namespace Dune

#endif
