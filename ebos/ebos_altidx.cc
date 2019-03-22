// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief This is an ebos variant which uses alternative phase and component indices than
 *        the default variant.
 *
 * It is purely for testing purposes and is supposed to produce bitwise identical
 * results.
 */
#include "config.h"

#include "ebos.hh"

namespace Ewoms {
class EclAlternativeBlackOilIndexTraits
{
    typedef Opm::BlackOilDefaultIndexTraits DIT;

public:
    static const unsigned waterPhaseIdx = DIT::oilPhaseIdx;
    static const unsigned oilPhaseIdx = DIT::gasPhaseIdx;
    static const unsigned gasPhaseIdx = DIT::waterPhaseIdx;

    static const unsigned waterCompIdx = DIT::gasCompIdx;
    static const unsigned oilCompIdx = DIT::waterCompIdx;
    static const unsigned gasCompIdx = DIT::oilCompIdx;
};
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(EbosAltIdxTypeTag, INHERITS_FROM(EbosTypeTag));

// use a fluid system with different indices than the default
SET_PROP(EbosAltIdxTypeTag, FluidSystem)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::BlackOilFluidSystem<Scalar, Ewoms::EclAlternativeBlackOilIndexTraits> type;
};

END_PROPERTIES

int main(int argc, char **argv)
{
    typedef TTAG(EbosAltIdxTypeTag) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
