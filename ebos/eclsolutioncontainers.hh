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
 * \copydoc Opm::EclProblem
 */
#ifndef ECL_SOLUTION_CONTAINERS_HH
#define ECL_SOLUTION_CONTAINERS_HH

#include <vector>

namespace Opm {

//! \brief Struct holding polymer extension data.
template<class Scalar>
struct PolymerSolutionContainer {
    std::vector<Scalar> maxAdsorption;
    std::vector<Scalar> concentration;
    std::vector<Scalar> moleWeight; // polymer molecular weight

    static PolymerSolutionContainer serializationTestObject();

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(maxAdsorption);
        serializer(concentration);
        serializer(moleWeight);
    }

    bool operator==(const PolymerSolutionContainer& rhs) const;
};

//! \brief Struct holding MICP extension data.
template<class Scalar>
struct MICPSolutionContainer {
    std::vector<Scalar> microbialConcentration;
    std::vector<Scalar> oxygenConcentration;
    std::vector<Scalar> ureaConcentration;
    std::vector<Scalar> biofilmConcentration;
    std::vector<Scalar> calciteConcentration;

    static MICPSolutionContainer serializationTestObject();

    //! \brief Resize vectors and zero initialize.
    void resize(const unsigned numElems);

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(microbialConcentration);
        serializer(oxygenConcentration);
        serializer(ureaConcentration);
        serializer(biofilmConcentration);
        serializer(calciteConcentration);
    }

    bool operator==(const MICPSolutionContainer& rhs) const;
};

} // namespace Opm

#endif
