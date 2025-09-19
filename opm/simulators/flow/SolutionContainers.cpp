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
 * \copydoc Opm::FlowProblem
 */
#include <config.h>
#include <opm/simulators/flow/SolutionContainers.hpp>

namespace Opm {

template<class Scalar>
PolymerSolutionContainer<Scalar>
PolymerSolutionContainer<Scalar>::serializationTestObject()
{
    return {{3.0, 4.0, 5.0},
            {12.0},
            {13.0, 14.0}};
}

template<class Scalar>
bool PolymerSolutionContainer<Scalar>::
operator==(const PolymerSolutionContainer<Scalar>& rhs) const
{
    return this->maxAdsorption == rhs.maxAdsorption &&
           this->concentration == rhs.concentration &&
           this->moleWeight == rhs.moleWeight;
}

template<class Scalar>
BioeffectsSolutionContainer<Scalar>
BioeffectsSolutionContainer<Scalar>::serializationTestObject()
{
    return {{16.0},
            {17.0},
            {18.0},
            {19.0},
            {20.0}};
}

template<class Scalar>
void BioeffectsSolutionContainer<Scalar>::resize(const unsigned numElems)
{
    microbialConcentration.resize(numElems, 0.0);
    oxygenConcentration.resize(numElems, 0.0);
    ureaConcentration.resize(numElems, 0.0);
    biofilmVolumeFraction.resize(numElems, 0.0);
    calciteVolumeFraction.resize(numElems, 0.0);
}

template<class Scalar>
bool BioeffectsSolutionContainer<Scalar>::
operator==(const BioeffectsSolutionContainer<Scalar>& rhs) const
{
    return this->microbialConcentration == rhs.microbialConcentration &&
           this->oxygenConcentration == rhs.oxygenConcentration &&
           this->ureaConcentration == rhs.ureaConcentration &&
           this->biofilmVolumeFraction == rhs.biofilmVolumeFraction &&
           this->calciteVolumeFraction == rhs.calciteVolumeFraction;
}

#define INSTANTIATE_TYPE(T) \
    template struct PolymerSolutionContainer<T>; \
    template struct BioeffectsSolutionContainer<T>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
