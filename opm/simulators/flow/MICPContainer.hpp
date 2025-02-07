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
 * \copydoc Opm::OutputBlackOilModule
 */
#ifndef OPM_MICP_CONTAINER_HPP
#define OPM_MICP_CONTAINER_HPP

#include <vector>

namespace Opm {

namespace data { class Solution; }
template<class Scalar> class MICPSolutionContainer;

template<class Scalar>
class MICPContainer
{
    using ScalarBuffer = std::vector<Scalar>;

public:
    void allocate(const unsigned bufferSize);

    void assign(const unsigned globalDofIdx,
                 const Scalar microbialConcentration,
                 const Scalar oxygenConcentration,
                 const Scalar ureaConcentration,
                 const Scalar biofilmConcentration,
                 const Scalar calciteConcentration);

    MICPSolutionContainer<Scalar> getSolution() const;

    void outputRestart(data::Solution& sol);

    void readRestart(const unsigned globalDofIdx,
                     const unsigned elemIdx,
                     const data::Solution& sol);

    bool allocated() const
    { return allocated_; }

private:
    bool allocated_ = false;
    ScalarBuffer cMicrobes_;
    ScalarBuffer cOxygen_;
    ScalarBuffer cUrea_;
    ScalarBuffer cBiofilm_;
    ScalarBuffer cCalcite_;
};

} // namespace Opm

#endif // OPM_MICP_CONTAINER_HPP
