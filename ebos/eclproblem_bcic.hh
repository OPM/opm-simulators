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
#ifndef OPM_ECL_PROBLEM_BCIC_HH
#define OPM_ECL_PROBLEM_BCIC_HH

#include <ebos/eclequilinitializer.hh>

#include <vector>

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Handling of boundary- and initial conditions for EclProblem.
 */
template <class TypeTag>
class EclProblemBCIC
{
public:
    using InitialFluidState = typename EclEquilInitializer<TypeTag>::ScalarFluidState;

    //! \brief Returns whether or not boundary conditions are trivial.
    bool nonTrivialBoundaryConditions() const
    {
        return nonTrivialBoundaryConditions_;
    }

    //! \brief Returns a const reference to initial fluid state for an element.
    const InitialFluidState& initialFluidState(const unsigned idx) const
    {
        return initialFluidStates_[idx];
    }

    bool nonTrivialBoundaryConditions_ = false;
    std::vector<InitialFluidState> initialFluidStates_;
};

} // namespace Opm

#endif
