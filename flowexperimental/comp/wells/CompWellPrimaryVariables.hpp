/*
  Copyright 2024, SINTEF Digital

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_COMP_WELL_PRIMARY_VARIABLES_HPP
#define OPM_COMP_WELL_PRIMARY_VARIABLES_HPP

#include <opm/material/densead/Evaluation.hpp>

#include <vector>

namespace Opm {

template <typename FluidSystem>
class CompWellPrimaryVariables {
public:
    static constexpr int numWellConservationEq = FluidSystem::numComponents;
    static constexpr int numWellControlEq = 1;
    static constexpr int numWellEq = numWellConservationEq + numWellControlEq;
    // the indices for the primary variables
    // the first primary variable will be the total mass rate
    // the last primary variable will be the BHP
    // the one in the middle with will the mole fractions for the numWellEq - 1 components
    static constexpr int QTotalMass = 0;
    static constexpr int Bhp = numWellEq - numWellControlEq;

    using Scalar = typename FluidSystem::Scalar;
    // we use hard-coded Evaluation type for now
    // TODO: we can try to use DyanmicEvaluation here
    using EvalWell = DenseAd::Evaluation<Scalar, numWellEq>;

private:
    std::vector<Scalar> value_;
    std::vector<EvalWell> evaluation_;
};

} // end of namespace Opm

#endif // OPM_COMP_WELL_PRIMARY_VARIABLES_HPP
