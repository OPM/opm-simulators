/*
  Copyright 2026 NORCE Research AS

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

#ifndef OPM_COMPOSITIONAL_MODEL_PARAMETERS_HPP
#define OPM_COMPOSITIONAL_MODEL_PARAMETERS_HPP

#include <opm/simulators/flow/BlackoilModelParameters.hpp>

namespace Opm {

/// Model parameters with compositional defaults for the solution change and residual tolerances
template <class Scalar>
struct CompositionalModelParameters : public BlackoilModelParameters<Scalar>
{
    /// Register the model parameters and set the compositional defaults
    static void registerParameters();
};

} // namespace Opm

#endif // OPM_COMPOSITIONAL_MODEL_PARAMETERS_HPP
