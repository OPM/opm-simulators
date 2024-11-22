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

#ifndef OPM_COMP_WELL_HPP
#define OPM_COMP_WELL_HPP

#include <opm/models/utils/propertysystem.hh>

#include "CompWellEquations.hpp"
#include "CompWellInterface.hpp"
#include "CompWellPrimaryVariables.hpp"

#include <string>

namespace Opm
{

template <typename TypeTag>
class CompWell  : public CompWellInterface<TypeTag>
{
public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = CompWellPrimaryVariables<FluidSystem>;
    using WellEquations = CompWellEquations<Scalar, PrimaryVariables::numWellEq, Indices::numEq>;

    CompWell(const Well& well,
             int index_of_well);
private:

    // primary variables
    PrimaryVariables primary_variables_;
    WellEquations well_equations_;

};

} // end of namespace Opm

#include "CompWell_impl.hpp"

#endif // OPM_COMPOSITIONAL_WELL_MODEL_HPP