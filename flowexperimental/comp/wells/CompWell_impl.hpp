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

#include <opm/material/fluidstates/CompositionalFluidState.hpp>

namespace Opm
{

template <typename TypeTag>
CompWell<TypeTag>::
CompWell(const Well& well,
         int index_of_well,
         const std::vector<CompConnectionData<Scalar>>& well_connection_data)
  : CompWellInterface<TypeTag>(well, index_of_well, well_connection_data)
{
}

template <typename TypeTag>
void
CompWell<TypeTag>::
init() {
    Base::init();
    // primary_variables_.init();
    well_equations_.init(this->number_of_connection_);
}

template <typename TypeTag>
void
CompWell<TypeTag>::
calculateExplitQuantities(const Simulator& simulator,
                          const SingleCompWellState<Scalar>& well_state)
{
    updatePrimaryVariables(simulator, well_state);
    updatePrimaryVariableEvaluation();

    // flash calculation in the wellbore
    using FluidState = CompositionalFluidState<Scalar, FluidSystem>;
    FluidState fluid_state = this->primary_variables_.toFluidStateScalar();
    PTFlash<Scalar, FluidSystem>::flash_solve_scalar_(fluid_state, "ssi", 1.e-6, 3);
    // calculating the mass within the wellbore

    // We should do a flash calculation here
    // assert(false && " calculateExplitQuantities is not implemented yet");
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updatePrimaryVariables(const Simulator& /* simulator */,
                       const SingleCompWellState<Scalar>& well_state)
{
    this->primary_variables_.update(well_state);
}

template <typename TypeTag>
void
CompWell<TypeTag>::
updatePrimaryVariableEvaluation()
{
    this->primary_variables_.updateEvaluation();
}



} // end of namespace Opm
