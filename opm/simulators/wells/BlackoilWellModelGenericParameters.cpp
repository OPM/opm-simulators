/*
  Copyright 2026 Equinor ASA.

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

#include <config.h>
#include <opm/simulators/wells/BlackoilWellModelGenericParameters.hpp>

#include <opm/simulators/flow/BlackoilModelParameters.hpp>

#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

template<class Scalar>
BlackoilWellModelGenericParameters<Scalar>::BlackoilWellModelGenericParameters()
    : nupcol_group_rate_tolerance_(Parameters::Get<Parameters::NupcolGroupRateTolerance<Scalar>>())
    , max_number_of_group_switches_(Parameters::Get<Parameters::MaximumNumberOfGroupSwitches>())
    , use_multisegment_well_(Parameters::Get<Parameters::UseMultisegmentWell>())
{}

template struct BlackoilWellModelGenericParameters<double>;

#if FLOW_INSTANTIATE_FLOAT
template struct BlackoilWellModelGenericParameters<float>;
#endif

} // namespace Opm
