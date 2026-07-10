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

#ifndef OPM_BLACKOILWELLMODEL_GENERIC_PARAMETERS_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_GENERIC_PARAMETERS_HEADER_INCLUDED

namespace Opm {

/// Parameter bundle consumed by BlackoilWellModelGeneric.
///
/// Deliberately narrow subset of BlackoilModelParameters: holds only
/// the fields the generic well-model base reads, so the base does not
/// depend on parameters specific to a particular reservoir model
/// (e.g. NLDD-solver knobs).
///
/// The default constructor populates fields from the global parameter
/// system, which must already have been initialised (typically via
/// BlackoilModelParameters::registerParameters()).
template<class Scalar>
struct BlackoilWellModelGenericParameters
{
    BlackoilWellModelGenericParameters();

    Scalar nupcol_group_rate_tolerance_;
    int max_number_of_group_switches_;
    bool use_multisegment_well_;
    bool enable_group_tree_balancer_;
};

} // namespace Opm

#endif // OPM_BLACKOILWELLMODEL_GENERIC_PARAMETERS_HEADER_INCLUDED
