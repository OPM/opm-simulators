/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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

#ifndef OPM_BLACKOILWELLMODEL_NETWORK_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_NETWORK_HEADER_INCLUDED

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/utils/basicproperties.hh>

#include <opm/simulators/wells/BlackoilWellModelNetworkGeneric.hpp>

#include <map>
#include <string>
#include <tuple>

namespace Opm {
    class DeferredLogger;
    template<class TypeTag> class BlackoilWellModel;
}

namespace Opm {

/// Class for handling the blackoil well network model.
template<typename TypeTag>
class BlackoilWellModelNetwork :
    public BlackoilWellModelNetworkGeneric<GetPropType<TypeTag, Properties::Scalar>,
                                           typename GetPropType<TypeTag, Properties::FluidSystem>::IndexTraitsType>
{
    using BaseType =
        BlackoilWellModelNetworkGeneric<GetPropType<TypeTag, Properties::Scalar>,
                                        typename GetPropType<TypeTag, Properties::FluidSystem>::IndexTraitsType>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    BlackoilWellModelNetwork(BlackoilWellModel<TypeTag>& well_model);

    std::tuple<bool, Scalar>
    update(const bool mandatory_network_balance,
           DeferredLogger& deferred_logger,
           const bool relax_network_tolerance = false);

    // Pre-step network solve at static reservoir conditions (group and well states might be updated)
    void doPreStepRebalance(DeferredLogger& deferred_logger);

protected:
    /// This function is to be used for well groups in an extended network that act as a subsea manifold
    /// The wells of such group should have a common THP and total phase rate(s) obeying (if possible)
    /// the well group constraint set by GCONPROD
    bool computeWellGroupThp(const double dt, DeferredLogger& local_deferredLogger);

    BlackoilWellModel<TypeTag>& well_model_;
    std::map<std::string, Scalar> well_group_thp_calc_;
};

} // namespace Opm

#include "BlackoilWellModelNetwork_impl.hpp"

#endif // OPM_BLACKOILWELLMODEL_NETWORK_HEADER_INCLUDED
