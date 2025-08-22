/*
  Copyright 2020 Equinor ASA.

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

#ifndef OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED
#define OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/GasLiftGroupInfo.hpp>

#include <optional>

namespace Opm {

template<class TypeTag> class WellInterface;

template<class TypeTag>
class GasLiftSingleWell : public GasLiftSingleWellGeneric<GetPropType<TypeTag, Properties::Scalar>,
                                         typename GetPropType<TypeTag, Properties::FluidSystem>::IndexTraitsType>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using GLiftSyncGroups = typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::GLiftSyncGroups;
    using BasicRates = typename GasLiftSingleWellGeneric<Scalar, IndexTraits>::BasicRates;

public:
    GasLiftSingleWell(const WellInterface<TypeTag>& well,
                      const Simulator& simulator,
                      const SummaryState& summary_state,
                      DeferredLogger& deferred_logger,
                      WellState<Scalar, IndexTraits>& well_state,
                      const GroupState<Scalar>& group_state,
                      GasLiftGroupInfo<Scalar, IndexTraits>& group_info,
                      GLiftSyncGroups& sync_groups,
                      const Parallel::Communication& comm,
                      bool glift_debug);

    const WellInterfaceGeneric<Scalar, IndexTraits>& getWell() const override { return well_; }

private:
    std::optional<Scalar>
    computeBhpAtThpLimit_(Scalar alq,
                          bool debug_ouput = true) const override;

    BasicRates computeWellRates_(Scalar bhp,
                                 bool bhp_is_limited,
                                 bool debug_output = true) const override;

    void setAlqMaxRate_(const GasLiftWell& well);
    void setupPhaseVariables_();
    bool checkThpControl_() const override;

    const Simulator& simulator_;
    const WellInterface<TypeTag>& well_;
};

} // namespace Opm

#include "GasLiftSingleWell_impl.hpp"

#endif // OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED
