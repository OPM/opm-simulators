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
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/GasLiftGroupInfo.hpp>
#include <opm/simulators/wells/WellInterface.hpp>

#include <optional>
#include <vector>
#include <utility>

namespace Opm
{
    template<class TypeTag>
    class GasLiftSingleWell : public GasLiftSingleWellGeneric
    {
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using GLiftSyncGroups = typename GasLiftSingleWellGeneric::GLiftSyncGroups;

    public:
        GasLiftSingleWell(
            const WellInterface<TypeTag> &well,
            const Simulator &ebos_simulator,
            const SummaryState &summary_state,
            DeferredLogger &deferred_logger,
            WellState &well_state,
            const GroupState& group_state,
            GasLiftGroupInfo &group_info,
            GLiftSyncGroups &sync_groups,
            const Parallel::Communication& comm,
            bool glift_debug
        );
        const WellInterfaceGeneric &getWell() const override { return well_; }

    private:
        std::optional<double> computeBhpAtThpLimit_(double alq, bool debug_ouput=true) const override;
        BasicRates computeWellRates_(
            double bhp, bool bhp_is_limited, bool debug_output=true) const override;
        void setAlqMaxRate_(const GasLiftWell& well);
        void setupPhaseVariables_();
        bool checkThpControl_() const override;


        const Simulator &ebos_simulator_;
        const WellInterface<TypeTag> &well_;
    };

} // namespace Opm

#include "GasLiftSingleWell_impl.hpp"

#endif // OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED
