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

#ifndef OPM_BLACKOILWELLMODEL_NETWORK_GENERIC_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_NETWORK_GENERIC_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>

#include <opm/output/data/Groups.hpp>

#include <opm/simulators/flow/NewtonIterationContext.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <map>
#include <optional>
#include <string>

namespace Opm {
    class Schedule;
    class UnitSystem;
    template<class Scalar, class IndexTraits> class BlackoilWellModelGeneric;
    template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
    template<typename Scalar> class VFPProdProperties;
}

namespace Opm {

/// Class for handling the blackoil well network model.
template<typename Scalar, typename IndexTraits>
class BlackoilWellModelNetworkGeneric
{
public:
    BlackoilWellModelNetworkGeneric(BlackoilWellModelGeneric<Scalar,IndexTraits>& well_model);

    virtual ~BlackoilWellModelNetworkGeneric() = default;

    /// return true if network is active (at least one network well in prediction mode)
    bool active() const
    { return active_; }

    const std::map<std::string, Scalar>&
    nodePressures() const { return node_pressures_; }

    // do not use, only needed for serialization testing
    void setNodePressures(const std::map<std::string, Scalar>& values)
    { node_pressures_ = values; }

    void setFromRestart(const std::optional<std::map<std::string, double>>& restart_pressures);

    //! \brief Initialize wells according to network configuration.
    void initialize(const int report_step);

    //! \brief Initialize a single well according to network configuration.
    void initializeWell(WellInterfaceGeneric<Scalar,IndexTraits>& well);

    /// Checks if network is active (at least one network well on prediction).
    void updateActiveState(const int report_step);

    /// Checks if there are reasons to perform a pre-step network re-balance.
    /// (Currently, the only reasons are network well status changes.)
    /// (TODO: Consider if adding network change events would be helpful.)
    bool needPreStepRebalance(const int report_step) const;

    bool shouldBalance(const int reportStepIndex,
                       const NewtonIterationContext& iterCtx) const;

    Scalar updatePressures(const int reportStepIdx,
                           const Scalar damping_factor,
                           const Scalar update_upper_bound);

    void assignNodeValues(std::map<std::string, data::NodeData>& nodevalues,
                          const int reportStepIdx) const;

    void commitState()
    { this->last_valid_node_pressures_ = this->node_pressures_; }

    void resetState()
    { this->node_pressures_ = this->last_valid_node_pressures_; }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(node_pressures_);
        serializer(last_valid_node_pressures_);
    }

    bool operator==(const BlackoilWellModelNetworkGeneric<Scalar,IndexTraits>& rhs) const;

protected:
    std::map<std::string, Scalar>
    computePressures(const Network::ExtNetwork& network,
                     const VFPProdProperties<Scalar>& vfp_prod_props,
                     const UnitSystem& unit_system,
                     const int reportStepIdx,
                     const Parallel::Communication& comm) const;


    bool active_{false};
    BlackoilWellModelGeneric<Scalar,IndexTraits>& well_model_;

    // Network pressures for output and initialization
    std::map<std::string, Scalar> node_pressures_;
    // Valid network pressures for output and initialization for safe restart after failed iterations
    std::map<std::string, Scalar> last_valid_node_pressures_;
};

} // namespace Opm

#endif
