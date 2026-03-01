/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_RESCOUP_RECEIVE_GROUP_TARGETS_HPP
#define OPM_RESCOUP_RECEIVE_GROUP_TARGETS_HPP
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlave.hpp>
#include <opm/simulators/wells/GroupStateHelper.hpp>
#include <opm/simulators/wells/GuideRateHandler.hpp>

namespace Opm {

template<class Scalar, class IndexTraits>
class RescoupReceiveGroupTargets {
public:
    RescoupReceiveGroupTargets(
        GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler,
        GroupStateHelper<Scalar, IndexTraits>& group_state_helper
    );
    void receiveGroupTargetsFromMaster();
private:
    GuideRateHandler<Scalar, IndexTraits>& guide_rate_handler_;
    GroupStateHelper<Scalar, IndexTraits>& group_state_helper_;
    ReservoirCouplingSlave<Scalar>& reservoir_coupling_slave_;
};

} // namespace Opm

#endif // OPM_RESCOUP_RECEIVE_GROUP_TARGETS_HPP
