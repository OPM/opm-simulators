// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 Equinor AS

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

#include <ebos/eclinterregflows.hh>

#include <algorithm>
#include <stdexcept>

Opm::EclInterRegFlowMapSingleFIP::
EclInterRegFlowMapSingleFIP(const std::vector<int>& region)
    : region_(region.size(), 0)
{
    if (! region.empty()) {
        this->maxLocalRegionID_= this->maxGlobalRegionID_ =
            *std::max_element(region.begin(), region.end());
    }

    std::transform(region.begin(), region.end(),
                   this->region_.begin(),
                   [](const int regID) { return regID - 1; });
}

void
Opm::EclInterRegFlowMapSingleFIP::
addConnection(const Cell& source,
              const Cell& destination,
              const data::InterRegFlowMap::FlowRates& rates)
{
    if (this->isReadFromStream_) {
        throw std::logic_error {
            "Cannot add new connection to deserialised object"
        };
    }

    if (! (source.isInterior || destination.isInterior)) {
        // Connection between two cells not on this process.  Unlikely, but
        // nothing to do here.
        return;
    }

    const auto r1 = this->region_[ source.activeIndex ];
    const auto r2 = this->region_[ destination.activeIndex ];

    if (r1 == r2) {
        // Connection is internal to a region.  Nothing to do.
        return;
    }

    if ((source.isInterior && destination.isInterior) ||
        (source.isInterior && (r1 < r2)) ||
        (destination.isInterior && (r2 < r1)))
    {
        // Inter-region connection internal to an MPI rank or this rank owns
        // the flow rate across this connection.
        this->iregFlow_.addConnection(r1, r2, rates);
    }
}

void Opm::EclInterRegFlowMapSingleFIP::compress()
{
    this->iregFlow_.compress(this->maxGlobalRegionID_);
}

void Opm::EclInterRegFlowMapSingleFIP::clear()
{
    this->iregFlow_.clear();
    this->isReadFromStream_ = false;
}

const Opm::data::InterRegFlowMap&
Opm::EclInterRegFlowMapSingleFIP::getInterRegFlows() const
{
    return this->iregFlow_;
}

std::size_t Opm::EclInterRegFlowMapSingleFIP::getLocalMaxRegionID() const
{
    return this->maxLocalRegionID_;
}

bool
Opm::EclInterRegFlowMapSingleFIP::
assignGlobalMaxRegionID(const std::size_t regID)
{
    if (regID < this->maxLocalRegionID_) {
        return false;
    }

    this->maxGlobalRegionID_ = regID;
    return true;
}
