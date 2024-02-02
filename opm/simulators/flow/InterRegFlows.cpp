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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/simulators/flow/InterRegFlows.hpp>

#include <algorithm>
#include <stdexcept>

Opm::InterRegFlowMapSingleFIP::
InterRegFlowMapSingleFIP(const std::vector<int>& region)
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
Opm::InterRegFlowMapSingleFIP::
addConnection(const Cell& source,
              const Cell& destination,
              const data::InterRegFlowMap::FlowRates& rates)
{
    if (this->isReadFromStream_) {
        throw std::logic_error {
            "Cannot add new connection to deserialised object"
        };
    }

    if (! source.isInterior ||
        (source.cartesianIndex > destination.cartesianIndex))
    {
        // Connection handled in different call.  Don't double-count
        // contributions.
        return;
    }

    const auto r1 = this->region_[ source.activeIndex ];
    const auto r2 = this->region_[ destination.activeIndex ];

    if (r1 == r2) {
        // Connection is internal to a region.  Nothing to do.
        return;
    }

    // Inter-region connection internal to an MPI rank or this rank owns
    // the flow rate across this connection.
    this->iregFlow_.addConnection(r1, r2, rates);
}

void Opm::InterRegFlowMapSingleFIP::compress()
{
    this->iregFlow_.compress(this->maxGlobalRegionID_);
}

void Opm::InterRegFlowMapSingleFIP::clear()
{
    this->iregFlow_.clear();
    this->isReadFromStream_ = false;
}

const Opm::data::InterRegFlowMap&
Opm::InterRegFlowMapSingleFIP::getInterRegFlows() const
{
    return this->iregFlow_;
}

std::size_t Opm::InterRegFlowMapSingleFIP::getLocalMaxRegionID() const
{
    return this->maxLocalRegionID_;
}

bool
Opm::InterRegFlowMapSingleFIP::
assignGlobalMaxRegionID(const std::size_t regID)
{
    if (regID < this->maxLocalRegionID_) {
        return false;
    }

    if (regID > this->maxGlobalRegionID_) {
        this->maxGlobalRegionID_ = regID;
    }

    return true;
}

// =====================================================================
//
// Implementation of EclInterRegFlowMap (wrapper for multiple arrays)
//
// =====================================================================

Opm::InterRegFlowMap
Opm::InterRegFlowMap::createMapFromNames(std::vector<std::string> names)
{
    auto map = InterRegFlowMap{};

    map.names_ = std::move(names);
    map.regionMaps_.resize(map.names_.size(), InterRegFlowMapSingleFIP{});

    return map;
}

Opm::InterRegFlowMap::
InterRegFlowMap(const std::size_t                numCells,
                const std::vector<SingleRegion>& regions,
                const std::size_t                declaredMaxRegID)
{
    this->regionMaps_.reserve(regions.size());
    this->names_.reserve(regions.size());

    this->numCells_ = numCells;

    for (const auto& region : regions) {
        this->regionMaps_.emplace_back(region.definition);
        this->names_.push_back(region.name);
    }

    if (declaredMaxRegID > std::size_t{0}) {
        for (auto& regionMap : this->regionMaps_) {
            regionMap.assignGlobalMaxRegionID(declaredMaxRegID);
        }
    }
}

void
Opm::InterRegFlowMap::
addConnection(const Cell& source,
              const Cell& destination,
              const data::InterRegFlowMap::FlowRates& rates)
{
    for (auto& regionMap : this->regionMaps_) {
        regionMap.addConnection(source, destination, rates);
    }
}

void Opm::InterRegFlowMap::compress()
{
    for (auto& regionMap : this->regionMaps_) {
        regionMap.compress();
    }
}

void Opm::InterRegFlowMap::clear()
{
    for (auto& regionMap : this->regionMaps_) {
        regionMap.clear();
    }

    this->readIsConsistent_ = true;
}

const std::vector<std::string>&
Opm::InterRegFlowMap::names() const
{
    return this->names_;
}

std::vector<Opm::data::InterRegFlowMap>
Opm::InterRegFlowMap::getInterRegFlows() const
{
    auto maps = std::vector<data::InterRegFlowMap>{};
    maps.reserve(this->regionMaps_.size());

    for (auto& regionMap : this->regionMaps_) {
        maps.push_back(regionMap.getInterRegFlows());
    }

    return maps;
}

std::vector<std::size_t>
Opm::InterRegFlowMap::getLocalMaxRegionID() const
{
    auto maxLocalRegionID = std::vector<std::size_t>{};
    maxLocalRegionID.reserve(this->regionMaps_.size());

    for (auto& regionMap : this->regionMaps_) {
        maxLocalRegionID.push_back(regionMap.getLocalMaxRegionID());
    }

    return maxLocalRegionID;
}

bool
Opm::InterRegFlowMap::
assignGlobalMaxRegionID(const std::vector<std::size_t>& regID)
{
    if (regID.size() != this->regionMaps_.size()) {
        return false;
    }

    auto assignmentOK = true;

    const auto numReg = regID.size();
    for (auto region = 0*numReg; region < numReg; ++region) {
        const auto ok = this->regionMaps_[region]
            .assignGlobalMaxRegionID(regID[region]);

        assignmentOK = assignmentOK && ok;
    }

    return assignmentOK;
}

bool Opm::InterRegFlowMap::readIsConsistent() const
{
    return this->readIsConsistent_;
}
