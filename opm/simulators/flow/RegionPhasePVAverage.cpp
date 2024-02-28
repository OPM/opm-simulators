// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 Equinor AS

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

#include <opm/simulators/flow/RegionPhasePVAverage.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FieldProps.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

namespace {
    std::vector<std::string> fipRegionNames(const std::vector<std::string>& regionNames)
    {
        auto regs = regionNames;

        std::sort(regs.begin(), regs.end());

        return { regs.begin(), std::unique(regs.begin(), regs.end()) };
    }

    std::vector<std::vector<double>::size_type>
    regionStartPointers(const std::vector<std::string>&                            regionNames,
                        const Opm::RegionPhasePoreVolAverage::RegionArrayAccessor& getRegionArray,
                        const Opm::Parallel::Communication&                        comm)
    {
        // All elements get an initial value of 1 to account for the maximum
        // region ID.  There should be "max ID + 1" elements for each region
        // set.  Unused IDs--e.g., 0 if IDs start at 1--will effectively be
        // wasted in this scheme.  If that becomes a problem we might
        // consider some kind of renumbering approach.
        auto start = std::vector<std::vector<double>::size_type>(regionNames.size() + 1, 1);

        const auto nset = regionNames.size();
        for (auto rset = 0*nset; rset < nset; ++rset) {
            const auto& reg = getRegionArray(regionNames[rset]);

            auto m = std::max_element(reg.begin(), reg.end());
            if (m == reg.end()) { // reg.empty()
                continue;
            }

            start[rset + 1] += *m;
        }

        comm.max(start.data(), start.size());

        std::partial_sum(start.begin(), start.end(), start.begin());

        return start;
    }
} // Anonymous namespace

Opm::RegionPhasePoreVolAverage::
RegionPhasePoreVolAverage(const Parallel::Communication&  comm,
                          const std::size_t               numPhases,
                          const std::vector<std::string>& regionNames,
                          RegionArrayAccessor             getRegionArray)
    : comm_           { std::cref(comm) }
    , np_             { numPhases }
    , rsetNames_      { fipRegionNames(regionNames) }
    , getRegionArray_ { std::move(getRegionArray) }
    , rsStart_        { regionStartPointers(rsetNames_, getRegionArray_, comm_) }
    , x_              (rsStart_.back() * numPhases * AvgType::NumTypes * Element::NumElem)
{}

double Opm::RegionPhasePoreVolAverage::fieldValue(const Phase& p) const
{
    return this->averageValueWithFallback(this->fieldStartIx(p.ix));
}

double
Opm::RegionPhasePoreVolAverage::
value(std::string_view rset, const Phase& p, const Region& r) const
{
    auto rsetPos = std::lower_bound(this->rsetNames_.begin(),
                                    this->rsetNames_.end(), rset);

    if ((rsetPos == this->rsetNames_.end()) || (*rsetPos != rset)) {
        // rset is not a known region set name.
        return 0.0;             // Maybe nullopt or throw here...
    }

    const auto rsetIx = std::distance(this->rsetNames_.begin(), rsetPos);
    return this->averageValueWithFallback(this->rsetStartIx(rsetIx, r.ix, p.ix));
}

void Opm::RegionPhasePoreVolAverage::prepareAccumulation()
{
    std::fill(this->x_.begin(), this->x_.end(), 0.0);
}

void Opm::RegionPhasePoreVolAverage::
addCell(const std::size_t activeCell,
        const Phase&      p,
        const CellValue&  cv)
{
    this->add(this->fieldStartIx(p.ix), cv);

    for (auto rset = 0*this->rsetNames_.size(); rset < this->rsetNames_.size(); ++rset) {
        this->add(this->rsetStartIx(rset, this->regionIndex(rset, activeCell), p.ix), cv);
    }
}

void Opm::RegionPhasePoreVolAverage::accumulateParallel()
{
    this->comm_.get().sum(this->x_.data(), this->x_.size());
}

double Opm::RegionPhasePoreVolAverage::averageValueWithFallback(const Ix start) const
{
    const auto spv_w = this->weight(start, AvgType::SatPV);

    return (spv_w > 0.0)
        ? this->averageValue(start, AvgType::SatPV)
        : this->averageValue(start, AvgType::PV);
}

double Opm::RegionPhasePoreVolAverage::averageValue(const Ix      start,
                                                    const AvgType type) const
{
    return this->value(start, type) / this->weight(start, type);
}

Opm::RegionPhasePoreVolAverage::Ix
Opm::RegionPhasePoreVolAverage::fieldStartIx(const unsigned int phase) const
{
    return this->startIx(0, phase);
}

Opm::RegionPhasePoreVolAverage::Ix
Opm::RegionPhasePoreVolAverage::rsetStartIx(const std::size_t  rset,
                                            const int          region,
                                            const unsigned int phase) const
{
    return this->startIx(this->rsStart_[rset] + region, phase);
}

Opm::RegionPhasePoreVolAverage::Ix
Opm::RegionPhasePoreVolAverage::startIx(const std::size_t offset,
                                        const unsigned int phase) const
{
    return (offset*this->np_ + phase) * AvgType::NumTypes * Element::NumElem;
}

int Opm::RegionPhasePoreVolAverage::regionIndex(const std::size_t rset,
                                                const std::size_t activeCell) const
{
    return this->getRegionArray_(this->rsetNames_[rset])[activeCell];
}

void Opm::RegionPhasePoreVolAverage::add(const Ix start, const CellValue& cv)
{
    this->add(start, AvgType::SatPV, cv.value, cv.sat * cv.porv);
    this->add(start, AvgType::PV   , cv.value,          cv.porv);
}

void Opm::RegionPhasePoreVolAverage::add(const Ix      start,
                                         const AvgType type,
                                         const double  x,
                                         const double  w)
{
    this->value (start, type) += w * x;
    this->weight(start, type) += w;
}

double& Opm::RegionPhasePoreVolAverage::value(const Ix start, const AvgType type)
{
    return this->x_[ this->valueArrayIndex(start, type, Element::Value) ];
}

double& Opm::RegionPhasePoreVolAverage::weight(const Ix start, const AvgType type)
{
    return this->x_[ this->valueArrayIndex(start, type, Element::Weight) ];
}

double Opm::RegionPhasePoreVolAverage::value(const Ix start, const AvgType type) const
{
    return this->x_[ this->valueArrayIndex(start, type, Element::Value) ];
}

double Opm::RegionPhasePoreVolAverage::weight(const Ix start, const AvgType type) const
{
    return this->x_[ this->valueArrayIndex(start, type, Element::Weight) ];
}

Opm::RegionPhasePoreVolAverage::Ix
Opm::RegionPhasePoreVolAverage::valueArrayIndex(const Ix      start,
                                                const AvgType type,
                                                const Element element) const
{
    return start + type*Element::NumElem + element;
}
