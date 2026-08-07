/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include <opm/simulators/utils/ParallelRegionVariableValues.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <memory>

Opm::ParallelRegionVariableValues::
ParallelRegionVariableValues(const Parallel::Communication& comm)
    : data::RegionVariableValues {}
    , comm_ { comm }
{}

std::unique_ptr<Opm::data::RegionVariableValues>
Opm::ParallelRegionVariableValues::clone() const
{
    return std::make_unique<ParallelRegionVariableValues>(*this);
}

void Opm::ParallelRegionVariableValues::communicateIncrement()
{
    this->comm_.sum(this->increment_.data(), this->increment_.size());
}
