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

#include <opm/simulators/utils/ParallelRegionsetVariableDescriptor.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <cassert>
#include <memory>

Opm::ParallelRegionsetVariableDescriptor::
ParallelRegionsetVariableDescriptor(const Parallel::Communication& comm)
    : data::RegionsetVariableDescriptor {}
    , comm_ { comm }
{}

std::unique_ptr<Opm::data::RegionsetVariableDescriptor>
Opm::ParallelRegionsetVariableDescriptor::clone() const
{
    return std::make_unique<ParallelRegionsetVariableDescriptor>(*this);
}

void Opm::ParallelRegionsetVariableDescriptor::communicateGlobalRegsetMaxIDs()
{
    assert (this->regsetMaxID_.has_value());

    this->comm_.max(this->regsetMaxID_->data(),
                    this->regsetMaxID_->size());
}
