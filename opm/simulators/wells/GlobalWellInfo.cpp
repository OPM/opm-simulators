/*
  Copyright 2021 Equinor

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
#include <opm/simulators/wells/GlobalWellInfo.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <stdexcept>

namespace Opm {



GlobalWellInfo::GlobalWellInfo(const Schedule& sched, std::size_t report_step, const std::vector<Well>& local_wells)
{
    auto num_wells = sched.numWells(report_step);
    this->m_in_injecting_group.resize(num_wells);
    this->m_in_producing_group.resize(num_wells);
    for (const auto& wname : sched.wellNames(report_step)) {
        const auto& well = sched.getWell(wname, report_step);
        auto global_well_index = well.seqIndex();
        this->name_map.emplace( well.name(), global_well_index );
    }

    for (const auto& well : local_wells)
        this->local_map.push_back( well.seqIndex() );
}


bool GlobalWellInfo::in_injecting_group(const std::string& wname) const {
    auto global_well_index = this->name_map.at(wname);
    return this->m_in_injecting_group[global_well_index];
}


bool GlobalWellInfo::in_producing_group(const std::string& wname) const {
    auto global_well_index = this->name_map.at(wname);
    return this->m_in_producing_group[global_well_index];
}


void GlobalWellInfo::update_injector(std::size_t well_index, Well::Status well_status, Well::InjectorCMode injection_cmode) {
    if (well_status == Well::Status::OPEN && injection_cmode == Well::InjectorCMode::GRUP)
        this->m_in_injecting_group[this->local_map[well_index]] = 1;
}

void GlobalWellInfo::update_producer(std::size_t well_index, Well::Status well_status, Well::ProducerCMode production_cmode) {
    if (well_status == Well::Status::OPEN && production_cmode == Well::ProducerCMode::GRUP)
        this->m_in_producing_group[this->local_map[well_index]] = 1;
}

void GlobalWellInfo::clear() {
    this->m_in_injecting_group.assign(this->name_map.size(), 0);
    this->m_in_producing_group.assign(this->name_map.size(), 0);
}


std::size_t GlobalWellInfo::well_index(const std::string& wname) const {
    return this->name_map.at(wname);
}

const std::string& GlobalWellInfo::well_name(std::size_t well_index) const {
    for (const auto& [name, index] : this->name_map) {
        if (index == well_index)
            return name;
    }
    throw std::logic_error("No well with index: " + std::to_string(well_index));
}
}
