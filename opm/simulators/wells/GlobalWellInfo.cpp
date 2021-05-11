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
#include <stdexcept>

#include <opm/simulators/wells/GlobalWellInfo.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>


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

    for (const auto& well : local_wells) {
        this->local_map.push_back( well.seqIndex() );
        this->is_injector.push_back( well.isInjector() );
    }
}


bool GlobalWellInfo::in_injecting_group(const std::string& wname) const {
    auto global_well_index = this->name_map.at(wname);
    return this->m_in_injecting_group[global_well_index];
}


bool GlobalWellInfo::in_producing_group(const std::string& wname) const {
    auto global_well_index = this->name_map.at(wname);
    return this->m_in_producing_group[global_well_index];
}


void GlobalWellInfo::update_group(const WellContainer<Well::Status>& well_status, const std::vector<Well::InjectorCMode>& injection_cmode, const std::vector<Well::ProducerCMode>& production_cmode) {
    if (well_status.size() != this->local_map.size())
        throw std::logic_error("Size mismatch");

    this->m_in_injecting_group.assign(this->name_map.size(), 0);
    this->m_in_producing_group.assign(this->name_map.size(), 0);
    for (std::size_t well_index = 0; well_index < well_status.size(); well_index++) {
        if (well_status[well_index] == Well::Status::OPEN) {
            if (this->is_injector[well_index]) {
                if (injection_cmode[well_index] == Well::InjectorCMode::GRUP)
                    this->m_in_injecting_group[this->local_map[well_index]] = 1;
            } else {
                if (production_cmode[well_index] == Well::ProducerCMode::GRUP)
                    this->m_in_producing_group[this->local_map[well_index]] = 1;
            }
        }
    }
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
