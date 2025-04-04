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

template<class Scalar>
GlobalWellInfo<Scalar>::
GlobalWellInfo(const Schedule& sched,
               std::size_t report_step,
               const std::vector<Well>& local_wells)
{
    auto num_wells = sched.numWells(report_step);
    this->m_in_injecting_group.resize(num_wells);
    this->m_in_producing_group.resize(num_wells);
    this->m_is_open.resize(num_wells);
    this->m_efficiency_scaling_factors.resize(num_wells);
    for (const auto& wname : sched.wellNames(report_step)) {
        const auto& well = sched.getWell(wname, report_step);
        auto global_well_index = well.seqIndex();
        this->name_map.emplace( well.name(), global_well_index );
    }

    for (const auto& well : local_wells)
        this->local_map.push_back( well.seqIndex() );
}

template<class Scalar>
bool GlobalWellInfo<Scalar>::
in_injecting_group(const std::string& wname) const
{
    auto global_well_index = this->name_map.at(wname);
    return this->m_in_injecting_group[global_well_index];
}

template<class Scalar>
bool GlobalWellInfo<Scalar>::
in_producing_group(const std::string& wname) const
{
    auto global_well_index = this->name_map.at(wname);
    return this->m_in_producing_group[global_well_index];
}

template<class Scalar>
bool GlobalWellInfo<Scalar>::
is_open(const std::string& wname) const
{
    auto global_well_index = this->name_map.at(wname);
    return this->m_is_open[global_well_index];
}

template<class Scalar>
void GlobalWellInfo<Scalar>::
update_injector(std::size_t well_index,
                Well::Status well_status,
                Well::InjectorCMode injection_cmode)
{
    if (well_status == Well::Status::OPEN) {
        this->m_is_open[this->local_map[well_index]] = 1;
        if (injection_cmode == Well::InjectorCMode::GRUP) {
            this->m_in_injecting_group[this->local_map[well_index]] = 1;
        } else {
            this->m_in_injecting_group[this->local_map[well_index]] = 0;
        }
    } else {
        this->m_is_open[this->local_map[well_index]] = 0;
    }
}

template<class Scalar>
void GlobalWellInfo<Scalar>::
update_producer(std::size_t well_index,
                Well::Status well_status,
                Well::ProducerCMode production_cmode)
{
    if (well_status == Well::Status::OPEN) {
        this->m_is_open[this->local_map[well_index]] = 1;
        if (production_cmode == Well::ProducerCMode::GRUP) {
            this->m_in_producing_group[this->local_map[well_index]] = 1;
        } else {
            this->m_in_producing_group[this->local_map[well_index]] = 0;
        }
    } else {
        this->m_is_open[this->local_map[well_index]] = 0;
    }
}

template<class Scalar>
void GlobalWellInfo<Scalar>::
update_efficiency_scaling_factor(std::size_t well_index,
                                 const Scalar efficiency_scaling_factor)
{
    this->m_efficiency_scaling_factors[this->local_map[well_index]] = efficiency_scaling_factor;
}

template<class Scalar>
void GlobalWellInfo<Scalar>::clear()
{
    this->m_in_injecting_group.assign(this->name_map.size(), 0);
    this->m_in_producing_group.assign(this->name_map.size(), 0);
    this->m_is_open.assign(this->name_map.size(), 0);
    this->m_efficiency_scaling_factors.assign(this->name_map.size(), 1.0);
}

template<class Scalar>
bool GlobalWellInfo<Scalar>::isRank0() const
{
    return is_rank0_;
}

template<class Scalar>
Scalar GlobalWellInfo<Scalar>::
efficiency_scaling_factor(const std::string& wname) const
{
    auto global_well_index = this->name_map.at(wname);
    return this->m_efficiency_scaling_factors[global_well_index];
}

template<class Scalar>
std::size_t GlobalWellInfo<Scalar>::
well_index(const std::string& wname) const
{
    return this->name_map.at(wname);
}

template<class Scalar>
const std::string& GlobalWellInfo<Scalar>::
well_name(std::size_t well_index) const
{
    for (const auto& [name, index] : this->name_map) {
        if (index == well_index)
            return name;
    }
    throw std::logic_error("No well with index: " + std::to_string(well_index));
}

template class GlobalWellInfo<double>;

#if FLOW_INSTANTIATE_FLOAT
template class GlobalWellInfo<float>;
#endif

}
