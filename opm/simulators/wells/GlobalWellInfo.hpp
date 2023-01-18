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

#ifndef OPM_GLOBAL_WELL_INFO_HEADER_INCLUDED
#define OPM_GLOBAL_WELL_INFO_HEADER_INCLUDED

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace Opm {

class Schedule;
class Well;
enum class WellInjectorCMode;
enum class WellProducerCMode;
enum class WellStatus;


/*
  The context of the GlobalWellInfo class is the situation where the wells are
  distributed among different processors. Most well processing only considers
  the wells defined on the local process, but in some cases we need global
  information about all the wells. This class maintains the following:

  - Mapping between global well index and well name.

  - Mapping between local well index and global index (only used internally in
    class).

  - Functionality to query well whether it is currently injecting or producing
    under group control.
*/


class GlobalWellInfo {
public:


    /*
      Will sum the m_in_injecting_group and m_in_producing_group vectors across
      all processes, so that all processes can query for an arbitrary well.
    */
    template <typename Comm>
    void communicate(const Comm& comm) {
        auto size = this->m_in_injecting_group.size();
        comm.sum( this->m_in_injecting_group.data(), size);
        comm.sum( this->m_in_producing_group.data(), size);
    }



    GlobalWellInfo(const Schedule& sched, std::size_t report_step, const std::vector<Well>& local_wells);
    bool in_producing_group(const std::string& wname) const;
    bool in_injecting_group(const std::string& wname) const;
    std::size_t well_index(const std::string& wname) const;
    const std::string& well_name(std::size_t well_index) const;
    void update_injector(std::size_t well_index, WellStatus well_status, WellInjectorCMode injection_cmode);
    void update_producer(std::size_t well_index, WellStatus well_status, WellProducerCMode production_cmode);
    void clear();

private:
    std::vector<std::size_t> local_map;    // local_index -> global_index

    std::map<std::string, std::size_t> name_map; // string -> global_index
    std::vector<int> m_in_injecting_group;       // global_index -> int/bool
    std::vector<int> m_in_producing_group;       // global_index -> int/bool
};


}

#endif

