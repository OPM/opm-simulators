/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.

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
#ifndef OPM_DYNAMICLISTECONLIMITED_HPP
#define OPM_DYNAMICLISTECONLIMITED_HPP

#include <vector>
#include <string>
#include <map>

#include <cassert>

namespace Opm
{

    /// to handle the wells and connections violating economic limits.
    class DynamicListEconLimited
    {
    public:

        DynamicListEconLimited() {
        }

        bool wellShutEconLimited(const std::string& well_name) const {
            return std::find(m_shut_wells.begin(), m_shut_wells.end(), well_name) != m_shut_wells.end();
        }

        void addShutWell(const std::string& well_name) {
            assert( !wellShutEconLimited(well_name) );
            assert( !wellStoppedEconLimited(well_name) );

            m_shut_wells.push_back(well_name);
        }

        bool wellStoppedEconLimited(const std::string& well_name) const {
            return std::find(m_stopped_wells.begin(), m_stopped_wells.end(), well_name) != m_stopped_wells.end();
        }

        void addStoppedWell(const std::string& well_name) {
            assert( !wellShutEconLimited(well_name) );
            assert( !wellStoppedEconLimited(well_name) );

            m_stopped_wells.push_back(well_name);
        }


        // TODO: maybe completion better here
        bool anyConnectionClosedForWell(const std::string& well_name) const {
            return (m_cells_closed_connections.find(well_name) != m_cells_closed_connections.end());
        }

        const std::vector<int>& getClosedConnectionsForWell(const std::string& well_name) const {
            return (m_cells_closed_connections.find(well_name)->second);
        }

        void addClosedConnectionsForWell(const std::string& well_name,
                                         const int cell_closed_connection) {
            if (!anyConnectionClosedForWell(well_name)) {
                // first time adding a connection for the well
                std::vector<int> vector_cells = {cell_closed_connection};
                m_cells_closed_connections[well_name] = vector_cells;
            } else {
                std::vector<int>& closed_connections = m_cells_closed_connections.find(well_name)->second;
                closed_connections.push_back(cell_closed_connection);
            }
        }

    private:
        std::vector <std::string> m_shut_wells;
        std::vector <std::string> m_stopped_wells;
        // using grid cell number to indicate the location of the connections
        std::map<std::string, std::vector<int>> m_cells_closed_connections;
    };

} // namespace Opm
#endif  /* OPM_DYNAMICLISTECONLIMITED_HPP */

