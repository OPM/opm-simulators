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
#define	OPM_DYNAMICLISTECONLIMITED_HPP

#include <vector>
#include <string>

#include <cassert>

namespace Opm
{

    /// to handle the wells and connections voilating economic limits.
    class DynamicListEconLimited
    {
    public:
        bool anyWellEconLimited() const {
            return !(m_shut_wells.empty());
        };

        bool wellEconLimited(const std::string& well_name) const {
            return std::find(m_shut_wells.begin(), m_shut_wells.end(), well_name) != m_shut_wells.end();
        };

        void addShuttedWell(const std::string& well_name) {
            // the well should not be in the list
            // TODO: not sure wheter a shutted well can
            // still be running through some other mechanism.
            assert( !wellEconLimited(well_name) );

            m_shut_wells.push_back(well_name);
        };

    private:
        std::vector <std::string> m_shut_wells;
    };

} // namespace Opm
#endif	/* OPM_DYNAMICLISTECONLIMITED_HPP */

