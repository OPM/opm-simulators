/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_AUTODIFF_VFPPROPERTIES_HPP_
#define OPM_AUTODIFF_VFPPROPERTIES_HPP_

#include <opm/simulators/wells/VFPInjProperties.hpp>
#include <opm/simulators/wells/VFPProdProperties.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>

#include <cstddef>
#include <map>

namespace Opm {

class VFPInjTable;
class VFPProdTable;

/**
 * A thin wrapper class that holds one VFPProdProperties and one
 * VFPInjProperties object.
 */
class VFPProperties {
public:
    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_tables A map of different VFPINJ tables.
     * @param prod_tables A map of different VFPPROD tables.
     */

    VFPProperties(const std::vector<std::reference_wrapper<const VFPInjTable>>& inj_tables,
                  const std::vector<std::reference_wrapper<const VFPProdTable>>& prod_tables,
                  const WellState& well_state)
                  :well_state_(well_state)
    {
        for (const auto& vfpinj : inj_tables)
            this->m_inj.addTable( vfpinj );

        for (const auto& vfpprod : prod_tables)
            this->m_prod.addTable( vfpprod );
    };

    /**
     * Returns the VFP properties for injection wells
     */
    const VFPInjProperties* getInj() const {
        return &m_inj;
    }

    /**
     * Returns the VFP properties for production wells
     */
    const VFPProdProperties* getProd() const {
        return &m_prod;
    }

    double getExplicitWFR(const int table_id, const std::size_t well_index) const {
        const auto& rates = well_state_.well(well_index).prev_surface_rates;
        const auto& pu = well_state_.phaseUsage();
        const auto& aqua = pu.phase_used[BlackoilPhases::Aqua]? rates[pu.phase_pos[BlackoilPhases::Aqua]]:0.0;
        const auto& liquid = pu.phase_used[BlackoilPhases::Liquid]? rates[pu.phase_pos[BlackoilPhases::Liquid]]:0.0;
        const auto& vapour = pu.phase_used[BlackoilPhases::Vapour]? rates[pu.phase_pos[BlackoilPhases::Vapour]]:0.0;
        const VFPProdTable& table = this->m_prod.getTable(table_id);
        return detail::getWFR(table, aqua, liquid, vapour);
    }

    double getExplicitGFR(const int table_id, const std::size_t well_index) const {
        const auto& rates = well_state_.well(well_index).prev_surface_rates;
        const auto& pu = well_state_.phaseUsage();
        const auto& aqua = pu.phase_used[BlackoilPhases::Aqua]? rates[pu.phase_pos[BlackoilPhases::Aqua]]:0.0;
        const auto& liquid = pu.phase_used[BlackoilPhases::Liquid]? rates[pu.phase_pos[BlackoilPhases::Liquid]]:0.0;
        const auto& vapour = pu.phase_used[BlackoilPhases::Vapour]? rates[pu.phase_pos[BlackoilPhases::Vapour]]:0.0;
        const VFPProdTable& table = this->m_prod.getTable(table_id);
        return detail::getGFR(table, aqua, liquid, vapour);
    }

    std::pair<double, double> 
    getFloIPR(const int table_id, const std::size_t well_index) const {
        std::pair<double,double>retval(0.0, 0.0);
        const VFPProdTable& table = this->m_prod.getTable(table_id);
        const auto& pu = well_state_.phaseUsage();
        const auto& ipr_a= well_state_.well(well_index).ipr_a;
        const auto& aqua_a = pu.phase_used[BlackoilPhases::Aqua]? ipr_a[pu.phase_pos[BlackoilPhases::Aqua]]:0.0;
        const auto& liquid_a = pu.phase_used[BlackoilPhases::Liquid]? ipr_a[pu.phase_pos[BlackoilPhases::Liquid]]:0.0;
        const auto& vapour_a = pu.phase_used[BlackoilPhases::Vapour]? ipr_a[pu.phase_pos[BlackoilPhases::Vapour]]:0.0;
        const auto& ipr_b= well_state_.well(well_index).ipr_b;
        const auto& aqua_b = pu.phase_used[BlackoilPhases::Aqua]? ipr_b[pu.phase_pos[BlackoilPhases::Aqua]]:0.0;
        const auto& liquid_b = pu.phase_used[BlackoilPhases::Liquid]? ipr_b[pu.phase_pos[BlackoilPhases::Liquid]]:0.0;
        const auto& vapour_b = pu.phase_used[BlackoilPhases::Vapour]? ipr_b[pu.phase_pos[BlackoilPhases::Vapour]]:0.0;
        return std::make_pair(detail::getFlo(table, aqua_a, liquid_a, vapour_a), 
                                detail::getFlo(table, aqua_b, liquid_b, vapour_b));
    }

private:
    VFPInjProperties m_inj;
    VFPProdProperties m_prod;
    const WellState& well_state_;

};


} //Namespace

#endif /* OPM_AUTODIFF_VFPPROPERTIES_HPP_ */
