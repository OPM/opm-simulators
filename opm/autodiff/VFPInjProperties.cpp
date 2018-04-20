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





#include "config.h"

#include <opm/autodiff/VFPInjProperties.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>

#include <opm/autodiff/VFPHelpers.hpp>



namespace Opm {


namespace detail {

} //Namespace detail









VFPInjProperties::VFPInjProperties() {

}



VFPInjProperties::VFPInjProperties(const VFPInjTable* table){
    m_tables[table->getTableNum()] = table;
}




VFPInjProperties::VFPInjProperties(const std::map<int, std::shared_ptr<const VFPInjTable> >& tables) {
    for (const auto& table : tables) {
        m_tables[table.first] = table.second.get();
    }
}

VFPInjProperties::ADB VFPInjProperties::bhp(const std::vector<int>& table_id,
                                            const Wells& wells,
                                            const ADB& qs,
                                            const ADB& thp_val) const {
    const int nw = wells.number_of_wells;

    //Short-hands for water / oil / gas phases
    //TODO enable support for two-phase.
    assert(wells.number_of_phases == 3);
    const ADB& w = subset(qs, Span(nw, 1, BlackoilPhases::Aqua*nw));
    const ADB& o = subset(qs, Span(nw, 1, BlackoilPhases::Liquid*nw));
    const ADB& g = subset(qs, Span(nw, 1, BlackoilPhases::Vapour*nw));

    return bhp(table_id, w, o, g, thp_val);
}

VFPInjProperties::ADB VFPInjProperties::bhp(const std::vector<int>& table_id,
                                            const ADB& aqua,
                                            const ADB& liquid,
                                            const ADB& vapour,
                                            const ADB& thp_arg) const {
    const int nw = thp_arg.size();

    std::vector<int> block_pattern = detail::commonBlockPattern(aqua, liquid, vapour, thp_arg);

    assert(static_cast<int>(table_id.size()) == nw);
    assert(aqua.size()     == nw);
    assert(liquid.size()   == nw);
    assert(vapour.size()   == nw);
    assert(thp_arg.size()      == nw);

    //Allocate data for bhp's and partial derivatives
    ADB::V value = ADB::V::Zero(nw);
    ADB::V dthp = ADB::V::Zero(nw);
    ADB::V dflo = ADB::V::Zero(nw);

    //Get the table for each well
    std::vector<const VFPInjTable*> well_tables(nw, nullptr);
    for (int i=0; i<nw; ++i) {
        if (table_id[i] > 0) {
            well_tables[i] = detail::getTable(m_tables, table_id[i]);
        }
    }

    //Get the right FLO variable for each well as a single ADB
    const ADB flo = detail::combineADBVars<VFPInjTable::FLO_TYPE>(well_tables, aqua, liquid, vapour);

    //Compute the BHP for each well independently
    for (int i=0; i<nw; ++i) {
        const VFPInjTable* table = well_tables[i];
        if (table != nullptr) {
            //First, find the values to interpolate between
            auto flo_i = detail::findInterpData(flo.value()[i], table->getFloAxis());
            auto thp_i = detail::findInterpData(thp_arg.value()[i], table->getTHPAxis());

            detail::VFPEvaluation bhp_val = detail::interpolate(table->getTable(), flo_i, thp_i);

            value[i] = bhp_val.value;
            dthp[i] = bhp_val.dthp;
            dflo[i] = bhp_val.dflo;
        }
        else {
            value[i] = -1e100; //Signal that this value has not been calculated properly, due to "missing" table
        }
    }

    //Create diagonal matrices from ADB::Vs
    ADB::M dthp_diag(dthp.matrix().asDiagonal());
    ADB::M dflo_diag(dflo.matrix().asDiagonal());

    //Calculate the Jacobians
    const int num_blocks = block_pattern.size();
    std::vector<ADB::M> jacs(num_blocks);
    for (int block = 0; block < num_blocks; ++block) {
        //Could have used fastSparseProduct and temporary variables
        //but may not save too much on that.
        jacs[block] = ADB::M(nw, block_pattern[block]);

        if (!thp_arg.derivative().empty()) {
            jacs[block] += dthp_diag * thp_arg.derivative()[block];
        }
        if (!flo.derivative().empty()) {
            jacs[block] += dflo_diag * flo.derivative()[block];
        }
    }

    ADB retval = ADB::function(std::move(value), std::move(jacs));
    return retval;
}





double VFPInjProperties::bhp(int table_id,
        const double& aqua,
        const double& liquid,
        const double& vapour,
        const double& thp_arg) const {
    const VFPInjTable* table = detail::getTable(m_tables, table_id);

    detail::VFPEvaluation retval = detail::bhp(table, aqua, liquid, vapour, thp_arg);
    return retval.value;
}








double VFPInjProperties::thp(int table_id,
        const double& aqua,
        const double& liquid,
        const double& vapour,
        const double& bhp_arg) const {
    const VFPInjTable* table = detail::getTable(m_tables, table_id);
    const VFPInjTable::array_type& data = table->getTable();

    //Find interpolation variables
    double flo = detail::getFlo(aqua, liquid, vapour, table->getFloType());

    const std::vector<double> thp_array = table->getTHPAxis();
    int nthp = thp_array.size();

    /**
     * Find the function bhp_array(thp) by creating a 1D view of the data
     * by interpolating for every value of thp. This might be somewhat
     * expensive, but let us assome that nthp is small
     */
    auto flo_i = detail::findInterpData(flo, table->getFloAxis());
    std::vector<double> bhp_array(nthp);
    for (int i=0; i<nthp; ++i) {
        auto thp_i = detail::findInterpData(thp_array[i], thp_array);
        bhp_array[i] = detail::interpolate(data, flo_i, thp_i).value;
    }

    double retval = detail::findTHP(bhp_array, thp_array, bhp_arg);
    return retval;
}






const VFPInjTable* VFPInjProperties::getTable(const int table_id) const {
    return detail::getTable(m_tables, table_id);
}







} //Namespace Opm
