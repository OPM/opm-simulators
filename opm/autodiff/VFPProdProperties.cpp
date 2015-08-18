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

#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/VFPProdTable.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>

#include <opm/autodiff/VFPHelpers.hpp>



namespace Opm {




VFPProdProperties::VFPProdProperties() {

}



VFPProdProperties::VFPProdProperties(const VFPProdTable* table){
    m_tables[table->getTableNum()] = table;
}




VFPProdProperties::VFPProdProperties(const std::map<int, VFPProdTable>& tables) {
    for (const auto& table : tables) {
        m_tables[table.first] = &table.second;
    }
}




VFPProdProperties::ADB VFPProdProperties::bhp(const std::vector<int>& table_id,
        const Wells& wells,
        const ADB& qs,
        const ADB& thp,
        const ADB& alq) const {
    const int np = wells.number_of_phases;
    const int nw = wells.number_of_wells;

    //Short-hands for water / oil / gas phases
    //TODO enable support for two-phase.
    assert(np == 3);
    const ADB& w = subset(qs, Span(nw, 1, BlackoilPhases::Aqua*nw));
    const ADB& o = subset(qs, Span(nw, 1, BlackoilPhases::Liquid*nw));
    const ADB& g = subset(qs, Span(nw, 1, BlackoilPhases::Vapour*nw));

    return bhp(table_id, w, o, g, thp, alq);
}






VFPProdProperties::ADB VFPProdProperties::bhp(const std::vector<int>& table_id,
        const ADB& aqua,
        const ADB& liquid,
        const ADB& vapour,
        const ADB& thp,
        const ADB& alq) const {
    const int nw = thp.size();

    std::vector<int> block_pattern = detail::commonBlockPattern(aqua, liquid, vapour, thp, alq);

    assert(static_cast<int>(table_id.size()) == nw);
    assert(aqua.size()     == nw);
    assert(liquid.size()   == nw);
    assert(vapour.size()   == nw);
    assert(thp.size()      == nw);
    assert(alq.size()      == nw);

    //Allocate data for bhp's and partial derivatives
    ADB::V value = ADB::V::Zero(nw);
    ADB::V dthp = ADB::V::Zero(nw);
    ADB::V dwfr = ADB::V::Zero(nw);
    ADB::V dgfr = ADB::V::Zero(nw);
    ADB::V dalq = ADB::V::Zero(nw);
    ADB::V dflo = ADB::V::Zero(nw);

    //Get the table for each well
    std::vector<const VFPProdTable*> well_tables(nw, NULL);
    for (int i=0; i<nw; ++i) {
        if (table_id[i] >= 0) {
            well_tables[i] = detail::getTable(m_tables, table_id[i]);
        }
    }

    //Get the right FLO/GFR/WFR variable for each well as a single ADB
    const ADB flo = detail::gather_vars<VFPProdTable::FLO_TYPE>(well_tables, aqua, liquid, vapour);
    const ADB wfr = detail::gather_vars<VFPProdTable::WFR_TYPE>(well_tables, aqua, liquid, vapour);
    const ADB gfr = detail::gather_vars<VFPProdTable::GFR_TYPE>(well_tables, aqua, liquid, vapour);

    //Compute the BHP for each well independently
    for (int i=0; i<nw; ++i) {
        const VFPProdTable* table = well_tables[i];
        if (table != NULL) {
            //First, find the values to interpolate between
            auto flo_i = detail::findInterpData(flo.value()[i], table->getFloAxis());
            auto thp_i = detail::findInterpData(thp.value()[i], table->getTHPAxis());
            auto wfr_i = detail::findInterpData(wfr.value()[i], table->getWFRAxis());
            auto gfr_i = detail::findInterpData(gfr.value()[i], table->getGFRAxis());
            auto alq_i = detail::findInterpData(alq.value()[i], table->getALQAxis());

            detail::adb_like bhp_val = detail::interpolate(table->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);

            /*
            static const int N=40;
            std::cout << "bhp=" << bhp_val.value << ";" << std::endl;
            std::cout << "flo=" << flo.value()[i] << ";" << std::endl;
            std::cout << "thp=" << thp.value()[i] << ";" << std::endl;
            std::cout << "wfr=" << wfr.value()[i] << ";" << std::endl;
            std::cout << "gfr=" << gfr.value()[i] << ";" << std::endl;
            std::cout << "alq=" << alq.value()[i] << ";" << std::endl;
            std::cout << "bhp_vfp=[" << std::endl;
            for (int j=0; j<N; ++j) {
                const double start = table->getFloAxis().front();
                const double end = table->getFloAxis().back();
                const double dist = end - start;
                double flo_d = start + (j/static_cast<double>(N-1)) * dist;
                auto flo_i = detail::findInterpData(flo_d, table->getFloAxis());
                detail::adb_like bhp_val = detail::interpolate(table->getTable(), flo_i, thp_i, wfr_i, gfr_i, alq_i);
                std::cout << bhp_val.value << ",";
            }
            std::cout << "];" << std::endl;

            std::cout << "flo_vfp=[" << std::endl;
            for (int j=0; j<N; ++j) {
                const double start = table->getFloAxis().front();
                const double end = table->getFloAxis().back();
                const double dist = end - start;
                double flo_d = start + (j/static_cast<double>(N-1)) * dist;
                std::cout << flo_d << ",";
            }
            std::cout << "];" << std::endl;
            */


            value[i] = bhp_val.value;
            dthp[i] = bhp_val.dthp;
            dwfr[i] = bhp_val.dwfr;
            dgfr[i] = bhp_val.dgfr;
            dalq[i] = bhp_val.dalq;
            dflo[i] = bhp_val.dflo;
        }
        else {
            value[i] = -1e100; //Signal that this value has not been calculated properly, due to "missing" table
        }
    }

    //Create diagonal matrices from ADB::Vs
    ADB::M dthp_diag = spdiag(dthp);
    ADB::M dwfr_diag = spdiag(dwfr);
    ADB::M dgfr_diag = spdiag(dgfr);
    ADB::M dalq_diag = spdiag(dalq);
    ADB::M dflo_diag = spdiag(dflo);

    //Calculate the Jacobians
    const int num_blocks = block_pattern.size();
    std::vector<ADB::M> jacs(num_blocks);
    for (int block = 0; block < num_blocks; ++block) {
        //Could have used fastSparseProduct and temporary variables
        //but may not save too much on that.
        jacs[block] = ADB::M(nw, block_pattern[block]);

        if (!thp.derivative().empty()) {
            jacs[block] += dthp_diag * thp.derivative()[block];
        }
        if (!wfr.derivative().empty()) {
            jacs[block] += dwfr_diag * wfr.derivative()[block];
        }
        if (!gfr.derivative().empty()) {
            jacs[block] += dgfr_diag * gfr.derivative()[block];
        }
        if (!alq.derivative().empty()) {
            jacs[block] += dalq_diag * alq.derivative()[block];
        }
        if (!flo.derivative().empty()) {
            jacs[block] -= dflo_diag * flo.derivative()[block];
        }
    }

    ADB retval = ADB::function(std::move(value), std::move(jacs));
    return retval;
}



double VFPProdProperties::bhp(int table_id,
        const double& aqua,
        const double& liquid,
        const double& vapour,
        const double& thp,
        const double& alq) const {
    const VFPProdTable* table = detail::getTable(m_tables, table_id);

    detail::adb_like retval = detail::bhp(table, aqua, liquid, vapour, thp, alq);
    return retval.value;
}



double VFPProdProperties::thp(int table_id,
        const double& aqua,
        const double& liquid,
        const double& vapour,
        const double& bhp,
        const double& alq) const {
    const VFPProdTable* table = detail::getTable(m_tables, table_id);
    const VFPProdTable::array_type& data = table->getTable();

    //Find interpolation variables
    double flo = detail::getFlo(aqua, liquid, vapour, table->getFloType());
    double wfr = detail::getWFR(aqua, liquid, vapour, table->getWFRType());
    double gfr = detail::getGFR(aqua, liquid, vapour, table->getGFRType());

    const std::vector<double> thp_array = table->getTHPAxis();
    int nthp = thp_array.size();

    /**
     * Find the function bhp_array(thp) by creating a 1D view of the data
     * by interpolating for every value of thp. This might be somewhat
     * expensive, but let us assome that nthp is small
     */
    auto flo_i = detail::findInterpData(flo, table->getFloAxis());
    auto wfr_i = detail::findInterpData(wfr, table->getWFRAxis());
    auto gfr_i = detail::findInterpData(gfr, table->getGFRAxis());
    auto alq_i = detail::findInterpData(alq, table->getALQAxis());
    std::vector<double> bhp_array(nthp);
    for (int i=0; i<nthp; ++i) {
        auto thp_i = detail::findInterpData(thp_array[i], thp_array);
        bhp_array[i] = detail::interpolate(data, flo_i, thp_i, wfr_i, gfr_i, alq_i).value;
    }

    double thp = detail::findTHP(bhp_array, thp_array, bhp);
    return thp;
}






const VFPProdTable* VFPProdProperties::getTable(const int table_id) const {
    return detail::getTable(m_tables, table_id);
}







}
