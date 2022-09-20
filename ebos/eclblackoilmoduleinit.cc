/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <ebos/eclblackoilmoduleinit.hh>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvtwsaltTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PermfactTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SaltSolubilityTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>

#include <opm/models/blackoil/blackoilbrineparams.hh>

#include <stdexcept>

namespace Opm {

template<bool enableSaltPrecipitation, class Scalar>
BlackOilBrineParams<Scalar>
setupBrineParams(bool enableBrine,
                 const EclipseState& eclState)
{
    BlackOilBrineParams<Scalar> params;

    // some sanity checks: if brine are enabled, the BRINE keyword must be
    // present, if brine are disabled the keyword must not be present.
    if (enableBrine && !eclState.runspec().phases().active(Phase::BRINE)) {
        throw std::runtime_error("Non-trivial brine treatment requested at compile time, but "
                                 "the deck does not contain the BRINE keyword");
    }
    else if (!enableBrine && eclState.runspec().phases().active(Phase::BRINE)) {
        throw std::runtime_error("Brine treatment disabled at compile time, but the deck "
                                 "contains the BRINE keyword");
    }

    if (!eclState.runspec().phases().active(Phase::BRINE))
        return params; // brine treatment is supposed to be disabled

    const auto& tableManager = eclState.getTableManager();

    unsigned numPvtRegions = tableManager.getTabdims().getNumPVTTables();
    params.referencePressure_.resize(numPvtRegions);

    const auto& pvtwsaltTables = tableManager.getPvtwSaltTables();

    // initialize the objects which deal with the BDENSITY keyword
    const auto& bdensityTables = tableManager.getBrineDensityTables();
    if (!bdensityTables.empty()) {
        params.bdensityTable_.resize(numPvtRegions);
        assert(numPvtRegions == bdensityTables.size());
        for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
            const auto& bdensityTable = bdensityTables[pvtRegionIdx];
            const auto& pvtwsaltTable = pvtwsaltTables[pvtRegionIdx];
            const auto& c = pvtwsaltTable.getSaltConcentrationColumn();
            params.bdensityTable_[pvtRegionIdx].setXYContainers(c, bdensityTable);
        }
    }

    if constexpr (enableSaltPrecipitation) {
        const TableContainer& permfactTables = tableManager.getPermfactTables();
        params.permfactTable_.resize(numPvtRegions);
        for (size_t i = 0; i < permfactTables.size(); ++i) {
            const PermfactTable& permfactTable = permfactTables.getTable<PermfactTable>(i);
            params.permfactTable_[i].setXYContainers(permfactTable.getPorosityChangeColumn(), permfactTable.getPermeabilityMultiplierColumn());
        }

        const TableContainer& saltsolTables = tableManager.getSaltsolTables();
        if (!saltsolTables.empty()) {
            params.saltsolTable_.resize(numPvtRegions);
            assert(numPvtRegions == saltsolTables.size());
            for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
                const SaltsolTable& saltsolTable = saltsolTables.getTable<SaltsolTable>(pvtRegionIdx );
                params.saltsolTable_[pvtRegionIdx] = saltsolTable.getSaltsolColumn().front();
            }
        }
    }

    return params;
}

template BlackOilBrineParams<double>
setupBrineParams<false,double>(bool, const EclipseState&);
template BlackOilBrineParams<double>
setupBrineParams<true,double>(bool, const EclipseState&);

}
