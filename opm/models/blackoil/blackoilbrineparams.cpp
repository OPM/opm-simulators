// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#include <opm/models/blackoil/blackoilbrineparams.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvtwsaltTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PcfactTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PermfactTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SaltSolubilityTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>
#endif

#include <cassert>
#include <cstddef>
#include <stdexcept>

namespace Opm {

#if HAVE_ECL_INPUT
template<class Scalar>
template<bool enableBrine, bool enableSaltPrecipitation>
void BlackOilBrineParams<Scalar>::
initFromState(const EclipseState& eclState)
{
    // some sanity checks: if brine are enabled, the BRINE keyword must be
    // present, if brine are disabled the keyword must not be present.
    if constexpr (enableBrine) {
        if (!eclState.runspec().phases().active(Phase::BRINE)) {
            throw std::runtime_error("Non-trivial brine treatment requested at compile time, but "
                                     "the deck does not contain the BRINE keyword");
        }
    }
    else {
        if (eclState.runspec().phases().active(Phase::BRINE)) {
            throw std::runtime_error("Brine treatment disabled at compile time, but the deck "
                                     "contains the BRINE keyword");
        }
    }

    if (!eclState.runspec().phases().active(Phase::BRINE)) {
        return; // brine treatment is supposed to be disabled
    }

    const auto& tableManager = eclState.getTableManager();

    unsigned numPvtRegions = tableManager.getTabdims().getNumPVTTables();
    referencePressure_.resize(numPvtRegions);

    const auto& pvtwsaltTables = tableManager.getPvtwSaltTables();

    // initialize the objects which deal with the BDENSITY keyword
    const auto& bdensityTables = tableManager.getBrineDensityTables();
    if (!bdensityTables.empty()) {
        bdensityTable_.resize(numPvtRegions);
        assert(numPvtRegions == bdensityTables.size());
        for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++pvtRegionIdx) {
            const auto& bdensityTable = bdensityTables[pvtRegionIdx];
            const auto& pvtwsaltTable = pvtwsaltTables[pvtRegionIdx];
            const auto& c = pvtwsaltTable.getSaltConcentrationColumn();
            bdensityTable_[pvtRegionIdx].setXYContainers(c, bdensityTable);
        }
    }

    if constexpr (enableSaltPrecipitation) {
        const TableContainer& permfactTables = tableManager.getPermfactTables();
        unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
        permfactTable_.resize(numSatRegions);
        for (std::size_t i = 0; i < permfactTables.size(); ++i) {
            const PermfactTable& permfactTable = permfactTables.getTable<PermfactTable>(i);
            permfactTable_[i].setXYContainers(permfactTable.getPorosityChangeColumn(), permfactTable.getPermeabilityMultiplierColumn());
        }

        const TableContainer& saltsolTables = tableManager.getSaltsolTables();
        if (!saltsolTables.empty()) {
            saltsolTable_.resize(numPvtRegions);
            saltdenTable_.resize(numPvtRegions);
            assert(numPvtRegions == saltsolTables.size());
            for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++pvtRegionIdx) {
                const SaltsolTable& saltsolTable = saltsolTables.getTable<SaltsolTable>(pvtRegionIdx );
                saltsolTable_[pvtRegionIdx] = saltsolTable.getSaltsolColumn().front();
                saltdenTable_[pvtRegionIdx] = saltsolTable.getSaltdenColumn().front();
            }
        }

        const TableContainer& pcfactTables = tableManager.getPcfactTables();
        if (!pcfactTables.empty()) {
            pcfactTable_.resize(numSatRegions);
            for (std::size_t i = 0; i < pcfactTables.size(); ++i) {
                const PcfactTable& pcfactTable = pcfactTables.getTable<PcfactTable>(i);
                pcfactTable_[i].setXYContainers(pcfactTable.getPorosityChangeColumn(), pcfactTable.getPcMultiplierColumn());
            }
        }
   }
}
#endif

#define INSTANTIATE_TYPE(T)                                                                \
    template struct BlackOilBrineParams<T>;                                                \
    template void BlackOilBrineParams<T>::initFromState<false,false>(const EclipseState&); \
    template void BlackOilBrineParams<T>::initFromState<true,false>(const EclipseState&);  \
    template void BlackOilBrineParams<T>::initFromState<false,true>(const EclipseState&);  \
    template void BlackOilBrineParams<T>::initFromState<true,true>(const EclipseState&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
