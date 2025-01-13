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
#include <opm/models/blackoil/blackoilfoamparams.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/FoamadsTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/FoammobTable.hpp>
#endif

#include <cstddef>
#include <stdexcept>

namespace Opm {

#if HAVE_ECL_INPUT
template<class Scalar>
template<bool enableFoam>
void BlackOilFoamParams<Scalar>::
initFromState(const EclipseState& eclState)
{
    // some sanity checks: if foam is enabled, the FOAM keyword must be
    // present, if foam is disabled the keyword must not be present.
    if constexpr (enableFoam) {
        if (!eclState.runspec().phases().active(Phase::FOAM)) {
            throw std::runtime_error("Non-trivial foam treatment requested at compile time, but "
                                     "the deck does not contain the FOAM keyword");
        }
    }
    else {
        if (eclState.runspec().phases().active(Phase::FOAM)) {
            throw std::runtime_error("Foam treatment disabled at compile time, but the deck "
                                     "contains the FOAM keyword");
        }
    }

    if (!eclState.runspec().phases().active(Phase::FOAM)) {
        return; // foam treatment is supposed to be disabled
    }

    transport_phase_ = eclState.getInitConfig().getFoamConfig().getTransportPhase();

    if (eclState.getInitConfig().getFoamConfig().getMobilityModel() != FoamConfig::MobilityModel::TAB) {
        throw std::runtime_error("In FOAMOPTS, only TAB is allowed for the gas mobility factor reduction model.");
    }

    const auto& tableManager = eclState.getTableManager();
    const unsigned int numSatRegions = tableManager.getTabdims().getNumSatTables();
    setNumSatRegions(numSatRegions);
    const unsigned int numPvtRegions = tableManager.getTabdims().getNumPVTTables();
    gasMobilityMultiplierTable_.resize(numPvtRegions);

    // Get and check FOAMROCK data.
    const FoamConfig& foamConf = eclState.getInitConfig().getFoamConfig();
    if (numSatRegions != foamConf.size()) {
        throw std::runtime_error("Inconsistent sizes, number of saturation regions differ from the number of elements "
                                 "in FoamConfig, which typically corresponds to the number of records in FOAMROCK.");
    }

    // Get and check FOAMADS data.
    const auto& foamadsTables = tableManager.getFoamadsTables();
    if (foamadsTables.empty()) {
        throw std::runtime_error("FOAMADS must be specified in FOAM runs");
    }
    if (numSatRegions != foamadsTables.size()) {
        throw std::runtime_error("Inconsistent sizes, number of saturation regions differ from the "
                                 "number of FOAMADS tables.");
    }

    // Set data that vary with saturation region.
    for (std::size_t satReg = 0; satReg < numSatRegions; ++satReg) {
        const auto& rec = foamConf.getRecord(satReg);
        foamCoefficients_[satReg] = FoamCoefficients();
        foamCoefficients_[satReg].fm_min = rec.minimumSurfactantConcentration();
        foamCoefficients_[satReg].fm_surf = rec.referenceSurfactantConcentration();
        foamCoefficients_[satReg].ep_surf = rec.exponent();
        foamRockDensity_[satReg] = rec.rockDensity();
        foamAllowDesorption_[satReg] = rec.allowDesorption();
        const auto& foamadsTable = foamadsTables.template getTable<FoamadsTable>(satReg);
        const auto& conc = foamadsTable.getFoamConcentrationColumn();
        const auto& ads = foamadsTable.getAdsorbedFoamColumn();
        adsorbedFoamTable_[satReg].setXYContainers(conc, ads);
    }

    // Get and check FOAMMOB data.
    const auto& foammobTables = tableManager.getFoammobTables();
    if (foammobTables.empty()) {
        // When in the future adding support for the functional
        // model, FOAMMOB will not be required anymore (functional
        // family of keywords can be used instead, FOAMFSC etc.).
        throw std::runtime_error("FOAMMOB must be specified in FOAM runs");
    }
    if (numPvtRegions != foammobTables.size()) {
        throw std::runtime_error("Inconsistent sizes, number of PVT regions differ from the "
                                 "number of FOAMMOB tables.");
    }

    // Set data that vary with PVT region.
    for (std::size_t pvtReg = 0; pvtReg < numPvtRegions; ++pvtReg) {
        const auto& foammobTable = foammobTables.template getTable<FoammobTable>(pvtReg);
        const auto& conc = foammobTable.getFoamConcentrationColumn();
        const auto& mobMult = foammobTable.getMobilityMultiplierColumn();
        gasMobilityMultiplierTable_[pvtReg].setXYContainers(conc, mobMult);
    }
}
#endif

template<class Scalar>
void BlackOilFoamParams<Scalar>::setNumSatRegions(unsigned numRegions)
{
    foamCoefficients_.resize(numRegions);
    foamRockDensity_.resize(numRegions);
    foamAllowDesorption_.resize(numRegions);
    adsorbedFoamTable_.resize(numRegions);
}

#define INSTANTIATE_TYPE(T)                                                         \
    template struct BlackOilFoamParams<T>;                                          \
    template void BlackOilFoamParams<T>::initFromState<false>(const EclipseState&); \
    template void BlackOilFoamParams<T>::initFromState<true>(const EclipseState&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
