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
#include <opm/models/blackoil/blackoilmicpparams.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/BiofilmTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/DiffMICPTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PermfactTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

#include <algorithm>
#include <stdexcept>
#include <type_traits>

namespace Opm {

#if HAVE_ECL_INPUT
template<class Scalar>
template<bool enableMICP>
void BlackOilMICPParams<Scalar>::
initFromState(const EclipseState& eclState)
{
    // some sanity checks: if MICP is enabled, the MICP keyword must be
    // present, if MICP is disabled the keyword must not be present.
    if constexpr (enableMICP) {
        if (!eclState.runspec().micp()) {
            throw std::runtime_error("Non-trivial MICP treatment requested at compile time, but "
                                     "the deck does not contain the MICP keyword");
        }
    }
    else {
        if (eclState.runspec().micp()) {
            throw std::runtime_error("MICP treatment disabled at compile time, but the deck "
                                     "contains the MICP keyword");
        }
    }

    if (!eclState.runspec().micp())
        return; // MICP treatment is supposed to be disabled*/

    const auto& tableManager = eclState.getTableManager();
    unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
    unsigned numPvtRegions = tableManager.getTabdims().getNumPVTTables();

    // initialize the objects which deal with the MICP parameters
    const TableContainer& biofilmTables = tableManager.getBiofilmTables();
    if (biofilmTables.empty()) {
        throw std::runtime_error("MICP requires the BIOFPARA keyword");
    }
    densityBiofilm_.resize(numSatRegions);
    densityCalcite_.resize(numSatRegions);
    detachmentRate_.resize(numSatRegions);
    detachmentExponent_.resize(numSatRegions);
    halfVelocityOxygen_.resize(numSatRegions);
    halfVelocityUrea_.resize(numSatRegions);
    maximumGrowthRate_.resize(numSatRegions);
    maximumUreaUtilization_.resize(numSatRegions);
    microbialAttachmentRate_.resize(numSatRegions);
    microbialDeathRate_.resize(numSatRegions);
    oxygenConsumptionFactor_.resize(numSatRegions);
    yieldGrowthCoefficient_.resize(numSatRegions);
    yieldUreaToCalciteCoefficient_.resize(numSatRegions);
    for (unsigned stnRegionIdx = 0; stnRegionIdx < numSatRegions; ++stnRegionIdx) {
        const BiofilmTable& biofilmTable = biofilmTables.getTable<BiofilmTable>(stnRegionIdx);
        densityBiofilm_[stnRegionIdx] = biofilmTable.getDensityBiofilm().front();
        densityCalcite_[stnRegionIdx] = biofilmTable.getDensityCalcite().front();
        detachmentRate_[stnRegionIdx] = biofilmTable.getDetachmentRate().front();
        detachmentExponent_[stnRegionIdx] = biofilmTable.getDetachmentExponent().front();
        halfVelocityOxygen_[stnRegionIdx] = biofilmTable.getHalfVelocityOxygen().front();
        halfVelocityUrea_[stnRegionIdx] = biofilmTable.getHalfVelocityUrea().front();
        maximumGrowthRate_[stnRegionIdx] = biofilmTable.getMaximumGrowthRate().front();
        maximumUreaUtilization_[stnRegionIdx] = biofilmTable.getMaximumUreaUtilization().front();
        microbialAttachmentRate_[stnRegionIdx] = biofilmTable.getMicrobialAttachmentRate().front();
        microbialDeathRate_[stnRegionIdx] = biofilmTable.getMicrobialDeathRate().front();
        oxygenConsumptionFactor_[stnRegionIdx] = biofilmTable.getOxygenConsumptionFactor().front();
        yieldGrowthCoefficient_[stnRegionIdx] = biofilmTable.getYieldGrowthCoefficient().front();
        yieldUreaToCalciteCoefficient_[stnRegionIdx] = biofilmTable.getYieldUreaToCalciteCoefficient().front();
    }

    const TableContainer& diffMICPTables = tableManager.getDiffMICPTables();
    microbialDiffusion_.resize(numPvtRegions);
    oxygenDiffusion_.resize(numPvtRegions);
    ureaDiffusion_.resize(numPvtRegions);
    for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++pvtRegionIdx) {
        if (!diffMICPTables.empty()) {
            const DiffMICPTable& diffMICPTable = diffMICPTables.getTable<DiffMICPTable>(pvtRegionIdx);
            microbialDiffusion_[pvtRegionIdx] = diffMICPTable.getMicrobialDiffusion().front();
            oxygenDiffusion_[pvtRegionIdx] = diffMICPTable.getOxygenDiffusion().front();
            ureaDiffusion_[pvtRegionIdx] = diffMICPTable.getUreaDiffusion().front();
        }
        else {
            microbialDiffusion_[pvtRegionIdx] = 0.0;
            oxygenDiffusion_[pvtRegionIdx] = 0.0;
            ureaDiffusion_[pvtRegionIdx] = 0.0;
        }
    }

    const TableContainer& permfactTables = tableManager.getPermfactTables();
    if (permfactTables.empty()) {
        throw std::runtime_error("MICP requires the PERMFACT keyword");
    }
    permfactTable_.resize(numSatRegions);
    for (std::size_t i = 0; i < permfactTables.size(); ++i) {
        const PermfactTable& permfactTable = permfactTables.getTable<PermfactTable>(i);
        permfactTable_[i].setXYContainers(permfactTable.getPorosityChangeColumn(), permfactTable.getPermeabilityMultiplierColumn());
    }
}
#endif

#define INSTANTIATE_TYPE(T)                                                         \
    template struct BlackOilMICPParams<T>;                                          \
    template void BlackOilMICPParams<T>::initFromState<false>(const EclipseState&); \
    template void BlackOilMICPParams<T>::initFromState<true>(const EclipseState&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
