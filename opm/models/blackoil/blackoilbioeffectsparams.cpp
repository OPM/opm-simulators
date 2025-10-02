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
#include <opm/models/blackoil/blackoilbioeffectsparams.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/BiofilmTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/DiffMICPTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PcfactTable.hpp>
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
template<bool enableBioeffects , bool enableMICP>
void BlackOilBioeffectsParams<Scalar>::
initFromState(const EclipseState& eclState)
{
    // some sanity checks:
    // if biofilm is enabled, the BIOFILM keyword must be present,
    // if biofilm is disabled the keyword must not be present,
    // if MICP is enabled, the MICP keyword must be present, and
    // if MICP is disabled the keyword must not be present
    if constexpr (enableBioeffects) {
        if constexpr (enableMICP) {
            if (!eclState.runspec().micp()) {
                throw std::runtime_error("Non-trivial MICP treatment requested at compile time, but "
                                         "the deck does not contain the MICP keyword");
            }
        }
        else {
            if (!eclState.runspec().biof()) {
                throw std::runtime_error("Non-trivial Biofilm effects requested at compile time, but "
                                         "the deck does not contain the BIOFILM keyword");
            }
        }
    }
    else {
        if (eclState.runspec().biof()) {
            throw std::runtime_error("Biofilm effects disabled at compile time, but the deck "
                                     "contains the BIOFILM keyword");
        }
        else if (eclState.runspec().micp()) {
            throw std::runtime_error("MICP treatment disabled at compile time, but the deck "
                                     "contains the MICP keyword");
        }
    }

    if (!eclState.runspec().micp() && !eclState.runspec().biof())
        return; // bioeffects are supposed to be disabled

    const auto& tableManager = eclState.getTableManager();
    unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
    unsigned numPvtRegions = tableManager.getTabdims().getNumPVTTables();

    // initialize the objects which deal with the bioparameters
    const TableContainer& biofilmTables = tableManager.getBiofilmTables();
    if (biofilmTables.empty()) {
        if constexpr (enableMICP) {
            throw std::runtime_error("MICP requires the BIOFPARA keyword");
        }
        else {
            throw std::runtime_error("BIOFILM requires the BIOFPARA keyword");
        }
    }
    densityBiofilm_.resize(numSatRegions);
    microbialDeathRate_.resize(numSatRegions);
    maximumGrowthRate_.resize(numSatRegions);
    halfVelocityGrowth_.resize(numSatRegions);
    yieldGrowthCoefficient_.resize(numSatRegions);
    detachmentRate_.resize(numSatRegions);
    detachmentExponent_.resize(numSatRegions);
    microbialAttachmentRate_.resize(numSatRegions);
    if constexpr (enableMICP) {
        densityCalcite_.resize(numSatRegions);
        halfVelocityUrea_.resize(numSatRegions);
        maximumUreaUtilization_.resize(numSatRegions);
        oxygenConsumptionFactor_.resize(numSatRegions);
        yieldUreaToCalciteCoefficient_.resize(numSatRegions);
    }
    for (unsigned stnRegionIdx = 0; stnRegionIdx < numSatRegions; ++stnRegionIdx) {
        const BiofilmTable& biofilmTable = biofilmTables.getTable<BiofilmTable>(stnRegionIdx);
        densityBiofilm_[stnRegionIdx] = biofilmTable.getDensityBiofilm().front();
        microbialDeathRate_[stnRegionIdx] = biofilmTable.getMicrobialDeathRate().front();
        maximumGrowthRate_[stnRegionIdx] = biofilmTable.getMaximumGrowthRate().front();
        halfVelocityGrowth_[stnRegionIdx] = biofilmTable.getHalfVelocityGrowth().front();
        yieldGrowthCoefficient_[stnRegionIdx] = biofilmTable.getYieldGrowthCoefficient().front();
        detachmentRate_[stnRegionIdx] = biofilmTable.getDetachmentRate().front();
        detachmentExponent_[stnRegionIdx] = biofilmTable.getDetachmentExponent().front();
        microbialAttachmentRate_[stnRegionIdx] = biofilmTable.getMicrobialAttachmentRate().front();
        if constexpr (enableMICP) {
            densityCalcite_[stnRegionIdx] = biofilmTable.getDensityCalcite().front();
            halfVelocityUrea_[stnRegionIdx] = biofilmTable.getHalfVelocityUrea().front();
            maximumUreaUtilization_[stnRegionIdx] = biofilmTable.getMaximumUreaUtilization().front();
            oxygenConsumptionFactor_[stnRegionIdx] = biofilmTable.getOxygenConsumptionFactor().front();
            yieldUreaToCalciteCoefficient_[stnRegionIdx] = biofilmTable.getYieldUreaToCalciteCoefficient().front();
        }
    }

    const TableContainer& diffMICPTables = tableManager.getDiffMICPTables();
    bioDiffCoefficient_.resize(numPvtRegions);
    for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++pvtRegionIdx) {
        bioDiffCoefficient_[pvtRegionIdx].resize(BlackOilBioeffectsParams::numDiffCoef, 0.0);
        if (!diffMICPTables.empty()) {
            const DiffMICPTable& diffMICPTable = diffMICPTables.getTable<DiffMICPTable>(pvtRegionIdx);
            bioDiffCoefficient_[pvtRegionIdx][BlackOilBioeffectsParams::micrDiffIdx] = diffMICPTable.getMicrobialDiffusion().front();
            bioDiffCoefficient_[pvtRegionIdx][BlackOilBioeffectsParams::oxygDiffIdx] = diffMICPTable.getOxygenDiffusion().front();
            bioDiffCoefficient_[pvtRegionIdx][BlackOilBioeffectsParams::ureaDiffIdx] = diffMICPTable.getUreaDiffusion().front();
        }
    }

    if constexpr (!enableMICP) {
        const TableContainer& pcfactTables = tableManager.getPcfactTables();
        if (!pcfactTables.empty()) {
            pcfactTable_.resize(numSatRegions);
            for (size_t i = 0; i < pcfactTables.size(); ++i) {
                const PcfactTable& pcfactTable = pcfactTables.getTable<PcfactTable>(i);
                pcfactTable_[i].setXYContainers(pcfactTable.getPorosityChangeColumn(), pcfactTable.getPcMultiplierColumn());
            }
        }
    }

    const TableContainer& permfactTables = tableManager.getPermfactTables();
    if (permfactTables.empty()) {
        if constexpr (enableMICP) {
            throw std::runtime_error("MICP requires the PERMFACT keyword");
        }
        else {
            throw std::runtime_error("BIOFILM requires the PERMFACT keyword");
        }
    }
    permfactTable_.resize(numSatRegions);
    for (std::size_t i = 0; i < permfactTables.size(); ++i) {
        const PermfactTable& permfactTable = permfactTables.getTable<PermfactTable>(i);
        permfactTable_[i].setXYContainers(permfactTable.getPorosityChangeColumn(), permfactTable.getPermeabilityMultiplierColumn());
    }
}
#endif

#define INSTANTIATE_TYPE(T)                                                               \
    template struct BlackOilBioeffectsParams<T>;                                                \
    template void BlackOilBioeffectsParams<T>::initFromState<false,false>(const EclipseState&); \
    template void BlackOilBioeffectsParams<T>::initFromState<true,false>(const EclipseState&);  \
    template void BlackOilBioeffectsParams<T>::initFromState<true,true>(const EclipseState&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
