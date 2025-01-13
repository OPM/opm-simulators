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
#include <opm/models/blackoil/blackoilsolventparams.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SsfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Sof2Table.hpp>
#include <opm/input/eclipse/EclipseState/Tables/MsfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PmiscTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/MiscTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SorwmisTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SgcwmisTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TlpmixpaTable.hpp>
#endif

namespace Opm {

#if HAVE_ECL_INPUT
template<class Scalar>
template<bool enableSolvent>
void BlackOilSolventParams<Scalar>::
initFromState(const EclipseState& eclState, const Schedule& schedule)
{
    // some sanity checks: if solvents are enabled, the SOLVENT keyword must be
    // present, if solvents are disabled the keyword must not be present.
    if constexpr (enableSolvent) {
        if (!eclState.runspec().phases().active(Phase::SOLVENT)) {
            throw std::runtime_error("Non-trivial solvent treatment requested at compile "
                                     "time, but the deck does not contain the SOLVENT keyword");
        }
    } else {
        if (eclState.runspec().phases().active(Phase::SOLVENT)) {
            throw std::runtime_error("Solvent treatment disabled at compile time, but the deck "
                                     "contains the SOLVENT keyword");
        }
    }

    if (!eclState.runspec().phases().active(Phase::SOLVENT)) {
        return; // solvent treatment is supposed to be disabled
    }

    co2sol_ = eclState.runspec().co2Sol();
    h2sol_ = eclState.runspec().h2Sol();

    if (co2sol_ && h2sol_) {
        throw std::runtime_error("CO2SOL and H2SOL can not be used together");
    }

    if (co2sol_ || h2sol_) {
        if (co2sol_) {
            co2GasPvt_.initFromState(eclState, schedule);
            brineCo2Pvt_.initFromState(eclState, schedule);
        } else {
            h2GasPvt_.initFromState(eclState, schedule);
            brineH2Pvt_.initFromState(eclState, schedule);
        }
        if (eclState.getSimulationConfig().hasDISGASW()) {
            rsSolw_active_ = true;
        }
    } else
        solventPvt_.initFromState(eclState, schedule);

    const auto& tableManager = eclState.getTableManager();
    // initialize the objects which deal with the SSFN keyword
    const auto& ssfnTables = tableManager.getSsfnTables();
    unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
    setNumSatRegions(numSatRegions);
    for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
        const auto& ssfnTable = ssfnTables.template getTable<SsfnTable>(satRegionIdx);
        ssfnKrg_[satRegionIdx].setXYContainers(ssfnTable.getSolventFractionColumn(),
                                               ssfnTable.getGasRelPermMultiplierColumn(),
                                               /*sortInput=*/true);
        ssfnKrs_[satRegionIdx].setXYContainers(ssfnTable.getSolventFractionColumn(),
                                               ssfnTable.getSolventRelPermMultiplierColumn(),
                                               /*sortInput=*/true);
    }

    // initialize the objects needed for miscible solvent and oil simulations
    isMiscible_ = false;
    if (!eclState.getTableManager().getMiscTables().empty()) {
        isMiscible_ = true;

        unsigned numMiscRegions = 1;

        // misicible hydrocabon relative permeability wrt water
        const auto& sof2Tables = tableManager.getSof2Tables();
        if (!sof2Tables.empty()) {
            // resize the attributes of the object
            sof2Krn_.resize(numSatRegions);
            for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
                const auto& sof2Table = sof2Tables.template getTable<Sof2Table>(satRegionIdx);
                sof2Krn_[satRegionIdx].setXYContainers(sof2Table.getSoColumn(),
                                                       sof2Table.getKroColumn(),
                                                       /*sortInput=*/true);
            }
        }
        else if (eclState.runspec().phases().active(Phase::OIL)) {
            throw std::runtime_error("SOF2 must be specified in MISCIBLE (SOLVENT and OIL) runs\n");
        }

        const auto& miscTables = tableManager.getMiscTables();
        if (!miscTables.empty()) {
            assert(numMiscRegions == miscTables.size());

            // resize the attributes of the object
            misc_.resize(numMiscRegions);
            for (unsigned miscRegionIdx = 0; miscRegionIdx < numMiscRegions; ++miscRegionIdx) {
                const auto& miscTable = miscTables.template getTable<MiscTable>(miscRegionIdx);

                // solventFraction = Ss / (Ss + Sg);
                const auto& solventFraction = miscTable.getSolventFractionColumn();
                const auto& misc = miscTable.getMiscibilityColumn();
                misc_[miscRegionIdx].setXYContainers(solventFraction, misc);
            }
        }
        else {
            throw std::runtime_error("MISC must be specified in MISCIBLE (SOLVENT) runs\n");
        }

        // resize the attributes of the object
        pmisc_.resize(numMiscRegions);
        const auto& pmiscTables = tableManager.getPmiscTables();
        if (!pmiscTables.empty()) {
            assert(numMiscRegions == pmiscTables.size());

            for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                const auto& pmiscTable = pmiscTables.template getTable<PmiscTable>(regionIdx);

                // Copy data
                const auto& po = pmiscTable.getOilPhasePressureColumn();
                const auto& pmisc = pmiscTable.getMiscibilityColumn();

                pmisc_[regionIdx].setXYContainers(po, pmisc);
            }
        }
        else {
            std::vector<double> x = {0.0,1.0e20};
            std::vector<double> y = {1.0,1.0};
            TabulatedFunction constant = TabulatedFunction(2, x, y);
            for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                pmisc_[regionIdx] = constant;
            }
        }

        // miscible relative permeability multipleiers
        msfnKrsg_.resize(numSatRegions);
        msfnKro_.resize(numSatRegions);
        const auto& msfnTables = tableManager.getMsfnTables();
        if (!msfnTables.empty()) {
            assert(numSatRegions == msfnTables.size());

            for (unsigned regionIdx = 0; regionIdx < numSatRegions; ++regionIdx) {
                const MsfnTable& msfnTable = msfnTables.template getTable<MsfnTable>(regionIdx);

                // Copy data
                // Ssg = Ss + Sg;
                const auto& Ssg = msfnTable.getGasPhaseFractionColumn();
                const auto& krsg = msfnTable.getGasSolventRelpermMultiplierColumn();
                const auto& kro = msfnTable.getOilRelpermMultiplierColumn();

                msfnKrsg_[regionIdx].setXYContainers(Ssg, krsg);
                msfnKro_[regionIdx].setXYContainers(Ssg, kro);
            }
        }
        else {
            std::vector<double> x = {0.0,1.0};
            std::vector<double> y = {1.0,0.0};
            TabulatedFunction unit = TabulatedFunction(2, x, x);
            TabulatedFunction invUnit = TabulatedFunction(2, x, y);

            for (unsigned regionIdx = 0; regionIdx < numSatRegions; ++regionIdx) {
                setMsfn(regionIdx, unit, invUnit);
            }
        }
        // resize the attributes of the object
        sorwmis_.resize(numMiscRegions);
        const auto& sorwmisTables = tableManager.getSorwmisTables();
        if (!sorwmisTables.empty()) {
            assert(numMiscRegions == sorwmisTables.size());

            for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                const auto& sorwmisTable = sorwmisTables.template getTable<SorwmisTable>(regionIdx);

                // Copy data
                const auto& sw = sorwmisTable.getWaterSaturationColumn();
                const auto& sorwmis = sorwmisTable.getMiscibleResidualOilColumn();

                sorwmis_[regionIdx].setXYContainers(sw, sorwmis);
            }
        }
        else {
            // default
            std::vector<double> x = {0.0,1.0};
            std::vector<double> y = {0.0,0.0};
            TabulatedFunction zero = TabulatedFunction(2, x, y);
            for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                sorwmis_[regionIdx] = zero;
            }
        }

        // resize the attributes of the object
        sgcwmis_.resize(numMiscRegions);
        const auto& sgcwmisTables = tableManager.getSgcwmisTables();
        if (!sgcwmisTables.empty()) {
            assert(numMiscRegions == sgcwmisTables.size());

            for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                const auto& sgcwmisTable = sgcwmisTables.template getTable<SgcwmisTable>(regionIdx);

                // Copy data
                const auto& sw = sgcwmisTable.getWaterSaturationColumn();
                const auto& sgcwmis = sgcwmisTable.getMiscibleResidualGasColumn();

                sgcwmis_[regionIdx].setXYContainers(sw, sgcwmis);
            }
        }
        else {
            // default
            std::vector<double> x = {0.0,1.0};
            std::vector<double> y = {0.0,0.0};
            TabulatedFunction zero = TabulatedFunction(2, x, y);
            for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx)
                sgcwmis_[regionIdx] = zero;
        }

        const auto& tlmixpar = eclState.getTableManager().getTLMixpar();
        if (!tlmixpar.empty()) {
            // resize the attributes of the object
            tlMixParamViscosity_.resize(numMiscRegions);
            tlMixParamDensity_.resize(numMiscRegions);

            assert(numMiscRegions == tlmixpar.size());
            for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                const auto& tlp = tlmixpar[regionIdx];
                tlMixParamViscosity_[regionIdx] = tlp.viscosity_parameter;
                tlMixParamDensity_[regionIdx] = tlp.density_parameter;
            }
        }
        else
            throw std::runtime_error("TLMIXPAR must be specified in MISCIBLE (SOLVENT) runs\n");

        // resize the attributes of the object
        tlPMixTable_.resize(numMiscRegions);
        if (!eclState.getTableManager().getTlpmixpaTables().empty()) {
            const auto& tlpmixparTables = tableManager.getTlpmixpaTables();
            if (!tlpmixparTables.empty()) {
                assert(numMiscRegions == tlpmixparTables.size());
                for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx) {
                    const auto& tlpmixparTable = tlpmixparTables.template getTable<TlpmixpaTable>(regionIdx);

                    // Copy data
                    const auto& po = tlpmixparTable.getOilPhasePressureColumn();
                    const auto& tlpmixpa = tlpmixparTable.getMiscibilityColumn();

                    tlPMixTable_[regionIdx].setXYContainers(po, tlpmixpa);
                }
            }
            else {
                // if empty keyword. Try to use the pmisc table as default.
                if (!pmisc_.empty()) {
                    tlPMixTable_ = pmisc_;
                }
                else {
                    throw std::invalid_argument("If the pressure dependent TL values in "
                                                "TLPMIXPA is defaulted (no entries), then "
                                                "the PMISC tables must be specified.");
                }
            }
        }
        else {
            // default
            std::vector<double> x = {0.0,1.0e20};
            std::vector<double> y = {1.0,1.0};
            TabulatedFunction ones = TabulatedFunction(2, x, y);
            for (unsigned regionIdx = 0; regionIdx < numMiscRegions; ++regionIdx)
                tlPMixTable_[regionIdx] = ones;
        }
    }
}
#endif

template<class Scalar>
void BlackOilSolventParams<Scalar>::
setNumSatRegions(unsigned numRegions)
{
    ssfnKrg_.resize(numRegions);
    ssfnKrs_.resize(numRegions);
}

template<class Scalar>
void BlackOilSolventParams<Scalar>::
setMsfn(unsigned satRegionIdx,
        const TabulatedFunction& msfnKrsg,
        const TabulatedFunction& msfnKro)
{
    msfnKrsg_[satRegionIdx] = msfnKrsg;
    msfnKro_[satRegionIdx] = msfnKro;
}

#define INSTANTIATE_TYPE(T)                                                                             \
    template struct BlackOilSolventParams<T>;                                                           \
    template void BlackOilSolventParams<T>::initFromState<false>(const EclipseState&, const Schedule&); \
    template void BlackOilSolventParams<T>::initFromState<true>(const EclipseState&, const Schedule&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
