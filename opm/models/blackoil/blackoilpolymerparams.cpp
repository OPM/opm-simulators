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
#include <opm/models/blackoil/blackoilpolymerparams.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PlyadsTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PlymaxTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PlyrockTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PlyshlogTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PlyviscTable.hpp>
#endif

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

namespace {

template<class To, class From>
std::vector<std::vector<To>>
convertVecToVec(const std::vector<std::vector<From>>& input)
{
    std::vector<std::vector<To>> output;
    output.reserve(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        output.emplace_back(input[i].begin(), input[i].end());
    }

    return output;
}

}

namespace Opm {

#if HAVE_ECL_INPUT
template<class Scalar>
template<bool enablePolymer, bool enablePolymerMolarWeight>
void BlackOilPolymerParams<Scalar>::
initFromState(const EclipseState& eclState)
{
    // some sanity checks: if polymers are enabled, the POLYMER keyword must be
    // present, if polymers are disabled the keyword must not be present.
    if constexpr (enablePolymer) {
        if (!eclState.runspec().phases().active(Phase::POLYMER)) {
            throw std::runtime_error("Non-trivial polymer treatment requested at compile time, but "
                                     "the deck does not contain the POLYMER keyword");
        }
    }
    else {
        if (eclState.runspec().phases().active(Phase::POLYMER)) {
            throw std::runtime_error("Polymer treatment disabled at compile time, but the deck "
                                     "contains the POLYMER keyword");
        }
    }

    if constexpr (enablePolymerMolarWeight) {
        if (!eclState.runspec().phases().active(Phase::POLYMW)) {
            throw std::runtime_error("Polymer molecular weight tracking is enabled at compile time, but "
                                     "the deck does not contain the POLYMW keyword");
        }
    }
    else {
        if (eclState.runspec().phases().active(Phase::POLYMW)) {
            throw std::runtime_error("Polymer molecular weight tracking is disabled at compile time, but the deck "
                                     "contains the POLYMW keyword");
        }
    }

    if constexpr (enablePolymerMolarWeight && !enablePolymer) {
        throw std::runtime_error("Polymer molecular weight tracking is enabled while polymer treatment "
                                 "is disabled at compile time");
    }

    if (!eclState.runspec().phases().active(Phase::POLYMER)) {
        return; // polymer treatment is supposed to be disabled
    }

    const auto& tableManager = eclState.getTableManager();

    unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
    setNumSatRegions(numSatRegions);

           // initialize the objects which deal with the PLYROCK keyword
    const auto& plyrockTables = tableManager.getPlyrockTables();
    if (!plyrockTables.empty()) {
        assert(numSatRegions == plyrockTables.size());
        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
            const auto& plyrockTable = plyrockTables.template getTable<PlyrockTable>(satRegionIdx);
            setPlyrock(satRegionIdx,
                       plyrockTable.getDeadPoreVolumeColumn()[0],
                       plyrockTable.getResidualResistanceFactorColumn()[0],
                       plyrockTable.getRockDensityFactorColumn()[0],
                       static_cast<AdsorptionBehaviour>(plyrockTable.getAdsorbtionIndexColumn()[0]),
                       plyrockTable.getMaxAdsorbtionColumn()[0]);
        }
    }
    else {
        throw std::runtime_error("PLYROCK must be specified in POLYMER runs\n");
    }

           // initialize the objects which deal with the PLYADS keyword
    const auto& plyadsTables = tableManager.getPlyadsTables();
    if (!plyadsTables.empty()) {
        assert(numSatRegions == plyadsTables.size());
        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
            const auto& plyadsTable = plyadsTables.template getTable<PlyadsTable>(satRegionIdx);
            // Copy data
            const auto& c = plyadsTable.getPolymerConcentrationColumn();
            const auto& ads = plyadsTable.getAdsorbedPolymerColumn();
            plyadsAdsorbedPolymer_[satRegionIdx].setXYContainers(c, ads);
        }
    }
    else {
        throw std::runtime_error("PLYADS must be specified in POLYMER runs\n");
    }


    unsigned numPvtRegions = tableManager.getTabdims().getNumPVTTables();
    plyviscViscosityMultiplierTable_.resize(numPvtRegions);

           // initialize the objects which deal with the PLYVISC keyword
    const auto& plyviscTables = tableManager.getPlyviscTables();
    if (!plyviscTables.empty()) {
        // different viscosity model is used for POLYMW
        if constexpr (enablePolymerMolarWeight) {
            OpmLog::warning("PLYVISC should not be used in POLYMW runs, "
                            "it will have no effect. A viscosity model based on PLYVMH is used instead.\n");
        }
        else {
            assert(numPvtRegions == plyviscTables.size());
            for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++pvtRegionIdx) {
                const auto& plyadsTable = plyviscTables.template getTable<PlyviscTable>(pvtRegionIdx);
                // Copy data
                const auto& c = plyadsTable.getPolymerConcentrationColumn();
                const auto& visc = plyadsTable.getViscosityMultiplierColumn();
                plyviscViscosityMultiplierTable_[pvtRegionIdx].setXYContainers(c, visc);
            }
        }
    }
    else if constexpr (!enablePolymerMolarWeight) {
        throw std::runtime_error("PLYVISC must be specified in POLYMER runs\n");
    }

           // initialize the objects which deal with the PLYMAX keyword
    const auto& plymaxTables = tableManager.getPlymaxTables();
    const unsigned numMixRegions = plymaxTables.size();
    setNumMixRegions(numMixRegions, enablePolymerMolarWeight);
    if (!plymaxTables.empty()) {
        for (unsigned mixRegionIdx = 0; mixRegionIdx < numMixRegions; ++mixRegionIdx) {
            const auto& plymaxTable = plymaxTables.template getTable<PlymaxTable>(mixRegionIdx);
            plymaxMaxConcentration_[mixRegionIdx] = plymaxTable.getPolymerConcentrationColumn()[0];
        }
    }
    else {
        throw std::runtime_error("PLYMAX must be specified in POLYMER runs\n");
    }

    if (!eclState.getTableManager().getPlmixparTable().empty()) {
        if constexpr (enablePolymerMolarWeight) {
            OpmLog::warning("PLMIXPAR should not be used in POLYMW runs, it will have no effect.\n");
        }
        else {
            const auto& plmixparTable = eclState.getTableManager().getPlmixparTable();
            // initialize the objects which deal with the PLMIXPAR keyword
            for (unsigned mixRegionIdx = 0; mixRegionIdx < numMixRegions; ++mixRegionIdx) {
                plymixparToddLongstaff_[mixRegionIdx] = plmixparTable[mixRegionIdx].todd_langstaff;
            }
        }
    }
    else if constexpr (!enablePolymerMolarWeight) {
        throw std::runtime_error("PLMIXPAR must be specified in POLYMER runs\n");
    }

    hasPlyshlog_ = eclState.getTableManager().hasTables("PLYSHLOG");
    hasShrate_ = eclState.getTableManager().useShrate();

    if ((hasPlyshlog_ || hasShrate_) && enablePolymerMolarWeight) {
        OpmLog::warning("PLYSHLOG and SHRATE should not be used in POLYMW runs, they will have no effect.\n");
    }

    if (hasPlyshlog_ && !enablePolymerMolarWeight) {
        const auto& plyshlogTables = tableManager.getPlyshlogTables();
        assert(numPvtRegions == plyshlogTables.size());
        plyshlogShearEffectRefMultiplier_.resize(numPvtRegions);
        plyshlogShearEffectRefLogVelocity_.resize(numPvtRegions);
        for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++pvtRegionIdx) {
            const auto& plyshlogTable = plyshlogTables.template getTable<PlyshlogTable>(pvtRegionIdx);

            Scalar plyshlogRefPolymerConcentration = plyshlogTable.getRefPolymerConcentration();
            auto waterVelocity = plyshlogTable.getWaterVelocityColumn().vectorCopy();
            auto shearMultiplier = plyshlogTable.getShearMultiplierColumn().vectorCopy();

                   // do the unit version here for the waterVelocity
            UnitSystem unitSystem = eclState.getDeckUnitSystem();
            double siFactor = hasShrate_? unitSystem.parse("1/Time").getSIScaling() : unitSystem.parse("Length/Time").getSIScaling();
            for (std::size_t i = 0; i < waterVelocity.size(); ++i) {
                waterVelocity[i] *= siFactor;
                // for plyshlog the input must be stored as logarithms
                // the interpolation is then done the log-space.
                waterVelocity[i] = std::log(waterVelocity[i]);
            }

            Scalar refViscMult = plyviscViscosityMultiplierTable_[pvtRegionIdx].eval(plyshlogRefPolymerConcentration, /*extrapolate=*/true);
            // convert the table using referece conditions
            for (std::size_t i = 0; i < waterVelocity.size(); ++i) {
                shearMultiplier[i] *= refViscMult;
                shearMultiplier[i] -= 1;
                shearMultiplier[i] /= (refViscMult - 1);
            }
            plyshlogShearEffectRefMultiplier_[pvtRegionIdx].resize(waterVelocity.size());
            plyshlogShearEffectRefLogVelocity_[pvtRegionIdx].resize(waterVelocity.size());

            for (std::size_t i = 0; i < waterVelocity.size(); ++i) {
                plyshlogShearEffectRefMultiplier_[pvtRegionIdx][i] = shearMultiplier[i];
                plyshlogShearEffectRefLogVelocity_[pvtRegionIdx][i] = waterVelocity[i];
            }
        }
    }

    if (hasShrate_ && !enablePolymerMolarWeight) {
        if (!hasPlyshlog_) {
            throw std::runtime_error("PLYSHLOG must be specified if SHRATE is used in POLYMER runs\n");
        }
        const auto& shrateTable = eclState.getTableManager().getShrateTable();
        shrate_.resize(numPvtRegions);
        for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++pvtRegionIdx) {
            if (shrateTable.empty()) {
                shrate_[pvtRegionIdx] = 4.8; //default;
            }
            else if (shrateTable.size() == numPvtRegions) {
                shrate_[pvtRegionIdx] = shrateTable[pvtRegionIdx].rate;
            }
            else {
                throw std::runtime_error("SHRATE must either have 0 or number of NUMPVT entries\n");
            }
        }
    }

    if constexpr (enablePolymerMolarWeight) {
        const auto& plyvmhTable = eclState.getTableManager().getPlyvmhTable();
        if (!plyvmhTable.empty()) {
            assert(plyvmhTable.size() == numMixRegions);
            for (std::size_t regionIdx = 0; regionIdx < numMixRegions; ++regionIdx) {
                plyvmhCoefficients_[regionIdx].k_mh = plyvmhTable[regionIdx].k_mh;
                plyvmhCoefficients_[regionIdx].a_mh = plyvmhTable[regionIdx].a_mh;
                plyvmhCoefficients_[regionIdx].gamma = plyvmhTable[regionIdx].gamma;
                plyvmhCoefficients_[regionIdx].kappa = plyvmhTable[regionIdx].kappa;
            }
        }
        else {
            throw std::runtime_error("PLYVMH keyword must be specified in POLYMW rus \n");
        }

               // handling PLYMWINJ keyword
        const auto& plymwinjTables = tableManager.getPlymwinjTables();
        for (const auto& table : plymwinjTables) {
            const int tableNumber = table.first;
            const auto& plymwinjtable = table.second;
            const std::vector<double>& throughput = plymwinjtable.getThroughputs();
            const std::vector<double>& watervelocity = plymwinjtable.getVelocities();
            const std::vector<std::vector<double>>& molecularweight = plymwinjtable.getMoleWeights();
            if constexpr (std::is_same_v<Scalar, float>) {
                const std::vector<Scalar> tp(throughput.begin(), throughput.end());
                const std::vector<Scalar> wv(watervelocity.begin(), watervelocity.end());
                const auto mw = convertVecToVec<float>(molecularweight);
                TabulatedTwoDFunction tablefunc(tp, wv, mw, true, false);
                plymwinjTables_[tableNumber] = std::move(tablefunc);
            } else {
                TabulatedTwoDFunction tablefunc(throughput, watervelocity, molecularweight, true, false);
                plymwinjTables_[tableNumber] = std::move(tablefunc);
            }
        }

               // handling SKPRWAT keyword
        const auto& skprwatTables = tableManager.getSkprwatTables();
        for (const auto& table : skprwatTables) {
            const int tableNumber = table.first;
            const auto& skprwattable = table.second;
            const std::vector<double>& throughput = skprwattable.getThroughputs();
            const std::vector<double>& watervelocity = skprwattable.getVelocities();
            const std::vector<std::vector<double>>& skinpressure = skprwattable.getSkinPressures();
            if constexpr (std::is_same_v<Scalar, float>) {
                const std::vector<Scalar> tp(throughput.begin(), throughput.end());
                const std::vector<Scalar> wv(watervelocity.begin(), watervelocity.end());
                const auto sp = convertVecToVec<float>(skinpressure);
                TabulatedTwoDFunction tablefunc(tp, wv, sp, true, false);
                skprwatTables_[tableNumber] = std::move(tablefunc);
            } else {
                TabulatedTwoDFunction tablefunc(throughput, watervelocity, skinpressure, true, false);
                skprwatTables_[tableNumber] = std::move(tablefunc);
            }
        }

               // handling SKPRPOLY keyword
        const auto& skprpolyTables = tableManager.getSkprpolyTables();
        for (const auto& table : skprpolyTables) {
            const int tableNumber = table.first;
            const auto& skprpolytable = table.second;
            const std::vector<double>& throughput = skprpolytable.getThroughputs();
            const std::vector<double>& watervelocity = skprpolytable.getVelocities();
            const std::vector<std::vector<double>>& skinpressure = skprpolytable.getSkinPressures();
            const double refPolymerConcentration = skprpolytable.referenceConcentration();
            if constexpr (std::is_same_v<Scalar, float>) {
                const std::vector<Scalar> tp(throughput.begin(), throughput.end());
                const std::vector<Scalar> wv(watervelocity.begin(), watervelocity.end());
                const auto sp = convertVecToVec<float>(skinpressure);
                SkprpolyTable tablefunc {
                    refPolymerConcentration,
                    TabulatedTwoDFunction(tp, wv, sp, true, false)
                };
                skprpolyTables_[tableNumber] = std::move(tablefunc);
            } else {
                SkprpolyTable tablefunc {
                    refPolymerConcentration,
                    TabulatedTwoDFunction(throughput, watervelocity, skinpressure, true, false)
                };
                skprpolyTables_[tableNumber] = std::move(tablefunc);
            }
        }
    }
}
#endif // HAVE_ECL_INPUT

template<class Scalar>
void BlackOilPolymerParams<Scalar>::
setNumSatRegions(unsigned numRegions)
{
    plyrockDeadPoreVolume_.resize(numRegions);
    plyrockResidualResistanceFactor_.resize(numRegions);
    plyrockRockDensityFactor_.resize(numRegions);
    plyrockAdsorbtionIndex_.resize(numRegions);
    plyrockMaxAdsorbtion_.resize(numRegions);
    plyadsAdsorbedPolymer_.resize(numRegions);
}

template<class Scalar>
void BlackOilPolymerParams<Scalar>::
setNumMixRegions(unsigned numRegions, bool enablePolymerMolarWeight)
{
    plymaxMaxConcentration_.resize(numRegions);
    plymixparToddLongstaff_.resize(numRegions);

    if (enablePolymerMolarWeight) {
        plyvmhCoefficients_.resize(numRegions);
    }
}

template<class Scalar>
void BlackOilPolymerParams<Scalar>::
setPlyrock(unsigned satRegionIdx,
           const Scalar& plyrockDeadPoreVolume,
           const Scalar& plyrockResidualResistanceFactor,
           const Scalar& plyrockRockDensityFactor,
           const Scalar& plyrockAdsorbtionIndex,
           const Scalar& plyrockMaxAdsorbtion)
{
    plyrockDeadPoreVolume_[satRegionIdx] = plyrockDeadPoreVolume;
    plyrockResidualResistanceFactor_[satRegionIdx] = plyrockResidualResistanceFactor;
    plyrockRockDensityFactor_[satRegionIdx] = plyrockRockDensityFactor;
    plyrockAdsorbtionIndex_[satRegionIdx] = plyrockAdsorbtionIndex;
    plyrockMaxAdsorbtion_[satRegionIdx] = plyrockMaxAdsorbtion;
}

#define INSTANTIATE_TYPE(T)                                                                  \
    template struct BlackOilPolymerParams<T>;                                                \
    template void BlackOilPolymerParams<T>::initFromState<false,false>(const EclipseState&); \
    template void BlackOilPolymerParams<T>::initFromState<true,false>(const EclipseState&);  \
    template void BlackOilPolymerParams<T>::initFromState<true,true>(const EclipseState&);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
