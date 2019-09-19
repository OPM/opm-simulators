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
/*!
 * \file
 *
 * \brief Contains the classes required to extend the black-oil model by polymer.
 */
#ifndef EWOMS_BLACK_OIL_POLYMER_MODULE_HH
#define EWOMS_BLACK_OIL_POLYMER_MODULE_HH

#include "blackoilproperties.hh"
#include <opm/models/io/vtkblackoilpolymermodule.hh>
#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/IntervalTabulated2DFunction.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyadsTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlymaxTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyrockTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyshlogTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyviscTable.hpp>
#endif

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <string>

namespace Opm {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by polymer.
 */
template <class TypeTag, bool enablePolymerV = GET_PROP_VALUE(TypeTag, EnablePolymer)>
class BlackOilPolymerModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

    typedef typename Opm::Tabulated1DFunction<Scalar> TabulatedFunction;
    typedef typename Opm::IntervalTabulated2DFunction<Scalar> TabulatedTwoDFunction;

    static constexpr unsigned polymerConcentrationIdx = Indices::polymerConcentrationIdx;
    static constexpr unsigned polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;
    static constexpr unsigned contiPolymerEqIdx = Indices::contiPolymerEqIdx;
    static constexpr unsigned contiPolymerMolarWeightEqIdx = Indices::contiPolymerMWEqIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;


    static constexpr unsigned enablePolymer = enablePolymerV;
    static constexpr bool enablePolymerMolarWeight = GET_PROP_VALUE(TypeTag, EnablePolymerMW);

    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr unsigned numPhases = FluidSystem::numPhases;

    struct SkprpolyTable {
        double refConcentration;
        TabulatedTwoDFunction table_func;
    };

public:
    enum AdsorptionBehaviour { Desorption = 1, NoDesorption = 2 };

    // a struct containing the constants to calculate polymer viscosity
    // based on Mark-Houwink equation and Huggins equation, the constants are provided
    // by the keyword PLYVMH
    struct PlyvmhCoefficients {
        Scalar k_mh;
        Scalar a_mh;
        Scalar gamma;
        Scalar kappa;
    };

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the polymer module
     */
    static void initFromDeck(const Opm::Deck& deck, const Opm::EclipseState& eclState)
    {
        // some sanity checks: if polymers are enabled, the POLYMER keyword must be
        // present, if polymers are disabled the keyword must not be present.
        if (enablePolymer && !deck.hasKeyword("POLYMER")) {
            throw std::runtime_error("Non-trivial polymer treatment requested at compile time, but "
                                     "the deck does not contain the POLYMER keyword");
        }
        else if (!enablePolymer && deck.hasKeyword("POLYMER")) {
            throw std::runtime_error("Polymer treatment disabled at compile time, but the deck "
                                     "contains the POLYMER keyword");
        }

        if (enablePolymerMolarWeight && !deck.hasKeyword("POLYMW")) {
            throw std::runtime_error("Polymer molecular weight tracking is enabled at compile time, but "
                                     "the deck does not contain the POLYMW keyword");
        }
        else if (!enablePolymerMolarWeight && deck.hasKeyword("POLYMW")) {
            throw std::runtime_error("Polymer molecular weight tracking is disabled at compile time, but the deck "
                                     "contains the POLYMW keyword");
        }

        if (enablePolymerMolarWeight && !enablePolymer) {
            throw std::runtime_error("Polymer molecular weight tracking is enabled while polymer treatment "
                                     "is disabled at compile time");
        }

        if (!deck.hasKeyword("POLYMER"))
            return; // polymer treatment is supposed to be disabled

        const auto& tableManager = eclState.getTableManager();

        unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
        setNumSatRegions(numSatRegions);

        // initialize the objects which deal with the PLYROCK keyword
        const auto& plyrockTables = tableManager.getPlyrockTables();
        if (!plyrockTables.empty()) {
            assert(numSatRegions == plyrockTables.size());
            for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++ satRegionIdx) {
                const auto& plyrockTable = plyrockTables.template getTable<Opm::PlyrockTable>(satRegionIdx);
                setPlyrock(satRegionIdx,
                           plyrockTable.getDeadPoreVolumeColumn()[satRegionIdx],
                           plyrockTable.getResidualResistanceFactorColumn()[satRegionIdx],
                           plyrockTable.getRockDensityFactorColumn()[satRegionIdx],
                           static_cast<AdsorptionBehaviour>(plyrockTable.getAdsorbtionIndexColumn()[satRegionIdx]),
                           plyrockTable.getMaxAdsorbtionColumn()[satRegionIdx]);
            }
        }
        else {
            throw std::runtime_error("PLYROCK must be specified in POLYMER runs\n");
        }

        // initialize the objects which deal with the PLYADS keyword
        const auto& plyadsTables = tableManager.getPlyadsTables();
        if (!plyadsTables.empty()) {
            assert(numSatRegions == plyadsTables.size());
            for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++ satRegionIdx) {
                const auto& plyadsTable = plyadsTables.template getTable<Opm::PlyadsTable>(satRegionIdx);
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
        setNumPvtRegions(numPvtRegions);

        // initialize the objects which deal with the PLYVISC keyword
        const auto& plyviscTables = tableManager.getPlyviscTables();
        if (!plyviscTables.empty()) {
            // different viscosity model is used for POLYMW
            if (enablePolymerMolarWeight) {
                Opm::OpmLog::warning("PLYVISC should not be used in POLYMW runs, "
                                     "it will have no effect. A viscosity model based on PLYVMH is used instead.\n");
            }
            else {

                assert(numPvtRegions == plyviscTables.size());
                for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
                    const auto& plyadsTable = plyviscTables.template getTable<Opm::PlyviscTable>(pvtRegionIdx);
                    // Copy data
                    const auto& c = plyadsTable.getPolymerConcentrationColumn();
                    const auto& visc = plyadsTable.getViscosityMultiplierColumn();
                    plyviscViscosityMultiplierTable_[pvtRegionIdx].setXYContainers(c, visc);
                }
            }
        }
        else if (!enablePolymerMolarWeight) {
            throw std::runtime_error("PLYVISC must be specified in POLYMER runs\n");
        }

        // initialize the objects which deal with the PLYMAX keyword
        const auto& plymaxTables = tableManager.getPlymaxTables();
        const unsigned numMixRegions = plymaxTables.size();
        setNumMixRegions(numMixRegions);
        if (!plymaxTables.empty()) {
            for (unsigned mixRegionIdx = 0; mixRegionIdx < numMixRegions; ++ mixRegionIdx) {
                const auto& plymaxTable = plymaxTables.template getTable<Opm::PlymaxTable>(mixRegionIdx);
                setPlymax(mixRegionIdx, plymaxTable.getPolymerConcentrationColumn()[mixRegionIdx]);
            }
        }
        else {
            throw std::runtime_error("PLYMAX must be specified in POLYMER runs\n");
        }

        if (deck.hasKeyword("PLMIXPAR")) {
            if (enablePolymerMolarWeight) {
                Opm::OpmLog::warning("PLMIXPAR should not be used in POLYMW runs, it will have no effect.\n");
            }
            else {
                // initialize the objects which deal with the PLMIXPAR keyword
                for (unsigned mixRegionIdx = 0; mixRegionIdx < numMixRegions; ++ mixRegionIdx) {
                    const auto& plmixparRecord = deck.getKeyword("PLMIXPAR").getRecord(mixRegionIdx);
                    setPlmixpar(mixRegionIdx, plmixparRecord.getItem("TODD_LONGSTAFF").getSIDouble(0));
                }
            }
        }
        else if (!enablePolymerMolarWeight) {
            throw std::runtime_error("PLMIXPAR must be specified in POLYMER runs\n");
        }

        hasPlyshlog_ = deck.hasKeyword("PLYSHLOG");
        hasShrate_ = deck.hasKeyword("SHRATE");

        if ((hasPlyshlog_ || hasShrate_) && enablePolymerMolarWeight) {
            Opm::OpmLog::warning("PLYSHLOG and SHRATE should not be used in POLYMW runs, they will have no effect.\n");
        }

        if (hasPlyshlog_ && !enablePolymerMolarWeight) {
            const auto& plyshlogTables = tableManager.getPlyshlogTables();
            assert(numPvtRegions == plyshlogTables.size());
            plyshlogShearEffectRefMultiplier_.resize(numPvtRegions);
            plyshlogShearEffectRefLogVelocity_.resize(numPvtRegions);
            for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
                const auto& plyshlogTable = plyshlogTables.template getTable<Opm::PlyshlogTable>(pvtRegionIdx);

                Scalar plyshlogRefPolymerConcentration = plyshlogTable.getRefPolymerConcentration();
                auto waterVelocity = plyshlogTable.getWaterVelocityColumn().vectorCopy();
                auto shearMultiplier = plyshlogTable.getShearMultiplierColumn().vectorCopy();

                // do the unit version here for the waterVelocity
                Opm::UnitSystem unitSystem = deck.getActiveUnitSystem();
                double siFactor = hasShrate_? unitSystem.parse("1/Time").getSIScaling() : unitSystem.parse("Length/Time").getSIScaling();
                for (size_t i = 0; i < waterVelocity.size(); ++i) {
                    waterVelocity[i] *= siFactor;
                    // for plyshlog the input must be stored as logarithms
                    // the interpolation is then done the log-space.
                    waterVelocity[i] = std::log(waterVelocity[i]);
                }

                Scalar refViscMult = plyviscViscosityMultiplierTable_[pvtRegionIdx].eval(plyshlogRefPolymerConcentration, /*extrapolate=*/true);
                // convert the table using referece conditions
                for (size_t i = 0; i < waterVelocity.size(); ++i) {
                    shearMultiplier[i] *= refViscMult;
                    shearMultiplier[i] -= 1;
                    shearMultiplier[i] /= (refViscMult - 1);
                    shearMultiplier[i] = shearMultiplier[i];
                }
                plyshlogShearEffectRefMultiplier_[pvtRegionIdx].resize(waterVelocity.size());
                plyshlogShearEffectRefLogVelocity_[pvtRegionIdx].resize(waterVelocity.size());

                for (size_t i = 0; i < waterVelocity.size(); ++i) {
                    plyshlogShearEffectRefMultiplier_[pvtRegionIdx][i] = shearMultiplier[i];
                    plyshlogShearEffectRefLogVelocity_[pvtRegionIdx][i] = waterVelocity[i];
                }
            }
        }

        if (hasShrate_ && !enablePolymerMolarWeight) {
            if(!hasPlyshlog_) {
                throw std::runtime_error("PLYSHLOG must be specified if SHRATE is used in POLYMER runs\n");
            }
            const auto& shrateKeyword = deck.getKeyword("SHRATE");
            const std::vector<double>& shrateFromDeck = shrateKeyword.getSIDoubleData();
            shrate_.resize(numPvtRegions);
            for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
                if (shrateFromDeck.empty()) {
                    shrate_[pvtRegionIdx] = 4.8; //default;
                }
                else if (shrateFromDeck.size() == numPvtRegions) {
                    shrate_[pvtRegionIdx] = shrateKeyword.getSIDoubleData()[pvtRegionIdx];
                }
                else {
                    throw std::runtime_error("SHRATE must either have 0 or number of NUMPVT entries\n");
                }
            }
        }

        if (enablePolymerMolarWeight) {
            const Opm::DeckKeyword& plyvmhKeyword = deck.getKeyword("PLYVMH");
            assert(plyvmhKeyword.size() == numMixRegions);
            if (plyvmhKeyword.size() > 0) {
                for (size_t regionIdx = 0; regionIdx < plyvmhKeyword.size(); ++regionIdx) {
                    const Opm::DeckRecord& record = plyvmhKeyword.getRecord(regionIdx);
                    plyvmhCoefficients_[regionIdx].k_mh = record.getItem("K_MH").getSIDouble(0);
                    plyvmhCoefficients_[regionIdx].a_mh = record.getItem("A_MH").getSIDouble(0);
                    plyvmhCoefficients_[regionIdx].gamma = record.getItem("GAMMA").getSIDouble(0);
                    plyvmhCoefficients_[regionIdx].kappa = record.getItem("KAPPA").getSIDouble(0);
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
                TabulatedTwoDFunction tablefunc(throughput, watervelocity, molecularweight, true, false);
                plymwinjTables_[tableNumber] = std::move(tablefunc);
            }

            // handling SKPRWAT keyword
            const auto& skprwatTables = tableManager.getSkprwatTables();
            for (const auto& table : skprwatTables) {
                const int tableNumber = table.first;
                const auto& skprwattable = table.second;
                const std::vector<double>& throughput = skprwattable.getThroughputs();
                const std::vector<double>& watervelocity = skprwattable.getVelocities();
                const std::vector<std::vector<double>>& skinpressure = skprwattable.getSkinPressures();
                TabulatedTwoDFunction tablefunc(throughput, watervelocity, skinpressure, true, false);
                skprwatTables_[tableNumber] = std::move(tablefunc);
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
                SkprpolyTable tablefunc = {refPolymerConcentration,
                                           TabulatedTwoDFunction(throughput, watervelocity, skinpressure, true, false)};
                skprpolyTables_[tableNumber] = std::move(tablefunc);
            }
        }
    }
#endif

    /*!
     * \brief Specify the number of satuation regions.
     *
     * This must be called before setting the PLYROCK and PLYADS of any region.
     */
    static void setNumSatRegions(unsigned numRegions)
    {
        plyrockDeadPoreVolume_.resize(numRegions);
        plyrockResidualResistanceFactor_.resize(numRegions);
        plyrockRockDensityFactor_.resize(numRegions);
        plyrockAdsorbtionIndex_.resize(numRegions);
        plyrockMaxAdsorbtion_.resize(numRegions);
        plyadsAdsorbedPolymer_.resize(numRegions);
    }

    /*!
     * \brief Specify the polymer rock properties a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    static void setPlyrock(unsigned satRegionIdx,
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

    /*!
     * \brief Specify the number of pvt regions.
     *
     * This must be called before setting the PLYVISC of any region.
     */
    static void setNumPvtRegions(unsigned numRegions)
    {
        plyviscViscosityMultiplierTable_.resize(numRegions);
    }

    /*!
     * \brief Specify the polymer viscosity a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    static void setPlyvisc(unsigned satRegionIdx,
                           const TabulatedFunction& plyviscViscosityMultiplierTable)
    {
        plyviscViscosityMultiplierTable_[satRegionIdx] = plyviscViscosityMultiplierTable;
    }

    /*!
     * \brief Specify the number of mix regions.
     *
     * This must be called before setting the PLYMAC and PLMIXPAR of any region.
     */
    static void setNumMixRegions(unsigned numRegions)
    {
        plymaxMaxConcentration_.resize(numRegions);
        plymixparToddLongstaff_.resize(numRegions);

        if (enablePolymerMolarWeight) {
            plyvmhCoefficients_.resize(numRegions);
        }
    }

    /*!
     * \brief Specify the maximum polymer concentration a single region.
     *
     * The index of specified here must be in range [0, numMixRegionIdx)
     */
    static void setPlymax(unsigned mixRegionIdx,
                          const Scalar& plymaxMaxConcentration)
    {
        plymaxMaxConcentration_[mixRegionIdx] = plymaxMaxConcentration;
    }

    /*!
     * \brief Specify the maximum polymer concentration a single region.
     *
     * The index of specified here must be in range [0, numMixRegionIdx)
     */
    static void setPlmixpar(unsigned mixRegionIdx,
                            const Scalar& plymixparToddLongstaff)
    {
        plymixparToddLongstaff_[mixRegionIdx] = plymixparToddLongstaff;
    }

    /*!
    * \brief get the PLYMWINJ table
    */
    static TabulatedTwoDFunction& getPlymwinjTable(const int tableNumber)
    {
        const auto iterTable = plymwinjTables_.find(tableNumber);
        if (iterTable != plymwinjTables_.end()) {
            return iterTable->second;
        }
        else {
            throw std::runtime_error(" the PLYMWINJ table " + std::to_string(tableNumber) + " does not exist\n");
        }
    }

    /*!
    * \brief get the SKPRWAT table
    */
    static TabulatedTwoDFunction& getSkprwatTable(const int tableNumber)
    {
        const auto iterTable = skprwatTables_.find(tableNumber);
        if (iterTable != skprwatTables_.end()) {
            return iterTable->second;
        }
        else {
            throw std::runtime_error(" the SKPRWAT table " + std::to_string(tableNumber) + " does not exist\n");
        }
    }

    /*!
    * \brief get the SKPRPOLY table
    */
    static SkprpolyTable& getSkprpolyTable(const int tableNumber)
    {
        const auto iterTable = skprpolyTables_.find(tableNumber);
        if (iterTable != skprpolyTables_.end()) {
            return iterTable->second;
        }
        else {
            throw std::runtime_error(" the SKPRPOLY table " + std::to_string(tableNumber) + " does not exist\n");
        }
    }

    /*!
     * \brief Register all run-time parameters for the black-oil polymer module.
     */
    static void registerParameters()
    {
        if (!enablePolymer)
            // polymers have been disabled at compile time
            return;

        Opm::VtkBlackOilPolymerModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all polymer specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enablePolymer)
            // polymers have been disabled at compile time
            return;

        model.addOutputModule(new Opm::VtkBlackOilPolymerModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enablePolymer)
            // polymers have been disabled at compile time
            return false;

        if (!enablePolymerMolarWeight)
           return pvIdx == polymerConcentrationIdx;

        // both enablePolymer and enablePolymerMolarWeight are true here
        return pvIdx == polymerConcentrationIdx || pvIdx == polymerMoleWeightIdx;
    }

    static std::string primaryVarName(unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        if (pvIdx == polymerConcentrationIdx) {
            return "polymer_waterconcentration";
        }
        else {
            return "polymer_molecularweight";
        }
    }

    static Scalar primaryVarWeight(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enablePolymer)
            return false;

        if (!enablePolymerMolarWeight)
            return eqIdx == contiPolymerEqIdx;

        // both enablePolymer and enablePolymerMolarWeight are true here
        return eqIdx == contiPolymerEqIdx || eqIdx == contiPolymerMolarWeightEqIdx;
    }

    static std::string eqName(unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        if (eqIdx == contiPolymerEqIdx)
            return "conti^polymer";
        else
            return "conti^polymer_molecularweight";
    }

    static Scalar eqWeight(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    // must be called after water storage is computed
    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enablePolymer)
            return;

        const auto& fs = intQuants.fluidState();

        LhsEval surfaceVolumeWater =
                Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx))
                * Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.porosity());

        // avoid singular matrix if no water is present.
        surfaceVolumeWater = Opm::max(surfaceVolumeWater, 1e-10);

        // polymer in water phase
        const LhsEval massPolymer = surfaceVolumeWater
                * Toolbox::template decay<LhsEval>(intQuants.polymerConcentration())
                * (1.0 - Toolbox::template decay<LhsEval>(intQuants.polymerDeadPoreVolume()));

        // polymer in solid phase
        const LhsEval adsorptionPolymer =
                Toolbox::template decay<LhsEval>(1.0 - intQuants.porosity())
                * Toolbox::template decay<LhsEval>(intQuants.polymerRockDensity())
                * Toolbox::template decay<LhsEval>(intQuants.polymerAdsorption());

        LhsEval accumulationPolymer = massPolymer + adsorptionPolymer;

        storage[contiPolymerEqIdx] += accumulationPolymer;

        // tracking the polymer molecular weight
        if (enablePolymerMolarWeight) {
            accumulationPolymer = Opm::max(accumulationPolymer, 1e-10);

            storage[contiPolymerMolarWeightEqIdx]  += accumulationPolymer
                                         * Toolbox::template decay<LhsEval> (intQuants.polymerMoleWeight());
        }
    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)

    {
        if (!enablePolymer)
            return;

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        const unsigned upIdx = extQuants.upstreamIndex(FluidSystem::waterPhaseIdx);
        const unsigned inIdx = extQuants.interiorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const unsigned contiWaterEqIdx = Indices::conti0EqIdx + Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);


        if (upIdx == inIdx) {
            flux[contiPolymerEqIdx] =
                    extQuants.volumeFlux(waterPhaseIdx)
                    *up.fluidState().invB(waterPhaseIdx)
                    *up.polymerViscosityCorrection()
                    /extQuants.polymerShearFactor()
                    *up.polymerConcentration();

            // modify water
            flux[contiWaterEqIdx] /=
                    extQuants.waterShearFactor();
        }
        else {
            flux[contiPolymerEqIdx] =
                    extQuants.volumeFlux(waterPhaseIdx)
                    *Opm::decay<Scalar>(up.fluidState().invB(waterPhaseIdx))
                    *Opm::decay<Scalar>(up.polymerViscosityCorrection())
                    /Opm::decay<Scalar>(extQuants.polymerShearFactor())
                    *Opm::decay<Scalar>(up.polymerConcentration());

            // modify water
            flux[contiWaterEqIdx] /=
                    Opm::decay<Scalar>(extQuants.waterShearFactor());
        }

        // flux related to transport of polymer molecular weight
        if (enablePolymerMolarWeight) {
            if (upIdx == inIdx)
                flux[contiPolymerMolarWeightEqIdx] =
                    flux[contiPolymerEqIdx]*up.polymerMoleWeight();
            else
                flux[contiPolymerMolarWeightEqIdx] =
                    flux[contiPolymerEqIdx]*Opm::decay<Scalar>(up.polymerMoleWeight());
        }

    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
    {
        // do not consider consider the cange of polymer primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enablePolymer)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        outstream << priVars[polymerConcentrationIdx];
        outstream << priVars[polymerMoleWeightIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enablePolymer)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        instream >> priVars0[polymerConcentrationIdx];
        instream >> priVars0[polymerMoleWeightIdx];

        // set the primary variables for the beginning of the current time step.
        priVars1[polymerConcentrationIdx] = priVars0[polymerConcentrationIdx];
        priVars1[polymerMoleWeightIdx] = priVars0[polymerMoleWeightIdx];
    }

    static const Scalar plyrockDeadPoreVolume(const ElementContext& elemCtx,
                                              unsigned scvIdx,
                                              unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockDeadPoreVolume_[satnumRegionIdx];
    }

    static const Scalar plyrockResidualResistanceFactor(const ElementContext& elemCtx,
                                                        unsigned scvIdx,
                                                        unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockResidualResistanceFactor_[satnumRegionIdx];
    }

    static const Scalar plyrockRockDensityFactor(const ElementContext& elemCtx,
                                                 unsigned scvIdx,
                                                 unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockRockDensityFactor_[satnumRegionIdx];
    }

    static const Scalar plyrockAdsorbtionIndex(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockAdsorbtionIndex_[satnumRegionIdx];
    }

    static const Scalar plyrockMaxAdsorbtion(const ElementContext& elemCtx,
                                             unsigned scvIdx,
                                             unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyrockMaxAdsorbtion_[satnumRegionIdx];
    }

    static const TabulatedFunction& plyadsAdsorbedPolymer(const ElementContext& elemCtx,
                                                          unsigned scvIdx,
                                                          unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyadsAdsorbedPolymer_[satnumRegionIdx];
    }

    static const TabulatedFunction& plyviscViscosityMultiplierTable(const ElementContext& elemCtx,
                                                                    unsigned scvIdx,
                                                                    unsigned timeIdx)
    {
        unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
    }

    static const TabulatedFunction& plyviscViscosityMultiplierTable(unsigned pvtnumRegionIdx)
    {
        return plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
    }

    static const Scalar plymaxMaxConcentration(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    {
        unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plymaxMaxConcentration_[polymerMixRegionIdx];
    }

    static const Scalar plymixparToddLongstaff(const ElementContext& elemCtx,
                                               unsigned scvIdx,
                                               unsigned timeIdx)
    {
        unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plymixparToddLongstaff_[polymerMixRegionIdx];
    }

    static const PlyvmhCoefficients& plyvmhCoefficients(const ElementContext& elemCtx,
                                                        const unsigned scvIdx,
                                                        const unsigned timeIdx)
    {
        const unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return plyvmhCoefficients_[polymerMixRegionIdx];
    }

    static bool hasPlyshlog()
    {
        return hasPlyshlog_;
    }

    static bool hasShrate()
    {
        return hasShrate_;
    }

    static const Scalar shrate(unsigned pvtnumRegionIdx)
    {
        return shrate_[pvtnumRegionIdx];
    }

    /*!
     * \brief Computes the shear factor
     *
     * Input is polymer concentration and either the water velocity or the shrate if hasShrate_ is true.
     * The pvtnumRegionIdx is needed to make sure the right table is used.
     */
    template <class Evaluation>
    static Evaluation computeShearFactor(const Evaluation& polymerConcentration,
                                         unsigned pvtnumRegionIdx,
                                         const Evaluation& v0)
    {
        typedef Opm::MathToolbox<Evaluation> ToolboxLocal;

        const auto& viscosityMultiplierTable = plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
        Scalar viscosityMultiplier = viscosityMultiplierTable.eval(Opm::scalarValue(polymerConcentration), /*extrapolate=*/true);

        const Scalar eps = 1e-14;
        // return 1.0 if the polymer has no effect on the water.
        if (std::abs((viscosityMultiplier - 1.0)) < eps)
            return ToolboxLocal::createConstant(v0, 1.0);

        const std::vector<Scalar>& shearEffectRefLogVelocity = plyshlogShearEffectRefLogVelocity_[pvtnumRegionIdx];
        auto v0AbsLog = Opm::log(Opm::abs(v0));
        // return 1.0 if the velocity /sharte is smaller than the first velocity entry.
        if (v0AbsLog < shearEffectRefLogVelocity[0])
            return ToolboxLocal::createConstant(v0, 1.0);

        // compute shear factor from input
        // Z = (1 + (P - 1) * M(v)) / P
        // where M(v) is computed from user input
        // and P = viscosityMultiplier
        const std::vector<Scalar>& shearEffectRefMultiplier = plyshlogShearEffectRefMultiplier_[pvtnumRegionIdx];
        size_t numTableEntries = shearEffectRefLogVelocity.size();
        assert(shearEffectRefMultiplier.size() == numTableEntries);

        std::vector<Scalar> shearEffectMultiplier(numTableEntries, 1.0);
        for (size_t i = 0; i < numTableEntries; ++i) {
            shearEffectMultiplier[i] = (1.0 + (viscosityMultiplier - 1.0)*shearEffectRefMultiplier[i]) / viscosityMultiplier;
            shearEffectMultiplier[i] = Opm::log(shearEffectMultiplier[i]);
        }
        // store the logarithmic velocity and logarithmic multipliers in a table for easy look up and
        // linear interpolation in the logarithmic space.
        TabulatedFunction logShearEffectMultiplier = TabulatedFunction(numTableEntries, shearEffectRefLogVelocity, shearEffectMultiplier, /*bool sortInputs =*/ false);

        // Find sheared velocity (v) that satisfies
        // F = log(v) + log (Z) - log(v0) = 0;

        // Set up the function
        // u = log(v)
        auto F = [&logShearEffectMultiplier, &v0AbsLog](const Evaluation& u) {
            return u + logShearEffectMultiplier.eval(u, true) - v0AbsLog;
        };
        // and its derivative
        auto dF = [&logShearEffectMultiplier](const Evaluation& u) {
            return 1 + logShearEffectMultiplier.evalDerivative(u, true);
        };

        // Solve F = 0 using Newton
        // Use log(v0) as initial value for u
        auto u = v0AbsLog;
        bool converged = false;
        for (int i = 0; i < 20; ++i) {
            auto f = F(u);
            auto df = dF(u);
            u -= f/df;
            if (std::abs(Opm::scalarValue(f)) < 1e-12) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            throw std::runtime_error("Not able to compute shear velocity. \n");
        }

        // return the shear factor
        return Opm::exp(logShearEffectMultiplier.eval(u, /*extrapolate=*/true));

    }

    const Scalar molarMass() const
    {
        return 0.25; // kg/mol
    }




private:
    static std::vector<Scalar> plyrockDeadPoreVolume_;
    static std::vector<Scalar> plyrockResidualResistanceFactor_;
    static std::vector<Scalar> plyrockRockDensityFactor_;
    static std::vector<Scalar> plyrockAdsorbtionIndex_;
    static std::vector<Scalar> plyrockMaxAdsorbtion_;
    static std::vector<TabulatedFunction> plyadsAdsorbedPolymer_;
    static std::vector<TabulatedFunction> plyviscViscosityMultiplierTable_;
    static std::vector<Scalar> plymaxMaxConcentration_;
    static std::vector<Scalar> plymixparToddLongstaff_;
    static std::vector<std::vector<Scalar>> plyshlogShearEffectRefMultiplier_;
    static std::vector<std::vector<Scalar>> plyshlogShearEffectRefLogVelocity_;
    static std::vector<Scalar> shrate_;
    static bool hasShrate_;
    static bool hasPlyshlog_;

    static std::vector<PlyvmhCoefficients> plyvmhCoefficients_;
    static std::map<int, TabulatedTwoDFunction> plymwinjTables_;
    static std::map<int, TabulatedTwoDFunction> skprwatTables_;

    static std::map<int, SkprpolyTable> skprpolyTables_;
};



template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockDeadPoreVolume_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockResidualResistanceFactor_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockRockDensityFactor_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockAdsorbtionIndex_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockMaxAdsorbtion_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedFunction>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyadsAdsorbedPolymer_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedFunction>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyviscViscosityMultiplierTable_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plymaxMaxConcentration_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plymixparToddLongstaff_;

template <class TypeTag, bool enablePolymerV>
std::vector<std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyshlogShearEffectRefMultiplier_;

template <class TypeTag, bool enablePolymerV>
std::vector<std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyshlogShearEffectRefLogVelocity_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
BlackOilPolymerModule<TypeTag, enablePolymerV>::shrate_;

template <class TypeTag, bool enablePolymerV>
bool
BlackOilPolymerModule<TypeTag, enablePolymerV>::hasShrate_;

template <class TypeTag, bool enablePolymerV>
bool
BlackOilPolymerModule<TypeTag, enablePolymerV>::hasPlyshlog_;

template <class TypeTag, bool enablePolymerV>
std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::PlyvmhCoefficients>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plyvmhCoefficients_;

template <class TypeTag, bool enablePolymerV>
std::map<int, typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedTwoDFunction>
BlackOilPolymerModule<TypeTag, enablePolymerV>::plymwinjTables_;

template <class TypeTag, bool enablePolymerV>
std::map<int, typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedTwoDFunction>
BlackOilPolymerModule<TypeTag, enablePolymerV>::skprwatTables_;

template <class TypeTag, bool enablePolymerV>
std::map<int, typename BlackOilPolymerModule<TypeTag, enablePolymerV>::SkprpolyTable>
BlackOilPolymerModule<TypeTag, enablePolymerV>::skprpolyTables_;

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilPolymerIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        polymers extension of the black-oil model.
 */
template <class TypeTag, bool enablePolymerV = GET_PROP_VALUE(TypeTag, EnablePolymer)>
class BlackOilPolymerIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef BlackOilPolymerModule<TypeTag> PolymerModule;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    static constexpr int polymerConcentrationIdx = Indices::polymerConcentrationIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr bool enablePolymerMolarWeight = GET_PROP_VALUE(TypeTag, EnablePolymerMW);
    static constexpr int polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;


public:

    /*!
     * \brief Update the intensive properties needed to handle polymers from the
     *        primary variables
     *
     */
    void polymerPropertiesUpdate_(const ElementContext& elemCtx,
                                  unsigned dofIdx,
                                  unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        polymerConcentration_ = priVars.makeEvaluation(polymerConcentrationIdx, timeIdx);
        if (enablePolymerMolarWeight) {
            polymerMoleWeight_ = priVars.makeEvaluation(polymerMoleWeightIdx, timeIdx);
        }
        const Scalar cmax = PolymerModule::plymaxMaxConcentration(elemCtx, dofIdx, timeIdx);

        // permeability reduction due to polymer
        const Scalar& maxAdsorbtion = PolymerModule::plyrockMaxAdsorbtion(elemCtx, dofIdx, timeIdx);
        const auto& plyadsAdsorbedPolymer = PolymerModule::plyadsAdsorbedPolymer(elemCtx, dofIdx, timeIdx);
        polymerAdsorption_ = plyadsAdsorbedPolymer.eval(polymerConcentration_, /*extrapolate=*/true);
        if (PolymerModule::plyrockAdsorbtionIndex(elemCtx, dofIdx, timeIdx) == PolymerModule::NoDesorption) {
            const Scalar& maxPolymerAdsorption = elemCtx.problem().maxPolymerAdsorption(elemCtx, dofIdx, timeIdx);
            polymerAdsorption_ = std::max(Evaluation(maxPolymerAdsorption) , polymerAdsorption_);
        }

        // compute resitanceFactor
        const Scalar& residualResistanceFactor = PolymerModule::plyrockResidualResistanceFactor(elemCtx, dofIdx, timeIdx);
        const Evaluation resistanceFactor = 1.0 + (residualResistanceFactor - 1.0) * polymerAdsorption_ / maxAdsorbtion;

        // compute effective viscosities
        if (!enablePolymerMolarWeight) {
            const auto& fs = asImp_().fluidState_;
            const Evaluation& muWater = fs.viscosity(waterPhaseIdx);
            const auto& viscosityMultiplier = PolymerModule::plyviscViscosityMultiplierTable(elemCtx, dofIdx, timeIdx);
            const Evaluation viscosityMixture = viscosityMultiplier.eval(polymerConcentration_, /*extrapolate=*/true) * muWater;

            // Do the Todd-Longstaff mixing
            const Scalar plymixparToddLongstaff = PolymerModule::plymixparToddLongstaff(elemCtx, dofIdx, timeIdx);
            const Evaluation viscosityPolymer = viscosityMultiplier.eval(cmax, /*extrapolate=*/true) * muWater;
            const Evaluation viscosityPolymerEffective = pow(viscosityMixture, plymixparToddLongstaff) * pow(viscosityPolymer, 1.0 - plymixparToddLongstaff);
            const Evaluation viscosityWaterEffective = pow(viscosityMixture, plymixparToddLongstaff) * pow(muWater, 1.0 - plymixparToddLongstaff);

            const Evaluation cbar = polymerConcentration_ / cmax;
            // waterViscosity / effectiveWaterViscosity
            waterViscosityCorrection_ = muWater * ((1.0 - cbar) / viscosityWaterEffective + cbar / viscosityPolymerEffective);
            // effectiveWaterViscosity / effectivePolymerViscosity
            polymerViscosityCorrection_ =  (muWater / waterViscosityCorrection_) / viscosityPolymerEffective;
        }
        else { // based on PLYVMH
            const auto& plyvmhCoefficients = PolymerModule::plyvmhCoefficients(elemCtx, dofIdx, timeIdx);
            const Scalar k_mh = plyvmhCoefficients.k_mh;
            const Scalar a_mh = plyvmhCoefficients.a_mh;
            const Scalar gamma = plyvmhCoefficients.gamma;
            const Scalar kappa = plyvmhCoefficients.kappa;

            // viscosity model based on Mark-Houwink equation and Huggins equation
            // 1000 is a emperical constant, most likely related to unit conversion
            const Evaluation intrinsicViscosity = k_mh * pow(polymerMoleWeight_ * 1000., a_mh);
            const Evaluation x = polymerConcentration_ * intrinsicViscosity;
            waterViscosityCorrection_ = 1.0 / (1.0 + gamma * (x + kappa * x * x));
            polymerViscosityCorrection_ = 1.0;
        }

        // adjust water mobility
        asImp_().mobility_[waterPhaseIdx] *= waterViscosityCorrection_ / resistanceFactor;

        // update rock properties
        polymerDeadPoreVolume_ = PolymerModule::plyrockDeadPoreVolume(elemCtx, dofIdx, timeIdx);
        polymerRockDensity_ = PolymerModule::plyrockRockDensityFactor(elemCtx, dofIdx, timeIdx);
    }

    const Evaluation& polymerConcentration() const
    { return polymerConcentration_; }

    const Evaluation& polymerMoleWeight() const
    {
        if (!enablePolymerMolarWeight)
            throw std::logic_error("polymerMoleWeight() is called but polymer milecular weight is disabled");

        return polymerMoleWeight_;
    }

    const Scalar& polymerDeadPoreVolume() const
    { return polymerDeadPoreVolume_; }

    const Evaluation& polymerAdsorption() const
    { return polymerAdsorption_; }

    const Scalar& polymerRockDensity() const
    { return polymerRockDensity_; }

    // effectiveWaterViscosity / effectivePolymerViscosity
    const Evaluation& polymerViscosityCorrection() const
    { return polymerViscosityCorrection_; }

    // waterViscosity / effectiveWaterViscosity
    const Evaluation& waterViscosityCorrection() const
    { return waterViscosityCorrection_; }


protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation polymerConcentration_;
    // polymer molecular weight
    Evaluation polymerMoleWeight_;
    Scalar polymerDeadPoreVolume_;
    Scalar polymerRockDensity_;
    Evaluation polymerAdsorption_;
    Evaluation polymerViscosityCorrection_;
    Evaluation waterViscosityCorrection_;


};

template <class TypeTag>
class BlackOilPolymerIntensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    void polymerPropertiesUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                                  unsigned scvIdx OPM_UNUSED,
                                  unsigned timeIdx OPM_UNUSED)
    { }

    const Evaluation& polymerMoleWeight() const
    { throw std::logic_error("polymerMoleWeight() called but polymer molecular weight is disabled"); }

    const Evaluation& polymerConcentration() const
    { throw std::runtime_error("polymerConcentration() called but polymers are disabled"); }

    const Evaluation& polymerDeadPoreVolume() const
    { throw std::runtime_error("polymerDeadPoreVolume() called but polymers are disabled"); }

    const Evaluation& polymerAdsorption() const
    { throw std::runtime_error("polymerAdsorption() called but polymers are disabled"); }

    const Evaluation& polymerRockDensity() const
    { throw std::runtime_error("polymerRockDensity() called but polymers are disabled"); }

    const Evaluation& polymerViscosityCorrection() const
    { throw std::runtime_error("polymerViscosityCorrection() called but polymers are disabled"); }

    const Evaluation& waterViscosityCorrection() const
    { throw std::runtime_error("waterViscosityCorrection() called but polymers are disabled"); }
};


/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilPolymerExtensiveQuantities
 *
 * \brief Provides the polymer specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enablePolymerV = GET_PROP_VALUE(TypeTag, EnablePolymer)>
class BlackOilPolymerExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr unsigned waterPhaseIdx =  FluidSystem::waterPhaseIdx;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef BlackOilPolymerModule<TypeTag> PolymerModule;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> DimEvalVector;

public:
    /*!
     * \brief Method which calculates the shear factor based on flow velocity
     *
     * This is the variant of the method which assumes that the problem is specified
     * using permeabilities, i.e., *not* via transmissibilities.
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
                                  // compiler errors if it is not instantiated!
    void updateShearMultipliersPerm(const ElementContext& elemCtx OPM_UNUSED,
                                    unsigned scvfIdx OPM_UNUSED,
                                    unsigned timeIdx OPM_UNUSED)
    {
        throw std::runtime_error("The extension of the blackoil model for polymers is not yet "
                                 "implemented for problems specified using permeabilities.");
    }

    /*!
     * \brief Method which calculates the shear factor based on flow velocity
     *
     * This is the variant of the method which assumes that the problem is specified
     * using transmissibilities, i.e., *not* via permeabilities.
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
    // compiler errors if it is not instantiated!
    void updateShearMultipliers(const ElementContext& elemCtx,
                                unsigned scvfIdx,
                                unsigned timeIdx)
    {

        waterShearFactor_ = 1.0;
        polymerShearFactor_ = 1.0;

        if (!PolymerModule::hasPlyshlog())
            return;

        const ExtensiveQuantities& extQuants = asImp_();
        unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        unsigned interiorDofIdx = extQuants.interiorIndex();
        unsigned exteriorDofIdx = extQuants.exteriorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);

        // compute water velocity from flux
        Evaluation poroAvg = intQuantsIn.porosity()*0.5 + intQuantsEx.porosity()*0.5;
        unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvfIdx, timeIdx);
        const Evaluation& Sw = up.fluidState().saturation(waterPhaseIdx);
        unsigned cellIdx = elemCtx.globalSpaceIndex(scvfIdx, timeIdx);
        const auto& materialLawManager = elemCtx.problem().materialLawManager();
        const auto& scaledDrainageInfo =
                materialLawManager->oilWaterScaledEpsInfoDrainage(cellIdx);
        const Scalar& Swcr = scaledDrainageInfo.Swcr;

        // guard against zero porosity and no mobile water
        Evaluation denom = Opm::max(poroAvg * (Sw - Swcr), 1e-12);
        Evaluation waterVolumeVelocity = extQuants.volumeFlux(waterPhaseIdx) / denom;

        // if shrate is specified. Compute shrate based on the water velocity
        if (PolymerModule::hasShrate()) {
            const Evaluation& relWater = up.relativePermeability(waterPhaseIdx);
            Scalar trans = elemCtx.problem().transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
            if (trans > 0.0) {
                Scalar faceArea = elemCtx.stencil(timeIdx).interiorFace(scvfIdx).area();
                auto dist = elemCtx.pos(interiorDofIdx, timeIdx) -  elemCtx.pos(exteriorDofIdx, timeIdx);
                // compute permeability from transmissibility.
                Scalar absPerm = trans / faceArea * dist.two_norm();
                waterVolumeVelocity *=
                    PolymerModule::shrate(pvtnumRegionIdx)*Opm::sqrt(poroAvg*Sw / (relWater*absPerm));
                assert(Opm::isfinite(waterVolumeVelocity));
            }
        }

        // compute share factors for water and polymer
        waterShearFactor_ =
            PolymerModule::computeShearFactor(up.polymerConcentration(),
                                              pvtnumRegionIdx,
                                              waterVolumeVelocity);
        polymerShearFactor_ =
            PolymerModule::computeShearFactor(up.polymerConcentration(),
                                              pvtnumRegionIdx,
                                              waterVolumeVelocity*up.polymerViscosityCorrection());

    }

    const Evaluation& polymerShearFactor() const
    { return polymerShearFactor_; }

    const Evaluation& waterShearFactor() const
    { return waterShearFactor_; }


private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation polymerShearFactor_;
    Evaluation waterShearFactor_;

};

template <class TypeTag>
class BlackOilPolymerExtensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

public:
    void updateShearMultipliers(const ElementContext& elemCtx OPM_UNUSED,
                                unsigned scvfIdx OPM_UNUSED,
                                unsigned timeIdx OPM_UNUSED)
    { }

    void updateShearMultipliersPerm(const ElementContext& elemCtx OPM_UNUSED,
                                    unsigned scvfIdx OPM_UNUSED,
                                    unsigned timeIdx OPM_UNUSED)
    { }

    const Evaluation& polymerShearFactor() const
    { throw std::runtime_error("polymerShearFactor() called but polymers are disabled"); }

    const Evaluation& waterShearFactor() const
    { throw std::runtime_error("waterShearFactor() called but polymers are disabled"); }
};


} // namespace Opm

#endif
