// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2015 by Andreas Lauser

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
*/
/*!
 * \file
 * \copydoc Opm::FluidSystems::BlackOil
 */
#ifndef OPM_BLACK_OIL_FLUID_SYSTEM_HPP
#define OPM_BLACK_OIL_FLUID_SYSTEM_HPP

#include "blackoilpvt/OilPvtMultiplexer.hpp"
#include "blackoilpvt/GasPvtMultiplexer.hpp"
#include "blackoilpvt/WaterPvtMultiplexer.hpp"

#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/Constants.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <memory>
#include <vector>
#include <array>

namespace Opm {
namespace FluidSystems {
/*!
 * \brief A fluid system which uses the black-oil model assumptions to calculate
 *        termodynamically meaningful quantities.
 *
 * \tparam Scalar The type used for scalar floating point values
 */
template <class Scalar>
class BlackOil : public BaseFluidSystem<Scalar, BlackOil<Scalar> >
{
public:
    typedef Opm::GasPvtMultiplexer<Scalar> GasPvt;
    typedef Opm::OilPvtMultiplexer<Scalar> OilPvt;
    typedef Opm::WaterPvtMultiplexer<Scalar> WaterPvt;

    //! \copydoc BaseFluidSystem::ParameterCache
    class ParameterCache : public Opm::NullParameterCache
    {
    public:
        ParameterCache(int /*regionIdx*/=0)
        { regionIdx_ = 0; }

        /*!
         * \brief Return the index of the region which should be used to determine the
         *        thermodynamic properties
         *
         * This is only required because "oil" and "gas" are pseudo-components, i.e. for
         * more comprehensive equations of state there would only be one "region".
         */
        unsigned regionIndex() const
        { return regionIdx_; }

        /*!
         * \brief Set the index of the region which should be used to determine the
         *        thermodynamic properties
         *
         * This is only required because "oil" and "gas" are pseudo-components, i.e. for
         * more comprehensive equations of state there would only be one "region".
         */
        void setRegionIndex(unsigned val)
        { regionIdx_ = val; }

    private:
        unsigned regionIdx_;
    };

    /****************************************
     * Initialization
     ****************************************/
#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize the fluid system using an ECL deck object
     */
    static void initFromDeck(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        auto densityKeyword = deck->getKeyword("DENSITY");
        size_t numRegions = densityKeyword->size();
        initBegin(numRegions);

        setEnableDissolvedGas(deck->hasKeyword("DISGAS"));
        setEnableVaporizedOil(deck->hasKeyword("VAPOIL"));

        // set the reference densities of all PVT regions
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            Opm::DeckRecordConstPtr densityRecord = densityKeyword->getRecord(regionIdx);
            setReferenceDensities(densityRecord->getItem("OIL")->getSIDouble(0),
                                  densityRecord->getItem("WATER")->getSIDouble(0),
                                  densityRecord->getItem("GAS")->getSIDouble(0),
                                  regionIdx);
        }

        gasPvt_ = std::make_shared<GasPvt>();
        gasPvt_->initFromDeck(deck, eclState);

        oilPvt_ = std::make_shared<OilPvt>();
        oilPvt_->initFromDeck(deck, eclState);

        waterPvt_ = std::make_shared<WaterPvt>();
        waterPvt_->initFromDeck(deck, eclState);

        initEnd();
    }
#endif // HAVE_OPM_PARSER

    /*!
     * \brief Begin the initialization of the black oil fluid system.
     *
     * After calling this method the reference densities, all dissolution and formation
     * volume factors, the oil bubble pressure, all viscosities and the water
     * compressibility must be set. Before the fluid system can be used, initEnd() must
     * be called to finalize the initialization.
     */
    static void initBegin(size_t numPvtRegions)
    {
        enableDissolvedGas_ = true;
        enableVaporizedOil_ = false;

        resizeArrays_(numPvtRegions);
    }

    /*!
     * \brief Specify whether the fluid system should consider that the gas component can
     *        dissolve in the oil phase
     *
     * By default, dissolved gas is considered.
     */
    static void setEnableDissolvedGas(bool yesno)
    { enableDissolvedGas_ = yesno; }

    /*!
     * \brief Specify whether the fluid system should consider that the oil component can
     *        dissolve in the gas phase
     *
     * By default, vaporized oil is not considered.
     */
    static void setEnableVaporizedOil(bool yesno)
    { enableVaporizedOil_ = yesno; }

    /*!
     * \brief Set the pressure-volume-saturation (PVT) relations for the gas phase.
     */
    static void setGasPvt(std::shared_ptr<GasPvt> pvtObj)
    { gasPvt_ = pvtObj; }

    /*!
     * \brief Set the pressure-volume-saturation (PVT) relations for the oil phase.
     */
    static void setOilPvt(std::shared_ptr<OilPvt> pvtObj)
    { oilPvt_ = pvtObj; }

    /*!
     * \brief Set the pressure-volume-saturation (PVT) relations for the water phase.
     */
    static void setWaterPvt(std::shared_ptr<WaterPvt> pvtObj)
    { waterPvt_ = pvtObj; }

    /*!
     * \brief Initialize the values of the reference densities
     *
     * \param rhoOil The reference density of (gas saturated) oil phase.
     * \param rhoWater The reference density of the water phase.
     * \param rhoGas The reference density of the gas phase.
     */
    static void setReferenceDensities(Scalar rhoOil,
                                      Scalar rhoWater,
                                      Scalar rhoGas,
                                      unsigned regionIdx)
    {
        referenceDensity_[regionIdx][oilPhaseIdx] = rhoOil;
        referenceDensity_[regionIdx][waterPhaseIdx] = rhoWater;
        referenceDensity_[regionIdx][gasPhaseIdx] = rhoGas;
    }

    /*!
     * \brief Finish initializing the black oil fluid system.
     */
    static void initEnd()
    {
        // calculate the final 2D functions which are used for interpolation.
        size_t numRegions = molarMass_.size();
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate molar masses

            // water is simple: 18 g/mol
            molarMass_[regionIdx][waterCompIdx] = 18e-3;

            // for gas, we take the density at standard conditions and assume it to be ideal
            Scalar p = surfacePressure;
            Scalar T = surfaceTemperature;
            Scalar rho_g = referenceDensity_[/*regionIdx=*/0][gasPhaseIdx];
            molarMass_[regionIdx][gasCompIdx] = Opm::Constants<Scalar>::R*T*rho_g / p;

            // finally, for oil phase, we take the molar mass from the
            // spe9 paper
            molarMass_[regionIdx][oilCompIdx] = 175e-3; // kg/mol
        }
    }

    /****************************************
     * Generic phase properties
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 3;

    //! Index of the water phase
    static const int waterPhaseIdx = 0;
    //! Index of the oil phase
    static const int oilPhaseIdx = 1;
    //! Index of the gas phase
    static const int gasPhaseIdx = 2;

    //! The pressure at the surface
    static const Scalar surfacePressure;

    //! The temperature at the surface
    static const Scalar surfaceTemperature;

    //! \copydoc BaseFluidSystem::phaseName
    static const char *phaseName(unsigned phaseIdx)
    {
        static const char *name[] = { "water", "oil", "gas" };

        assert(0 <= phaseIdx && phaseIdx < numPhases + 1);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gasPhaseIdx;
    }

    /****************************************
     * Generic component related properties
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 3;

    //! Index of the oil component
    static const int oilCompIdx = 0;
    //! Index of the water component
    static const int waterCompIdx = 1;
    //! Index of the gas component
    static const int gasCompIdx = 2;

    //! \copydoc BaseFluidSystem::componentName
    static const char *componentName(unsigned compIdx)
    {
        static const char *name[] = { "Oil", "Water", "Gas" };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx, unsigned regionIdx = 0)
    { return molarMass_[regionIdx][compIdx]; }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        // fugacity coefficients are only pressure dependent -> we
        // have an ideal mixture
        return true;
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(unsigned /*phaseIdx*/)
    { return true; /* all phases are compressible */ }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(unsigned /*phaseIdx*/)
    { return false; }


    /****************************************
     * Black-oil specific properties
     ****************************************/
    /*!
     * \brief Returns the number of PVT regions which are considered.
     *
     * By default, this is 1.
     */
    static size_t numRegions()
    { return molarMass_.size(); }

    /*!
     * \brief Returns whether the fluid system should consider that the gas component can
     *        dissolve in the oil phase
     *
     * By default, dissolved gas is considered.
     */
    static bool enableDissolvedGas()
    { return enableDissolvedGas_; }

    /*!
     * \brief Returns whether the fluid system should consider that the oil component can
     *        dissolve in the gas phase
     *
     * By default, vaporized oil is not considered.
     */
    static bool enableVaporizedOil()
    { return enableVaporizedOil_; }

    /*!
     * \brief Returns the density of a fluid phase at surface pressure [kg/m^3]
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    static Scalar referenceDensity(unsigned phaseIdx, unsigned regionIdx)
    { return referenceDensity_[regionIdx][phaseIdx]; }

    /****************************************
     * thermodynamic quantities (generic version, only isothermal)
     ****************************************/
    //! \copydoc BaseFluidSystem::density
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval density(const FluidState &fluidState,
                           ParameterCache &paramCache,
                           unsigned phaseIdx)
    { return density<FluidState, LhsEval>(fluidState, phaseIdx, paramCache.regionIndex()); }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval fugacityCoefficient(const FluidState &fluidState,
                                       const ParameterCache &paramCache,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    { return fugacityCoefficient<FluidState, LhsEval>(fluidState, phaseIdx, compIdx, paramCache.regionIndex()); }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval viscosity(const FluidState &fluidState,
                             const ParameterCache &paramCache,
                             unsigned phaseIdx)
    { return viscosity<FluidState, LhsEval>(fluidState, phaseIdx, paramCache.regionIndex()); }


    /****************************************
     * thermodynamic quantities (black-oil specific version: Note that the PVT region
     * index is explicitly passed instead of a parameter cache object)
     ****************************************/
    //! \copydoc BaseFluidSystem::density
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval density(const FluidState &fluidState,
                           unsigned phaseIdx,
                           unsigned regionIdx)
    {
        assert(0 <= phaseIdx && phaseIdx <= numPhases);
        assert(0 <= regionIdx && regionIdx <= numRegions());

        switch (phaseIdx) {
        case oilPhaseIdx: {
            if (!enableDissolvedGas()) {
                // immiscible oil
                const auto& bo = inverseFormationVolumeFactor<FluidState, LhsEval>(fluidState, oilPhaseIdx, regionIdx);
                return referenceDensity(phaseIdx, regionIdx)*bo;
            }

            // miscible oil
            const auto& bo = inverseFormationVolumeFactor<FluidState, LhsEval>(fluidState, oilPhaseIdx, regionIdx);
            const auto& Rs = getRs_<LhsEval, FluidState>(fluidState, regionIdx);

            return
                bo*referenceDensity(oilPhaseIdx, regionIdx)
                + Rs*bo*referenceDensity(gasPhaseIdx, regionIdx);
        }

        case gasPhaseIdx: {
            if (!enableVaporizedOil()) {
                // immiscible gas
                const auto& bg = inverseFormationVolumeFactor<FluidState, LhsEval>(fluidState, gasPhaseIdx, regionIdx);
                return bg*referenceDensity(phaseIdx, regionIdx);
            }

            // miscible gas
            const auto& bg = inverseFormationVolumeFactor<FluidState, LhsEval>(fluidState, gasPhaseIdx, regionIdx);
            const auto& Rv = getRv_<LhsEval, FluidState>(fluidState, regionIdx);

            return
                bg*referenceDensity(gasPhaseIdx, regionIdx)
                + Rv*bg*referenceDensity(oilPhaseIdx, regionIdx);
        }

        case waterPhaseIdx:
            return
                referenceDensity(waterPhaseIdx, regionIdx)
                *inverseFormationVolumeFactor<FluidState, LhsEval>(fluidState, waterPhaseIdx, regionIdx);
        }

        OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
    }

    /*!
     * \brief Compute the density of a saturated fluid phase.
     *
     * This means the density of the given fluid phase if the dissolved component (gas
     * for the oil phase and oil for the gas phase) is at the thermodynamically possible
     * maximum. For the water phase, there's no difference to the density() method.
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval saturatedDensity(const FluidState &fluidState,
                                    unsigned phaseIdx,
                                    unsigned regionIdx)
    {
        assert(0 <= phaseIdx && phaseIdx <= numPhases);
        assert(0 <= regionIdx && regionIdx <= numRegions());

        switch (phaseIdx) {
        case oilPhaseIdx: {
            if (!enableDissolvedGas()) {
                // immiscible oil
                const auto& bo = inverseFormationVolumeFactor<FluidState, LhsEval>(fluidState, oilPhaseIdx, regionIdx);
                return referenceDensity(phaseIdx, regionIdx)*bo;
            }

            // miscible oil
            const auto& bo = saturatedInverseFormationVolumeFactor<FluidState, LhsEval>(fluidState, oilPhaseIdx, regionIdx);
            const auto& Rs = saturatedDissolutionFactor<FluidState, LhsEval>(fluidState, oilPhaseIdx, regionIdx);

            return
                bo*referenceDensity(oilPhaseIdx, regionIdx)
                + Rs*bo*referenceDensity(gasPhaseIdx, regionIdx);
        }

        case gasPhaseIdx: {
            if (!enableVaporizedOil()) {
                // immiscible gas
                const auto& bg = inverseFormationVolumeFactor<FluidState, LhsEval>(fluidState, gasPhaseIdx, regionIdx);
                return referenceDensity(phaseIdx, regionIdx)*bg;
            }

            // miscible gas
            const auto& bg = saturatedInverseFormationVolumeFactor<FluidState, LhsEval>(fluidState, gasPhaseIdx, regionIdx);
            const auto& Rv = saturatedDissolutionFactor<FluidState, LhsEval>(fluidState, gasPhaseIdx, regionIdx);

            return
                bg*referenceDensity(gasPhaseIdx, regionIdx)
                + Rv*bg*referenceDensity(oilPhaseIdx, regionIdx);
        }

        case waterPhaseIdx:
            return
                referenceDensity(waterPhaseIdx, regionIdx)
                *saturatedInverseFormationVolumeFactor<FluidState, LhsEval>(fluidState, waterPhaseIdx, regionIdx);
        }

        OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the formation volume factor \f$B_\alpha\f$ of an "undersaturated"
     *        fluid phase
     *
     * For the oil (gas) phase, "undersaturated" means that the concentration of the gas
     * (oil) component is not assumed to be at the thermodynamically possible maximum at
     * the given temperature and pressure.
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval inverseFormationVolumeFactor(const FluidState& fluidState,
                                                unsigned phaseIdx,
                                                unsigned regionIdx)
    {
        assert(0 <= phaseIdx && phaseIdx <= numPhases);
        assert(0 <= regionIdx && regionIdx <= numRegions());

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));

        switch (phaseIdx) {
        case oilPhaseIdx: {
            if (enableDissolvedGas()) {
                if (fluidState.saturation(gasPhaseIdx) > 0.0) {
                    if (fluidState.saturation(gasPhaseIdx) < 1e-4) {
                        // here comes the relatively expensive case: first calculate and then
                        // interpolate between the saturated and undersaturated quantities to
                        // avoid a discontinuity
                        const auto& Rs = getRs_<LhsEval, FluidState>(fluidState, regionIdx);
                        const auto& alpha = FsToolbox::template toLhs<LhsEval>(fluidState.saturation(gasPhaseIdx))/1e-4;
                        const auto& bSat = oilPvt_->saturatedInverseFormationVolumeFactor(regionIdx, T, p);
                        const auto& bUndersat = oilPvt_->inverseFormationVolumeFactor(regionIdx, T, p, Rs);
                        return alpha*bSat + (1.0 - alpha)*bUndersat;
                    }

                    return oilPvt_->saturatedInverseFormationVolumeFactor(regionIdx, T, p);
                }

                const auto& Rs = getRs_<LhsEval, FluidState>(fluidState, regionIdx);
                return oilPvt_->inverseFormationVolumeFactor(regionIdx, T, p, Rs);
            }

            const LhsEval Rs(0.0);
            return oilPvt_->inverseFormationVolumeFactor(regionIdx, T, p, Rs);
        }
        case gasPhaseIdx: {
            if (enableVaporizedOil()) {
                if (fluidState.saturation(oilPhaseIdx) > 0.0) {
                    if (fluidState.saturation(oilPhaseIdx) < 1e-4) {
                        // here comes the relatively expensive case: first calculate and then
                        // interpolate between the saturated and undersaturated quantities to
                        // avoid a discontinuity
                        const auto& Rv = getRv_<LhsEval, FluidState>(fluidState, regionIdx);
                        const auto& alpha = FsToolbox::template toLhs<LhsEval>(fluidState.saturation(oilPhaseIdx))/1e-4;
                        const auto& bSat = gasPvt_->saturatedInverseFormationVolumeFactor(regionIdx, T, p);
                        const auto& bUndersat = gasPvt_->inverseFormationVolumeFactor(regionIdx, T, p, Rv);
                        return alpha*bSat + (1.0 - alpha)*bUndersat;
                    }

                    return gasPvt_->saturatedInverseFormationVolumeFactor(regionIdx, T, p);
                }

                const auto& Rv = getRv_<LhsEval, FluidState>(fluidState, regionIdx);
                return gasPvt_->inverseFormationVolumeFactor(regionIdx, T, p, Rv);
            }

            const LhsEval Rv(0.0);
            return gasPvt_->inverseFormationVolumeFactor(regionIdx, T, p, Rv);
        }
        case waterPhaseIdx:
            return waterPvt_->inverseFormationVolumeFactor(regionIdx, T, p);
        default: OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
        }
    }

    /*!
     * \brief Returns the formation volume factor \f$B_\alpha\f$ of a "saturated" fluid
     *        phase
     *
     * For the oil phase, this means that it is gas saturated, the gas phase is oil
     * saturated and for the water phase, there is no difference to formationVolumeFactor()
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval saturatedInverseFormationVolumeFactor(const FluidState& fluidState,
                                                         unsigned phaseIdx,
                                                         unsigned regionIdx)
    {
        assert(0 <= phaseIdx && phaseIdx <= numPhases);
        assert(0 <= regionIdx && regionIdx <= numRegions());

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));

        switch (phaseIdx) {
        case oilPhaseIdx: return oilPvt_->saturatedInverseFormationVolumeFactor(regionIdx, T, p);
        case gasPhaseIdx: return gasPvt_->saturatedInverseFormationVolumeFactor(regionIdx, T, p);
        case waterPhaseIdx: return waterPvt_->inverseFormationVolumeFactor(regionIdx, T, p);
        default: OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
        }
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval fugacityCoefficient(const FluidState &fluidState,
                                       unsigned phaseIdx,
                                       unsigned compIdx,
                                       unsigned regionIdx)
    {
        assert(0 <= phaseIdx && phaseIdx <= numPhases);
        assert(0 <= compIdx && compIdx <= numComponents);
        assert(0 <= regionIdx && regionIdx <= numRegions());

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));

        // for the fugacity coefficient of the oil component in the oil phase, we use
        // some pseudo-realistic value for the vapor pressure to ease physical
        // interpretation of the results
        const LhsEval phi_oO = 20e3/p;

        // for the gas component in the gas phase, assume it to be an ideal gas
        const Scalar phi_gG = 1.0;

        // for the fugacity coefficient of the water component in the water phase, we use
        // the same approach as for the oil component in the oil phase
        const LhsEval phi_wW = 30e3/p;

        switch (phaseIdx) {
        case gasPhaseIdx: // fugacity coefficients for all components in the gas phase
            switch (compIdx) {
            case gasCompIdx:
                return phi_gG;

            // for the oil component, we calculate the Rv value for saturated gas and Rs
            // for saturated oil, and compute the fugacity coefficients at the
            // equilibrium. for this, we assume that in equilibrium the fugacities of the
            // oil component is the same in both phases.
            case oilCompIdx: {
                if (!enableVaporizedOil())
                    // if there's no vaporized oil, the gas phase is assumed to be
                    // immiscible with the oil component
                    return phi_gG*1e6;

                const auto& R_vSat = gasPvt_->saturatedOilVaporizationFactor(regionIdx, T, p);
                const auto& X_gOSat = convertRvToXgO(R_vSat, regionIdx);
                const auto& x_gOSat = convertXgOToxgO(X_gOSat, regionIdx);

                const auto& R_sSat = oilPvt_->saturatedGasDissolutionFactor(regionIdx, T, p);
                const auto& X_oGSat = convertRsToXoG(R_sSat, regionIdx);
                const auto& x_oGSat = convertXoGToxoG(X_oGSat, regionIdx);
                const auto& x_oOSat = 1.0 - x_oGSat;

                const auto& p_o = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(oilPhaseIdx));
                const auto& p_g = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(gasPhaseIdx));

                return phi_oO*p_o*x_oOSat / (p_g*x_gOSat);
            }

            case waterCompIdx:
                // the water component is assumed to be never miscible with the gas phase
                return phi_gG*1e6;

            default:
                OPM_THROW(std::logic_error,
                          "Invalid component index " << compIdx);
            }

        case oilPhaseIdx: // fugacity coefficients for all components in the oil phase
            switch (compIdx) {
            case oilCompIdx:
                return phi_oO;

            // for the oil and water components, we proceed analogous to the gas and
            // water components in the gas phase
            case gasCompIdx: {
                if (!enableDissolvedGas())
                    // if there's no dissolved gas, the oil phase is assumed to be
                    // immiscible with the gas component
                    return phi_oO*1e6;

                const auto& R_vSat = gasPvt_->saturatedOilVaporizationFactor(regionIdx, T, p);
                const auto& X_gOSat = convertRvToXgO(R_vSat, regionIdx);
                const auto& x_gOSat = convertXgOToxgO(X_gOSat, regionIdx);
                const auto& x_gGSat = 1.0 - x_gOSat;

                const auto& R_sSat = oilPvt_->saturatedGasDissolutionFactor(regionIdx, T, p);
                const auto& X_oGSat = convertRsToXoG(R_sSat, regionIdx);
                const auto& x_oGSat = convertXoGToxoG(X_oGSat, regionIdx);

                const auto& p_o = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(oilPhaseIdx));
                const auto& p_g = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(gasPhaseIdx));

                return phi_gG*p_g*x_gGSat / (p_o*x_oGSat);
            }

            case waterCompIdx:
                return phi_oO*1e6;

            default:
                OPM_THROW(std::logic_error,
                          "Invalid component index " << compIdx);
            }

        case waterPhaseIdx: // fugacity coefficients for all components in the water phase
            // the water phase fugacity coefficients are pretty simple: because the water
            // phase is assumed to consist entirely from the water component, we just
            // need to make sure that the fugacity coefficients for the other components
            // are a few orders of magnitude larger than the one of the water
            // component. (i.e., the affinity of the gas and oil components to the water
            // phase is lower by a few orders of magnitude)
            switch (compIdx) {
            case waterCompIdx: return phi_wW;
            case oilCompIdx: return 1.1e6*phi_wW;
            case gasCompIdx: return 1e6*phi_wW;
            default:
                OPM_THROW(std::logic_error,
                          "Invalid component index " << compIdx);
            }

        default:
            OPM_THROW(std::logic_error,
                      "Invalid phase index " << phaseIdx);
        }

        OPM_THROW(std::logic_error, "Unhandled phase or component index");
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval viscosity(const FluidState &fluidState,
                             unsigned phaseIdx,
                             unsigned regionIdx)
    {
        assert(0 <= phaseIdx && phaseIdx <= numPhases);
        assert(0 <= regionIdx && regionIdx <= numRegions());

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));

        switch (phaseIdx) {
        case oilPhaseIdx: {
            if (enableDissolvedGas()) {
                if (fluidState.saturation(gasPhaseIdx) > 0.0) {
                    if (fluidState.saturation(gasPhaseIdx) < 1e-4) {
                        // here comes the relatively expensive case: first calculate and then
                        // interpolate between the saturated and undersaturated quantities to
                        // avoid a discontinuity
                        const auto& Rs = getRs_<LhsEval, FluidState>(fluidState, regionIdx);
                        const auto& alpha = FsToolbox::template toLhs<LhsEval>(fluidState.saturation(gasPhaseIdx))/1e-4;
                        const auto& muSat = oilPvt_->saturatedViscosity(regionIdx, T, p);
                        const auto& muUndersat = oilPvt_->viscosity(regionIdx, T, p, Rs);
                        return alpha*muSat + (1.0 - alpha)*muUndersat;
                    }

                    return oilPvt_->saturatedViscosity(regionIdx, T, p);
                }

                const auto& Rs = getRs_<LhsEval, FluidState>(fluidState, regionIdx);
                return oilPvt_->viscosity(regionIdx, T, p, Rs);
            }

            const LhsEval Rs(0.0);
            return oilPvt_->viscosity(regionIdx, T, p, Rs);
        }

        case gasPhaseIdx: {
            if (enableVaporizedOil()) {
                if (fluidState.saturation(oilPhaseIdx) > 0.0) {
                    if (fluidState.saturation(oilPhaseIdx) < 1e-4) {
                        // here comes the relatively expensive case: first calculate and then
                        // interpolate between the saturated and undersaturated quantities to
                        // avoid a discontinuity
                        const auto& Rv = getRv_<LhsEval, FluidState>(fluidState, regionIdx);
                        const auto& alpha = FsToolbox::template toLhs<LhsEval>(fluidState.saturation(oilPhaseIdx))/1e-4;
                        const auto& muSat = gasPvt_->saturatedViscosity(regionIdx, T, p);
                        const auto& muUndersat = gasPvt_->viscosity(regionIdx, T, p, Rv);
                        return alpha*muSat + (1.0 - alpha)*muUndersat;
                    }

                    return gasPvt_->saturatedViscosity(regionIdx, T, p);
                }

                const auto& Rv = getRv_<LhsEval, FluidState>(fluidState, regionIdx);
                return gasPvt_->viscosity(regionIdx, T, p, Rv);
            }

            const LhsEval Rv(0.0);
            return gasPvt_->viscosity(regionIdx, T, p, Rv);
        }

        case waterPhaseIdx:
            // since water is always assumed to be immiscible in the black-oil model,
            // there is no "saturated water"
            return waterPvt_->viscosity(regionIdx, T, p);
        }

        OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the dissolution factor \f$R_\alpha\f$ of a saturated fluid phase
     *
     * For the oil (gas) phase, this means the R_S and R_V factors, for the water phase,
     * it is always 0.
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval saturatedDissolutionFactor(const FluidState& fluidState,
                                              unsigned phaseIdx,
                                              unsigned regionIdx)
    {
        assert(0 <= phaseIdx && phaseIdx <= numPhases);
        assert(0 <= regionIdx && regionIdx <= numRegions());

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));

        switch (phaseIdx) {
        case oilPhaseIdx: return oilPvt_->saturatedGasDissolutionFactor(regionIdx, T, p);
        case gasPhaseIdx: return gasPvt_->saturatedOilVaporizationFactor(regionIdx, T, p);
        case waterPhaseIdx: return 0.0;
        default: OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
        }
    }

    /*!
     * \brief Returns the saturation pressure of a given phase [Pa] depending on its
     *        composition.
     *
     * In the black-oil model, the saturation pressure it the pressure at which the fluid
     * phase is in equilibrium with the gas phase, i.e., it is the inverse of the
     * "dissolution factor". Note that a-priori this quantity is undefined for the water
     * phase (because water is assumed to be immiscible with everything else). This method
     * here just returns 0, though.
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval saturationPressure(const FluidState& fluidState,
                                      unsigned phaseIdx,
                                      unsigned regionIdx)
    {
        assert(0 <= phaseIdx && phaseIdx <= numPhases);
        assert(0 <= regionIdx && regionIdx <= numRegions());

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));

        switch (phaseIdx) {
        case oilPhaseIdx: return oilPvt_->saturationPressure(regionIdx, T, getRs_<LhsEval>(fluidState, regionIdx));
        case gasPhaseIdx: return gasPvt_->saturationPressure(regionIdx, T, getRv_<LhsEval>(fluidState, regionIdx));
        case waterPhaseIdx: return 0.0;
        default: OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
        }
    }

    /****************************************
     * Auxiliary and convenience methods for the black-oil model
     ****************************************/
    /*!
     * \brief Convert the mass fraction of the gas component in the oil phase to the
     *        corresponding gas dissolution factor.
     */
    template <class LhsEval>
    static LhsEval convertXoGToRs(const LhsEval& XoG, unsigned regionIdx)
    {
        Scalar rho_oRef = referenceDensity_[regionIdx][oilPhaseIdx];
        Scalar rho_gRef = referenceDensity_[regionIdx][gasPhaseIdx];

        return XoG/(1.0 - XoG)*(rho_oRef/rho_gRef);
    }

    /*!
     * \brief Convert the mass fraction of the oil component in the gas phase to the
     *        corresponding oil vaporization factor.
     */
    template <class LhsEval>
    static LhsEval convertXgOToRv(const LhsEval& XgO, unsigned regionIdx)
    {
        Scalar rho_oRef = referenceDensity_[regionIdx][oilPhaseIdx];
        Scalar rho_gRef = referenceDensity_[regionIdx][gasPhaseIdx];

        return XgO/(1.0 - XgO)*(rho_gRef/rho_oRef);
    }

    /*!
     * \brief Convert a gas dissolution factor to the the corresponding mass fraction
     *        of the gas component in the oil phase.
     */
    template <class LhsEval>
    static LhsEval convertRsToXoG(const LhsEval& Rs, unsigned regionIdx)
    {
        Scalar rho_oRef = referenceDensity_[regionIdx][oilPhaseIdx];
        Scalar rho_gRef = referenceDensity_[regionIdx][gasPhaseIdx];

        const LhsEval& rho_oG = Rs*rho_gRef;
        return rho_oG/(rho_oRef + rho_oG);
    }

    /*!
     * \brief Convert an oil vaporization factor to the corresponding mass fraction
     *        of the oil component in the gas phase.
     */
    template <class LhsEval>
    static LhsEval convertRvToXgO(const LhsEval& Rv, unsigned regionIdx)
    {
        Scalar rho_oRef = referenceDensity_[regionIdx][oilPhaseIdx];
        Scalar rho_gRef = referenceDensity_[regionIdx][gasPhaseIdx];

        const LhsEval& rho_gO = Rv*rho_oRef;
        return rho_gO/(rho_gRef + rho_gO);
    }

    /*!
     * \brief Convert a gas mass fraction in the oil phase the corresponding mole fraction.
     */
    template <class LhsEval>
    static LhsEval convertXoGToxoG(const LhsEval& XoG, unsigned regionIdx)
    {
        Scalar MO = molarMass_[regionIdx][oilCompIdx];
        Scalar MG = molarMass_[regionIdx][gasCompIdx];

        return XoG*MO / (MG*(1 - XoG) + XoG*MO);
    }

    /*!
     * \brief Convert a oil mass fraction in the gas phase the corresponding mole fraction.
     */
    template <class LhsEval>
    static LhsEval convertXgOToxgO(const LhsEval& XgO, unsigned regionIdx)
    {
        Scalar MO = molarMass_[regionIdx][oilCompIdx];
        Scalar MG = molarMass_[regionIdx][gasCompIdx];

        return XgO*MG / (MO*(1 - XgO) + XgO*MG);
    }

    /*!
     * \brief Return a reference to the low-level object which calculates the gas phase
     *        quantities.
     *
     * \note It is not recommended to use this method directly, but the black-oil
     *       specific methods of the fluid systems from above should be used instead.
     */
    static const GasPvt& gasPvt()
    { return *gasPvt_; }

    /*!
     * \brief Return a reference to the low-level object which calculates the oil phase
     *        quantities.
     *
     * \note It is not recommended to use this method directly, but the black-oil
     *       specific methods of the fluid systems from above should be used instead.
     */
    static const OilPvt& oilPvt()
    { return *oilPvt_; }

    /*!
     * \brief Return a reference to the low-level object which calculates the water phase
     *        quantities.
     *
     * \note It is not recommended to use this method directly, but the black-oil
     *       specific methods of the fluid systems from above should be used instead.
     */
    static const WaterPvt& waterPvt()
    { return *waterPvt_; }

private:
    static void resizeArrays_(size_t numRegions)
    {
        molarMass_.resize(numRegions);
        referenceDensity_.resize(numRegions);
    }

    template <class LhsEval, class FluidState>
    static LhsEval getRs_(const FluidState& fluidState, unsigned regionIdx)
    {
        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& XoG =
            FsToolbox::template toLhs<LhsEval>(fluidState.massFraction(oilPhaseIdx, gasCompIdx));
        return convertXoGToRs(XoG, regionIdx);
    }

    template <class LhsEval, class FluidState>
    static LhsEval getRv_(const FluidState& fluidState, unsigned regionIdx)
    {
        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& XgO =
            FsToolbox::template toLhs<LhsEval>(fluidState.massFraction(gasPhaseIdx, oilCompIdx));
        return convertXgOToRv(XgO, regionIdx);
    }

    static std::shared_ptr<GasPvt> gasPvt_;
    static std::shared_ptr<OilPvt> oilPvt_;
    static std::shared_ptr<WaterPvt> waterPvt_;

    static bool enableDissolvedGas_;
    static bool enableVaporizedOil_;

    // HACK for GCC 4.4: the array size has to be specified using the literal value '3'
    // here, because GCC 4.4 seems to be unable to determine the number of phases from
    // the BlackOil fluid system in the attribute declaration below...
    static std::vector<std::array<Scalar, /*numPhases=*/3> > referenceDensity_;
    static std::vector<std::array<Scalar, /*numComponents=*/3> > molarMass_;
};

template <class Scalar>
const Scalar
BlackOil<Scalar>::surfaceTemperature = 273.15 + 15.56; // [K]

template <class Scalar>
const Scalar
BlackOil<Scalar>::surfacePressure = 101325.0; // [Pa]

template <class Scalar>
bool BlackOil<Scalar>::enableDissolvedGas_;

template <class Scalar>
bool BlackOil<Scalar>::enableVaporizedOil_;

template <class Scalar>
std::shared_ptr<OilPvtMultiplexer<Scalar> >
BlackOil<Scalar>::oilPvt_;

template <class Scalar>
std::shared_ptr<Opm::GasPvtMultiplexer<Scalar> >
BlackOil<Scalar>::gasPvt_;

template <class Scalar>
std::shared_ptr<WaterPvtMultiplexer<Scalar> >
BlackOil<Scalar>::waterPvt_;

template <class Scalar>
std::vector<std::array<Scalar, 3> >
BlackOil<Scalar>::referenceDensity_;

template <class Scalar>
std::vector<std::array<Scalar, 3> >
BlackOil<Scalar>::molarMass_;
}} // namespace Opm, FluidSystems

#endif
