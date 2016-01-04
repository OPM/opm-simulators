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
 * \brief A fluid system which uses the black-oil parameters
 *        to calculate termodynamically meaningful quantities.
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
         */
        unsigned regionIndex() const
        { return regionIdx_; }

        /*!
         * \brief Set the index of the region which should be used to determine the
         *        thermodynamic properties
         */
        void setRegionIndex(unsigned val)
        { regionIdx_ = val; }

    private:
        unsigned regionIdx_;
    };

    /****************************************
     * Fluid phase parameters
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

        gasPvt_->initEnd();
        oilPvt_->initEnd();
        waterPvt_->initEnd();

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

    //! \copydoc BaseFluidSystem::phaseName
    static const char *phaseName(const unsigned phaseIdx)
    {
        static const char *name[] = { "water", "oil", "gas" };

        assert(0 <= phaseIdx && phaseIdx < numPhases + 1);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(const unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gasPhaseIdx;
    }

    /****************************************
     * Component related parameters
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
     * thermodynamic relations
     ****************************************/
    //! \copydoc BaseFluidSystem::density
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval density(const FluidState &fluidState,
                           ParameterCache &paramCache,
                           const unsigned phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        typedef typename FluidState::Scalar FsEval;
        typedef Opm::MathToolbox<FsEval> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));
        unsigned regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case waterPhaseIdx: return waterDensity<LhsEval>(T, p, regionIdx);
        case gasPhaseIdx: {
            if (fluidState.saturation(oilPhaseIdx) > 0.0)
                return gasPvt_->saturatedDensity(regionIdx, T, p);
            else {
                // undersaturated oil
                const auto& Rv = getRv_<LhsEval>(fluidState, regionIdx);
                return gasPvt_->density(regionIdx, T, p, Rv);
            }
        }
        case oilPhaseIdx: {
            if (fluidState.saturation(gasPhaseIdx) > 0.0)
                return oilPvt_->saturatedDensity(regionIdx, T, p);
            else {
                // undersaturated gas
                const auto& Rs = getRs_<LhsEval>(fluidState, regionIdx);
                return oilPvt_->density(regionIdx, T, p, Rs);
            }
        }
        }

        OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval fugacityCoefficient(const FluidState &fluidState,
                                       const ParameterCache &paramCache,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);
        assert(0 <= compIdx  && compIdx <= numComponents);

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));
        unsigned regionIdx = paramCache.regionIndex();

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

            // for the oil component, we calculate the Rv value for saturated oil,
            // convert that to a mole fraction and convert that to the corresponding
            // fugacity coefficient
            case oilCompIdx: {
                const auto& R_vSat = gasPvt_->saturatedOilVaporizationFactor(regionIdx, T, p);
                const auto& X_gOSat = convertRvToXgO(R_vSat, regionIdx);
                const auto& x_gOSat = convertXgOToxgO(X_gOSat, regionIdx);
                return phi_oO / x_gOSat;
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
                const auto& R_sSat = oilPvt_->saturatedGasDissolutionFactor(regionIdx, T, p);
                const auto& X_oGSat = convertRsToXoG(R_sSat, regionIdx);
                const auto& x_oGSat = convertXoGToxoG(X_oGSat, regionIdx);
                return phi_gG / x_oGSat;
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
                             const ParameterCache &paramCache,
                             unsigned phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));
        unsigned regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case oilPhaseIdx: {
            if (fluidState.saturation(gasPhaseIdx) > 0.0)
                return oilPvt_->saturatedViscosity(regionIdx, T, p);
            else {
                // undersaturated oil
                const auto& Rs = getRs_<LhsEval>(fluidState, regionIdx);
                return oilPvt_->viscosity(regionIdx, T, p, Rs);
            }
        }
        case waterPhaseIdx:
            return waterPvt_->viscosity(regionIdx, T, p);
        case gasPhaseIdx: {
            if (fluidState.saturation(oilPhaseIdx) > 0.0)
                return gasPvt_->saturatedViscosity(regionIdx, T, p);
            else {
                // undersaturated gas
                const auto& Rv = getRv_<LhsEval>(fluidState, regionIdx);
                return gasPvt_->viscosity(regionIdx, T, p, Rv);
            }
        }
        }

        OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
    }

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

    /*!
     * \brief Returns the oil formation volume factor \f$B_o\f$ of saturated oil for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    template <class LhsEval>
    static LhsEval saturatedOilFormationVolumeFactor(const LhsEval& temperature,
                                                     const LhsEval& pressure,
                                                     unsigned regionIdx)
    { return oilPvt_->saturatedFormationVolumeFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the gas formation volume factor \f$B_o\f$ of saturated gas for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    template <class LhsEval>
    static LhsEval saturatedGasFormationVolumeFactor(const LhsEval& temperature,
                                                     const LhsEval& pressure,
                                                     unsigned regionIdx)
    { return gasPvt_->saturatedFormationVolumeFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Return the formation volume factor of water.
     */
    template <class LhsEval>
    static LhsEval waterFormationVolumeFactor(const LhsEval& temperature,
                                              const LhsEval& pressure,
                                              unsigned regionIdx)
    { return waterPvt_->formationVolumeFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    template <class LhsEval>
    static LhsEval gasDissolutionFactor(const LhsEval& temperature,
                                        const LhsEval& pressure,
                                        unsigned regionIdx)
    { return oilPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    template <class LhsEval>
    static LhsEval oilVaporizationFactor(const LhsEval& temperature,
                                         const LhsEval& pressure,
                                         unsigned regionIdx)
    { return gasPvt_->saturatedOilVaporizationFactor(regionIdx, temperature, pressure); }

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
     * \brief Returns the saturation pressure of the oil phase [Pa] depending on its mass
     *        fraction of the gas component
     *
     * \param Rs The surface volume of gas component dissolved in what will yield one cubic meter of oil at the surface [-]
     */
    template <class LhsEval>
    static LhsEval oilSaturationPressure(const LhsEval& temperature,
                                         const LhsEval& Rs,
                                         unsigned regionIdx)
    { return oilPvt_->saturationPressure(regionIdx, temperature, Rs); }

    /*!
     * \brief The maximum mass fraction of the gas component in the oil phase.
     */
    template <class LhsEval>
    static LhsEval saturatedOilGasMassFraction(const LhsEval& temperature,
                                               const LhsEval& pressure,
                                               unsigned regionIdx)
    {
        const auto& Rs = oilPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure);
        return convertRsToXoG(Rs, regionIdx);
    }

    /*!
     * \brief Convert an oil vaporization factor to the corresponding mass fraction
     *        of the oil component in the gas phase.
     */
    template <class LhsEval>
    static LhsEval saturatedOilGasMoleFraction(const LhsEval& temperature,
                                               const LhsEval& pressure,
                                               unsigned regionIdx)
    {
        const auto& Rs = oilPvt_->saturatedGasDissolutionFactor(regionIdx, temperature, pressure);
        const auto& XoG = convertRsToXoG(Rs, regionIdx);
        return convertXoGToxoG(XoG, regionIdx);
    }

    /*!
     * \brief Return a reference to the low-level object which calculates the gas phase
     *        quantities.
     *
     * \param Rv The surface volume of oil component dissolved in what will yield one cubic meter of gas at the surface [-]
     */
    template <class LhsEval>
    static LhsEval gasSaturationPressure(const LhsEval& temperature,
                                         const LhsEval& Rv,
                                         unsigned regionIdx)
    { return gasPvt_->saturationPressure(regionIdx, temperature, Rv); }

    /*!
     * \brief The maximum mass fraction of the oil component in the gas phase.
     */
    template <class LhsEval>
    static LhsEval saturatedGasOilMassFraction(const LhsEval& temperature,
                                               const LhsEval& pressure,
                                               unsigned regionIdx)
    {
        const auto& Rv = gasPvt_->saturatedOilVaporizationFactor(regionIdx, temperature, pressure);
        return convertRvToXgO(Rv, regionIdx);
    }

    /*!
     * \brief The maximum mole fraction of the oil component in the gas phase.
     */
    template <class LhsEval>
    static LhsEval saturatedGasOilMoleFraction(const LhsEval& temperature,
                                               const LhsEval& pressure,
                                               unsigned regionIdx)
    {
        const auto& Rv = gasPvt_->saturatedOilVaporizationFactor(regionIdx, temperature, pressure);
        const auto& XgO = convertRvToXgO(Rv, regionIdx);
        return convertXgOToxgO(XgO, regionIdx);
    }

    /*!
     * \brief Return the normalized formation volume factor of (potentially)
     *        under-saturated oil.
     */
    template <class LhsEval>
    static LhsEval oilFormationVolumeFactor(const LhsEval& temperature,
                                            const LhsEval& pressure,
                                            const LhsEval& Rs,
                                            unsigned regionIdx)
    { return oilPvt_->formationVolumeFactor(regionIdx, temperature, pressure, Rs); }

    /*!
     * \brief Return the density of (potentially) under-saturated oil.
     */
    template <class LhsEval>
    static LhsEval oilDensity(const LhsEval& temperature,
                              const LhsEval& pressure,
                              const LhsEval& Rs,
                              unsigned regionIdx)
    { return oilPvt_->density(regionIdx, temperature, pressure, Rs); }

    /*!
     * \brief Return the density of gas-saturated oil.
     */
    template <class LhsEval>
    static LhsEval saturatedOilDensity(const LhsEval& temperature,
                                       const LhsEval& pressure,
                                       unsigned regionIdx)
    { return oilPvt_->saturatedDensity(regionIdx, temperature, pressure); }

    /*!
     * \brief Return the formation volume factor of gas.
     */
    template <class LhsEval>
    static LhsEval gasFormationVolumeFactor(const LhsEval& temperature,
                                            const LhsEval& pressure,
                                            const LhsEval& Rv,
                                            unsigned regionIdx)
    { return gasPvt_->formationVolumeFactor(regionIdx, temperature, pressure, Rv); }

    /*!
     * \brief Return the density of dry gas.
     */
    template <class LhsEval>
    static LhsEval gasDensity(const LhsEval& temperature,
                              const LhsEval& pressure,
                              const LhsEval& Rv,
                              unsigned regionIdx)
    { return gasPvt_->density(regionIdx, temperature, pressure, Rv); }

    /*!
     * \brief Return the density of gas-saturated oil.
     */
    template <class LhsEval>
    static LhsEval saturatedGasDensity(const LhsEval& temperature,
                                       const LhsEval& pressure,
                                       unsigned regionIdx)
    { return oilPvt_->saturatedDensity(regionIdx, temperature, pressure);  }

    /*!
     * \brief Return the density of water.
     */
    template <class LhsEval>
    static LhsEval waterDensity(const LhsEval& temperature,
                                const LhsEval& pressure,
                                unsigned regionIdx)
    { return waterPvt_->density(regionIdx, temperature, pressure); }

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
