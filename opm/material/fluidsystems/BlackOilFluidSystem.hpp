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

#include "blackoilpvt/OilPvtInterface.hpp"
#include "blackoilpvt/GasPvtInterface.hpp"
#include "blackoilpvt/WaterPvtInterface.hpp"

#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/Constants.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/ErrorMacros.hpp>

#include <memory>
#include <vector>
#include <array>

namespace Opm {
namespace FluidSystems {
/*!
 * \brief A fluid system which uses the black-oil parameters
 *        to calculate termodynamically meaningful quantities.
 */
template <class Scalar, class Evaluation = Scalar>
class BlackOil : public BaseFluidSystem<Scalar, BlackOil<Scalar, Evaluation> >
{
    typedef Opm::GasPvtInterface<Scalar, Evaluation> GasPvtInterface;
    typedef Opm::OilPvtInterface<Scalar, Evaluation> OilPvtInterface;
    typedef Opm::WaterPvtInterface<Scalar, Evaluation> WaterPvtInterface;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    class ParameterCache : public Opm::NullParameterCache
    {
    public:
        ParameterCache(int regionIdx=0)
        { regionIdx_ = 0; }

        /*!
         * \brief Return the index of the region which should be used to determine the
         *        thermodynamic properties
         */
        int regionIndex() const
        { return regionIdx_; }

        /*!
         * \brief Set the index of the region which should be used to determine the
         *        thermodynamic properties
         */
        void setRegionIndex(int val)
        { regionIdx_ = val; }

    private:
        int regionIdx_;
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

    /*!
     * \copydoc BaseFluidSystem::init
     *
     * \attention For this fluid system, this method just throws a
     *            <tt>std::logic_error</tt> as there is no
     *            way to generically calculate the required black oil
     *            parameters. Instead of this method, use
     * \code
     * FluidSystem::initBegin();
     * // set the black oil parameters
     * FluidSystem::initEnd();
     * \endcode
     */
    static void init()
    {
        OPM_THROW(std::logic_error,
                  "There is no generic init() method for this fluid system. The "
                  << "black-oil fluid system must be initialized using:\n"
                  << "    FluidSystem::initBegin()\n"
                  << "    // set black oil parameters\n"
                  << "    FluidSystem::initEnd()\n");
    }

    /*!
     * \brief Begin the initialization of the black oil fluid system.
     *
     * After calling this method the reference densities, all dissolution and formation
     * volume factors, the oil bubble pressure, all viscosities and the water
     * compressibility must be set. Before the fluid system can be used, initEnd() must
     * be called to finalize the initialization.
     */
    static void initBegin(int numPvtRegions)
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
    static void setGasPvt(std::shared_ptr<const GasPvtInterface> pvtObj)
    { gasPvt_ = pvtObj; }

    /*!
     * \brief Set the pressure-volume-saturation (PVT) relations for the oil phase.
     */
    static void setOilPvt(std::shared_ptr<const OilPvtInterface> pvtObj)
    { oilPvt_ = pvtObj; }

    /*!
     * \brief Set the pressure-volume-saturation (PVT) relations for the water phase.
     */
    static void setWaterPvt(std::shared_ptr<const WaterPvtInterface> pvtObj)
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
                                      int regionIdx)
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
        int numRegions = molarMass_.size();
        for (int regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
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
    static const char *phaseName(const int phaseIdx)
    {
        static const char *name[] = { "water", "oil", "gas" };

        assert(0 <= phaseIdx && phaseIdx < numPhases + 1);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(const int phaseIdx)
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
    static const char *componentName(int compIdx)
    {
        static const char *name[] = { "Oil", "Water", "Gas" };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(int compIdx, int regionIdx = 0)
    { return molarMass_[regionIdx][compIdx]; }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(int phaseIdx)
    {
        // fugacity coefficients are only pressure dependent -> we
        // have an ideal mixture
        return true;
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(int phaseIdx)
    { return true; /* all phases are compressible */ }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(int phaseIdx)
    { return false; }

    /****************************************
     * thermodynamic relations
     ****************************************/
    //! \copydoc BaseFluidSystem::density
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval density(const FluidState &fluidState,
                           ParameterCache &paramCache,
                           const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        typedef typename FluidState::Scalar FsEval;
        typedef Opm::MathToolbox<FsEval> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));
        int regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case waterPhaseIdx: return waterDensity<LhsEval>(T, p, regionIdx);
        case gasPhaseIdx: {
            const auto& XgO = FsToolbox::template toLhs<LhsEval>(fluidState.massFraction(gasPhaseIdx, oilCompIdx));
            return gasDensity<LhsEval>(T, p, XgO, regionIdx);
        }
        case oilPhaseIdx: {
            const auto& XoG = FsToolbox::template toLhs<LhsEval>(fluidState.massFraction(oilPhaseIdx, gasCompIdx));
            return oilDensity<LhsEval>(T, p, XoG, regionIdx);
        }
        }

        OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval fugacityCoefficient(const FluidState &fluidState,
                                       const ParameterCache &paramCache,
                                       int phaseIdx,
                                       int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);
        assert(0 <= compIdx  && compIdx <= numComponents);

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));
        int regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case waterPhaseIdx: return fugCoefficientInWater<LhsEval>(compIdx, T, p, regionIdx);
        case gasPhaseIdx: return fugCoefficientInGas<LhsEval>(compIdx, T, p, regionIdx);
        case oilPhaseIdx: return fugCoefficientInOil<LhsEval>(compIdx, T, p, regionIdx);
        }

        OPM_THROW(std::logic_error, "Unhandled phase or component index");
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState, class LhsEval = typename FluidState::Scalar>
    static LhsEval viscosity(const FluidState &fluidState,
                             const ParameterCache &paramCache,
                             int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const auto& p = FsToolbox::template toLhs<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& T = FsToolbox::template toLhs<LhsEval>(fluidState.temperature(phaseIdx));
        int regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case oilPhaseIdx: {
            const auto& XoG = FsToolbox::template toLhs<LhsEval>(fluidState.massFraction(oilPhaseIdx, gasCompIdx));
            return oilPvt_->viscosity(regionIdx, T, p, XoG);
        }
        case waterPhaseIdx:
            return waterPvt_->viscosity(regionIdx, T, p);
        case gasPhaseIdx: {
            const auto& XgO = FsToolbox::template toLhs<LhsEval>(fluidState.massFraction(gasPhaseIdx, oilCompIdx));
            return gasPvt_->viscosity(regionIdx, T, p, XgO);
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
    static Scalar referenceDensity(int phaseIdx, int regionIdx)
    { return referenceDensity_[regionIdx][phaseIdx]; }

    /*!
     * \brief Returns the oil formation volume factor \f$B_o\f$ of saturated oil for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    template <class LhsEval>
    static LhsEval saturatedOilFormationVolumeFactor(const LhsEval& temperature,
                                                     const LhsEval& pressure,
                                                     int regionIdx)
    {
        Valgrind::CheckDefined(pressure);

        // calculate the mass fractions of gas and oil
        const auto& XoG = saturatedOilGasMassFraction(temperature, pressure, regionIdx);

        return oilFormationVolumeFactor(temperature, pressure, XoG, regionIdx);
    }

    /*!
     * \brief Return the formation volume factor of water.
     */
    template <class LhsEval>
    static LhsEval waterFormationVolumeFactor(const LhsEval& temperature,
                                              const LhsEval& pressure,
                                              int regionIdx)
    { return waterPvt_->formationVolumeFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    template <class LhsEval>
    static LhsEval gasDissolutionFactor(const LhsEval& temperature,
                                        const LhsEval& pressure,
                                        int regionIdx)
    { return oilPvt_->gasDissolutionFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    template <class LhsEval>
    static LhsEval oilVaporizationFactor(const LhsEval& temperature,
                                         const LhsEval& pressure,
                                         int regionIdx)
    { return gasPvt_->oilVaporizationFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the water phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    template <class LhsEval>
    static LhsEval fugCoefficientInWater(int compIdx,
                                         const LhsEval& temperature,
                                         const LhsEval& pressure,
                                         int regionIdx)
    { return waterPvt_->fugacityCoefficient(regionIdx, temperature, pressure, compIdx); }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the gas phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    template <class LhsEval>
    static LhsEval fugCoefficientInGas(int compIdx,
                                       const LhsEval& temperature,
                                       const LhsEval& pressure,
                                       int regionIdx)
    { return gasPvt_->fugacityCoefficient(regionIdx, temperature, pressure, compIdx); }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the oil phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    template <class LhsEval>
    static LhsEval fugCoefficientInOil(int compIdx,
                                       const LhsEval& temperature,
                                       const LhsEval& pressure,
                                       int regionIdx)
    { return oilPvt_->fugacityCoefficient(regionIdx, temperature, pressure, compIdx); }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param XoG The mass fraction of the gas component in the oil phase [-]
     */
    template <class LhsEval>
    static LhsEval oilSaturationPressure(const LhsEval& temperature,
                                         const LhsEval& XoG,
                                         int regionIdx)
    { return oilPvt_->oilSaturationPressure(regionIdx, temperature, XoG); }

    /*!
     * \brief The maximum mass fraction of the gas component in the oil phase.
     */
    template <class LhsEval>
    static LhsEval saturatedOilGasMassFraction(const LhsEval& temperature,
                                               const LhsEval& pressure,
                                               int regionIdx)
    { return oilPvt_->saturatedOilGasMassFraction(regionIdx, temperature, pressure); }

    /*!
     * \brief The maximum mole fraction of the gas component in the oil phase.
     */
    template <class LhsEval>
    static LhsEval saturatedOilGasMoleFraction(const LhsEval& temperature,
                                               const LhsEval& pressure,
                                               int regionIdx)
    { return oilPvt_->saturatedOilGasMoleFraction(regionIdx, temperature, pressure); }

    /*!
     * \brief The maximum mass fraction of the oil component in the gas phase.
     */
    template <class LhsEval>
    static LhsEval saturatedGasOilMassFraction(const LhsEval& temperature,
                                               const LhsEval& pressure,
                                               int regionIdx)
    { return gasPvt_->saturatedGasOilMassFraction(regionIdx, temperature, pressure); }

    /*!
     * \brief The maximum mole fraction of the oil component in the gas phase.
     */
    template <class LhsEval>
    static LhsEval saturatedGasOilMoleFraction(const LhsEval& temperature,
                                               const LhsEval& pressure,
                                               int regionIdx)
    { return gasPvt_->saturatedGasOilMoleFraction(regionIdx, temperature, pressure); }

    /*!
     * \brief Return the normalized formation volume factor of (potentially)
     *        under-saturated oil.
     */
    template <class LhsEval>
    static LhsEval oilFormationVolumeFactor(const LhsEval& temperature,
                                            const LhsEval& pressure,
                                            const LhsEval& XoG,
                                            int regionIdx)
    { return oilPvt_->formationVolumeFactor(regionIdx, temperature, pressure, XoG); }

    /*!
     * \brief Return the density of (potentially) under-saturated oil.
     */
    template <class LhsEval>
    static LhsEval oilDensity(const LhsEval& temperature,
                              const LhsEval& pressure,
                              const LhsEval& XoG,
                              int regionIdx)
    { return oilPvt_->density(regionIdx, temperature, pressure, XoG); }

    /*!
     * \brief Return the density of gas-saturated oil.
     */
    template <class LhsEval>
    static LhsEval saturatedOilDensity(const LhsEval& temperature,
                                       const LhsEval& pressure,
                                       int regionIdx)
    {
        // mass fraction of gas-saturated oil
        const LhsEval& XoG = saturatedOilGasMassFraction(temperature, pressure, regionIdx);
        return oilPvt_->density(regionIdx, temperature, pressure, XoG);
    }

    /*!
     * \brief Return the formation volume factor of gas.
     */
    template <class LhsEval>
    static LhsEval gasFormationVolumeFactor(const LhsEval& temperature,
                                            const LhsEval& pressure,
                                            const LhsEval& XgO,
                                            int regionIdx)
    { return gasPvt_->formationVolumeFactor(regionIdx, temperature, pressure, XgO); }

    /*!
     * \brief Return the density of dry gas.
     */
    template <class LhsEval>
    static LhsEval gasDensity(const LhsEval& temperature,
                              const LhsEval& pressure,
                              const LhsEval& XgO,
                              int regionIdx)
    { return gasPvt_->density(regionIdx, temperature, pressure, XgO); }

    /*!
     * \brief Return the density of water.
     */
    template <class LhsEval>
    static LhsEval waterDensity(const LhsEval& temperature,
                                const LhsEval& pressure,
                                int regionIdx)
    { return waterPvt_->density(regionIdx, temperature, pressure); }

private:
    static void resizeArrays_(int numRegions)
    {
        molarMass_.resize(numRegions);
        referenceDensity_.resize(numRegions);
    }

    static std::shared_ptr<const Opm::GasPvtInterface<Scalar, Evaluation> > gasPvt_;
    static std::shared_ptr<const Opm::OilPvtInterface<Scalar, Evaluation> > oilPvt_;
    static std::shared_ptr<const Opm::WaterPvtInterface<Scalar, Evaluation> > waterPvt_;

    static bool enableDissolvedGas_;
    static bool enableVaporizedOil_;

    // HACK for GCC 4.4: the array size has to be specified using the literal value '3'
    // here, because GCC 4.4 seems to be unable to determine the number of phases from
    // the BlackOil fluid system in the attribute declaration below...
    static std::vector<std::array<Scalar, /*numPhases=*/3> > referenceDensity_;
    static std::vector<std::array<Scalar, /*numComponents=*/3> > molarMass_;
};

template <class Scalar, class Evaluation>
const Scalar
BlackOil<Scalar, Evaluation>::surfaceTemperature = 273.15 + 15.56; // [K]

template <class Scalar, class Evaluation>
const Scalar
BlackOil<Scalar, Evaluation>::surfacePressure = 101325.0; // [Pa]

template <class Scalar, class Evaluation>
bool BlackOil<Scalar, Evaluation>::enableDissolvedGas_;

template <class Scalar, class Evaluation>
bool BlackOil<Scalar, Evaluation>::enableVaporizedOil_;

template <class Scalar, class Evaluation>
std::shared_ptr<const OilPvtInterface<Scalar, Evaluation> >
BlackOil<Scalar, Evaluation>::oilPvt_;

template <class Scalar, class Evaluation>
std::shared_ptr<const GasPvtInterface<Scalar, Evaluation> >
BlackOil<Scalar, Evaluation>::gasPvt_;

template <class Scalar, class Evaluation>
std::shared_ptr<const WaterPvtInterface<Scalar, Evaluation> >
BlackOil<Scalar, Evaluation>::waterPvt_;

template <class Scalar, class Evaluation>
std::vector<std::array<Scalar, 3> >
BlackOil<Scalar, Evaluation>::referenceDensity_;

template <class Scalar, class Evaluation>
std::vector<std::array<Scalar, 3> >
BlackOil<Scalar, Evaluation>::molarMass_;
}} // namespace Opm, FluidSystems

#endif
