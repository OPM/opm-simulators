/*
  Copyright (C) 2011-2013 by Andreas Lauser

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

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <memory>
#include <vector>
#include <array>

namespace Opm {
template <class Scalar>
class OilPvtInterface;

template <class Scalar>
class GasPvtInterface;

template <class Scalar>
class WaterPvtInterface;

namespace FluidSystems {

/*!
 * \brief A fluid system which uses the black-oil parameters
 *        to calculate termodynamically meaningful quantities.
 */
template <class Scalar>
class BlackOil : public BaseFluidSystem<Scalar, BlackOil<Scalar> >
{
    typedef Opm::GasPvtInterface<Scalar> GasPvtInterface;
    typedef Opm::OilPvtInterface<Scalar> OilPvtInterface;
    typedef Opm::WaterPvtInterface<Scalar> WaterPvtInterface;

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

    //! Index of the oil phase
    static const int oilPhaseIdx = 0;
    //! Index of the water phase
    static const int waterPhaseIdx = 1;
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
        static const char *name[] = { "oil", "water", "gas" };

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
        static const char *name[] = { "O", "W", "G" };

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
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          ParameterCache &paramCache,
                          const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        Scalar p = fluidState.pressure(phaseIdx);
        Scalar T = fluidState.temperature(phaseIdx);
        int regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case waterPhaseIdx: return waterDensity(T, p, regionIdx);
        case gasPhaseIdx: {
            Scalar XgO = fluidState.massFraction(gasPhaseIdx, oilCompIdx);
            return gasDensity(T, p, XgO, regionIdx);
        }
        case oilPhaseIdx: {
            Scalar XoG = fluidState.massFraction(oilPhaseIdx, gasCompIdx);
            return oilDensity(T, p, XoG, regionIdx);
        }
        }

        OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);
        assert(0 <= compIdx  && compIdx <= numComponents);

        Scalar p = fluidState.pressure(phaseIdx);
        Scalar T = fluidState.temperature(phaseIdx);
        int regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case waterPhaseIdx: return fugCoefficientInWater(compIdx, T, p, regionIdx);
        case gasPhaseIdx: return fugCoefficientInGas(compIdx, T, p, regionIdx);
        case oilPhaseIdx: return fugCoefficientInOil(compIdx, T, p, regionIdx);
        }

        OPM_THROW(std::logic_error, "Unhandled phase or component index");
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        Scalar p = fluidState.pressure(phaseIdx);
        Scalar T = fluidState.temperature(phaseIdx);
        int regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case oilPhaseIdx: {
            Scalar XoG = fluidState.massFraction(oilPhaseIdx, gasCompIdx);
            return oilPvt_->viscosity(regionIdx, T, p, XoG);
        }
        case waterPhaseIdx:
            return waterPvt_->viscosity(regionIdx, T, p);
        case gasPhaseIdx: {
            Scalar XgO = fluidState.massFraction(gasPhaseIdx, oilCompIdx);
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
    static Scalar saturatedOilFormationVolumeFactor(Scalar temperature,
                                                    Scalar pressure,
                                                    int regionIdx)
    {
        Valgrind::CheckDefined(pressure);

        // calculate the mass fractions of gas and oil
        Scalar XoG = saturatedOilGasMassFraction(temperature, pressure, regionIdx);

        // ATTENTION: XoG is represented by the _first_ axis!
        return oilFormationVolumeFactor(temperature, pressure, XoG, regionIdx);
    }

    /*!
     * \brief Return the formation volume factor of water.
     */
    static Scalar waterFormationVolumeFactor(Scalar temperature, Scalar pressure, int regionIdx)
    { return waterPvt_->formationVolumeFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar gasDissolutionFactor(Scalar temperature, Scalar pressure, int regionIdx)
    { return oilPvt_->gasDissolutionFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the oil vaporization factor \f$R_v\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar oilVaporizationFactor(Scalar temperature, Scalar pressure, int regionIdx)
    { return gasPvt_->oilVaporizationFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the water phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar fugCoefficientInWater(int compIdx, Scalar temperature, Scalar pressure, int regionIdx)
    { return waterPvt_->fugacityCoefficient(regionIdx, temperature, pressure, compIdx); }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the gas phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar fugCoefficientInGas(int compIdx, Scalar temperature, Scalar pressure, int regionIdx)
    { return gasPvt_->fugacityCoefficient(regionIdx, temperature, pressure, compIdx); }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the oil phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar fugCoefficientInOil(int compIdx, Scalar temperature, Scalar pressure, int regionIdx)
    { return oilPvt_->fugacityCoefficient(regionIdx, temperature, pressure, compIdx); }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param XoG The mass fraction of the gas component in the oil phase [-]
     */
    static Scalar oilSaturationPressure(Scalar temperature, Scalar XoG, int regionIdx)
    { return oilPvt_->oilSaturationPressure(regionIdx, temperature, XoG); }

    /*!
     * \brief The maximum mass fraction of the gas component in the oil phase.
     */
    static Scalar saturatedOilGasMassFraction(Scalar temperature, Scalar pressure, int regionIdx)
    { return oilPvt_->saturatedOilGasMassFraction(regionIdx, temperature, pressure); }

    /*!
     * \brief The maximum mole fraction of the gas component in the oil phase.
     */
    static Scalar saturatedOilGasMoleFraction(Scalar temperature, Scalar pressure, int regionIdx)
    { return oilPvt_->saturatedOilGasMoleFraction(regionIdx, temperature, pressure); }

    /*!
     * \brief The maximum mass fraction of the oil component in the gas phase.
     */
    static Scalar saturatedGasOilMassFraction(Scalar temperature, Scalar pressure, int regionIdx)
    { return gasPvt_->saturatedGasOilMassFraction(regionIdx, temperature, pressure); }

    /*!
     * \brief The maximum mole fraction of the oil component in the gas phase.
     */
    static Scalar saturatedGasOilMoleFraction(Scalar temperature, Scalar pressure, int regionIdx)
    { return gasPvt_->saturatedGasOilMoleFraction(regionIdx, temperature, pressure); }

    /*!
     * \brief Return the normalized formation volume factor of (potentially)
     *        under-saturated oil.
     */
    static Scalar oilFormationVolumeFactor(Scalar temperature,
                                           Scalar pressure,
                                           Scalar XoG,
                                           int regionIdx)
    { return oilPvt_->formationVolumeFactor(regionIdx, temperature, pressure, XoG); }

    /*!
     * \brief Return the density of (potentially) under-saturated oil.
     */
    static Scalar oilDensity(Scalar temperature, Scalar pressure, Scalar XoG, int regionIdx)
    { return oilPvt_->density(regionIdx, temperature, pressure, XoG); }

    /*!
     * \brief Return the density of gas-saturated oil.
     */
    static Scalar saturatedOilDensity(Scalar temperature, Scalar pressure, int regionIdx)
    {
        // mass fraction of gas-saturated oil
        Scalar XoG = saturatedOilGasMassFraction(temperature, pressure, regionIdx);
        return oilPvt_->density(regionIdx, temperature, pressure, XoG);
    }

    /*!
     * \brief Return the formation volume factor of gas.
     */
    static Scalar gasFormationVolumeFactor(Scalar temperature, Scalar pressure, Scalar XgO, int regionIdx)
    { return gasPvt_->formationVolumeFactor(regionIdx, temperature, pressure, XgO); }

    /*!
     * \brief Return the density of dry gas.
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure, Scalar XgO, int regionIdx)
    { return gasPvt_->density(regionIdx, temperature, pressure, XgO); }

    /*!
     * \brief Return the density of water.
     */
    static Scalar waterDensity(Scalar temperature, Scalar pressure, int regionIdx)
    { return waterPvt_->density(regionIdx, temperature, pressure); }

private:
    static void resizeArrays_(int numRegions)
    {
        molarMass_.resize(numRegions);
        referenceDensity_.resize(numRegions);
    }

    static std::shared_ptr<const Opm::GasPvtInterface<Scalar> > gasPvt_;
    static std::shared_ptr<const Opm::OilPvtInterface<Scalar> > oilPvt_;
    static std::shared_ptr<const Opm::WaterPvtInterface<Scalar> > waterPvt_;

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
std::shared_ptr<const OilPvtInterface<Scalar> >
BlackOil<Scalar>::oilPvt_;

template <class Scalar>
std::shared_ptr<const GasPvtInterface<Scalar> >
BlackOil<Scalar>::gasPvt_;

template <class Scalar>
std::shared_ptr<const WaterPvtInterface<Scalar> >
BlackOil<Scalar>::waterPvt_;

template <class Scalar>
std::vector<std::array<Scalar, 3> >
BlackOil<Scalar>::referenceDensity_;

template <class Scalar>
std::vector<std::array<Scalar, 3> >
BlackOil<Scalar>::molarMass_;
}} // namespace Opm, FluidSystems

#endif
