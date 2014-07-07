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

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Spline.hpp>

#include <opm/material/IdealGas.hpp>
#include <opm/material/Constants.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/fluidsystems/NullParameterCache.hpp>

#include <opm/parser/eclipse/Utility/PvtoTable.hpp>
#include <opm/parser/eclipse/Utility/PvtwTable.hpp>
#include <opm/parser/eclipse/Utility/PvdgTable.hpp>

#include <array>
#include <vector>
#include <iostream>

namespace Opm {
namespace FluidSystems {

/*!
 * \brief A fluid system which uses the black-oil parameters
 *        to calculate termodynamically meaningful quantities.
 */
template <class Scalar>
class BlackOil
    : public BaseFluidSystem<Scalar, BlackOil<Scalar> >
{
    typedef Opm::Spline<Scalar> Spline;
    typedef std::vector<std::pair<Scalar, Scalar> > SplineSamplingPoints;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    typedef Opm::NullParameterCache ParameterCache;

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
     * After calling this method the surface densities, all formation
     * and volume factors, the oil bubble pressure, all viscosities
     * and the water compressibility must be set. Before the fluid
     * system can be used, initEnd() must be called to finalize the
     * initialization.
     */
    static void initBegin()
    {}

    /*!
     * \brief Sets the pressure-dependent oil viscosity, density and
     *        gas content using a table stemming from the Eclipse PVTO
     *        keyword.
     */
    static void setPvtoTable(const PvtoTable &pvtoTable)
    {
        const auto saturatedTable = pvtoTable.getOuterTable();

        int numSamples = saturatedTable->numRows();
        assert(numSamples > 1);

        gasDissolutionFactorSpline_.setXYArrays(numSamples,
                                                saturatedTable->getPressureColumn(),
                                                saturatedTable->getGasSolubilityColumn(),
                                                /*type=*/Spline::Monotonic);
        oilFormationVolumeFactorSpline_.setXYArrays(numSamples,
                                                    saturatedTable->getPressureColumn(),
                                                    saturatedTable->getOilFormationFactorColumn(),
                                                    /*type=*/Spline::Monotonic);

        oilViscositySpline_.setXYArrays(numSamples,
                                        saturatedTable->getPressureColumn(),
                                        saturatedTable->getOilViscosityColumn(),
                                        /*type=*/Spline::Monotonic);

        // we don't properly deal with undersaturated oil yet
    }

    /*!
     * \brief Sets the pressure-dependent water viscosity and density
     *        using a table stemming from the Eclipse PVTW keyword.
     *
     * This function also sets the surface viscosity and the surface
     * density of water, but these can be overwritten using setSurface*().
     */
    static void setPvtwTable(const PvtwTable &pvtwTable)
    {
        assert(pvtwTable.numRows() > 0);

        waterCompressibilityScalar_ = pvtwTable.getCompressibilityColumn()[0];
        waterViscosityScalar_ = pvtwTable.getViscosityColumn()[0];
    }

    /*!
     * \brief Sets the pressure-dependent viscosity and density of dry
     *        gas using a table stemming from the Eclipse PVDG
     *        keyword.
     *
     * This function also sets the surface viscosity and the surface
     * density of gas, but these can be overwritten using
     * setSurface*().
     */
    static void setPvdgTable(const PvdgTable &pvdgTable)
    {
        int numSamples = pvdgTable.numRows();

        assert(numSamples > 1);
        gasFormationVolumeFactorSpline_.setXYArrays(numSamples,
                                                    pvdgTable.getPressureColumn(),
                                                    pvdgTable.getFormationFactorColumn(),
                                                    /*type=*/Spline::Monotonic);

        gasViscositySpline_.setXYArrays(numSamples,
                                        pvdgTable.getPressureColumn(),
                                        pvdgTable.getFormationFactorColumn(),
                                        /*type=*/Spline::Monotonic);

    }

    /*!
     * \brief Initialize the values of the surface densities
     *
     * \param rhoOil The surface density of (gas saturated) oil phase.
     * \param rhoWater The surface density of the water phase.
     * \param rhoGas The surface density of the gas phase.
     */
    static void setSurfaceDensities(Scalar rhoOil,
                                    Scalar rhoWater,
                                    Scalar rhoGas)
    {
        surfaceDensity_[oilPhaseIdx] = rhoOil;
        surfaceDensity_[waterPhaseIdx] = rhoWater;
        surfaceDensity_[gasPhaseIdx] = rhoGas;
    }

    /*!
     * \brief Initialize the spline for the gas dissolution factor \f$R_s\f$
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setSaturatedOilGasDissolutionFactor(const SplineSamplingPoints &samplePoints)
    {
        gasDissolutionFactorSpline_.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);
        assert(gasDissolutionFactorSpline_.monotonic());
    }

    /*!
     * \brief Initialize the spline for the oil formation volume factor
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setSaturatedOilFormationVolumeFactor(const SplineSamplingPoints &samplePoints)
    {
        oilFormationVolumeFactorSpline_.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);
        assert(oilFormationVolumeFactorSpline_.monotonic());
    }

    /*!
     * \brief Initialize the spline for the viscosity of gas-saturated
     *        oil.
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setSaturatedOilViscosity(const SplineSamplingPoints &samplePoints)
    {
        oilViscositySpline_.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);
        assert(oilViscositySpline_.monotonic());
    }

    /*!
     * \brief Set the volume factor of a phase at reference pressure.
     *
     * This should usually be 1, but some codes like ECLiPSE also
     * allow different values for unknown reasons.
     *
     * \param phaseIdx The index of the fluid phase.
     * \param val The value of the reference volume factor.
     */
    static void setReferenceVolumeFactor(int phaseIdx, Scalar val)
    { refFormationVolumeFactor_[phaseIdx] = val; }

    /*!
     * \brief Initialize the spline for the formation volume factor of dry gas
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setGasFormationVolumeFactor(const SplineSamplingPoints &samplePoints)
    {
        gasFormationVolumeFactorSpline_.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);
        assert(gasFormationVolumeFactorSpline_.monotonic());
    }

    /*!
     * \brief Initialize the spline for the viscosity of dry gas
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setGasViscosity(const SplineSamplingPoints &samplePoints)
    {
        gasViscositySpline_.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);
        assert(gasViscositySpline_.monotonic());
    }

    /*!
     * \brief Set the water viscosity [Pa s]
     *
     * \param muWater The dynamic viscosity of the water phase.
     */
    static void setWaterViscosity(Scalar muWater)
    { waterViscosityScalar_ = muWater; }

    /*!
     * \brief Set the water compressibility [1 / Pa]
     *
     * \param cWater The compressibility of the water phase.
     */
    static void setWaterCompressibility(Scalar cWater)
    { waterCompressibilityScalar_ = cWater; }

    /*!
     * \brief Finish initializing the black oil fluid system.
     */
    static void initEnd()
    {
        // calculate molar masses

        // water is simple: 18 g/mol
        molarMass_[waterCompIdx] = 18e-3;

        // for gas, we take the density at standard pressure and
        // temperature and assume it to be ideal
        Scalar p = 1.0135e5;
        Scalar rho_g = surfaceDensity_[gasPhaseIdx];
        Scalar T = 297.15;
        molarMass_[gasCompIdx] = Opm::Constants<Scalar>::R*T*rho_g / p;

        // finally, for oil phase, we take the molar mass from the
        // spe9 paper
        molarMass_[oilCompIdx] = 175e-3; // kg/mol

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        int n = gasDissolutionFactorSpline_.numSamples()*5;
        int delta =
            (gasDissolutionFactorSpline_.xMax() - gasDissolutionFactorSpline_.xMin())/(n + 1);

        SplineSamplingPoints pSatSamplePoints;
        Scalar X_oG = 0;
        for (int i=0; i <= n; ++ i) {
            Scalar pSat = gasDissolutionFactorSpline_.xMin() + i*delta;
            X_oG = saturatedOilGasMassFraction(pSat);

            std::pair<Scalar, Scalar> val(X_oG, pSat);
            pSatSamplePoints.push_back(val);
        }
        saturationPressureSpline_.setContainerOfTuples(pSatSamplePoints, /*type=*/Spline::Monotonic);
    }

    //! \copydoc BaseFluidSystem::phaseName
    static const char *phaseName(const int phaseIdx)
    {
        static const char *name[] = { "o", "w", "g" };

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
    static const char *componentName(const int compIdx)
    {
        static const char *name[] = { "O", "W", "G" };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(const int compIdx)
    { return molarMass_[compIdx]; }

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

        switch (phaseIdx) {
        case waterPhaseIdx: return waterDensity_(p);
        case gasPhaseIdx: return gasDensity_(p);
        case oilPhaseIdx: return oilDensity(p, fluidState.massFraction(oilPhaseIdx, gasCompIdx));
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
        switch (phaseIdx) {
        case waterPhaseIdx: return fugCoefficientInWater(compIdx, p);
        case gasPhaseIdx: return fugCoefficientInGas(compIdx, p);
        case oilPhaseIdx: return fugCoefficientInOil(compIdx, p);
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

        switch (phaseIdx) {
        case oilPhaseIdx: return oilViscositySpline_.eval(p, /*extrapolate=*/true);
        case waterPhaseIdx: return waterViscosity_(p);
        case gasPhaseIdx: return gasViscosity_(p);
        }

        OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the density of a fluid phase at surface pressure [kg/m^3]
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    static Scalar surfaceDensity(int phaseIdx)
    { return surfaceDensity_[phaseIdx]; }

    /*!
     * \brief Returns the oil formation volume factor \f$B_o\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar oilFormationVolumeFactor(Scalar pressure)
    {
        Valgrind::CheckDefined(pressure);
        return oilFormationVolumeFactorSpline_.eval(pressure,
                                                    /*extrapolate=*/true);
    }

    /*!
     * \brief Returns the gas formation volume factor \f$B_g\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar gasFormationVolumeFactor(Scalar pressure)
    {
        if (pressure < gasFormationVolumeFactorSpline_.xMin()) {
            return (pressure - Bg0_.first)*mBg0_ + Bg0_.second;
        }
        return gasFormationVolumeFactorSpline_.eval(pressure,
                                                    /*extrapolate=*/true);
    }

    /*!
     * \brief Returns the gas formation factor \f$R_s\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar gasDissolutionFactor(Scalar pressure)
    {
        return gasDissolutionFactorSpline_.eval(pressure,
                                              /*extrapolate=*/true);
    }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the water phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar fugCoefficientInWater(int compIdx, Scalar pressure)
    {
        // set the affinity of the gas and oil components to the water
        // phase to be 6 orders of magnitute smaller than that of the
        // water component. for this we use a pseudo-realistic vapor
        // pressure of water as a starting point. (we just set it to
        // 30 kPa to ease interpreting the results.)
        const Scalar pvWater = 30e3;
        if (compIdx == oilCompIdx)
            return 1e3*pvWater / pressure;
        else if (compIdx == gasCompIdx)
            return 1e6*pvWater / pressure;

        return pvWater / pressure;
    }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the gas phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar fugCoefficientInGas(int compIdx, Scalar pressure)
    {
        // assume an ideal gas
        return 1.0;
    }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the oil phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar fugCoefficientInOil(int compIdx, Scalar pressure)
    {
        // set the oil component fugacity coefficient in oil phase
        // arbitrarily. we use some pseudo-realistic value for the vapor
        // pressure to ease physical interpretation of the results
        Scalar phi_oO = 20e3/pressure;

        if (compIdx == oilCompIdx)
            return phi_oO;
        else if (compIdx == waterCompIdx)
            // assume that the affinity of the water component to the
            // oil phase is one million times smaller than that of the
            // oil component
            return 1e6*phi_oO;

        /////////////
        // the rest of this method determines the fugacity coefficient
        // of the gas component:
        //
        // first, retrieve the mole fraction of gas a saturated oil
        // would exhibit at the given pressure
        Scalar x_oGf = saturatedOilGasMoleFraction(pressure);

        // then, scale the gas component's gas phase fugacity
        // coefficient, so that the oil phase ends up at the right
        // composition if we were doing a flash experiment
        Scalar phi_gG = fugCoefficientInGas(gasCompIdx, pressure);
        return phi_gG / x_oGf;
    }

    /*!
     * \brief Return the oil phase compressibility at _constant_ composition
     */
    static Scalar oilCompressibility()
    {
        return (1.1200 - 1.1189)/((5000 - 4000)*6894.76); // [kg/m^3 / Pa)]
    }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param X_oG The mass fraction of the gas component in the oil phase [-]
     */
    static Scalar oilSaturationPressure(Scalar X_oG)
    {
        // use the saturation pressure spline to get a pretty good
        // initial value
        Scalar pSat = saturationPressureSpline_.eval(X_oG, /*extrapolate=*/true);

        // Newton method to do the remaining work. If the initial
        // value is good, this should only take two to three
        // iterations...
        for (int i = 0; i < 20; ++i) {
            Scalar f = saturatedOilGasMassFraction(pSat) - X_oG;
            Scalar eps = pSat*1e-11;
            Scalar fPrime = ((saturatedOilGasMassFraction(pSat + eps) - X_oG) - f)/eps;

            Scalar delta = f/fPrime;
            pSat -= delta;

            if (std::abs(delta) < pSat * 1e-10)
                return pSat;
        }

        OPM_THROW(NumericalProblem, "Could find the oil saturation pressure for X_o^g = " << X_oG);
    }

    // the mass fraction of the gas component in the oil phase in a
    // flash experiment
    static Scalar saturatedOilGasMassFraction(Scalar pressure)
    {
        Scalar rho_gRef = surfaceDensity_[gasPhaseIdx];

        // calculate the mass of the gas component [kg/m^3] in the oil phase. This is
        // equivalent to the gas formation factor [m^3/m^3] at current pressure times the
        // gas density [kg/m^3] at standard pressure
        Scalar rho_oG = gasDissolutionFactor(pressure) * rho_gRef;

        // we now have the total density of saturated oil and the partial density of the
        // gas component within it. The gas mass fraction is the ratio of these two.
        return rho_oG/(surfaceDensity(oilPhaseIdx) + rho_oG);
    }

    // the mole fraction of the gas component of a gas-saturated oil phase
    static Scalar saturatedOilGasMoleFraction(Scalar pressure)
    {
        // calculate the mass fractions of gas and oil
        Scalar XoG = saturatedOilGasMassFraction(pressure);

        // which can be converted to mole fractions, given the
        // components' molar masses
        Scalar MG = molarMass(gasCompIdx);
        Scalar MO = molarMass(oilCompIdx);

        Scalar avgMolarMass = MO/(1 + XoG*(MO/MG - 1));
        Scalar xoG = XoG*avgMolarMass/MG;

        return xoG;
    }

    // density of gas-saturated oil
    static Scalar saturatedOilDensity(Scalar pressure)
    {
        // oil formation volume factor at reservoir pressure
        Scalar Bo = oilFormationVolumeFactor(pressure);

        // oil formation volume factor at standard pressure
        Scalar BoRef  = refFormationVolumeFactor_[oilPhaseIdx];

        // surface density of oil
        Scalar rhoRef = surfaceDensity_[oilPhaseIdx];

        // reservoir density is surface density scaled by the ratio of
        // the volume formation factors
        return rhoRef * BoRef/Bo;
    }


    // density of gas (potentially) under-saturated oil
    static Scalar oilDensity(Scalar oilPressure, Scalar massFractionGasInOil)
    {
        Scalar poSat = oilSaturationPressure(massFractionGasInOil);

        // retrieve the gas formation factor and the oil formation volume factor
        Scalar Rs = gasDissolutionFactorSpline_.eval(poSat, /*extrapolate=*/true);
        Scalar Bo = oilFormationVolumeFactor(poSat);

        // retrieve the derivatives of the oil formation volume
        // factor and the gas formation factor regarding pressure
        Scalar dBo_dp = oilFormationVolumeFactorSpline_.evalDerivative(poSat, /*extrapolate=*/true);
        Scalar dRs_dp = gasDissolutionFactorSpline_.evalDerivative(poSat, /*extrapolate=*/true);

        // define the derivatives of oil regarding oil component
        // mass fraction and pressure
        Scalar drhoo_dXoO =
            surfaceDensity_[oilPhaseIdx]
            * (1 + oilCompressibility()*(oilPressure - 1.0135e5));
        Scalar drhoo_dp = oilCompressibility();

        // Calculate the derivative of the density of saturated
        // oil regarding pressure
        Scalar drhoosat_dp = - surfaceDensity_[oilPhaseIdx]*dBo_dp / (Bo * Bo);

        // calculate the derivative of the gas component mass
        // fraction regarding pressure in saturated oil
        Scalar dXoOsat_dp =
            - surfaceDensity_[gasPhaseIdx]/surfaceDensity_[oilPhaseIdx]
            *(Bo * dRs_dp + Rs * dBo_dp);

        // Using the previous derivatives, define a derivative
        // for the oil density in regard to the gas mass fraction.
        Scalar drhoo_dXoG =
            drhoo_dXoO + (drhoo_dp - drhoosat_dp) / dXoOsat_dp;

        // calculate the composition of saturated oil.
        Scalar XoGsat = surfaceDensity_[gasPhaseIdx]/surfaceDensity_[oilPhaseIdx] * Rs * Bo;
        Scalar XoOsat = 1.0 - XoGsat;

        Scalar rhoo =
            surfaceDensity_[oilPhaseIdx]/Bo*(1 + drhoo_dp*(oilPressure - poSat))
            + (XoOsat - (1 - massFractionGasInOil))*drhoo_dXoO
            + (XoGsat - massFractionGasInOil)*drhoo_dXoG;

        return rhoo;
    }

private:
    static Scalar gasDensity_(Scalar pressure)
    {
        // gas formation volume factor at reservoir pressure
        Scalar Bg = gasFormationVolumeFactor(pressure);

        // gas formation volume factor at standard pressure
        Scalar BgRef = refFormationVolumeFactor_[gasPhaseIdx];

        // surface density of gas
        Scalar rhoRef = surfaceDensity_[gasPhaseIdx];

        // reservoir density is surface density scaled by the ratio of
        // the volume formation factors
        return rhoRef * BgRef/Bg;
    }

    static Scalar waterDensity_(Scalar pressure)
    {
        // compressibility of water times standard density
        Scalar rhoRef = surfaceDensity_[waterPhaseIdx];
        return rhoRef *
            (1 +
             waterCompressibilityScalar_
             * (pressure - 1.0135e5));
    }

    static Scalar gasViscosity_(Scalar pressure)
    {
        return gasViscositySpline_.eval(pressure,
                                        /*extrapolate=*/true);
    }

    static Scalar waterViscosity_(Scalar pressure)
    { return waterViscosityScalar_; }

    static Spline oilFormationVolumeFactorSpline_;
    static Spline oilViscositySpline_;

    static Spline gasDissolutionFactorSpline_;
    static Spline gasFormationVolumeFactorSpline_;
    static Spline saturationPressureSpline_;

    static Spline gasViscositySpline_;
    static std::pair<Scalar, Scalar> Bg0_;
    static Scalar mBg0_;

    static Scalar waterCompressibilityScalar_;
    static Scalar waterViscosityScalar_;

    static Scalar refFormationVolumeFactor_[numPhases];
    static Scalar surfaceDensity_[numPhases];
    static Scalar molarMass_[numComponents];
};

template <class Scalar>
typename BlackOil<Scalar>::Spline
BlackOil<Scalar>::oilFormationVolumeFactorSpline_;

template <class Scalar>
typename BlackOil<Scalar>::Spline
BlackOil<Scalar>::oilViscositySpline_;

template <class Scalar>
typename BlackOil<Scalar>::Spline
BlackOil<Scalar>::gasDissolutionFactorSpline_;

template <class Scalar>
typename BlackOil<Scalar>::Spline
BlackOil<Scalar>::gasFormationVolumeFactorSpline_;

template <class Scalar>
typename BlackOil<Scalar>::Spline
BlackOil<Scalar>::saturationPressureSpline_;

template <class Scalar>
typename BlackOil<Scalar>::Spline
BlackOil<Scalar>::gasViscositySpline_;

template <class Scalar>
std::pair<Scalar, Scalar>
BlackOil<Scalar>::Bg0_;

template <class Scalar>
Scalar BlackOil<Scalar>::mBg0_;

template <class Scalar>
Scalar
BlackOil<Scalar>::refFormationVolumeFactor_[BlackOil<Scalar>::numPhases];

template <class Scalar>
Scalar
BlackOil<Scalar>::surfaceDensity_[BlackOil<Scalar>::numPhases];

template <class Scalar>
Scalar
BlackOil<Scalar>::molarMass_[BlackOil<Scalar>::numComponents];

template <class Scalar>
Scalar
BlackOil<Scalar>::waterCompressibilityScalar_;

template <class Scalar>
Scalar
BlackOil<Scalar>::waterViscosityScalar_;
} // namespace FluidSystems
} // namespace Opm

#endif
