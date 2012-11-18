// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Ewoms::FluidSystems::BlackOil
 */
#ifndef EWOMS_BLACK_OIL_FLUID_SYSTEM_HH
#define EWOMS_BLACK_OIL_FLUID_SYSTEM_HH

#include <ewoms/common/exceptions.hh>
#include <ewoms/common/spline.hh>

#include <ewoms/material/constants.hh>
#include <ewoms/material/components/h2o.hh>
#include <ewoms/material/fluidsystems/basefluidsystem.hh>
#include <ewoms/material/fluidsystems/nullparametercache.hh>

#include <array>

namespace Ewoms {
namespace FluidSystems {

/*!
 * \brief A fluid system which uses the black-oil parameters
 *        to calculate termodynamically meaningful quantities.
 */
template <class Scalar>
class BlackOil
    : public BaseFluidSystem<Scalar, BlackOil<Scalar> >
{
    typedef Ewoms::Spline<Scalar, -1> Spline;
    typedef std::vector<std::pair<Scalar, Scalar> > SplineSamplingPoints;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    typedef Ewoms::NullParameterCache ParameterCache;

    /****************************************
     * Fluid phase parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 3;

    //! Index of the oil phase
    static const int oPhaseIdx = 0;
    //! Index of the water phase
    static const int wPhaseIdx = 1;
    //! Index of the gas phase
    static const int gPhaseIdx = 2;

    /*!
     * \copydoc BaseFluidSystem::init
     *
     * \attention For this fluid system, this method just throws a
     *            <tt>Dune::InvalidStateException</tt> as there is no
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
        DUNE_THROW(Dune::InvalidStateException, "No generic init() method for this fluid system. The black-oil fluid system must be initialized with:\n"
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
        surfaceDensity_[oPhaseIdx] = rhoOil;
        surfaceDensity_[wPhaseIdx] = rhoWater;
        surfaceDensity_[gPhaseIdx] = rhoGas;
    }

    /*!
     * \brief Initialize the spline for the gas formation volume factor
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setGasFormationFactor(const SplineSamplingPoints &samplePoints)
    {
        // we discard the last sample point because that is for
        // undersaturated oil
        //SplineSamplingPoints tmp(samplePoints.begin(), --samplePoints.end());
        SplineSamplingPoints tmp(samplePoints.begin(), samplePoints.end());
        gasFormationFactorSpline_.setContainerOfTuples(tmp);
    }

    /*!
     * \brief Initialize the spline for the oil formation volume factor
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setOilFormationVolumeFactor(const SplineSamplingPoints &samplePoints)
    {
        oilFormationVolumeFactorSpline_.setContainerOfTuples(samplePoints);
        assert(oilFormationVolumeFactorSpline_.monotonic());
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
    static void setReferenceVolumeFactor(int phaseIdx,
                                         Scalar val)
    {
        refFormationVolumeFactor_[phaseIdx] = val;
    }

    /*!
     * \brief Set the bubble pressure of the oil in the reservoir [Pa]
     *
     * (That's the pressure below which gas starts to appear.)
     *
     * \param val The value of the bubble pressure [Pa].
     */
    static void setBubblePressure(Scalar val)
    { bubblePressure_ = val; }

    /*!
     * \brief Initialize the spline for the oil viscosity
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setOilViscosity(const SplineSamplingPoints &samplePoints)
    {
        oilViscositySpline_.setContainerOfTuples(samplePoints);
        assert(oilViscositySpline_.monotonic());
    }

    /*!
     * \brief Initialize the spline for the gas formation volume factor
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setGasFormationVolumeFactor(const SplineSamplingPoints &samplePoints)
    {
        // we use linear interpolation between the first and the
        // second sampling point to avoid a non-monotonic spline...
        Bg0_ = samplePoints[0];
        mBg0_ =
            (samplePoints[1].second - samplePoints[0].second) /
            (samplePoints[1].first - samplePoints[0].first);

        SplineSamplingPoints tmp(++samplePoints.begin(), samplePoints.end());
        gasFormationVolumeFactorSpline_.setContainerOfTuples(tmp);
        assert(gasFormationVolumeFactorSpline_.monotonic());
    }

    /*!
     * \brief Initialize the spline for the gas viscosity
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setGasViscosity(const SplineSamplingPoints &samplePoints)
    {
        gasViscositySpline_.setContainerOfTuples(samplePoints);
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
        molarMass_[wCompIdx] = 18e-3;

        // for gas, we take the density at surface pressure and assume
        // it to be ideal
        Scalar p = 1.0135e5;
        Scalar rho_g = gasDensity_(p);
        Scalar T = 311;
        molarMass_[gCompIdx] = Ewoms::Constants<Scalar>::R*T*rho_g / p;

        // finally, for oil phase, we take the molar mass from the
        // spe9 paper
        molarMass_[oCompIdx] = 175e-3; // kg/mol

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        int n = gasFormationFactorSpline_.numSamples()*5;
        int delta =
            (gasFormationFactorSpline_.xMax() - gasFormationFactorSpline_.xMin())
            /(n + 1);

        SplineSamplingPoints pSatSamplePoints;
        Scalar X_oG = 0;
        for (int i=0; i <= n && X_oG < 0.7; ++ i) {
            Scalar pSat = i*delta + gasFormationFactorSpline_.xMin();
            X_oG = flashGasMassFracInOil_(pSat);

            std::pair<Scalar, Scalar> val(X_oG, pSat);
            pSatSamplePoints.push_back(val);;
        };
        saturationPressureSpline_.setContainerOfTuples(pSatSamplePoints);

#if 0
        {
            int n = 1000;
            Scalar pMin = 10e6;
            Scalar pMax = 35e6;
            for (int i = 0; i < n; ++i) {
                Scalar p = Scalar(i)/n*(pMax - pMin) + pMin;

                Scalar Bg = gasFormationVolumeFactor(p);
                Scalar Bo = oilFormationVolumeFactor(p);

                std::cerr.precision(16);
                std::cerr << p << " "
                          << Bo << " "
                          << Bg << " "
                          << "\n";
            };
            exit(1);
        }
#endif
    }

    //! \copydoc BaseFluidSystem::phaseName
    static const char *phaseName(const int phaseIdx)
    {
        static const char *name[] = { "o", "w", "g" };

        assert(0 <= phaseIdx && phaseIdx < numPhases + 1);
        return name[phaseIdx];
    }

    /*!
     * \brief Returns the pressure [Pa] at which the oil started to
     *        degas in the flash experiment which determined the
     *        black-oil parameters.
     */
    static Scalar bubblePressure()
    { return bubblePressure_; }

    //! \copydoc BaseFluidSystem::isLiquid
    static constexpr bool isLiquid(const int phaseIdx)
    {
        // assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gPhaseIdx;
    }

    /****************************************
     * Component related parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 3;

    //! Index of the oil component
    static const int oCompIdx = 0;
    //! Index of the water component
    static const int wCompIdx = 1;
    //! Index of the gas component
    static const int gCompIdx = 2;

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
    static constexpr bool isIdealMixture(int phaseIdx)
    {
        // fugacity coefficients are only pressure dependent -> we
        // have an ideal mixture
        return true;
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static constexpr bool isCompressible(int phaseIdx)
    { return true; /* all phases are compressible */ }

    //! \copydoc BaseFluidSystem::isIdealGas
    static constexpr bool isIdealGas(int phaseIdx)
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
        case wPhaseIdx: return waterDensity_(p);
        case gPhaseIdx: return gasDensity_(p);
        case oPhaseIdx: {
            Scalar X_oG = fluidState.massFraction(oPhaseIdx, gCompIdx);
            Scalar X_oO = fluidState.massFraction(oPhaseIdx, oCompIdx);
            Scalar sumX =
                std::max(1e-20,
                         fluidState.massFraction(oPhaseIdx, oCompIdx)
                         + fluidState.massFraction(oPhaseIdx, wCompIdx)
                         + fluidState.massFraction(oPhaseIdx, gCompIdx));

            Scalar pAlpha = oilSaturationPressure(X_oG/sumX);
            Scalar X_oGAlpha = fluidState.massFraction(oPhaseIdx, gCompIdx) / sumX;
            //Scalar X_oOAlpha = 1 - X_oGAlpha;
            Scalar rho_oAlpha =
                flashOilDensity_(pAlpha)
                +
                oilCompressibility()
                *(p - pAlpha);

            Scalar pBeta = p;
            Scalar X_oGBeta = flashGasMassFracInOil_(pBeta);
            Scalar X_oOBeta = 1 - X_oGBeta;
            Scalar rho_oBeta = flashOilDensity_(pBeta);

            Scalar rho_oPure =
                surfaceDensity_[oPhaseIdx]
                + oilCompressibility()*(p - 1.0135e5);

            Scalar drho_dXoO =
                (1.0*rho_oPure - X_oOBeta*rho_oBeta)
                / (1.0 - X_oOBeta);
            Scalar drho_dXoG =
                (X_oGAlpha*rho_oAlpha - X_oGBeta*rho_oBeta)
                / (X_oGAlpha - X_oOBeta);

            Scalar rho_o =
                rho_oBeta*(1 - X_oGBeta)
                + drho_dXoG*X_oG
                + drho_dXoO*(X_oO - X_oOBeta);

            return std::max(250.0, std::min(1250.0, rho_o));
        }
        }

        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
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
        case wPhaseIdx: return fugCoefficientInWater(compIdx, p);
        case gPhaseIdx: return fugCoefficientInGas(compIdx, p);
        case oPhaseIdx: return fugCoefficientInOil(compIdx, p);
        }

        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase or component index");
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
        case oPhaseIdx: return oilViscositySpline_.eval(p, /*extrapolate=*/true);
        case wPhaseIdx: return waterViscosity_(p);
        case gPhaseIdx: return gasViscosity_(p);
        }

        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
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
    static Scalar gasFormationFactor(Scalar pressure)
    {
        return gasFormationFactorSpline_.eval(pressure,
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
        if (compIdx == oCompIdx ||
            compIdx == gCompIdx)
        {
            return 1e6*pvWater / pressure;
        }

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

        if (compIdx == oCompIdx)
            return phi_oO;
        else if (compIdx == wCompIdx)
            // assume that the affinity of the water component to the
            // oil phase is 1million times smaller than that of the
            // oil component
            return 1e6*phi_oO;

        /////////////
        // the rest of this method determines the fugacity coefficient
        // of the gas component:
        //
        // first, retrieve the mole fraction of gas a "flashed" oil
        // would exhibit at the given pressure
        Scalar x_oGf = flashGasMoleFracInOil_(pressure);

        // then, scale the gas component's gas phase fugacity
        // coefficient, so that the oil phase ends up at the right
        // composition if we were doing a flash experiment
        Scalar phi_gG = fugCoefficientInGas(gCompIdx, pressure);
        return phi_gG / x_oGf;
    }

    /*!
     * \brief Return the oil phase compressibility at _constant_ composition
     */
    static Scalar oilCompressibility()
    {
        return
            (1.1200 - 1.1189)/((5000 - 4000)*6894.76); // [kg/m^3 / Pa)]
    }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param X_oG The mass fraction of the gas component in the oil phase [-]
     */
    static Scalar oilSaturationPressure(Scalar X_oG)
    { return saturationPressureSpline_.eval(X_oG, /*extrapolate=*/true); }

private:
    static Scalar gasDensity_(Scalar pressure)
    {
        // gas formation volume factor at reservoir pressure
        Scalar Bg = gasFormationVolumeFactor(pressure);

        // gas formation volume factor at standard pressure
        Scalar BgRef = refFormationVolumeFactor_[gPhaseIdx];

        // surface density of gas
        Scalar rhoRef = surfaceDensity_[gPhaseIdx];

        // reservoir density is surface density scaled by the ratio of
        // the volume formation factors
        return rhoRef * BgRef/Bg;
    }

    // density of oil in a flash experiment
    static Scalar flashOilDensity_(Scalar pressure)
    {
        // oil formation volume factor at reservoir pressure
        Scalar Bo = oilFormationVolumeFactor(pressure);

        // oil formation volume factor at standard pressure
        Scalar BoRef  = refFormationVolumeFactor_[oPhaseIdx];

        // surface density of oil
        Scalar rhoRef = surfaceDensity_[oPhaseIdx];

        // reservoir density is surface density scaled by the ratio of
        // the volume formation factors
        return rhoRef * BoRef/Bo;
    }

    // the mass fraction of the gas component in the oil phase in a
    // flash experiment
    static Scalar flashGasMassFracInOil_(Scalar pressure)
    {
        // first, we calculate the total reservoir oil phase density
        // [kg/m^3]
        Scalar rho_gRef = surfaceDensity_[gPhaseIdx];

        Scalar B_o = oilFormationVolumeFactor(pressure);
        Scalar rho_o = surfaceDensity(oPhaseIdx)/B_o;


        // then, we calculate the mass of the gas component [kg/m^3]
        // in the oil phase. This is equivalent to the gas formation
        // factor [m^3/m^3] at current pressure times the gas density
        // [kg/m^3] at standard pressure
        Scalar rho_oG = gasFormationFactor(pressure) * rho_gRef;

        // with these quantities it is pretty easy to calculate to the
        // composition of the oil phase in terms of mass fractions
        return rho_oG/rho_o;
    }

    // the mole fraction of the gas component in the oil phase in a
    // flash experiment
    static Scalar flashGasMoleFracInOil_(Scalar pressure)
    {
        // calculate the mass fractions of gas and oil
        Scalar X_oG = flashGasMassFracInOil_(pressure);
        //Scalar X_oO = 1 - X_oG;

        // which can be converted to mole fractions, given the
        // components' molar masses
        Scalar M_G = molarMass(gCompIdx);
        Scalar M_O = molarMass(oCompIdx);
        Scalar x_oO = M_G*(1 - X_oG) / (M_G + X_oG*(M_O - M_G));

        return 1.0 - x_oO;
    }

    static Scalar waterDensity_(Scalar pressure)
    {
        // compressibility of water times standard density
        Scalar rhoRef = surfaceDensity_[wPhaseIdx];
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

    static Spline gasFormationFactorSpline_;
    static Spline gasFormationVolumeFactorSpline_;
    static Spline saturationPressureSpline_;
    static Scalar bubblePressure_;

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
BlackOil<Scalar>::gasFormationFactorSpline_;

template <class Scalar>
typename BlackOil<Scalar>::Spline
BlackOil<Scalar>::gasFormationVolumeFactorSpline_;

template <class Scalar>
typename BlackOil<Scalar>::Spline
BlackOil<Scalar>::saturationPressureSpline_;

template <class Scalar>
Scalar BlackOil<Scalar>::bubblePressure_;

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
} // end namepace
} // end namepace

#endif
