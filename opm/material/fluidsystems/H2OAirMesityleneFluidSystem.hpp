/*
  Copyright (C) 2009-2013 by Andreas Lauser
  Copyright (C) 2011 by Benjamin Faigle

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
 * \copydoc Opm::FluidSystems::H2OAirMesitylene
 */
#ifndef OPM_H2O_AIR_MESITYLENE_FLUID_SYSTEM_HPP
#define OPM_H2O_AIR_MESITYLENE_FLUID_SYSTEM_HPP

#include "BaseFluidSystem.hpp"
#include "NullParameterCache.hpp"

#include <opm/material/IdealGas.hpp>
#include <opm/material/components/N2.hpp>
#include <opm/material/components/Air.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Mesitylene.hpp>
#include <opm/material/components/TabulatedComponent.hpp>
#include <opm/material/binarycoefficients/H2O_Air.hpp>
#include <opm/material/binarycoefficients/H2O_N2.hpp>
#include <opm/material/binarycoefficients/H2O_Mesitylene.hpp>
#include <opm/material/binarycoefficients/Air_Mesitylene.hpp>

#include <iostream>

namespace Opm {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A fluid system with water, gas and NAPL as phases and
 *        water, air and mesitylene (DNAPL) as components.
 */
template <class Scalar>
class H2OAirMesitylene
    : public BaseFluidSystem<Scalar, H2OAirMesitylene<Scalar> >
{
    typedef H2OAirMesitylene<Scalar> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

    typedef Opm::H2O<Scalar> IapwsH2O;
    typedef Opm::TabulatedComponent<Scalar, IapwsH2O, /*alongVaporPressure=*/false> TabulatedH2O;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    typedef NullParameterCache ParameterCache;

    //! The type of the mesithylene/napl component
    typedef Opm::Mesitylene<Scalar> NAPL;

    //! The type of the air component
    typedef Opm::Air<Scalar> Air;

    //! The type of the water component
    //typedef SimpleH2O H2O;
    typedef TabulatedH2O H2O;
    //typedef IapwsH2O H2O;

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 3;
    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 3;

    //! The index of the water phase
    static const int waterPhaseIdx = 0;
    //! The index of the NAPL phase
    static const int naplPhaseIdx = 1;
    //! The index of the gas phase
    static const int gasPhaseIdx = 2;

    //! The index of the water component
    static const int H2OIdx = 0;
    //! The index of the NAPL component
    static const int NAPLIdx = 1;
    //! The index of the air pseudo-component
    static const int airIdx = 2;

    //! \copydoc BaseFluidSystem::init
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/100,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water [K]
     * \param tempMax The maximum temperature used for tabulation of water [K]
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water [Pa]
     * \param pressMax The maximum pressure used for tabulation of water [Pa]
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (H2O::isTabulated) {
            TabulatedH2O::init(tempMin, tempMax, nTemp,
                               pressMin, pressMax, nPress);
        }
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gasPhaseIdx;
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(int phaseIdx)
    { return phaseIdx == gasPhaseIdx && H2O::gasIsIdeal() && Air::gasIsIdeal() && NAPL::gasIsIdeal(); }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // gases are always compressible
        return (phaseIdx == gasPhaseIdx)
            ? true
            : (phaseIdx == waterPhaseIdx)
            ? H2O::liquidIsCompressible()
            : NAPL::liquidIsCompressible();
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    //! \copydoc BaseFluidSystem::phaseName
    static const char *phaseName(int phaseIdx)
    {
        switch (phaseIdx) {
        case waterPhaseIdx: return "w";
        case naplPhaseIdx: return "n";
        case gasPhaseIdx: return "g";;
        };
        OPM_THROW(std::logic_error, "Invalid phase index " << phaseIdx);
    }

    //! \copydoc BaseFluidSystem::componentName
    static const char *componentName(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::name();
        case airIdx: return Air::name();
        case NAPLIdx: return NAPL::name();
        };
        OPM_THROW(std::logic_error, "Invalid component index " << compIdx);
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(int compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::molarMass()
            : (compIdx == airIdx)
            ? Air::molarMass()
            : (compIdx == NAPLIdx)
            ? NAPL::molarMass()
            : -1e10;
        //OPM_THROW(std::logic_error, "Invalid component index " << compIdx);
    }

    //! \copydoc BaseFluidSystem::density
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx) ;

        if (phaseIdx == waterPhaseIdx) {
            // See: Ochs 2008
            Scalar p = H2O::liquidIsCompressible()?fluidState.pressure(phaseIdx):1e100;
            Scalar rholH2O = H2O::liquidDensity(T, p);
            Scalar clH2O = rholH2O/H2O::molarMass();

            // this assumes each dissolved molecule displaces exactly one
            // water molecule in the liquid
            return
                clH2O*(H2O::molarMass()*fluidState.moleFraction(waterPhaseIdx, H2OIdx)
                       +
                       Air::molarMass()*fluidState.moleFraction(waterPhaseIdx, airIdx)
                       +
                       NAPL::molarMass()*fluidState.moleFraction(waterPhaseIdx, NAPLIdx));
        }
        else if (phaseIdx == naplPhaseIdx) {
            // assume pure NAPL for the NAPL phase
            Scalar p = NAPL::liquidIsCompressible()?fluidState.pressure(phaseIdx):1e100;
            return NAPL::liquidDensity(T, p);
        }

        assert (phaseIdx == gasPhaseIdx);
        Scalar pH2O =
            fluidState.moleFraction(gasPhaseIdx, H2OIdx)  *
            fluidState.pressure(gasPhaseIdx);
        Scalar pAir =
            fluidState.moleFraction(gasPhaseIdx, airIdx)  *
            fluidState.pressure(gasPhaseIdx);
        Scalar pNAPL =
            fluidState.moleFraction(gasPhaseIdx, NAPLIdx)  *
            fluidState.pressure(gasPhaseIdx);
        return
            H2O::gasDensity(T, pH2O) +
            Air::gasDensity(T, pAir) +
            NAPL::gasDensity(T, pNAPL);
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == waterPhaseIdx) {
            // assume pure water viscosity

            return H2O::liquidViscosity(T,
                                        p);
        }
        else if (phaseIdx == naplPhaseIdx) {
            // assume pure NAPL viscosity
            return NAPL::liquidViscosity(T, p);
        }

        assert (phaseIdx == gasPhaseIdx);

        /* Wilke method. See:
         *
         * See: R. Reid, et al.: The Properties of Gases and Liquids,
         * 4th edition, McGraw-Hill, 1987, 407-410
         * 5th edition, McGraw-Hill, 20001, p. 9.21/22
         *
         * in this case, we use a simplified version in order to avoid
         * computationally costly evaluation of sqrt and pow functions and
         * divisions
         * -- compare e.g. with Promo Class p. 32/33
         */
        Scalar muResult;
        const Scalar mu[numComponents] = {
            H2O::gasViscosity(T, H2O::vaporPressure(T)),
            Air::gasViscosity(T, p),
            NAPL::gasViscosity(T, NAPL::vaporPressure(T))
        };
        // molar masses
        const Scalar M[numComponents] = {
            H2O::molarMass(),
            Air::molarMass(),
            NAPL::molarMass()
        };

        Scalar muAW = mu[airIdx]*fluidState.moleFraction(gasPhaseIdx, airIdx)
            + mu[H2OIdx]*fluidState.moleFraction(gasPhaseIdx, H2OIdx)
            / (fluidState.moleFraction(gasPhaseIdx, airIdx)
               + fluidState.moleFraction(gasPhaseIdx, H2OIdx));
        Scalar xAW = fluidState.moleFraction(gasPhaseIdx, airIdx)
            + fluidState.moleFraction(gasPhaseIdx, H2OIdx);

        Scalar MAW = (fluidState.moleFraction(gasPhaseIdx, airIdx)*Air::molarMass()
                      + fluidState.moleFraction(gasPhaseIdx, H2OIdx)*H2O::molarMass())
            / xAW;

        Scalar phiCAW = 0.3; // simplification for this particular system
        /* actually like this
         * Scalar phiCAW = std::pow(1.+std::sqrt(mu[NAPLIdx]/muAW)*std::pow(MAW/M[NAPLIdx],0.25),2)
         *                 / std::sqrt(8.*(1.+M[NAPLIdx]/MAW));
         */
        Scalar phiAWC = phiCAW * muAW*M[NAPLIdx]/(mu[NAPLIdx]*MAW);

        muResult = (xAW*muAW)/(xAW+fluidState.moleFraction(gasPhaseIdx, NAPLIdx)*phiAWC)
            + (fluidState.moleFraction(gasPhaseIdx, NAPLIdx) * mu[NAPLIdx])
            / (fluidState.moleFraction(gasPhaseIdx, NAPLIdx) + xAW*phiCAW);
        return muResult;
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       const ParameterCache &paramCache,
                                       int phaseIdx,
                                       int compIdx)
    {
        return 0;
#if 0
        Scalar T = fluidState.temperature(phaseIdx) ;
        Scalar p = fluidState.pressure(phaseIdx);
        Scalar diffCont;

        if (phaseIdx==gasPhaseIdx) {
            Scalar diffAC = Opm::BinaryCoeff::Air_Mesitylene::gasDiffCoeff(T, p);
            Scalar diffWC = Opm::BinaryCoeff::H2O_Mesitylene::gasDiffCoeff(T, p);
            Scalar diffAW = Opm::BinaryCoeff::H2O_Air::gasDiffCoeff(T, p);

            const Scalar xga = fluidState.moleFraction(gasPhaseIdx, airIdx);
            const Scalar xgw = fluidState.moleFraction(gasPhaseIdx, H2OIdx);
            const Scalar xgc = fluidState.moleFraction(gasPhaseIdx, NAPLIdx);

            if (compIdx==NAPLIdx) return (1 - xgw)/(xga/diffAW + xgc/diffWC);
            else if (compIdx==H2OIdx) return (1 - xgc)/(xgw/diffWC + xga/diffAC);
            else if (compIdx==airIdx) OPM_THROW(std::logic_error,
                                                 "Diffusivity of air in the gas phase "
                                                 "is constraint by sum of diffusive fluxes = 0 !\n");
        }
        else if (phaseIdx==waterPhaseIdx){
            Scalar diffACl = 1.e-9; // BinaryCoeff::Air_Mesitylene::liquidDiffCoeff(temperature, pressure);
            Scalar diffWCl = 1.e-9; // BinaryCoeff::H2O_Mesitylene::liquidDiffCoeff(temperature, pressure);
            Scalar diffAWl = 1.e-9; // BinaryCoeff::H2O_Air::liquidDiffCoeff(temperature, pressure);

            Scalar xwa = fluidState.moleFraction(waterPhaseIdx, airIdx);
            Scalar xww = fluidState.moleFraction(waterPhaseIdx, H2OIdx);
            Scalar xwc = fluidState.moleFraction(waterPhaseIdx, NAPLIdx);

            switch (compIdx) {
            case NAPLIdx:
                diffCont = (1.- xww)/(xwa/diffAWl + xwc/diffWCl);
                return diffCont;
            case airIdx:
                diffCont = (1.- xwc)/(xww/diffWCl + xwa/diffACl);
                return diffCont;
            case H2OIdx:
                OPM_THROW(std::logic_error,
                           "Diffusivity of water in the water phase "
                           "is constraint by sum of diffusive fluxes = 0 !\n");
            };
        }
        else if (phaseIdx==naplPhaseIdx) {
            OPM_THROW(std::logic_error,
                       "Diffusion coefficients of "
                       "substances in liquid phase are undefined!\n");
        }
        return 0;
#endif
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(p);

        if (phaseIdx == waterPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            else if (compIdx == airIdx)
                return Opm::BinaryCoeff::H2O_N2::henry(T)/p;
            else if (compIdx == NAPLIdx)
                return Opm::BinaryCoeff::H2O_Mesitylene::henry(T)/p;
            assert(false);
        }
        // for the NAPL phase, we assume currently that nothing is
        // dissolved. this means that the affinity of the NAPL
        // component to the NAPL phase is much higher than for the
        // other components, i.e. the fugacity cofficient is much
        // smaller.
        else if (phaseIdx == naplPhaseIdx) {
            Scalar phiNapl = NAPL::vaporPressure(T)/p;
            if (compIdx == NAPLIdx)
                return phiNapl;
            else if (compIdx == airIdx)
                return 1e6*phiNapl;
            else if (compIdx == H2OIdx)
                return 1e6*phiNapl;
            assert(false);
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        assert(phaseIdx == gasPhaseIdx);
        return 1.0;
    }


    //! \copydoc BaseFluidSystem::enthalpy
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const ParameterCache &paramCache,
                           int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx) ;
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == waterPhaseIdx) {
            return H2O::liquidEnthalpy(T, p);
        }
        else if (phaseIdx == naplPhaseIdx) {
            return NAPL::liquidEnthalpy(T, p);
        }
        else if (phaseIdx == gasPhaseIdx) {
            // gas phase enthalpy depends strongly on composition
            Scalar result = 0;
            result += H2O::gasEnthalpy(T, p) * fluidState.massFraction(gasPhaseIdx, H2OIdx);
            result += NAPL::gasEnthalpy(T, p) * fluidState.massFraction(gasPhaseIdx, airIdx);
            result += Air::gasEnthalpy(T, p) * fluidState.massFraction(gasPhaseIdx, NAPLIdx);

            return result;
        }
        OPM_THROW(std::logic_error, "Invalid phase index " << phaseIdx);
    }

    //! \copydoc BaseFluidSystem::thermalConductivity
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx) ;
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == waterPhaseIdx){ // water phase
            return H2O::liquidThermalConductivity(T, p);
        }
        else if (phaseIdx == gasPhaseIdx) { // gas phase
            Scalar lambdaDryAir = Air::gasThermalConductivity(T, p);
            return lambdaDryAir;
        }

        assert(phaseIdx == naplPhaseIdx);

        // Taken from:
        //
        // D. K. H. Briggs: "Thermal Conductivity of Liquids",
        // Ind. Eng. Chem., 1957, 49 (3), pp 418–421
        //
        // Convertion to SI units:
        // 344e-6 cal/(s cm K) = 0.0143964 J/(s m K)
        return 0.0143964;
    }
};
} // namespace FluidSystems
} // namespace Opm

#endif
