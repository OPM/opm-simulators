// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Klaus Mosthaf                                *
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::FluidSystems::H2OAir
 */
#ifndef EWOMS_H2O_AIR_SYSTEM_HH
#define EWOMS_H2O_AIR_SYSTEM_HH

#include "basefluidsystem.hh"
#include "nullparametercache.hh"

#include <ewoms/material/idealgas.hh>
#include <ewoms/material/binarycoefficients/h2o_air.hh>
#include <ewoms/material/components/air.hh>
#include <ewoms/material/components/h2o.hh>
#include <ewoms/material/components/tabulatedcomponent.hh>
#include <ewoms/common/valgrind.hh>
#include <ewoms/common/exceptions.hh>

#include <iostream>
#include <cassert>

namespace Ewoms {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 *
 * \brief A fluid system with a liquid and a gaseous phase and water and air
 *        as components.
 *
 *  This fluidsystem is applied by default with the tabulated version of
 *  water of the IAPWS-formulation.
 */
template <class Scalar,
          //class H2Otype = Ewoms::SimpleH2O<Scalar>,
          class H2Otype = Ewoms::TabulatedComponent<Scalar, Ewoms::H2O<Scalar> >,
          bool useComplexRelations = true>
class H2OAir
    : public BaseFluidSystem<Scalar, H2OAir<Scalar, H2Otype, useComplexRelations> >
{
    typedef H2OAir<Scalar,H2Otype, useComplexRelations > ThisType;
    typedef BaseFluidSystem <Scalar, ThisType> Base;
    typedef Ewoms::IdealGas<Scalar> IdealGas;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    typedef NullParameterCache ParameterCache;

    //! The type of the water component used for this fluid system
    typedef H2Otype H2O;
    //! The type of the air component used for this fluid system
    typedef Ewoms::Air<Scalar> Air;

    //! \copydoc BaseFluidSystem::numPhases
    static constexpr int numPhases = 2;

    //! The index of the liquid phase
    static constexpr int lPhaseIdx = 0;
    //! The index of the gas phase
    static constexpr int gPhaseIdx = 1;

    //! \copydoc BaseFluidSystem::phaseName
    static const char *phaseName(int phaseIdx)
    {
        switch (phaseIdx) {
        case lPhaseIdx: return "liquid";
        case gPhaseIdx: return "gas";
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static constexpr bool isLiquid(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gPhaseIdx;
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static constexpr bool isCompressible(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return (phaseIdx == gPhaseIdx)
            // ideal gases are always compressible
            ? true
            :
            // the water component decides for the liquid phase...
            H2O::liquidIsCompressible();
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static constexpr bool isIdealGas(int phaseIdx)
    {
        return
            (phaseIdx == gPhaseIdx)
            ? H2O::gasIsIdeal() && Air::gasIsIdeal()
            : false;
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static constexpr bool isIdealMixture(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static constexpr int numComponents = 2;

    //! The index of the water component
    static constexpr int H2OIdx = 0;
    //! The index of the air component
    static constexpr int AirIdx = 1;

    //! \copydoc BaseFluidSystem::componentName
    static const char *componentName(int compIdx)
    {
        switch (compIdx)
        {
        case H2OIdx: return H2O::name();
        case AirIdx: return Air::name();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    //! \copydoc BaseFluidSystem::molarMass
    static constexpr Scalar molarMass(int compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::molarMass()
            : (compIdx == AirIdx)
            ? Air::molarMass()
            : 1e100;
   }


    /*!
     * \brief Critical temperature of a component [K].
     *
     * \param compIdx The index of the component to consider
     */
    static constexpr Scalar criticalTemperature(int compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::criticalTemperature()
            : (compIdx == AirIdx)
            ? Air::criticalTemperature()
            : 1e100;
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \param compIdx The index of the component to consider
     */
    static constexpr Scalar criticalPressure(int compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::criticalPressure()
            : (compIdx == AirIdx)
            ? Air::criticalPressure()
            : 1e100;
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \param compIdx The index of the component to consider
     */
    static constexpr Scalar acentricFactor(int compIdx)
    {
        return
            (compIdx == H2OIdx)
            ? H2O::acentricFactor()
            : (compIdx == AirIdx)
            ? Air::acentricFactor()
            : 1e100;
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \copydoc BaseFluidSystem::init
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        if (H2O::isTabulated)
            init(/*tempMin=*/273.15,
                 /*tempMax=*/623.15,
                 /*numTemp=*/100,
                 /*pMin=*/-10,
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
        if (useComplexRelations)
            std::cout << "Using complex H2O-Air fluid system\n";
        else
            std::cout << "Using fast H2O-Air fluid system\n";

        if (H2O::isTabulated) {
            std::cout << "Initializing tables for the H2O fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            H2O::init(tempMin, tempMax, nTemp,
                               pressMin, pressMax, nPress);
        }
    }

    //! \copydoc BaseFluidSystem::density
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p;
        if (isCompressible(phaseIdx))
            p = fluidState.pressure(phaseIdx);
        else {
            // random value which will hopefully cause things to blow
            // up if it is used in a calculation!
            p = - 1e100;
            Valgrind::SetUndefined(p);
        }


        Scalar sumMoleFrac = 0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            sumMoleFrac += fluidState.moleFraction(phaseIdx, compIdx);

        if (phaseIdx == lPhaseIdx)
        {
            if (!useComplexRelations)
                // assume pure water
                return H2O::liquidDensity(T, p);
            else
            {
                // See: Ochs 2008 (2.6)
                Scalar rholH2O = H2O::liquidDensity(T, p);
                Scalar clH2O = rholH2O/H2O::molarMass();

                return
                    clH2O
                    * (H2O::molarMass()*fluidState.moleFraction(lPhaseIdx, H2OIdx)
                           +
                           Air::molarMass()*fluidState.moleFraction(lPhaseIdx, AirIdx))
                   / sumMoleFrac;
            }
        }
        else if (phaseIdx == gPhaseIdx)
        {
            if (!useComplexRelations)
                // for the gas phase assume an ideal gas
                return
                    IdealGas::molarDensity(T, p)
                    * fluidState.averageMolarMass(gPhaseIdx)
                    / std::max(1e-5, sumMoleFrac);

            Scalar partialPressureH2O =
                fluidState.moleFraction(gPhaseIdx, H2OIdx)  *
                fluidState.pressure(gPhaseIdx);

            Scalar partialPressureAir =
                fluidState.moleFraction(gPhaseIdx, AirIdx)  *
                fluidState.pressure(gPhaseIdx);

            return
                H2O::gasDensity(T, partialPressureH2O) +
                Air::gasDensity(T, partialPressureAir);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == lPhaseIdx)
        {
            // assume pure water for the liquid phase
            // TODO: viscosity of mixture
            // couldn't find a way to solve the mixture problem
            return H2O::liquidViscosity(T, p);
        }
        else if (phaseIdx == gPhaseIdx)
        {
            if(!useComplexRelations){
                return Air::gasViscosity(T, p);
            }
            else //using a complicated version of this fluid system
            {
                /* Wilke method. See:
                 *
                 * See: R. Reid, et al.: The Properties of Gases and Liquids,
                 * 4th edition, McGraw-Hill, 1987, 407-410 or
                 * 5th edition, McGraw-Hill, 2000, p. 9.21/22
                 *
                 */

                Scalar muResult = 0;
                const Scalar mu[numComponents] = {
                    H2O::gasViscosity(T,
                                      H2O::vaporPressure(T)),
                    Air::gasViscosity(T, p)
                };

                // molar masses
                const Scalar M[numComponents] =  {
                    H2O::molarMass(),
                    Air::molarMass()
                };

                for (int i = 0; i < numComponents; ++i)
                {
                    Scalar divisor = 0;
                    for (int j = 0; j < numComponents; ++j)
                    {
                        Scalar phiIJ = 1 + std::sqrt(mu[i]/mu[j]) * // 1 + (mu[i]/mu[j]^1/2
                            std::pow(M[j]/M[i], 1./4.0);   // (M[i]/M[j])^1/4

                        phiIJ *= phiIJ;
                        phiIJ /= std::sqrt(8*(1 + M[i]/M[j]));
                        divisor += fluidState.moleFraction(phaseIdx, j)*phiIJ;
                    }
                    muResult += fluidState.moleFraction(phaseIdx, i)*mu[i] / divisor;
                }
                return muResult;
            }
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
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

        if (phaseIdx == lPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            return Ewoms::BinaryCoeff::H2O_Air::henry(T)/p;
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             const ParameterCache &paramCache,
                                             int phaseIdx,
                                             int compIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == lPhaseIdx)
            return BinaryCoeff::H2O_Air::liquidDiffCoeff(T, p);

        assert(phaseIdx == gPhaseIdx);
        return BinaryCoeff::H2O_Air::gasDiffCoeff(T, p);
    }

    //! \copydoc BaseFluidSystem::enthalpy
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const ParameterCache &paramCache,
                           int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(p);

        if (phaseIdx == lPhaseIdx)
        {
            // TODO: correct way to deal with the solutes???
            return H2O::liquidEnthalpy(T, p);
        }

        else if (phaseIdx == gPhaseIdx)
        {
            Scalar result = 0.0;
            result +=
                H2O::gasEnthalpy(T, p) *
                fluidState.massFraction(gPhaseIdx, H2OIdx);

            result +=
                Air::gasEnthalpy(T, p) *
                fluidState.massFraction(gPhaseIdx, AirIdx);
            return result;
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    //! \copydoc BaseFluidSystem::thermalConductivity
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature  = fluidState.temperature(phaseIdx) ;
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == lPhaseIdx){// liquid phase
            return H2O::liquidThermalConductivity(temperature, pressure);
        }
        else{// gas phase
            Scalar lambdaDryAir = Air::gasThermalConductivity(temperature, pressure);

            if (useComplexRelations){
                Scalar xAir = fluidState.moleFraction(phaseIdx, AirIdx);
                Scalar xH2O = fluidState.moleFraction(phaseIdx, H2OIdx);
                Scalar lambdaAir = xAir * lambdaDryAir;

                // Assuming Raoult's, Daltons law and ideal gas
                // in order to obtain the partial density of water in the air phase
                Scalar partialPressure  = pressure * xH2O;

                Scalar lambdaH2O =
                    xH2O
                    * H2O::gasThermalConductivity(temperature, partialPressure);
                return lambdaAir + lambdaH2O;
            }
            else
                return lambdaDryAir; // conductivity of Nitrogen [W / (m K ) ]
        }
    }
};

} // end namepace FluidSystems

} // end namepace

#endif
