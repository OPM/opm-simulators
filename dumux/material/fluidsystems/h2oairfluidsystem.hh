// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2012 by Christoph Grueninger                              *
 *   Copyright (C) 2011-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2011-2012 by Klaus Mosthaf                                *
 *   Copyright (C) 2011 by Philipp Nuske                                     *
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
 * \ingroup Fluidsystems
 *
 * \brief A fluid system with a liquid and a gaseous phase and \f$H_2O\f$ and \f$Air\f$
 *        as components.
 */
#ifndef DUMUX_H2O_AIR_SYSTEM_HH
#define DUMUX_H2O_AIR_SYSTEM_HH

#include "basefluidsystem.hh"
#include "nullparametercache.hh"

#include <dumux/material/idealgas.hh>
#include <dumux/material/binarycoefficients/h2o_air.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/exceptions.hh>

#include <iostream>
#include <cassert>

namespace Dumux
{
namespace FluidSystems
{

/*!
 * \ingroup Fluidsystems
 *
 * \brief A compositional twophase fluid system with water and air as
 *        components in both, the liquid and the gas phase.
 *
 *  This fluidsystem is applied by default with the tabulated version of
 *  water of the IAPWS-formulation.
 */
template <class Scalar,
          //class H2Otype = Dumux::SimpleH2O<Scalar>,
          class H2Otype = Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar> >,
          bool useComplexRelations = true>
class H2OAir
    : public BaseFluidSystem<Scalar, H2OAir<Scalar, H2Otype, useComplexRelations> >
{
    typedef H2OAir<Scalar,H2Otype, useComplexRelations > ThisType;
    typedef BaseFluidSystem <Scalar, ThisType> Base;
    typedef Dumux::IdealGas<Scalar> IdealGas;

public:
    typedef NullParameterCache ParameterCache;

    typedef H2Otype H2O;
    typedef Dumux::Air<Scalar> Air;

    static constexpr int numPhases = 2;

    static constexpr int lPhaseIdx = 0; // index of the liquid phase
    static constexpr int gPhaseIdx = 1; // index of the gas phase

    /*!
     * \brief Return the human readable name of a phase
     *
     * \param phaseIdx index of the phase
     */
    static const char *phaseName(int phaseIdx)
    {
        switch (phaseIdx) {
        case lPhaseIdx: return "liquid";
        case gPhaseIdx: return "gas";
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return whether a phase is liquid
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isLiquid(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gPhaseIdx;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are indepent on the fluid composition. This assumption is true
     * if Henry's law and Rault's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealMixture(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
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

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealGas(int phaseIdx)
    {
        return
            (phaseIdx == gPhaseIdx)
            ? H2O::gasIsIdeal() && Air::gasIsIdeal()
            : false;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 2;

    static constexpr int H2OIdx = 0;
    static constexpr int AirIdx = 1;

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx index of the component
     */
    static const char *componentName(int compIdx)
    {
        switch (compIdx)
        {
        case H2OIdx: return H2O::name();
        case AirIdx: return Air::name();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component [kg/mol].
     *
     * \param compIdx index of the component
     */
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
     * \brief Initialize the fluid system's static parameters generically
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
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

    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density [kg/m^3].
     *
     * Formula (2.6)
     * in
     * S.O.Ochs: "Development of a multiphase multicomponent
     * model for PEMFC - Technical report: IRTG-NUPUS",
     * University of Stuttgart, 2008
     *
     *
     * \param phaseIdx index of the phase
     * \param temperature phase temperature in [K]
     * \param pressure phase pressure in [Pa]
     * \param fluidState the fluid state
     *
     * \tparam FluidState the fluid state class
     */
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


    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase [Pa*s]
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
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

    /*!
     * \brief Returns the fugacity coefficient [-] of a component in a
     *        phase.
     *
     * The fugacity coefficient \f$\phi^\kappa_\alpha\f$ of
     * component \f$\kappa\f$ in phase \f$\alpha\f$ is connected to
     * the fugacity \f$f^\kappa_\alpha\f$ and the component's mole
     * fraction \f$x^\kappa_\alpha\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     \f]
     * where \f$p_\alpha\f$ is the pressure of the fluid phase.
     *
     * For liquids with very low miscibility this boils down to the
     * inverse Henry constant for the solutes and the saturated vapor pressure
     * both divided by phase pressure.
     */
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
            return Dumux::BinaryCoeff::H2O_Air::henry(T)/p;
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the first component to consider
     * \param compJIdx The index of the second component to consider
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             const ParameterCache &paramCache,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        if (compIIdx > compJIdx)
            std::swap(compIIdx, compJIdx);

#ifndef NDEBUG
        if (compIIdx == compJIdx ||
            phaseIdx > numPhases - 1 ||
            compJIdx > numComponents - 1)
        {
            DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx
                       << " in phase " << phaseIdx << " is undefined!\n");
        }
#endif

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        switch (phaseIdx)
        {
        case lPhaseIdx:
            switch (compIIdx) {
            case H2OIdx:
                switch (compJIdx) {
                case AirIdx:
                    return BinaryCoeff::H2O_Air::liquidDiffCoeff(T,
                                                                 p);
                }
            default:
                DUNE_THROW(Dune::InvalidStateException,
                           "Binary diffusion coefficients of trace "
                           "substances in liquid phase is undefined!\n");
            }
        case gPhaseIdx:
            switch (compIIdx){
            case H2OIdx:
                switch (compJIdx){
                case AirIdx:
                    return BinaryCoeff::H2O_Air::gasDiffCoeff(T,
                                                              p);
                }
            }
        }

        DUNE_THROW(Dune::InvalidStateException,
                   "Binary diffusion coefficient of components "
                   << compIIdx << " and " << compJIdx
                   << " in phase " << phaseIdx << " is undefined!\n");
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific enthalpy [J/kg].
     *
     * See:
     * Class Class 2000
     * Theorie und numerische Modellierung nichtisothermer Mehrphasenprozesse in NAPL-kontaminierten porösen Medien
     * Chapter 2.1.13 Innere Energie, Wäremekapazität, Enthalpie
     *
     * Formula (2.42):
     * the specifiv enthalpy of a gasphase result from the sum of (enthalpies*mass fraction) of the components
     */
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

    /*!
     * \brief Thermal conductivity of a fluid phase [W/(m K)].
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
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
