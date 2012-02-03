// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Benjamin Faigle                                   *
 *   Copyright (C) 2011 by Klaus Mosthaf                                     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 *
 * \ingroup Fluidsystems
 * \brief A fluid system with a liquid and a gaseous phase and \f$H_2O\f$ and \f$Air\f$
 *        as components.
 */
#ifndef DUMUX_H2O_AIR_SYSTEM_HH
#define DUMUX_H2O_AIR_SYSTEM_HH

#include <dumux/material/idealgas.hh>

#include <dumux/material/fluidsystems/basefluidsystem.hh>

#include <dumux/material/binarycoefficients/h2o_air.hh>
#include <dumux/material/fluidsystems/defaultcomponents.hh>
#include <dumux/material/components/air.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/common/exceptions.hh>

#include <assert.h>

#ifdef DUMUX_PROPERTIES_HH
#include <dumux/common/propertysystem.hh>
#include <dumux/common/basicproperties.hh>
#endif

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
 *
 *  To change the component formulation (i.e. to use nontabulated or
 *  incompressible water), or to switch on verbosity of tabulation,
 *  specify the water formulation via template arguments or via the property
 *  system, as described in the TypeTag Adapter at the end of the file.
 *
            // Select fluid system
            SET_PROP(TestDecTwoPTwoCProblem, FluidSystem)
            {
                typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
                typedef Dumux::FluidSystems::H2OAir<Scalar, Dumux::SimpleH2O<Scalar> > type;
            };

 *   Also remember to initialize tabulated components (FluidSystem::init()), while this
 *   is not necessary for non-tabularized ones.
 *
 * This FluidSystem can be used without the PropertySystem that is applied in Dumux,
 * as all Parameters are defined via template parameters. Hence it is in an
 * additional namespace Dumux::FluidSystem::.
 * An adapter class using Dumux::FluidSystem<TypeTag> is also provided
 * at the end of this file.
 */
template <class Scalar,
          class H2Otype = Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar> >,
          bool useComplexRelations = true>
class H2OAir
: public BaseFluidSystem<Scalar, H2OAir<Scalar, H2Otype, useComplexRelations> >
{
    typedef H2OAir<Scalar,H2Otype, useComplexRelations > ThisType;
    typedef BaseFluidSystem <Scalar, ThisType> Base;

    typedef Dumux::IdealGas<Scalar> IdealGas;

public:
    typedef H2Otype H2O;
    typedef Dumux::Air<Scalar> Air;

    static constexpr int numPhases = 2;

    static constexpr int lPhaseIdx = 0; // index of the liquid phase
    static constexpr int gPhaseIdx = 1; // index of the gas phase

    static constexpr int wPhaseIdx = lPhaseIdx; // index of the wetting phase
    static constexpr int nPhaseIdx = gPhaseIdx; // index of the non-wetting phase

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
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
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
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
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
    static bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // ideal gases are always compressible
        if (phaseIdx == gPhaseIdx)
            return true;
        // the water component decides for the liquid phase...
        return H2O::liquidIsCompressible();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluids decide
        if (phaseIdx == gPhaseIdx)
            return H2O::gasIsIdeal() && Air::gasIsIdeal();
        return false; // not a gas
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
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx)
        {
        case H2OIdx: return H2O::molarMass();
        case AirIdx: return Air::molarMass();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }


    /*!
     * \brief Critical temperature of a component [K].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(int compIdx)
    {
        static const Scalar Tcrit[] = {
            H2O::criticalTemperature(),
            Air::criticalTemperature()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return Tcrit[compIdx];
    };

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(int compIdx)
    {
        static const Scalar pcrit[] = {
            H2O::criticalPressure(),
            Air::criticalPressure()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return pcrit[compIdx];
    };

    /*!
     * \brief Molar volume of a component at the critical point [m^3/mol].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "H2OAirFluidSystem::criticalMolarVolume()");
        return 0;
    };

    /*!
     * \brief The acentric factor of a component [].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(int compIdx)
    {
        static const Scalar accFac[] = {
            H2O::acentricFactor(), // H2O (from Reid, et al.)
            Air::acentricFactor()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return accFac[compIdx];
    };

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
    using Base::density;
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
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
            //return 1000.;
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
            //return 1.189;
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
    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == lPhaseIdx)
        {
            return 0.0013;
            // assume pure water for the liquid phase
            // TODO: viscosity of mixture
            // couldn't find a way to solve the mixture problem
            return H2O::liquidViscosity(T, p);
        }
        else if (phaseIdx == gPhaseIdx)
        {
            return 1.8e-5;
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
                        Scalar phiIJ = 1 + sqrt(mu[i]/mu[j]) * // 1 + (mu[i]/mu[j]^1/2
                            pow(M[j]/M[i], 1./4.0);   // (M[i]/M[j])^1/4

                        phiIJ *= phiIJ;
                        phiIJ /= sqrt(8*(1 + M[i]/M[j]));
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
    using Base::fugacityCoefficient;
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
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

    using Base::diffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OAir::diffusionCoefficient()");
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
    using Base::binaryDiffusionCoefficient;
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        static Scalar undefined(1e10);
        Valgrind::SetUndefined(undefined);

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
    /*!
     *  \todo This system neglects the contribution of gas-molecules in the liquid phase.
     *        This contribution is probably not big. Somebody would have to find out the enthalpy of solution for this system. ...
     */
    using Base::enthalpy;
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
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
     * Use the conductivity of air and water as a first approximation.
     * Source:
     * http://en.wikipedia.org/wiki/List_of_thermal_conductivities
     */
    using Base::thermalConductivity;
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OAir::thermalConductivity()");
    }

    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/kg]}\f$.
     *
     * \param params    mutable parameters
     * \param phaseIdx  for which phase to give back the heat capacity
     */
    using Base::heatCapacity;
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OAir::heatCapacity()");
    }
};

} // end namepace FluidSystems

#ifdef DUMUX_PROPERTIES_HH
// forward defintions of the property tags
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Components);
}

/*!
 * \brief A twophase fluid system with water and air as components.
 *
 * This is an adapter to use Dumux::H2OAirFluidSystem<TypeTag>, as is
 * done with most other classes in Dumux.
 *  This fluidsystem is applied by default with the tabulated version of
 *  water of the IAPWS-formulation.
 *
 *  To change the component formulation (ie to use nontabulated or
 *  incompressible water), or to switch on verbosity of tabulation,
 *  use the property system and the property "Components":
 *
        // Select desired version of the component
        SET_PROP(myApplicationProperty, Components) : public GET_PROP(TypeTag, DefaultComponents)
        {
            typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

        // Do not use the defaults !
        //    typedef Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar> > H2O;

        // Apply e.g. untabulated water:
        typedef Dumux::H2O<Scalar> H2O;
        };

 *   Also remember to initialize tabulated components (FluidSystem::init()), while this
 *   is not necessary for non-tabularized ones.
 */
template<class TypeTag>
class H2OAirFluidSystem
: public FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar),
                              typename GET_PROP(TypeTag, Components)::H2O,
                             GET_PROP_VALUE(TypeTag, EnableComplicatedFluidSystem)>
{};
#endif

} // end namepace

#endif
