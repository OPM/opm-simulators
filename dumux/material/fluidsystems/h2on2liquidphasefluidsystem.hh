// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Katherina Baber
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
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
 * \ingroup Fluidsystems
 *
 * \brief A one-phase (water-phase) fluid system with water and nitrogen as components.
 */
#ifndef DUMUX_H2O_N2_LIQUIDPHASE_FLUID_SYSTEM_HH
#define DUMUX_H2O_N2_LIQUIDPHASE_FLUID_SYSTEM_HH

#include "basefluidsystem.hh"
#include "nullparametercache.hh"

#include <dumux/material/idealgas.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/binarycoefficients/h2o_n2.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/exceptions.hh>

#include <iostream>
#include <assert.h>

namespace Dumux
{
namespace FluidSystems
{

/*!
 * \ingroup Fluidsystems
 *
 * \brief A one-phase (water-phase) fluid system with water and nitrogen as components.
 *
 * This FluidSystem can be used without the PropertySystem that is applied in Dumux,
 * as all Parameters are defined via template parameters. Hence it is in an
 * additional namespace Dumux::FluidSystems::.
 */
template <class Scalar, bool useComplexRelations = true>
class H2ON2LiquidPhase
    : public BaseFluidSystem<Scalar, H2ON2LiquidPhase<Scalar, useComplexRelations> >
{
    typedef H2ON2LiquidPhase<Scalar, useComplexRelations> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

    // convenience typedefs
    typedef Dumux::H2O<Scalar> IapwsH2O;
    typedef Dumux::TabulatedComponent<Scalar, IapwsH2O > TabulatedH2O;
    typedef Dumux::N2<Scalar> SimpleN2;

public:
    typedef NullParameterCache ParameterCache;
 
    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! Number of phases in the fluid system
    static constexpr int numPhases = 1;

    //! Index of the liquid phase
    static constexpr int lPhaseIdx = 0;

    /*!
     * \brief Return the human readable name of a fluid phase.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static const char *phaseName(int phaseIdx)
    {
        static const char *name[] = {
            "l"
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    /*!
     * \brief Return whether a phase is liquid.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isLiquid(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return true; //only water phase present
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
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
        // the water component decides for the liquid phase...
        return H2O::liquidIsCompressible();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealGas(int phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return false; // not a gas (only liquid phase present)
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 2;

    static constexpr int H2OIdx = 0;
    static constexpr int N2Idx = 1;

    //! The components for pure water
    typedef TabulatedH2O H2O;
    //typedef SimpleH2O H2O;
    //typedef IapwsH2O H2O;

    //! The components for pure nitrogen
    typedef SimpleN2 N2;

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static const char *componentName(int compIdx)
    {
        static const char *name[] = {
            H2O::name(),
            N2::name()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     *
     * \param compIdx The index of the component to consider
     */
    static constexpr Scalar molarMass(int compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx == H2OIdx)
            ? H2O::molarMass()
            : (compIdx == N2Idx)
            ? N2::molarMass()
            : 1e100;
    }

    /*!
     * \brief Critical temperature of a component [K].
     *
     * \param compIdx The index of the component to consider
     */
    static constexpr Scalar criticalTemperature(int compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx == H2OIdx)
            ? H2O::criticalTemperature()
            : (compIdx == N2Idx)
            ? N2::criticalTemperature()
            : 1e100;
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \param compIdx The index of the component to consider
     */
    static constexpr Scalar criticalPressure(int compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx == H2OIdx)
            ? H2O::criticalPressure()
            : (compIdx == N2Idx)
            ? N2::criticalPressure()
            : 1e100;
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \param compIdx The index of the component to consider
     */
    static constexpr Scalar acentricFactor(int compIdx)
    {
        //assert(0 <= compIdx && compIdx < numComponents);
        return (compIdx == H2OIdx)
            ? H2O::acentricFactor()
            : (compIdx == N2Idx)
            ? N2::acentricFactor()
            : 1e100;
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters generically.
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
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
     *        problem specific temperature and pressure ranges.
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
            std::cout << "Using complex H2O-N2 fluid system\n";
        else
            std::cout << "Using fast H2O-N2 fluid system\n";

        if (H2O::isTabulated) {
            std::cout << "Initializing tables for the H2O fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            TabulatedH2O::init(tempMin, tempMax, nTemp,
                               pressMin, pressMax, nPress);
        }
    }

    /*!
     * \brief Calculate the density [kg/m^3] of a fluid phase.
     *
     * If useComplexRelations == true, we apply
     * Formula (2.6) from S.O.Ochs:
     * "Development of a multi-phase multicomponent
     * model for PEMFC - Technical report: IRTG-NUPUS",
     * University of Stuttgart, 2008
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        Scalar sumMoleFrac = 0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            sumMoleFrac += fluidState.moleFraction(phaseIdx, compIdx);

        // liquid phase
        if (phaseIdx == lPhaseIdx) {
            if (!useComplexRelations)
                // assume pure water
                return H2O::liquidDensity(T, p);
            else
            {
                // See: Ochs 2008
                Scalar rholH2O = H2O::liquidDensity(T, p);
                Scalar clH2O = rholH2O/H2O::molarMass();

                // this assumes each nitrogen molecule displaces exactly one
                // water molecule in the liquid
                return
                    clH2O
                    * (H2O::molarMass()*fluidState.moleFraction(lPhaseIdx, H2OIdx)
                       +
                       N2::molarMass()*fluidState.moleFraction(lPhaseIdx, N2Idx))
                    / sumMoleFrac;
            }
        }
        else //gas phase index
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase [Pa*s].
     *
     * \param fluidState An arbitrary fluid state
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

        // liquid phase
        if (phaseIdx == lPhaseIdx) {
            // assume pure water for the liquid phase
            return H2O::liquidViscosity(T, p);
        }
       else //gas phase index
           DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Calculate the fugacity coefficient [Pa] of an individual
     *        component in a fluid phase.
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
     * The quantity "fugacity" itself is just an other way to express
     * the chemical potential \f$\zeta^\kappa_\alpha\f$ of the
     * component. It is defined via
     *
     * \f[
     f^\kappa_\alpha := \exp\left\{\frac{\zeta^\kappa_\alpha}{k_B T_\alpha} \right\}
     \f]
     * where \f$k_B = 1.380\cdot10^{-23}\;J/K\f$ is the Boltzmann constant.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
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

        // liquid phase
        if (phaseIdx == lPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            return Dumux::BinaryCoeff::H2O_N2::henry(T)/p;
        }
        else //gas phase index
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
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

        // liquid phase
        if (phaseIdx == lPhaseIdx) {
            if (compIIdx == H2OIdx && compJIdx == N2Idx)
                return BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p);
            return undefined;
        }
        else //gas phase index
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific enthalpy [J/kg].
     *
     *  \todo This fluid system neglects the contribution of
     *        gas-molecules in the liquid phase. This contribution is
     *        probably not big. Somebody would have to find out the
     *        enthalpy of solution for this system. ...
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
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

        // liquid phase
        if (phaseIdx == lPhaseIdx) {
            // TODO: correct way to deal with the solutes???
            return H2O::liquidEnthalpy(T, p);
        }
        else //gas phase index
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Thermal conductivity of a fluid phase [W/(m K)].
     *
     * Use the conductivity of water as a first approximation.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        if (phaseIdx == lPhaseIdx){// liquid phase
            if(useComplexRelations){
                Scalar temperature  = fluidState.temperature(phaseIdx) ;
                Scalar pressure = fluidState.pressure(phaseIdx);
                return H2O::liquidThermalConductivity(temperature, pressure);
            }
            else
                return  0.578078;   // conductivity of water[W / (m K ) ] IAPWS evaluated at p=.1 MPa, T=8C
        }
        else //gas phase index
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/kg]}\f$.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const ParameterCache &paramCache,
                               int phaseIdx)
    {
        if (phaseIdx == lPhaseIdx) {
            return H2O::liquidHeatCapacity(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
        }
        else //gas phase index
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }
};

} // end namepace FluidSystems

} // end namepace

#endif
