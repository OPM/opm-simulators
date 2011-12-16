/*****************************************************************************
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
 *   Copyright (C) 2009-2011 by Markus Wolff                                 *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief A fluid system for two-phase models assuming immiscibility and
 *        thermodynamic equilibrium
 *
 * The wetting and the non-wetting phase can be defined via their
 * individual components.
 */
#ifndef DUMUX_2P_IMMISCIBLE_FLUID_SYSTEM_HH
#define DUMUX_2P_IMMISCIBLE_FLUID_SYSTEM_HH

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/gasphase.hh>

#include <dune/common/exceptions.hh>

#include "basefluidsystem.hh"
#include "../MpNcfluidstates/immisciblefluidstate.hh"

namespace Dumux {

/*!
 * \ingroup Fluidsystems
 *
 * \brief A fluid system for two-phase models assuming immiscibility and
 *        thermodynamic equilibrium
 *
 * The wetting and the non-wetting phase can be defined individually
 * via Dumux::LiquidPhase<Component> and
 * Dumux::GasPhase<Component>. These phases consist of one pure
 * component. With the help of this adapter class, the phase
 * properties can be accessed. This is suitable for pure two-phase
 * systems without compositional effects.
 */
template <class Scalar, class WettingPhase, class NonWettingPhase>
class TwoPImmiscibleFluidSystem
: public BaseFluidSystem<Scalar, TwoPImmiscibleFluidSystem<Scalar, WettingPhase, NonWettingPhase> >
{
    // do not try to instanciate this class, it has only static members!
    TwoPImmiscibleFluidSystem()
    {}

    typedef TwoPImmiscibleFluidSystem<Scalar, WettingPhase, NonWettingPhase> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;
public:
    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! Number of phases in the fluid system
    static constexpr int numPhases = 2;

    //! Index of the wetting phase
    static constexpr int wPhaseIdx = 0;
    //! Index of the non-wetting phase
    static constexpr int nPhaseIdx = 1;

    /*!
     * \brief Return the human readable name of a fluid phase
     */
    static const char *phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        static const char *name[] = {
            "w",
            "n"
        };
        return name[phaseIdx];
    }

    /*!
     * \brief Return whether a phase is liquid
     */
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == wPhaseIdx)
            return WettingPhase::isLiquid();
        return NonWettingPhase::isLiquid();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are indepent on the fluid composition. This assumtion is true
     * if immiscibility is assumed. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     */
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // we assume immisibility
        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 2;

    //! Index of the wetting phase's component
    static constexpr int wCompIdx = 0;
    //! Index of the non-wetting phase's component
    static constexpr int nCompIdx = 1;

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx index of the component
     */
    static const char *componentName(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::name();
        return NonWettingPhase::name();
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     *
     * \param compIdx index of the component
     */
    static Scalar molarMass(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::molarMass();
        return NonWettingPhase::molarMass();
    }

    /*!
     * \brief Critical temperature of a component [K].
     */
    static Scalar criticalTemperature(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::criticalTemperature();
        return NonWettingPhase::criticalTemperature();
    };

    /*!
     * \brief Critical pressure of a component [Pa].
     */
    static Scalar criticalPressure(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::criticalPressure();
        return NonWettingPhase::criticalPressure();
    };

    /*!
     * \brief The acentric factor of a component [].
     */
    static Scalar acentricFactor(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (compIdx == wCompIdx)
            return WettingPhase::acentricFactor();
        return NonWettingPhase::acentricFactor();
    };

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static void init()
    {
        // two gaseous phases at once do not make sense physically!
        // (But two liquids are fine)
        assert(WettingPhase::isLiquid() || NonWettingPhase::isLiquid());
    }

    /*!
     * \brief Return the density of a phase [kg/m^3].
     *
     * \param fluidState The fluid state of the two-phase model
     * \param phaseIdx Index of the fluid phase
     *
     * \tparam FluidState the fluid state class of the two-phase model
     */
    using Base::density;
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::density(temperature, pressure);
        return NonWettingPhase::density(temperature, pressure);
    }

    /*!
     * \brief Return the viscosity of a phase [Pa*s].
     *
     */
    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::viscosity(temperature, pressure);
        return NonWettingPhase::viscosity(temperature, pressure);
    }

    /*!
     * \brief Calculate the fugacity coefficient [Pa] of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\phi_\kappa\f$ is connected to the
     * fugacity \f$f_\kappa\f$ and the component's molarity
     * \f$x_\kappa\f$ by means of the relation
     *
     * \f[ f_\kappa = \phi_\kappa * x_{\kappa} \f]
     */
    using Base::fugacityCoefficient;
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        if (phaseIdx == compIdx)
            // TODO (?): calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
            return 1.0;
        return std::numeric_limits<Scalar>::infinity();
    }

    /*!
     * \brief Calculate the binary molecular diffusion coefficient for
     *        a component in a fluid phase [mol^2 * s / (kg*m^3)]
     *
     * Molecular diffusion of a compoent \f$\kappa\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \grad mu_\kappa \f]
     *
     * where \f$\mu_\kappa\f$ is the component's chemical potential,
     * \f$D\f$ is the diffusion coefficient and \f$J\f$ is the
     * diffusive flux. \f$mu_\kappa\f$ is connected to the component's
     * fugacity \f$f_\kappa\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$p_\alpha\f$ and \f$T_\alpha\f$ are the fluid phase'
     * pressure and temperature.
     */
    using Base::diffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Diffusion coefficients of components are meaningless if"
                   " immiscibility is assumed");
    };

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     */
    using Base::binaryDiffusionCoefficient;
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    {
        DUNE_THROW(Dune::InvalidStateException,
                   "Binary diffusion coefficients of components are meaningless if"
                   " immiscibility is assumed");
    };

    /*!
     * \brief Return the specific enthalpy of a fluid phase [J/kg].
     */
    using Base::enthalpy;
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                                 int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::enthalpy(temperature, pressure);
        return NonWettingPhase::enthalpy(temperature, pressure);
    }

    /*!
     * \brief Thermal conductivity of a fluid phase [W/(m^2 K/m)].
     */
    using Base::thermalConductivity;
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::thermalConductivity(temperature, pressure);
        return NonWettingPhase::thermalConductivity(temperature, pressure);
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
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
            return WettingPhase::heatCapacity(temperature, pressure);
        return NonWettingPhase::heatCapacity(temperature, pressure);
    }
};

} // end namepace

#endif
