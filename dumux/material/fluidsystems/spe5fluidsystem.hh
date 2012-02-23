// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \brief The mixing rule for the oil and the gas phases of the SPE5 problem.
 *
 * This problem comprises \f$H_2O\f$, \f$C_1\f$, \f$C_3\f$, \f$C_6\f$,
 * \f$C_10\f$, \f$C_15\f$ and \f$C_20\f$ as components.
 *
 * See:
 *
 * J.E. Killough, et al.: Fifth Comparative Solution Project:
 * Evaluation of Miscible Flood Simulators, Ninth SPE Symposium on
 * Reservoir Simulation, 1987
 */
#ifndef DUMUX_SPE5_FLUID_SYSTEM_HH
#define DUMUX_SPE5_FLUID_SYSTEM_HH

#include "spe5parametercache.hh"

#include <dumux/common/spline.hh>
#include <dumux/material/eos/pengrobinsonmixture.hh>

namespace Dumux
{
namespace FluidSystems
{
/*!
 * \brief The fluid system for the SPE-5 benchmark problem.
 *
 * This problem comprises \f$H_2O\f$, \f$C_1\f$, \f$C_3\f$, \f$C_6\f$,
 * \f$C_10\f$, \f$C_15\f$ and \f$C_20\f$ as components.
 *
 * See:
 *
 * J.E. Killough, et al.: Fifth Comparative Solution Project:
 * Evaluation of Miscible Flood Simulators, Ninth SPE Symposium on
 * Reservoir Simulation, 1987
 */
template <class Scalar>
class Spe5
{
    typedef Dumux::FluidSystems::Spe5<Scalar> ThisType;

    typedef typename Dumux::PengRobinsonMixture<Scalar, ThisType> PengRobinsonMixture;
    typedef typename Dumux::PengRobinson<Scalar> PengRobinson;

public:
    typedef Dumux::Spe5ParameterCache<Scalar, ThisType> ParameterCache;

    /****************************************
     * Fluid phase parameters
     ****************************************/

    //! Number of phases in the fluid system
    static const int numPhases = 3;

    //! Index of the gas phase
    static const int gPhaseIdx = 0;
    //! Index of the water phase
    static const int wPhaseIdx = 1;
    //! Index of the oil phase
    static const int oPhaseIdx = 2;

    //! The component for pure water to be used
    typedef Dumux::H2O<Scalar> H2O;

    /*!
     * \brief Return the human readable name of a fluid phase
     */
    static const char *phaseName(int phaseIdx)
    {
        static const char *name[] = {
            "g",
            "w",
            "o",
        };

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    /*!
     * \brief Return whether a phase is liquid
     */
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gPhaseIdx;
    }

    /*!
     * \brief Return whether a phase is compressible
     *
     * In the SPE-5 problems all fluids are compressible...
     */
    static bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are indepent on the fluid composition. This assumtion is true
     * if Henry's law and Rault's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     */
    static bool isIdealMixture(int phaseIdx)
    {
        // always use the reference oil for the fugacity coefficents,
        // so they cannot be dependent on composition and they the
        // phases thus always an ideal mixture
        return phaseIdx == wPhaseIdx;
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
        return false; // gas is not ideal here!
    }

    /****************************************
     * Component related parameters
     ****************************************/

    //! Number of components in the fluid system
    static const int numComponents = 7;

    static const int H2OIdx = 0;
    static const int C1Idx = 1;
    static const int C3Idx = 2;
    static const int C6Idx = 3;
    static const int C10Idx = 4;
    static const int C15Idx = 5;
    static const int C20Idx = 6;

    /*!
     * \brief Return the human readable name of a component
     */
    static const char *componentName(int compIdx)
    {
        static const char *name[] = {
            H2O::name(),
            "C1",
            "C3",
            "C6",
            "C10",
            "C15",
            "C20"
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     */
    static Scalar molarMass(int compIdx)
    {
        static const Scalar M[] = {
            H2O::molarMass(),
            16.04e-3, // C1
            44.10e-3, // C3
            86.18e-3, // C6
            142.29e-3, // C10
            206.00e-3, // C15
            282.00e-3 // C20
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return M[compIdx];
    }

    /*!
     * \brief Critical temperature of a component [K].
     */
    static Scalar criticalTemperature(int compIdx)
    {
        static const Scalar Tcrit[] = {
            H2O::criticalTemperature(), // H2O
            343.0*5/9, // C1
            665.7*5/9, // C3
            913.4*5/9, // C6
            1111.8*5/9, // C10
            1270.0*5/9, // C15
            1380.0*5/9 // C20
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return Tcrit[compIdx];
    };

    /*!
     * \brief Critical pressure of a component [Pa].
     */
    static Scalar criticalPressure(int compIdx)
    {
        static const Scalar pcrit[] = {
            H2O::criticalPressure(), // H2O
            667.8 * 6894.7573, // C1
            616.3 * 6894.7573, // C3
            436.9 * 6894.7573, // C6
            304.0 * 6894.7573, // C10
            200.0 * 6894.7573, // C15
            162.0 * 6894.7573  // C20
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return pcrit[compIdx];
    };

    /*!
     * \brief Molar volume of a component at the critical point [m^3/mol].
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        static const Scalar R = 8.314472;
        static const Scalar vcrit[] = {
            H2O::criticalMolarVolume(), // H2O
            0.290*R*criticalTemperature(1)/criticalPressure(1), // C1
            0.277*R*criticalTemperature(2)/criticalPressure(2), // C3
            0.264*R*criticalTemperature(3)/criticalPressure(3), // C6
            0.257*R*criticalTemperature(4)/criticalPressure(4), // C10
            0.245*R*criticalTemperature(5)/criticalPressure(5), // C15
            0.235*R*criticalTemperature(6)/criticalPressure(6) // C20
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return vcrit[compIdx];
    };

    /*!
     * \brief The acentric factor of a component [].
     */
    static Scalar acentricFactor(int compIdx)
    {
        static const Scalar accFac[] = {
            0.344, // H2O (from Reid, et al.)
            0.0130, // C1
            0.1524, // C3
            0.3007, // C6
            0.4885, // C10
            0.6500, // C15
            0.8500  // C20
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return accFac[compIdx];
    };

    /*!
     * \brief Returns the interaction coefficient for two components.
     *
     * The values are given by the SPE5 paper.
     */
    static Scalar interactionCoefficient(int comp1Idx, int comp2Idx)
    {
        int i = std::min(comp1Idx, comp2Idx);
        int j = std::max(comp1Idx, comp2Idx);
        if (i == C1Idx && (j == C15Idx || j == C20Idx))
            return 0.05;
        else if (i == C3Idx && (j == C15Idx || j == C20Idx))
            return 0.005;
        return 0;
    }

    /****************************************
     * Methods which compute stuff
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static void init()
    {
        Dumux::PengRobinsonParamsMixture<Scalar, ThisType, gPhaseIdx, /*useSpe5=*/true> prParams;

        // find envelopes of the 'a' and 'b' parameters for the range
        // 273.15K <= T <= 373.15K and 10e3 Pa <= p <= 100e6 Pa. For
        // this we take advantage of the fact that 'a' and 'b' for
        // mixtures is just a convex combination of the attractive and
        // repulsive parameters of the pure components
        Scalar minT = 273.15, maxT = 373.15;
        Scalar minP = 1e4, maxP = 100e6;

        Scalar minA = 1e100, maxA = -1e100;
        Scalar minB = 1e100, maxB = -1e100;

        prParams.updatePure(minT, minP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(maxT, minP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(minT, maxP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(maxT, maxP);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        PengRobinson::init(/*aMin=*/minA, /*aMax=*/maxA, /*na=*/100,
                           /*bMin=*/minB, /*bMax=*/maxB, /*nb=*/200);
    }

    /*!
     * \brief Calculate the density [kg/m^3] of a fluid phase
     *
     *
     * \param fluidState An abitrary fluid state
     * \param paramCache Container for cache parameters
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        return fluidState.averageMolarMass(phaseIdx)/paramCache.molarVolume(phaseIdx);
    }

    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase [Pa*s]
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fs,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        if (phaseIdx == gPhaseIdx) {
            // given by SPE-5 in table on page 64. we use a constant
            // viscosity, though...
            return 0.0170e-2 * 0.1;
        }
        else if (phaseIdx == wPhaseIdx)
            // given by SPE-5: 0.7 centi-Poise  = 0.0007 Pa s
            return 0.7e-2 * 0.1;
        else {
            assert(phaseIdx == oPhaseIdx);
            // given by SPE-5 in table on page 64. we use a constant
            // viscosity, though...
            return 0.208e-2 * 0.1;
        }
    }

    /*!
     * \brief Calculate the fugacity coefficient [Pa] of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\phi^\kappa_\alpha\f$ is connected
     * to the fugacity \f$f^\kappa_\alpha\f$ and the component's mole
     * fraction in a phase \f$x^\kappa_\alpha\f$ by means of the
     * relation
     *
     * \f[ f^\kappa_\alpha = \phi^\kappa_\alpha \cdot x^\kappa_\alpha \cdot p_alpha \f]
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fs,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);
        assert(0 <= compIdx  && compIdx <= numComponents);

        if (phaseIdx == oPhaseIdx || phaseIdx == gPhaseIdx)
            return PengRobinsonMixture::computeFugacityCoefficient(fs,
                                                                   paramCache,
                                                                   phaseIdx,
                                                                   compIdx);
        else {
            assert(phaseIdx == wPhaseIdx);
            return
                henryCoeffWater_(compIdx, fs.temperature(wPhaseIdx))
                / fs.pressure(wPhaseIdx);
        }
    }


    /*!
     * \brief Calculate the binary molecular diffusion coefficient for
     *        a component in a fluid phase [mol^2 * s / (kg*m^3)]
     *
     * Molecular diffusion of a compoent \f$\kappa\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D grad \mu_\kappa \f]
     *
     * where \f$\mu_\kappa\f$ is the component's chemical potential,
     * \f$D\f$ is the diffusion coefficient and \f$J\f$ is the
     * diffusive flux. \f$\mu_\kappa\f$ is connected to the component's
     * fugacity \f$f_\kappa\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$p_\alpha\f$ and \f$T_\alpha\f$ are the fluid phase'
     * pressure and temperature.
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fs,
                                       const ParameterCache &paramCache,
                                       int phaseIdx,
                                       int compIdx)
    { DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients"); }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     *
     * \param fluidState An abitrary fluid state
     * \param paramCache Container for cache parameters
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
    { DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficients"); }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        calculate its specific enthalpy [J/kg].
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fs,
                           const ParameterCache &paramCache,
                           int phaseIdx)
    { DUNE_THROW(Dune::NotImplemented, "Enthalpies"); }

    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx)
    { DUNE_THROW(Dune::NotImplemented, "Thermal conductivities"); }

    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const ParameterCache &paramCache,
                               int phaseIdx)
    { DUNE_THROW(Dune::NotImplemented, "Heat capacities"); }


private:
    static Scalar henryCoeffWater_(int compIdx, Scalar temperature)
    {
        // use henry's law for the solutes and the vapor pressure for
        // the solvent.
        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(temperature);

            // the values of the Henry constant for the solutes have
            // are faked so far...
        case C1Idx: return 5.57601e+09;
        case C3Idx: return 1e10;
        case C6Idx: return 1e10;
        case C10Idx: return 1e10;
        case C15Idx: return 1e10;
        case C20Idx: return 1e10;
        default: DUNE_THROW(Dune::InvalidStateException, "Unknown component index " << compIdx);
        }
    };
};

} // end namepace
} // end namepace

#endif
