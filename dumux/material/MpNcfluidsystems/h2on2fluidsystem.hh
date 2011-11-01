/*****************************************************************************
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
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
 * \brief A twophase fluid system with water and nitrogen as components.
 */
#ifndef DUMUX_H2O_N2_FLUID_SYSTEM_HH
#define DUMUX_H2O_N2_FLUID_SYSTEM_HH

#include <dumux/material/MpNcfluidstates/genericfluidstate.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/idealgas.hh>

#include <dumux/material/binarycoefficients/h2o_n2.hh>
#include <dumux/material/MpNcfluidstates/genericfluidstate.hh>

#include <dumux/common/exceptions.hh>

#include "nullparametercache.hh"

namespace Dumux
{
/*!
 * \brief A twophase fluid system with water and nitrogen as components.
 */
template <class Scalar>
class H2ON2FluidSystem
{
    // convenience typedefs
    typedef Dumux::IdealGas<Scalar> IdealGas;

    typedef Dumux::H2O<Scalar> IapwsH2O;
    typedef Dumux::SimpleH2O<Scalar> SimpleH2O;

    typedef Dumux::TabulatedComponent<Scalar, IapwsH2O > TabulatedH2O;

    typedef Dumux::N2<Scalar> SimpleN2;

public:
    //! The type of parameter cache objects
    typedef Dumux::NullParameterCache ParameterCache;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! Number of phases in the fluid system
    static constexpr int numPhases = 2;

    //! Index of the liquid phase
    static constexpr int lPhaseIdx = 0;
    //! Index of the gas phase
    static constexpr int gPhaseIdx = 1;
    
    //! The components for pure water
    typedef TabulatedH2O H2O;
    //typedef SimpleH2O H2O;
    //typedef IapwsH2O H2O;

    //! The components for pure nitrogen
    typedef SimpleN2 N2;
 
    /*!
     * \brief Return the human readable name of a fluid phase
     */
    static const char *phaseName(int phaseIdx)
    {
        static const char *name[] = {
            "l",
            "g"          
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
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 2;
    
    static constexpr int H2OIdx = 0;
    static constexpr int N2Idx = 1;

    /*!
     * \brief Return the human readable name of a component
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
     */
    static Scalar molarMass(int compIdx)
    {
        static const Scalar M[] = {
            H2O::molarMass(),
            N2::molarMass(),
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
            N2::criticalTemperature(), // H2O
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
            H2O::criticalPressure(),
            N2::criticalPressure()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return pcrit[compIdx];
    };

    /*!
     * \brief Molar volume of a component at the critical point [m^3/mol].
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "H2ON2StaticParams::criticalMolarVolume()");
    };

    /*!
     * \brief The acentric factor of a component [].
     */
    static Scalar acentricFactor(int compIdx)
    {
        static const Scalar accFac[] = {
            H2O::acentricFactor(), // H2O (from Reid, et al.)
            N2::acentricFactor()
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
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (H2O::isTabulated) {
            std::cout << "Initializing tables for the H2O fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";
            
            TabulatedH2O::init(tempMin, tempMax, nTemp,
                               pressMin, pressMax, nPress);
        }
    }

    /*!
     * \brief Calculate the molar volume [m^3/mol] of a fluid phase
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const ParameterCache &paramCache,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        switch (phaseIdx) {
        case lPhaseIdx:
            // assume pure water where one water molecule gets
            // replaced by one nitrogen molecule
            return H2O::liquidDensity(T, p);
        case gPhaseIdx:
            // assume ideal gas
            return
                IdealGas::molarDensity(T, p) 
                * fluidState.averageMolarMass(gPhaseIdx);
        }
        
        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
    };

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
        switch (phaseIdx) {
        case lPhaseIdx: 
            switch (compIdx) {
            case H2OIdx: return H2O::vaporPressure(T)/p;
            case N2Idx: return BinaryCoeff::H2O_N2::henry(T)/p;
            };
        case gPhaseIdx:
            return 1.0; // ideal gas
        };

        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase or component index");
    }
    
    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase [Pa*s]
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        switch (phaseIdx) {
        case lPhaseIdx:
            // assume pure water for the liquid phase
            return H2O::liquidViscosity(T, p);
        case gPhaseIdx:
            // assume pure water for the gas phase
            return N2::gasViscosity(T, p);
        }
        
        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
    };

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
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       const ParameterCache &paramCache,
                                       int phaseIdx,
                                       int compIdx)
    {
        // TODO!
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    };

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient for components
     *        \f$i\f$ and \f$j\f$ in this phase.
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

        switch (phaseIdx) {
        case lPhaseIdx:
            switch (compIIdx) {
            case H2OIdx:
                switch (compJIdx) {
                case N2Idx: return BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p);
                }
            default:
                DUNE_THROW(Dune::InvalidStateException,
                           "Binary diffusion coefficients of trace "
                           "substances in liquid phase is undefined!\n");
            }
        case gPhaseIdx:
            switch (compIIdx) {
            case H2OIdx:
                switch (compJIdx) {
                case N2Idx: return BinaryCoeff::H2O_N2::gasDiffCoeff(T, p);
                }
            }
        }

        DUNE_THROW(Dune::InvalidStateException,
                   "Binary diffusion coefficient of components "
                   << compIIdx << " and " << compJIdx
                   << " in phase " << phaseIdx << " is undefined!\n");
    };

    /*!
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific enthalpy [J/kg].
     *
     *  \todo This fluid system neglects the contribution of
     *        gas-molecules in the liquid phase. This contribution is
     *        probably not big. Somebody would have to find out the
     *        enthalpy of solution for this system. ...
     */
    template <class FluidState>
    static Scalar internalEnergy(const FluidState &fluidState, 
                                 const ParameterCache &paramCache,
                                 int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(p);
        if (phaseIdx == lPhaseIdx) {
            // TODO: correct way to deal with the solutes???
            return 
                H2O::liquidInternalEnergy(T, p) ;
        }
        else {
            // assume ideal gas
            Scalar XH2O = fluidState.massFrac(gPhaseIdx, H2OIdx);
            Scalar XN2 = fluidState.massFrac(gPhaseIdx, N2Idx);          
            Scalar result = 0;
            result += XH2O*H2O::gasInternalEnergy(T, p);
            result += XN2*N2::gasInternalEnergy(T, p);
            return result;
        }
    }

    /*!
     * \brief Thermal conductivity of a fluid phase [W/(m^2 K/m)].
     *
     * Use the conductivity of air and water as a first approximation.
     * Source:
     * http://en.wikipedia.org/wiki/List_of_thermal_conductivities
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx)
    {
//        TODO thermal conductivity is a function of:
//        Scalar p = fluidState.pressure(phaseIdx);
//        Scalar T = fluidState.temperature(phaseIdx);
//        Scalar x = fluidState.moleFrac(phaseIdx,compIdx);
#warning: so far rough estimates from wikipedia
        switch (phaseIdx) {
        case lPhaseIdx: // use conductivity of pure water
            return  0.6;   // conductivity of water[W / (m K ) ]
        case gPhaseIdx:// use conductivity of pure air
            return 0.025; // conductivity of air [W / (m K ) ]
        }
        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
    }
    
    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/kg]}\f$.
     *
     * \param params    mutable parameters
     * \param phaseIdx  for which phase to give back the heat capacity
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const ParameterCache &paramCache,
                               int phaseIdx)
    {
//        http://en.wikipedia.org/wiki/Heat_capacity
#warning: so far rough estimates from wikipedia
//      TODO heatCapacity is a function of composition.
//        Scalar p = fluidState.pressure(phaseIdx);
//        Scalar T = fluidState.temperature(phaseIdx);
//        Scalar x = fluidState.moleFrac(phaseIdx,compIdx);
        switch (phaseIdx) {
        case lPhaseIdx: // use heat capacity of pure liquid water
            return  4181.3;  // @(25°C) !!!
            /* [J/(kg K)]*/ /* not working because ddgamma_ddtau is not defined*/ /* Dumux::H2O<Scalar>::liquidHeatCap_p(T,
                                             p); */
        case gPhaseIdx:
            return  1003.5 ; // @ (0°C) !!!
            /* [J/(kg K)]*/ /* not working because ddgamma_ddtau is not defined*/ /*Dumux::H2O<Scalar>::gasHeatCap_p(T,
                                          p) ;*/
        }
        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
    }
};

} // end namepace

#endif
