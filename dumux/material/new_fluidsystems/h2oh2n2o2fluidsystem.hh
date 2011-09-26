/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
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
#ifndef DUMUX_H2O_H2_N2_O2_FLUID_SYSTEM_HH
#define DUMUX_H2O_H2_N2_O2_FLUID_SYSTEM_HH

#include <dumux/material/fluidstates/genericfluidstate.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/h2.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/o2.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/idealgas.hh>

#include <dumux/material/binarycoefficients/h2o_h2.hh>
#include <dumux/material/binarycoefficients/h2o_n2.hh>
#include <dumux/material/binarycoefficients/h2o_o2.hh>
#include <dumux/material/binarycoefficients/h2_n2.hh>
#include <dumux/material/binarycoefficients/h2_o2.hh>
#include <dumux/material/binarycoefficients/n2_o2.hh>

#include <dumux/material/fluidstates/genericfluidstate.hh>

#include <dumux/common/exceptions.hh>

namespace Dumux
{


template <class Scalar>
struct H2OH2N2O2StaticParameters {
    /****************************************
     * Fluid phase parameters
     ****************************************/

    //! Number of phases in the fluid system
    static const int numPhases = 2;

    //! Index of the liquid phase
    static const int lPhaseIdx = 0;
    //! Index of the gas phase
    static const int gPhaseIdx = 1;
    
    //! The component for pure water 
    typedef Dumux::SimpleH2O<Scalar> SimpleH2O;
    typedef Dumux::H2O<Scalar> IapwsH2O;
    typedef Dumux::TabulatedComponent<Scalar, IapwsH2O> TabulatedH2O;
    
    //! The component for pure water
    //typedef IapwsH2O H2O;
    //typedef TabulatedH2O H2O;
    typedef SimpleH2O H2O;
    
    //! The component for pure hydrogen
    typedef Dumux::H2<Scalar> H2;

    //! The component for pure nitrogen
    typedef Dumux::N2<Scalar> N2;

    //! The component for pure oxygen
    typedef Dumux::O2<Scalar> O2;

    /*!
     * \brief Initialize the static parameters.
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

    /****************************************
     * Component related parameters
     ****************************************/

    //! Number of components in the fluid system
    static const int numComponents = 4;
    
    static const int H2OIdx = 0;
    static const int H2Idx = 1;
    static const int N2Idx = 2;
    static const int O2Idx = 3;

    /*!
     * \brief Return the human readable name of a component
     */
    static const char *componentName(int compIdx)
    {
        static const char *name[] = {
            H2O::name(),
            H2::name(),
            N2::name(),
            O2::name(),
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
            H2::molarMass(),
            N2::molarMass(),
            O2::molarMass(),
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
            H2::criticalTemperature(), // H2O
            N2::criticalTemperature(), // H2O
            O2::criticalTemperature(), // H2O
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
            H2::criticalPressure(),
            N2::criticalPressure(),
            O2::criticalPressure(),
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
                   "H2OH2N2O2StaticParams::criticalMolarVolume()");
    };

    /*!
     * \brief The acentric factor of a component [].
     */
    static Scalar acentricFactor(int compIdx)
    {
        static const Scalar accFac[] = {
            H2O::acentricFactor(), // H2O (from Reid, et al.)
            H2::acentricFactor(),
            N2::acentricFactor(),
            O2::acentricFactor(),
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return accFac[compIdx];
    };
};

/*!
 * \brief A twophase fluid system with water, hydrogen, nitrogen and oxygen as components.
 */
template <class Scalar>
class H2OH2N2O2FluidSystem : public H2OH2N2O2StaticParameters<Scalar>
{
public:
    typedef Dumux::H2OH2N2O2StaticParameters<Scalar> StaticParameters;
    typedef Dumux::GenericFluidState<Scalar, StaticParameters> MutableParameters;

private:
    // convenience typedefs
    typedef StaticParameters SP;
    typedef Dumux::IdealGas<Scalar> IdealGas;

public:

    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        SP::init(tempMin, tempMax, nTemp,
                 pressMin, pressMax, nPress);
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
        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /*!
     * \brief Calculate the molar volume [m^3/mol] of a fluid phase
     */
    static Scalar computeMolarVolume(MutableParameters &params, 
                                     int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= SP::numPhases);

        params.updateMeanMolarMass(phaseIdx);
        Scalar T = params.temperature(phaseIdx);
        Scalar p = params.pressure(phaseIdx);

        switch (phaseIdx) {
        case SP::lPhaseIdx:
            // assume pure water where one water molecule gets
            // replaced by one nitrogen molecule
            return SP::H2O::molarMass()/SP::H2O::liquidDensity(T, p);
        case SP::gPhaseIdx:
            // assume ideal gas
            return 1.0 / IdealGas::concentration(T, p);
        }
        
        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
    };

    /*!
     * \brief Calculate the fugacity of a component in a fluid phase
     *        [Pa]
     *
     * The components chemical \f$mu_\kappa\f$ potential is connected
     * to the component's fugacity \f$f_\kappa\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$p_\alpha\f$ and \f$T_\alpha\f$ are the fluid phase'
     * pressure and temperature.
     */
    static Scalar computeFugacity(MutableParameters &params,
                                  int phaseIdx, 
                                  int compIdx)
    {
        Scalar x = params.moleFrac(phaseIdx,compIdx);
        Scalar p = params.pressure(phaseIdx);
        return x*p*computeFugacityCoeff(params, phaseIdx, compIdx);
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
    static Scalar computeFugacityCoeff(MutableParameters &params,
                                       int phaseIdx, 
                                       int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= SP::numPhases);
        assert(0 <= compIdx  && compIdx <= SP::numComponents);

        params.updateMeanMolarMass(phaseIdx);
        Scalar T = params.temperature(phaseIdx);
        Scalar p = params.pressure(phaseIdx);
        switch (phaseIdx) {
        case SP::lPhaseIdx: 
            switch (compIdx) {
            case SP::H2OIdx: return SP::H2O::vaporPressure(T)/p;
            case SP::H2Idx: return BinaryCoeff::H2O_H2::henry(T)/p;
            case SP::N2Idx: return BinaryCoeff::H2O_N2::henry(T)/p;
            case SP::O2Idx: return BinaryCoeff::H2O_O2::henry(T)/p;
            };
        case SP::gPhaseIdx:
            return 1.0; // ideal gas
        };

        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase or component index");
    }
    
    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase [Pa*s]
     */
    static Scalar computeViscosity(MutableParameters &params,
                                   int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= SP::numPhases);

        params.updateMeanMolarMass(phaseIdx);
        Scalar T = params.temperature(phaseIdx);
        Scalar p = params.pressure(phaseIdx);
        switch (phaseIdx) {
        case SP::lPhaseIdx:
            // assume pure water for the liquid phase
            return SP::H2O::liquidViscosity(T, p);
        case SP::gPhaseIdx:
            // assume pure nitrogen for the gas phase
            return SP::N2::gasViscosity(T, p);
        }
        
        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
    };

    /*!
     * \brief Calculate the binary molecular diffusion coefficient for
     *        a component in a fluid phase [mol^2 * s / (kg*m^3)]
     *
     * Molecular diffusion of a compoent $\kappa$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \grad mu_\kappa \f] 
     *
     * where \f$\mu_\kappa\$ is the component's chemical potential,
     * \f$D\f$ is the diffusion coefficient and \f$J\f$ is the
     * diffusive flux. \f$mu_\kappa\f$ is connected to the component's
     * fugacity \f$f_\kappa\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$p_\alpha\f$ and \f$T_\alpha\f$ are the fluid phase'
     * pressure and temperature.
     */
    static Scalar computeDiffusionCoeff(MutableParameters &params, 
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
    static Scalar binaryDiffCoeff(MutableParameters &params, 
                                  int phaseIdx,
                                  int compIIdx,
                                  int compJIdx)
                                  
    {
        params.updateMeanMolarMass(phaseIdx);
        if (compIIdx > compJIdx)
            std::swap(compIIdx, compJIdx);

#ifndef NDEBUG
        if (compIIdx == compJIdx ||
            phaseIdx > SP::numPhases - 1 ||
            compJIdx > SP::numComponents - 1)
        {
            DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx
                       << " in phase " << phaseIdx << " is undefined!\n");
        }
#endif

        Scalar T = params.temperature(phaseIdx);
        Scalar p = params.pressure(phaseIdx);       

        switch (phaseIdx) {
        case SP::lPhaseIdx:
            switch (compIIdx) {
            case SP::H2OIdx:
                switch (compJIdx) {
                case SP::H2Idx: return BinaryCoeff::H2O_H2::liquidDiffCoeff(T, p);
                case SP::N2Idx: return BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p);
                case SP::O2Idx: return BinaryCoeff::H2O_O2::liquidDiffCoeff(T, p);
                }
            default:
                DUNE_THROW(Dune::InvalidStateException,
                           "Binary diffusion coefficients of trace "
                           "substances in liquid phase is undefined!\n");
            }
        case SP::gPhaseIdx:
            switch (compIIdx) {
            case SP::H2OIdx:
                switch (compJIdx) {
                case SP::H2Idx: return BinaryCoeff::H2O_H2::gasDiffCoeff(T, p);
                case SP::N2Idx: return BinaryCoeff::H2O_N2::gasDiffCoeff(T, p);
                case SP::O2Idx: return BinaryCoeff::H2O_O2::gasDiffCoeff(T, p);
                }
            
            case SP::H2Idx:
                switch (compJIdx) {
                case SP::N2Idx: return BinaryCoeff::H2_N2::gasDiffCoeff(T, p);
                case SP::O2Idx: return BinaryCoeff::H2_O2::gasDiffCoeff(T, p);
                }
                
            case SP::N2Idx:
                switch (compJIdx) {
                case SP::O2Idx: return BinaryCoeff::N2_O2::gasDiffCoeff(T, p);
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
     */
    static Scalar computeEnthalpy(MutableParameters &params, 
                                  int phaseIdx)
    {
        params.updateMeanMolarMass(phaseIdx);
        Scalar T = params.temperature(phaseIdx);
        Scalar p = params.pressure(phaseIdx);
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(p);
        if (phaseIdx == SP::lPhaseIdx) {
            Scalar cH2 = params.molarity(SP::lPhaseIdx, SP::H2Idx);
            Scalar cN2 = params.molarity(SP::lPhaseIdx, SP::N2Idx);
            Scalar cO2 = params.molarity(SP::lPhaseIdx, SP::O2Idx);
            Scalar pH2 = SP::H2::gasPressure(T, cN2*SP::H2::molarMass());
            Scalar pN2 = SP::N2::gasPressure(T, cN2*SP::N2::molarMass());
            Scalar pO2 = SP::O2::gasPressure(T, cN2*SP::O2::molarMass());

            Scalar XH2O = params.massFrac(SP::lPhaseIdx, SP::H2OIdx);
            Scalar XH2 = params.massFrac(SP::lPhaseIdx, SP::H2Idx);
            Scalar XN2 = params.massFrac(SP::lPhaseIdx, SP::N2Idx);
            Scalar XO2 = params.massFrac(SP::lPhaseIdx, SP::O2Idx);
            
            // TODO: correct way to deal with the solutes???
            return 
                (XH2O*SP::H2O::liquidEnthalpy(T, p) +
                 XH2*SP::N2::gasEnthalpy(T, pH2) + 
                 XN2*SP::N2::gasEnthalpy(T, pN2) + 
                 XO2*SP::N2::gasEnthalpy(T, pO2))
                /
                (XH2O + XN2);
        }
        else {
            Scalar cH2O = params.molarity(SP::gPhaseIdx, SP::H2OIdx);
            Scalar cH2 = params.molarity(SP::gPhaseIdx, SP::H2Idx);
            Scalar cN2 = params.molarity(SP::gPhaseIdx, SP::N2Idx);
            Scalar cO2 = params.molarity(SP::gPhaseIdx, SP::O2Idx);
            
            // assume ideal mixture for gas
            Scalar pH2O = SP::H2O::gasPressure(T, cH2O*SP::H2O::molarMass());
            Scalar pH2 = SP::H2::gasPressure(T, cH2*SP::H2::molarMass());
            Scalar pN2 = SP::N2::gasPressure(T, cN2*SP::N2::molarMass());
            Scalar pO2 = SP::O2::gasPressure(T, cO2*SP::O2::molarMass());

            Scalar XH2O = params.massFrac(SP::gPhaseIdx, SP::H2OIdx);
            Scalar XH2 = params.massFrac(SP::gPhaseIdx, SP::H2Idx);
            Scalar XN2 = params.massFrac(SP::gPhaseIdx, SP::N2Idx);
            Scalar XO2 = params.massFrac(SP::gPhaseIdx, SP::O2Idx);

            // add up all the "partial enthalpies"
            Scalar result = 0;
            result += SP::H2O::gasEnthalpy(T, pH2O);
            result += SP::H2::gasEnthalpy(T, pH2);
            result += SP::N2::gasEnthalpy(T, pN2);
            result += SP::O2::gasEnthalpy(T, pO2);

            // normalize the result to 100%
            result /= XH2O + XH2 + XN2 + XO2;

            return result;
        }
    };

    /*!
     * \brief Thermal conductivity of phases [W m / (m^2 K)]
     * 		Use the cunductivity of air and water as a first approximation.
     * 		Source: http://en.wikipedia.org/wiki/List_of_thermal_conductivities
     */
    static Scalar thermalConductivity(const MutableParameters &params,
    								  const int phaseIdx)
    {
//    	TODO thermal conductivity is a function of:
//        Scalar p = params.pressure(phaseIdx);
//        Scalar T = params.temperature(phaseIdx);
//        Scalar x = params.moleFrac(phaseIdx,compIdx);
        switch (phaseIdx) {
        case SP::lPhaseIdx: // use conductivity of pure water
            return  0.6;   // conductivity of water[W / (m K ) ]
        case SP::gPhaseIdx:// use conductivity of pure air
            return  0.025; // conductivity of air [W / (m K ) ]
        }
        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
    }
};

} // end namepace

#endif
