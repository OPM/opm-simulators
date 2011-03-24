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
 * \brief Base class for all fluid systems which do not require
 *        mutable parameters.
 */
#ifndef DUMUX_SIMPLE_FLUID_SYSTEM_HH
#define DUMUX_SIMPLE_FLUID_SYSTEM_HH

#include <dumux/material/fluidstates/genericfluidstate.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/idealgas.hh>

#include <dumux/material/binarycoefficients/h2o_n2.hh>

#include <dumux/common/exceptions.hh>

#include "basicmutableparameters.hh"

namespace Dumux
{
/*!
 * \brief Base class for all fluid systems which do not require
 *        mutable parameters.
 */
template <class Scalar, class StaticParametersT, class Implementation>
class SimpleFluidSystem : public StaticParametersT
{
public:
    typedef StaticParametersT StaticParameters;
    typedef SimpleMutableParameters<Scalar, StaticParameters> MutableParameters;
    typedef typename MutableParameters::FluidState FluidState;

public:
    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static void init()
    {}

    static Scalar computeMolarVolume(const MutableParameters &mutParams, int phaseIdx)
    { Implementation::computeMolarVolume_(mutParams.fluidState(), phaseIdx); }

    /*!
     * \brief Calculate the molar volume [m^3/mol] of a fluid phase.
     *
     * This is the method which does not require any mutable
     * parameters and can thus only gets a fluid state object as
     * argument.
     */
    static Scalar computeMolarVolume_(const FluidState &fs, int phaseIdx)
    {
       DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
    };

    static Scalar computePureFugacity(const MutableParameters &params,
                                      int phaseIdx, 
                                      int compIdx)
    {
        const FluidState &fs = params.fluidState();
        Scalar T = fs.temperature(phaseIdx);
        switch (phaseIdx) {
        case SP::lPhaseIdx: 
            switch (compIdx) {
            case SP::H2OIdx: return SP::H2O::vaporPressure(T);
            case SP::N2Idx: return BinaryCoeff::H2O_N2::henry(T);
            };
        case SP::gPhaseIdx: {
            switch (compIdx) {
            case SP::H2OIdx: 
                SP::H2O::vaporPressure(T);
            case SP::N2Idx:
                // ideal gas
                return fs.pressure(phaseIdx);
            };
        };

        default: DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
        }
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
    static Scalar computeFugacity(const MutableParameters &params,
                                  int phaseIdx, 
                                  int compIdx)
    {
        const FluidState &fs = params.fluidState();

        Scalar x = fs.moleFrac(phaseIdx,compIdx);
        Scalar p = fs.pressure(phaseIdx);
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
    static Scalar computeFugacityCoeff(const MutableParameters &params,
                                       int phaseIdx, 
                                       int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= SP::numPhases);
        assert(0 <= compIdx  && compIdx <= SP::numComponents);

        const FluidState &fs = params.fluidState();

        Scalar T = fs.temperature(phaseIdx);
        Scalar p = fs.pressure(phaseIdx);
        switch (phaseIdx) {
        case SP::lPhaseIdx: 
            switch (compIdx) {
            case SP::H2OIdx: return SP::H2O::vaporPressure(T)/p;
            case SP::N2Idx: return BinaryCoeff::H2O_N2::henry(T)/p;
            };
        case SP::gPhaseIdx:
            return 1.0; // ideal gas
        };

        DUNE_THROW(Dune::InvalidStateException, "Unhandled phase or component index");
    }
    
    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase [Pa*s]
     */
    static Scalar computeViscosity(const MutableParameters &params,
                                   int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= SP::numPhases);

        const FluidState &fs = params.fluidState();
        Scalar T = fs.temperature(phaseIdx);
        Scalar p = fs.pressure(phaseIdx);
        switch (phaseIdx) {
        case SP::lPhaseIdx:
            // assume pure water for the liquid phase
            return SP::H2O::liquidViscosity(T, p);
        case SP::gPhaseIdx:
            // assume pure water for the gas phase
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
    static Scalar computeDiffusionCoeff(const MutableParameters &params, 
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
    static Scalar binaryDiffCoeff(const MutableParameters &params, 
                                  int phaseIdx,
                                  int compIIdx,
                                  int compJIdx)
                                  
    {
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

        const FluidState &fs = params.fluidState();
        Scalar T = fs.temperature(phaseIdx);
        Scalar p = fs.pressure(phaseIdx);       

        switch (phaseIdx) {
        case SP::lPhaseIdx:
            switch (compIIdx) {
            case SP::H2OIdx:
                switch (compJIdx) {
                case SP::N2Idx: return BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p);
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
                case SP::N2Idx: return BinaryCoeff::H2O_N2::gasDiffCoeff(T, p);
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
    static Scalar computeEnthalpy(const MutableParameters &params, 
                                  int phaseIdx)
    {
        const FluidState &fs = params.fluidState();
        Scalar T = fs.temperature(phaseIdx);
        Scalar p = fs.pressure(phaseIdx);

        if (phaseIdx == SP::lPhaseIdx) {
#warning hack
            T = 300.0;

            Scalar cN2 = fs.molarity(SP::lPhaseIdx, SP::N2Idx);
            Scalar pN2 = SP::N2::gasPressure(T, cN2*SP::N2::molarMass());

            Scalar XH2O = fs.massFrac(SP::lPhaseIdx, SP::H2OIdx);
            Scalar XN2 = fs.massFrac(SP::lPhaseIdx, SP::N2Idx);
            // TODO: correct way to deal with the solutes???
            return
                (XH2O*SP::H2O::liquidEnthalpy(T, p)
                 +
                 XN2*SP::N2::gasEnthalpy(T, pN2))
                /
                (XH2O + XN2);
        }
        else {
            Scalar cH2O = fs.molarity(SP::gPhaseIdx, SP::H2OIdx);
            Scalar cN2 = fs.molarity(SP::gPhaseIdx, SP::N2Idx);
            
            // assume ideal gas
            Scalar pH2O = SP::H2O::gasPressure(T, cH2O*SP::H2O::molarMass());
            Scalar pN2 = SP::N2::gasPressure(T, cN2*SP::H2O::molarMass());

            Scalar XH2O = fs.massFrac(SP::lPhaseIdx, SP::H2OIdx);
            Scalar XN2 = fs.massFrac(SP::lPhaseIdx, SP::N2Idx);          
            Scalar result = 0;
            result +=
                SP::H2O::gasEnthalpy(T, pH2O) *
                fs.massFrac(SP::gPhaseIdx, SP::H2OIdx);
            result +=
                SP::N2::gasEnthalpy(T, pN2) *
                fs.massFrac(SP::gPhaseIdx, SP::N2Idx);
            result /= XH2O + XN2;

            return result;
        }
    };

    /*!
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific internal energy [J/kg].
     */
    static Scalar internalEnergy(const MutableParameters &params, 
                                 int phaseIdx)
    {
        const FluidState &fs = params.fluidState();
        Scalar T = fs.temperature(phaseIdx);
        Scalar p = fs.pressure(phaseIdx);
        Scalar rho = fs.density(phaseIdx);
        return
            computeEnthalpy(params, phaseIdx) -
            p/rho;

    };
};

} // end namepace

#endif
