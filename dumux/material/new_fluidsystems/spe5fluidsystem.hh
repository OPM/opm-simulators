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

#include "dumux/common/spline.hh"
#include "spe5/spe5staticparameters.hh"
#include "spe5/spe5mutableparameters.hh"
#include "spe5/spe5pengrobinsonparams.hh"
#include "../eos/pengrobinsonmixture.hh"

namespace Dumux
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
class Spe5FluidSystem
{
public:
    typedef Dumux::Spe5StaticParameters<Scalar> StaticParameters;
    typedef Dumux::Spe5MutableParameters<Scalar> MutableParameters;

private:
    typedef typename Dumux::PengRobinsonMixture<Scalar, StaticParameters> PengRobinsonMixture;
    typedef typename Dumux::PengRobinson<Scalar> PengRobinson;

public:
    // copy number of phases and components from the static parameters
    enum { numPhases = StaticParameters::numPhases };
    enum { numComponents = StaticParameters::numComponents };

    // copy phase indices from the static parameters
    enum {
        wPhaseIdx = StaticParameters::wPhaseIdx,
        oPhaseIdx = StaticParameters::oPhaseIdx,
        gPhaseIdx = StaticParameters::gPhaseIdx
    };

    // copy component indices from the static parameters
    enum {
        H2OIdx = StaticParameters::H2OIdx,
        C1Idx = StaticParameters::C1Idx,
        C3Idx = StaticParameters::C3Idx,
        C6Idx = StaticParameters::C6Idx,
        C10Idx = StaticParameters::C10Idx,
        C15Idx = StaticParameters::C15Idx,
        C20Idx = StaticParameters::C20Idx,
    };

    // copy the component types from the static parameters
    typedef typename StaticParameters::H2O H2O;

    /*!
     * \brief Initialize the fluid system's static parameters
     */
    static void init()
    {
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
     * \brief Return the human readable name of a component
     */
    static const char *componentName(int compIdx)
    {
        assert(0 <= compIdx  && compIdx <= numComponents);
        return StaticParameters::componentName(compIdx);
    }

    /*!
     * \brief Return the human readable name of a phase
     */
    static const char *phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);
        return StaticParameters::phaseName(phaseIdx);
    }

    /*!
     * \brief Return the molar mass [kg/mol] of a component
     */
    static Scalar molarMass(int compIdx)
    {
        assert(0 <= compIdx  && compIdx <= numComponents);
        return StaticParameters::molarMass(compIdx);
    }
    
    /*!
     * \brief Calculate the molar volume [m^3/mol] of a fluid phase
     */
    static Scalar computeMolarVolume(MutableParameters &params, 
                                     int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        params.update(phaseIdx);
        return params.molarVolume(phaseIdx);
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
        params.update(phaseIdx);
        return 
            params.moleFrac(phaseIdx, compIdx) *
            params.pressure(phaseIdx) *
            computeFugacityCoeff(params, phaseIdx, compIdx);
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
//        const Scalar phi_g[numComponents] = {0.315214, 0.86048, 0.378601, 0.12922, 0.0320002, 0.00813658, 0.00174178 };
//        const Scalar phi_o[numComponents] = {0.0633375, 1.73469, 0.19746, 0.0147604, 0.000630321, 2.48063e-05, 7.74427e-07 };
//        const Scalar pRefOil = 1.5e+07;

        assert(0 <= phaseIdx  && phaseIdx <= numPhases);
        assert(0 <= compIdx  && compIdx <= numComponents);

        params.update(phaseIdx);
        switch (phaseIdx) {
        case gPhaseIdx:
        case oPhaseIdx: {
            params.update(phaseIdx);             
            Scalar phi = PengRobinsonMixture::computeFugacityCoeff(params,
                                                                   phaseIdx, 
                                                                   compIdx);
            return phi;
        }
        case wPhaseIdx:
            return 
                henryCoeffWater_(compIdx, params.temperature(wPhaseIdx))
                /
                params.pressure(wPhaseIdx);
        default: DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
        }
    }
    
    /*
    static Scalar evalTable_(const Scalar *y, Scalar p)
    {
        Scalar tmp = pIdx_(p);
        int pIdx = static_cast<int>(tmp);

        Scalar p1 = (pIdx + 0) * (pMax_ - pMin_)/Scalar(numEntries_) + pMin_;
        Scalar p2 = (pIdx + 1) * (pMax_ - pMin_)/Scalar(numEntries_) + pMin_;
        
#if 0
        Scalar alpha = tmp - pIdx;
        return y[pIdx]*(1 - alpha) + y[pIdx + 1]*alpha; 
#else
        Scalar pPrimeLeft =  (y[pIdx + 1] - y[pIdx - 1])/( 2*(pMax_ - pMin_)/numEntries_ );
        Scalar pPrimeRight =  (y[pIdx + 2] - y[pIdx]    )/( 2*(pMax_ - pMin_)/numEntries_ );
        Spline<Scalar> sp(p1, p2,
                          y[pIdx], y[pIdx + 1],
                          pPrimeLeft, pPrimeRight);
        return sp.eval(p);
#endif
    }
    */

    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase [Pa*s]
     */
    static Scalar computeViscosity(MutableParameters &params,
                                   int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);
        params.update(phaseIdx);
        
        switch (phaseIdx) {
        case gPhaseIdx: {
            // given by SPE-5 in table on page 64. we use a constant
            // viscosity, though...
            return 0.0170e-2 * 0.1; 
        }
        case wPhaseIdx: 
            // given by SPE-5: 0.7 centi-Poise  = 0.0007 Pa s
            return 0.7e-2 * 0.1;
        case oPhaseIdx: {
            // given by SPE-5 in table on page 64. we use a constant
            // viscosity, though...
            return 0.208e-2 * 0.1; 
        }
        default: DUNE_THROW(Dune::InvalidStateException, "Unhandled phase index " << phaseIdx);
        }
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
     *        calculate its specific enthalpy [J/kg].
     */
    static Scalar computeEnthalpy(MutableParameters &params, 
                                  int phaseIdx)
    {
        // TODO!
        DUNE_THROW(Dune::NotImplemented, "Enthalpies");
    };

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        calculate its specific internal energy [J/kg].
     */
    static Scalar internalEnergy(MutableParameters &params, 
                                 int phaseIdx)
    {
        // TODO!
        DUNE_THROW(Dune::NotImplemented, "Enthalpies");
    };

private:   
    static Scalar henryCoeffWater_(int compIdx, Scalar temperature)
    {
        // use henry's law for the solutes and the vapor pressure for
        // the solvent.
        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(temperature);
            
            // the values of the Henry constant for the solutes have
            // been computed using the Peng-Robinson equation of state
            // (-> slope of the component's fugacity function at
            // almost 100% water content)
        case C1Idx: return 5.57601e+09;
        case C3Idx: return 1.89465e+10;
        case C6Idx: return 5.58969e+12;
        case C10Idx: return 4.31947e+17;
        case C15Idx: return 4.27283e+28;
        case C20Idx: return 3.39438e+36;
        default: DUNE_THROW(Dune::InvalidStateException, "Unknown component index " << compIdx);
        }
    };
};

} // end namepace

#endif
