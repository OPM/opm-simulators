// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
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
 * \ingroup Fluidsystems
 * \brief A fluid system with a liquid and a gaseous phase and \f$H_2O\f$ and \f$N_2\f$
 *        as components.
 */
#ifndef DUMUX_H2O_N2_SYSTEM_HH
#define DUMUX_H2O_N2_SYSTEM_HH

#include <dumux/material/idealgas.hh>
#include <dumux/material/settablephase.hh>

#include <dumux/material/binarycoefficients/h2o_n2.hh>

#include "defaultcomponents.hh"

#include <dumux/common/propertysystem.hh>

namespace Dumux
{
// forward defintions of the property tags
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Components);
};

/*!
 * \ingroup Fluidsystems
 *
 * \brief A compositional fluid with water and molecular nitrogen as
 *        components in both, the liquid and the gas phase.
 */
template <class TypeTag>
class H2O_N2_System
{
    typedef H2O_N2_System<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef Dumux::IdealGas<Scalar> IdealGas;
    typedef Dumux::SettablePhase<Scalar, ThisType> SettablePhase;

    typedef typename GET_PROP(TypeTag, PTAG(Components)) Components;

public:
    typedef typename Components::H2O H2O;
    typedef typename Components::N2 N2;

    static const int numComponents = 2;
    static const int numPhases = 2;

    static const int lPhaseIdx = 0; // index of the liquid phase
    static const int gPhaseIdx = 1; // index of the gas phase

    static const int wPhaseIdx = lPhaseIdx; // index of the wetting phase
    static const int nPhaseIdx = gPhaseIdx; // index of the non-wetting phase

    static const int H2OIdx = 0;
    static const int N2Idx = 1;

    static void init()
    { Components::init(); }

    /*!
     * \brief Return the human readable name of a phase
     *
     * \param phaseIdx index of the phase
     */
    static const char *phaseName(int phaseIdx)
    {
        switch (phaseIdx) {
        case lPhaseIdx: return "l";
        case gPhaseIdx: return "g";
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx index of the component
     */
    static const char *componentName(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::name();
        case N2Idx: return N2::name();
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
        switch (compIdx) {
        case H2OIdx: return H2O::molarMass();
        case N2Idx: return N2::molarMass();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density [kg/m^3].
     *
     * \param phaseIdx index of the phase
     * \param temperature phase temperature in [K]
     * \param pressure phase pressure in [Pa]
     * \param fluidState the fluid state
     *
     * \tparam FluidState the fluid state class
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &fluidState)
    {
        if (phaseIdx == lPhaseIdx) {
            // See: Ochs 2008
            // \todo: proper citation
            Scalar rholH2O = H2O::liquidDensity(temperature, pressure);
            Scalar clH2O = rholH2O/H2O::molarMass();

            // this assumes each nitrogen molecule displaces exactly one
            // water molecule in the liquid
            return
                clH2O*(H2O::molarMass()*fluidState.moleFrac(lPhaseIdx, H2OIdx)
                       +
                       N2::molarMass()*fluidState.moleFrac(lPhaseIdx, N2Idx));
        }
        else if (phaseIdx == gPhaseIdx) {
            Scalar fugH2O =
                fluidState.moleFrac(gPhaseIdx, H2OIdx)  *
                fluidState.phasePressure(gPhaseIdx);
            Scalar fugN2 =
                fluidState.moleFrac(gPhaseIdx, N2Idx)  *
                fluidState.phasePressure(gPhaseIdx);
            return
                H2O::gasDensity(temperature, fugH2O) +
                N2::gasDensity(temperature, fugN2);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its viscosity.
     */
    template <class FluidState>
    static Scalar phaseViscosity(int phaseIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 const FluidState &fluidState)
    {
        if (phaseIdx == lPhaseIdx) {
            // assume pure water for the liquid phase
            // TODO: viscosity of mixture
            return H2O::liquidViscosity(temperature,
                                        pressure);
        }
        else {
            return N2::gasViscosity(temperature,
                                    pressure);
            /* Wilke method. See:
             *
             * See: R. Reid, et al.: The Properties of Gases and Liquids,
             * 4th edition, McGraw-Hill, 1987, 407-410
             * 5th edition, McGraw-Hill, 20001, p. 9.21/22
             */
            Scalar muResult = 0;
            const Scalar mu[numComponents] = {
                H2O::gasViscosity(temperature,
                                  H2O::vaporPressure(temperature)),
                N2::gasViscosity(temperature,
                                 pressure)
            };
            // molar masses
            const Scalar M[numComponents] = {
                H2O::molarMass(),
                N2::molarMass()
            };

            for (int i = 0; i < numComponents; ++i) {
                Scalar divisor = 0;
                for (int j = 0; j < numComponents; ++j) {
                    Scalar phiIJ = 1 + sqrt(mu[i]/mu[j] *
                                            pow(M[j]/M[i], 1/4.0));
                    phiIJ *= phiIJ;
                    phiIJ /= sqrt(8*(1 + M[i]/M[j]));
                    divisor += fluidState.moleFrac(phaseIdx, j)*phiIJ;
                }
                muResult += fluidState.moleFrac(phaseIdx, i)*mu[i] / divisor;
            }

            return muResult;
        }
    }


    /*!
     * \brief Assuming the composition of a single phase and the
     *        pressure of all phases is known or that all phases are
     *        present, compute the thermodynamic equilibrium from the
     *        temperature and phase pressures. If the known phase
     *        index
     *
     */
    template <class FluidState>
    static void computeEquilibrium(FluidState &fluidState,
                                   int knownPhaseIdx = -1)
    {
        const Scalar T = fluidState.temperature();
        const Scalar pg = fluidState.phasePressure(gPhaseIdx);
        const Scalar pl = fluidState.phasePressure(lPhaseIdx);

        const Scalar betaH2O = H2O::vaporPressure(T);
        const Scalar betaN2 = BinaryCoeff::H2O_N2::henry(T);

        if (knownPhaseIdx < 0)
        {
            // we only have all phase pressures and temperature and
            // know that all phases are present
            Scalar xlH2O = (pg - betaN2)/(betaH2O - betaN2);
            Scalar xlN2 = 1 - xlH2O;
            Scalar rhol = liquidPhaseDensity_(T, pl, xlH2O, xlN2);

            Scalar cgH2O = H2O::gasDensity(T, betaH2O*xlH2O)/H2O::molarMass();
            Scalar cgN2 = N2::gasDensity(T, betaN2*xlN2)/N2::molarMass();

            Scalar xgH2O = cgH2O/(cgH2O + cgN2);
            Scalar xgN2 = cgN2/(cgH2O + cgN2);

            // set the liquid phase composition
            SettablePhase liquid;
            liquid.moleFrac_[H2OIdx] = xlH2O;
            liquid.moleFrac_[N2Idx] = xlN2;
            liquid.pressure_ = pl;
            liquid.density_ = rhol;
            liquid.xToX(); // compute mass fractions from mole fractions
            fluidState.assignPhase(lPhaseIdx, liquid);

            // set the gas phase composition
            SettablePhase gas;
            gas.moleFrac_[H2OIdx] = xgH2O;
            gas.moleFrac_[N2Idx] = xgN2;
            gas.pressure_ = pg;
            gas.density_ = cgH2O*H2O::molarMass() + cgN2*N2::molarMass();
            gas.xToX(); // compute mass fractions from mole fractions
            fluidState.assignPhase(gPhaseIdx, gas);
        }
        else if (knownPhaseIdx == lPhaseIdx) {
            // the composition of the liquid phase is given

            // retrieve the known mole fractions from the fluid state
            Scalar xlH2O = fluidState.moleFrac(lPhaseIdx, H2OIdx);
            Scalar xlN2 = fluidState.moleFrac(lPhaseIdx, N2Idx);

            // calculate the component contentrations in the gas phase
            Scalar pH2O = betaH2O*xlH2O; // fugacity of water
            Scalar pN2 = betaN2*xlN2; // fugacity of nitrogen
            Scalar cgH2O = H2O::gasDensity(T, pH2O)/H2O::molarMass();
            Scalar cgN2 = N2::gasDensity(T, pN2)/N2::molarMass();

            // convert concentrations to mole fractions
            Scalar xgH2O = cgH2O/(cgH2O + cgN2) * (pH2O + pN2)/pg;
            Scalar xgN2 = cgN2/(cgH2O + cgN2) * (pH2O + pN2)/pg;

            // set gas phase composition
            SettablePhase gas;
            gas.moleFrac_[H2OIdx] = xgH2O;
            gas.moleFrac_[N2Idx] = xgN2;
            gas.pressure_ = pg;
            gas.density_ = cgH2O*H2O::molarMass() + cgN2*N2::molarMass();
            gas.xToX(); // update mass fractions from mole fractions
            fluidState.assignPhase(gPhaseIdx, gas);
        }
        else if (knownPhaseIdx == gPhaseIdx) {
            // the composition of the gas phase is given

            Scalar xgH2O = fluidState.moleFrac(gPhaseIdx, H2OIdx);
            Scalar xgN2 = fluidState.moleFrac(gPhaseIdx, N2Idx);
            Scalar pgH2O = pg*xgH2O;
            Scalar pgN2 = pg*xgN2;

            Scalar xlH2O = pgH2O/betaH2O;
            Scalar xlN2 = pgN2/betaN2;

            SettablePhase liquid;
            liquid.moleFrac_[H2OIdx] = xlH2O;
            liquid.moleFrac_[N2Idx] = xlN2;
            liquid.pressure_ = pl;
            liquid.density_ = liquidPhaseDensity_(T, pl, xlH2O, xlN2);
            liquid.xToX(); // update mass fractions from mole fractions
            fluidState.assignPhase(lPhaseIdx, liquid);
        }
    }

    /*!
     * \brief Returns the activity coefficient of a component in a
     *        phase.
     *
     * We define the activity coefficent \f$\gamma_{\alpha,\kappa}\f$
     * of component \f$\kappa\f$ by the following equation:
     *  \f[ f_\kappa = p_\alpha \gamma_{\alpha,\kappa} \f]
     * where \f$f_\kappa\f$  is the component's fugacity and \f$p_\alpha\f$
     * is the phase' pressure. To yield a phase-independent name for this
     * method, it also treats gaseous phases although the term "activity"
     * might be misleading in that context.
     *
     * For liquids with very low miscibility this boils down to the
     * inverse Henry constant for the solutes and the partial pressure
     * for the solvent.
     *
     * For ideal gases this is equivalent to the gas pressure, in real
     * gases it is the gas pressure times the component's fugacity
     * coefficient.
     */
    template <class FluidState>
    static Scalar activityCoeff(int phaseIdx,
                                int compIdx,
                                Scalar temperature,
                                Scalar pressure,
                                const FluidState &state)
    {
        if (phaseIdx == gPhaseIdx) {
            return pressure;
            Scalar fugH2O = std::max(1e-3, state.fugacity(H2OIdx));
            Scalar fugN2 = std::max(1e-3, state.fugacity(N2Idx));
            Scalar cH2O = H2O::gasDensity(temperature, fugH2O) / H2O::molarMass();
            Scalar cN2 = N2::gasDensity(temperature, fugN2) / N2::molarMass();

            Scalar alpha = (fugH2O + fugN2)/pressure;

            if (compIdx == H2OIdx)
                return fugH2O/(alpha*cH2O/(cH2O + cN2));
            else if (compIdx == N2Idx)
                return fugN2/(alpha*cN2/(cH2O + cN2));

            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }

        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(temperature);
        case N2Idx: return BinaryCoeff::H2O_N2::henry(temperature);
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     *
     */
    template <class FluidState>
    static Scalar diffCoeff(int phaseIdx,
                            int compIIdx,
                            int compJIdx,
                            Scalar temperature,
                            Scalar pressure,
                            const FluidState &fluidState)
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


        switch (phaseIdx) {
        case lPhaseIdx:
            switch (compIIdx) {
            case H2OIdx:
                switch (compJIdx) {
                case N2Idx: return BinaryCoeff::H2O_N2::liquidDiffCoeff(temperature,
                                                                        pressure);
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
                case N2Idx: return BinaryCoeff::H2O_N2::gasDiffCoeff(temperature,
                                                                     pressure);
                }
            }
        }

        DUNE_THROW(Dune::InvalidStateException,
                   "Binary diffusion coefficient of components "
                   << compIIdx << " and " << compJIdx
                   << " in phase " << phaseIdx << " is undefined!\n");
    };

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific enthalpy [J/kg].
     */
    template <class FluidState>
    static Scalar phaseEnthalpy(int phaseIdx,
                                Scalar temperature,
                                Scalar pressure,
                                const FluidState &fluidState)
    {
        if (phaseIdx == lPhaseIdx) {
            Scalar cN2 = fluidState.concentration(lPhaseIdx, N2Idx);
            Scalar pN2 = N2::gasPressure(temperature, cN2*N2::molarMass());

            // TODO: correct way to deal with the solutes???
            return
                fluidState.massFrac(lPhaseIdx, H2OIdx)*
                H2O::liquidEnthalpy(temperature, pressure)
                +
                fluidState.massFrac(lPhaseIdx, N2Idx)*
                N2::gasEnthalpy(temperature, pN2);
        }
        else {
            Scalar cH2O = fluidState.concentration(gPhaseIdx, H2OIdx);
            Scalar cN2 = fluidState.concentration(gPhaseIdx, N2Idx);

            Scalar pH2O = H2O::gasPressure(temperature, cH2O*H2O::molarMass());
            Scalar pN2 = N2::gasPressure(temperature, cN2*N2::molarMass());

            Scalar result = 0;
            result +=
                H2O::gasEnthalpy(temperature, pH2O) *
                fluidState.massFrac(gPhaseIdx, H2OIdx);
            result +=
                N2::gasEnthalpy(temperature, pN2) *
                fluidState.massFrac(gPhaseIdx, N2Idx);

            return result;
        }
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific internal energy [J/kg].
     */
    template <class FluidState>
    static Scalar phaseInternalEnergy(int phaseIdx,
                                      Scalar temperature,
                                      Scalar pressure,
                                      const FluidState &fluidState)
    {
        return
            phaseEnthalpy(phaseIdx, temperature, pressure, fluidState) -
            pressure/phaseDensity(phaseIdx, temperature, pressure, fluidState);
    }

private:
    static Scalar liquidPhaseDensity_(Scalar T, Scalar pl, Scalar xlH2O, Scalar xlN2)
    {
        // See: Ochs 2008
        // Formula (2.6)
        // \todo: proper citation
        Scalar rholH2O = H2O::liquidDensity(T, pl);
        Scalar clH2O = rholH2O/H2O::molarMass();

        // this assumes each nitrogen molecule displaces exactly one
        // water molecule in the liquid
        return
            clH2O*(xlH2O*H2O::molarMass()
                   +
                   xlN2*N2::molarMass());
    }

    static Scalar gasPhaseDensity_(Scalar T, Scalar pg, Scalar xgH2O, Scalar xgN2)
    {
        return H2O::gasDensity(T, pg*xgH2O) + N2::gasDensity(T, pg*xgN2);
    };
};

} // end namepace

#endif
