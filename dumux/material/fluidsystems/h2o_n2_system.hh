// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A fluid system with water and gas as phases and \f$H_2O\f$ and \f$N_2\f$
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
/*!
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
     * \brief Return the human readable name of a component
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
     * \brief Given the gas phase's composition, temperature and
     *        pressure, compute the partial presures of all components
     *        in [Pa] and assign it to the FluidState.
     *
     * This is required for models which cannot calculate the the
     * partial pressures of the components in the gas phase from the
     * degasPressure(). To use this method, the FluidState must have a
     * setPartialPressure(componentIndex, pressure) method.
     */
    template <class FluidState>
    static void computePartialPressures(Scalar temperature,
                                        Scalar pg,
                                        FluidState &fluidState)
    {
        Scalar X1 = fluidState.massFrac(gPhaseIdx, H2OIdx);

        // We use the newton method for this. For the initial value we
        // assume all components to be an ideal gas
        Scalar pH2O = fluidState.moleFrac(gPhaseIdx, H2OIdx) * pg;
        Scalar eps = pg*1e-9;

        Scalar deltaP = pH2O;
        Valgrind::CheckDefined(pH2O);
        Valgrind::CheckDefined(deltaP);
        for (int i = 0; i < 5 && std::abs(deltaP/pg) > 1e-9; ++i) {
            Scalar f =
                H2O::gasDensity(temperature, pH2O)*(X1 - 1) +
                X1*N2::gasDensity(temperature, pg - pH2O);

            Scalar df_dp;
            df_dp  =
                H2O::gasDensity(temperature, pH2O + eps)*(X1 - 1) +
                X1*N2::gasDensity(temperature, pg - (pH2O + eps));
            df_dp -=
                H2O::gasDensity(temperature, pH2O - eps)*(X1 - 1) +
                X1*N2::gasDensity(temperature, pg - (pH2O - eps));
            df_dp /=
                2*eps;

            deltaP = - f/df_dp;

            pH2O += deltaP;
            Valgrind::CheckDefined(pH2O);
            Valgrind::CheckDefined(deltaP);
        }

        fluidState.setPartialPressure(H2OIdx, pH2O);
        fluidState.setPartialPressure(N2Idx, pg - pH2O);
    };

    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density [kg/m^3].
     */
    template <class FluidState>
    static Scalar phaseDensity(int phaseIdx,
                               Scalar temperature,
                               Scalar pressure,
                               const FluidState &fluidState)
    {
        if (phaseIdx == lPhaseIdx)
            return liquidPhaseDensity_(temperature,
                                       pressure,
                                       fluidState.moleFrac(lPhaseIdx, H2OIdx),
                                       fluidState.moleFrac(lPhaseIdx, N2Idx));
        else if (phaseIdx == gPhaseIdx)
            return gasPhaseDensity_(temperature,
                                    pressure,
                                    fluidState.moleFrac(gPhaseIdx, H2OIdx),
                                    fluidState.moleFrac(gPhaseIdx, N2Idx));
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
             * S.O.Ochs: "Development of a multiphase multicomponent
             * model for PEMFC - Technical report: IRTG-NUPUS",
             * University of Stuttgart, 2008
             *
             * and:
             *
             * See: R. Reid, et al.: The Properties of Gases and Liquids, 4th
             * edition, McGraw-Hill, 1987, 407-410
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
                                            pow(M[i]/M[j], 1/4.0));
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

            // check for consistency of the gasDensity_() method
            checkConsistentGasDensity_(gas.density_,
                                       fluidState.phasePressure(gPhaseIdx),
                                       fluidState);
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

            // check for consistency of the gasDensity_() method if
            // the gas phase is present
            if (std::abs((pN2 + pH2O - pg)/pg) < 1e-8)
                checkConsistentGasDensity_(gas.density_,
                                           betaH2O*xlH2O + betaN2*xlN2,
                                           fluidState);
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
     * \brief Returns the fugacity coefficient of a component in a
     *        phase.
     *
     * For solutions with only traces in a solvent this boils down to
     * the inverse Henry constant for the solutes and the partial
     * pressure for the solvent.
     */
    template <class FluidState>
    static Scalar fugacityCoeff(int phaseIdx,
                                int compIdx,
                                Scalar temperature,
                                Scalar pressure,
                                const FluidState &state)
    {
        if (phaseIdx != lPhaseIdx)
            DUNE_THROW(Dune::NotImplemented,
                       "Fugacities in the gas phase");

        switch (compIdx) {
        case H2OIdx: return H2O::vaporPressure(temperature);
        case N2Idx: return BinaryCoeff::H2O_N2::henry(temperature);
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given a component's pressure and temperature, return its
     *        density in a phase [kg/m^3].
     */
    static Scalar componentDensity(int phaseIdx,
                                   int compIdx,
                                   Scalar temperature,
                                   Scalar pressure)
    {
        if (phaseIdx == lPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::liquidDensity(temperature, pressure);
            else if (compIdx == N2Idx)
                return N2::liquidDensity(temperature, pressure);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
        else if (phaseIdx == gPhaseIdx) {
            if (compIdx == H2OIdx) {
                return H2O::gasDensity(temperature, pressure);
            }
            else if (compIdx == N2Idx)
                return N2::gasDensity(temperature, pressure);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    };

    /*!
     * \brief Given a component's density and temperature, return the
     *        corresponding pressure in a phase [Pa].
     */
    static Scalar componentPressure(int phaseIdx,
                                    int compIdx,
                                    Scalar temperature,
                                    Scalar density)
    {
        if (phaseIdx == lPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::liquidPressure(temperature, density);
            else if (compIdx == N2Idx)
                return N2::liquidPressure(temperature, density);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
        else if (phaseIdx == gPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::gasPressure(temperature, density);
            else if (compIdx == N2Idx)
                return N2::gasPressure(temperature, density);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    };

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficent for components
     *        \f$i\f$ and \f$j\f$ in this phase.
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
        if (phaseIdx == lPhaseIdx)  {
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
            Scalar pWater = fluidState.partialPressure(H2OIdx);
            Scalar pN2 = fluidState.partialPressure(N2Idx);

            Scalar result = 0;
            result +=
                H2O::gasEnthalpy(temperature, pWater) *
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

    // defect of gas density
    static inline Scalar gasDefect_(Scalar pH2O, Scalar T, Scalar pg, Scalar xgH2O, Scalar xgN2)
    {
        // this assumes that sum of the fugacities of the individual
        // components add up to the gas pressure
        return
            H2O::gasDensity(T, pH2O)/H2O::molarMass()*(xgH2O - 1) +
            xgH2O*N2::gasDensity(T, pg - pH2O)/N2::molarMass();
    }

    static Scalar gasPhaseDensity_(Scalar T, Scalar pg, Scalar xgH2O, Scalar xgN2)
    {
        // assume ideal gas for the initial condition
        Scalar pH2O = pg*xgH2O;
        Scalar delta = 1e100;

        // newton method. makes sure that the total pressure of the
        // gas phase is the sum of the individual gas phase fugacities
        // of the individual components
        for (int i = 0; i < 10; ++i) {
            Scalar f;
            Scalar df_dpH2O;

            f = gasDefect_(pH2O, T, pg, xgH2O, xgN2);
            Scalar eps = (std::abs(pg) + 1)*1e-9;
            df_dpH2O = gasDefect_(pH2O + eps, T, pg, xgH2O, xgN2);
            df_dpH2O -= gasDefect_(pH2O - eps, T, pg, xgH2O, xgN2);
            df_dpH2O /= 2*eps;

            delta = - f/df_dpH2O;
            pH2O += delta;

            if (std::abs(delta) < 1e-10*pg) {
                return H2O::gasDensity(T, pH2O) + N2::gasDensity(T, pg - pH2O);
            }
        };
        DUNE_THROW(NumericalProblem,
                   "Can not calculate the gas density at "
                   << " T=" << T
                   << " pg=" << pg
                   << " xgH2O=" << xgH2O
                   << " xgN2=" << xgN2
                   << "\n");
    };

    template <class FluidState>
    static void checkConsistentGasDensity_(Scalar rhoToTest,
                                           Scalar pg,
                                           const FluidState &fluidState)
    {
#ifndef NDEBUG
        Scalar rho = gasPhaseDensity_(fluidState.temperature(),
                                      pg,
                                      fluidState.moleFrac(gPhaseIdx, H2OIdx),
                                      fluidState.moleFrac(gPhaseIdx, N2Idx));
        static const Scalar eps = 1e-8;
        if (std::abs(rhoToTest/rho - 1.0) > eps) {
            DUNE_THROW(Dune::InvalidStateException,
                       "Density of gas phase is inconsistent: rho1/rho2 = "
                       << rhoToTest/rho);
        }
#endif
    }
    };

} // end namepace

#endif
