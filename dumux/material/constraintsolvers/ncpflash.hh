// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \brief Determines the phase compositions, pressures and saturations
 *        given the total mass of all components.
 */
#ifndef DUMUX_NCP_FLASH_HH
#define DUMUX_NCP_FLASH_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {

/*!
 * \brief Determines the phase compositions, pressures and saturations
 *        given the total mass of all components.
 */
template <class Scalar, class FluidSystem>
class NcpFlash
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    typedef typename FluidSystem::ParameterCache ParameterCache;

    // In a M-phase, N-component context, we have the following
    // unknowns:
    //
    // - M pressures
    // - M saturations
    // - M*N mole fractions
    //
    // This sums up to M*(N + 2). On the equations side of things,
    // we have:
    //
    // - (M - 1)*N equation stemming from the fact that the
    //   fugacity of any component is the same in all phases
    // - 1 equation from the closure condition of all saturations
    //   (they sum up to 1)
    // - M - 1 constraints from the capillary pressures
    //   (-> p_\beta = p_\alpha + p_c\alpha,\beta)
    // - N constraints from the fact that the total mass of each
    //   component is given (-> sum_\alpha rhoMolar_\alpha *
    //   x_\alpha^\kappa = const)
    // - M model constraints. Here we use the NCP constraints
    //   (-> 0 = min{S_\alpha, 1 - \sum_\kappa x_\alpha^\kappa})
    //
    // this also sums up to M*(N + 2).
    //
    // We use the following catches: Capillary pressures are taken
    // into accout expicitly, so that only the pressure of the first
    // phase is solved implicitly, also the closure condition for the
    // saturations is taken into account explicitly, which means, that
    // we don't need to implicitly solve for the last
    // saturation. These two measures reduce the number of unknowns to
    // M*(N + 1), namely:
    //
    // - 1 pressure
    // - M - 1 saturations
    // - M*N mole fractions
    static constexpr int numEq = numPhases*(numComponents + 1);

    typedef Dune::FieldMatrix<Scalar, numEq, numEq> Matrix;
    typedef Dune::FieldVector<Scalar, numEq> Vector;

public:
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    /*!
     * \brief Guess initial values for all quantities.
     */
    template <class FluidState>
    static void guessInitial(FluidState &fluidState,
                             ParameterCache &paramCache,
                             const ComponentVector &globalMolarities)
    {
        // the sum of all molarities
        Scalar sumMoles = 0;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            sumMoles += globalMolarities[compIdx];

        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // composition
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                fluidState.setMoleFraction(phaseIdx,
                                           compIdx,
                                           globalMolarities[compIdx]/sumMoles);

            // pressure. use atmospheric pressure as initial guess
            fluidState.setPressure(phaseIdx, 2e5);

            // saturation. assume all fluids to be equally distributed
            fluidState.setSaturation(phaseIdx, 1.0/numPhases);
        }

        // set the fugacity coefficients of all components in all phases
        paramCache.updateAll(fluidState);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }
    }

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     *
     * The phase's fugacities must already be set.
     */
    template <class MaterialLaw, class FluidState>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      const typename MaterialLaw::Params &matParams,
                      const ComponentVector &globalMolarities)
    {
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-25);

        /////////////////////////
        // Newton method
        /////////////////////////

        // Jacobian matrix
        Matrix J;
        // solution, i.e. phase composition
        Vector deltaX;
        // right hand side
        Vector b;

        Valgrind::SetUndefined(J);
        Valgrind::SetUndefined(deltaX);
        Valgrind::SetUndefined(b);

        // make the fluid state consistent with the fluid system.
        completeFluidState_<MaterialLaw>(fluidState,
                                         paramCache,
                                         matParams);

        /*
        std::cout << "--------------------\n";
        std::cout << "globalMolarities: ";
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
            std::cout << globalMolarities[compIdx] << " ";
        std::cout << "\n";
        */

        const int nMax = 50; // <- maximum number of newton iterations
        for (int nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate Jacobian matrix and right hand side
            linearize_<MaterialLaw>(J, b, fluidState, paramCache, matParams, globalMolarities);
            Valgrind::CheckDefined(J);
            Valgrind::CheckDefined(b);

            // Solve J*x = b
            deltaX = 0;

            try { J.solve(deltaX, b); }
            catch (Dune::FMatrixError e)
            {
                /*
                printFluidState_(fluidState);
                std::cout << "error: " << e << "\n";
                std::cout << "b: " << b << "\n";
                std::cout << "J: " << J << "\n";
                */

                throw Dumux::NumericalProblem(e.what());
            }
            Valgrind::CheckDefined(deltaX);

            /*
            printFluidState_(fluidState);
            std::cout << "J:\n";
            for (int i = 0; i < numEq; ++i) {
                for (int j = 0; j < numEq; ++j) {
                    std::ostringstream os;
                    os << J[i][j];

                    std::string s(os.str());
                    do {
                        s += " ";
                    } while (s.size() < 20);
                    std::cout << s;
                }
                std::cout << "\n";
            };

            std::cout << "deltaX: " << deltaX << "\n";
            std::cout << "---------------\n";
            */

            // update the fluid quantities.
            Scalar relError = update_<MaterialLaw>(fluidState, paramCache, matParams, deltaX);

            if (relError < 1e-9)
                return;
        }

        /*
        printFluidState_(fluidState);
        std::cout << "globalMolarities: ";
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
            std::cout << globalMolarities[compIdx] << " ";
        std::cout << "\n";
        */

        DUNE_THROW(NumericalProblem,
                   "Flash calculation failed."
                   " {c_alpha^kappa} = {" << globalMolarities << "}, T = "
                   << fluidState.temperature(/*phaseIdx=*/0));
    }


protected:
    template <class FluidState>
    static void printFluidState_(const FluidState &fs)
    {
        std::cout << "saturations: ";
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fs.saturation(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "pressures: ";
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fs.pressure(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "densities: ";
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fs.density(phaseIdx) << " ";
        std::cout << "\n";

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            std::cout << "composition " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                std::cout << fs.moleFraction(phaseIdx, compIdx) << " ";
            }
            std::cout << "\n";
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            std::cout << "fugacities " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                std::cout << fs.fugacity(phaseIdx, compIdx) << " ";
            }
            std::cout << "\n";
        }

        std::cout << "global component molarities: ";
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            Scalar sum = 0;
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                sum += fs.saturation(phaseIdx)*fs.molarity(phaseIdx, compIdx);
            }
            std::cout << sum << " ";
        }
        std::cout << "\n";
    }

    template <class MaterialLaw, class FluidState>
    static void linearize_(Matrix &J,
                           Vector &b,
                           FluidState &fluidState,
                           ParameterCache &paramCache,
                           const typename MaterialLaw::Params &matParams,
                           const ComponentVector &globalMolarities)
    {
        FluidState origFluidState(fluidState);
        ParameterCache origParamCache(paramCache);

        Vector tmp;

        // reset jacobian
        J = 0;

        Valgrind::SetUndefined(b);
        calculateDefect_(b, fluidState, fluidState, globalMolarities);
        Valgrind::CheckDefined(b);

        // assemble jacobian matrix
        for (int pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            ////////
            // approximately calculate partial derivatives of the
            // fugacity defect of all components in regard to the mole
            // fraction of the i-th component. This is done via
            // forward differences

            // deviate the mole fraction of the i-th component
            Scalar x_i = getQuantity_(fluidState, pvIdx);
            const Scalar eps = 1e-8/quantityWeight_(fluidState, pvIdx);
            setQuantity_<MaterialLaw>(fluidState, paramCache, matParams, pvIdx, x_i + eps);
            assert(getQuantity_(fluidState, pvIdx) == x_i + eps);

            // compute derivative of the defect
            calculateDefect_(tmp, origFluidState, fluidState, globalMolarities);
            tmp -= b;
            tmp /= eps;

            // store derivative in jacobian matrix
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                J[eqIdx][pvIdx] = tmp[eqIdx];

            // fluid state and parameter cache to their original values
            fluidState = origFluidState;
            paramCache = origParamCache;

            // end forward differences
            ////////
        }
    }

    template <class FluidState>
    static void calculateDefect_(Vector &b,
                                 const FluidState &fluidStateEval,
                                 const FluidState &fluidState,
                                 const ComponentVector &globalMolarities)
    {
        int eqIdx = 0;

        // fugacity of any component must be equal in all phases
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
                b[eqIdx] =
                    fluidState.fugacity(/*phaseIdx=*/0, compIdx) -
                    fluidState.fugacity(phaseIdx, compIdx);
                ++eqIdx;
            }
        }

        assert(eqIdx == numComponents*(numPhases - 1));

        // the fact saturations must sum up to 1 is included explicitly!

        // capillary pressures are explicitly included!

        // global molarities are given
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            b[eqIdx] = 0.0;
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                b[eqIdx] +=
                    fluidState.saturation(phaseIdx)
                    * fluidState.molarity(phaseIdx, compIdx);
            }

            b[eqIdx] -= globalMolarities[compIdx];
            ++eqIdx;
        }

        // model assumptions (-> non-linear complementarity functions)
        // must be adhered
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar sumMoleFracEval = 0.0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                sumMoleFracEval += fluidStateEval.moleFraction(phaseIdx, compIdx);

            if (1.0 - sumMoleFracEval > fluidStateEval.saturation(phaseIdx)) {
                b[eqIdx] = fluidState.saturation(phaseIdx);
            }
            else {
                Scalar sumMoleFrac = 0.0;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    sumMoleFrac += fluidState.moleFraction(phaseIdx, compIdx);
                b[eqIdx] = 1.0 - sumMoleFrac;
            }

            ++eqIdx;
        }
    }

    template <class MaterialLaw, class FluidState>
    static Scalar update_(FluidState &fluidState,
                          ParameterCache &paramCache,
                          const typename MaterialLaw::Params &matParams,
                          const Vector &deltaX)
    {
        Scalar relError = 0;
        for (int pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            Scalar tmp = getQuantity_(fluidState, pvIdx);
            Scalar delta = deltaX[pvIdx];

            relError = std::max(relError, std::abs(delta)*quantityWeight_(fluidState, pvIdx));

            if (isSaturationIdx_(pvIdx)) {
                // dampen to at most 20% change in saturation per
                // iteration
                delta = std::min(0.2, std::max(-0.2, delta));
            }
            else if (isMoleFracIdx_(pvIdx)) {
                // dampen to at most 15% change in mole fraction per
                // iteration
                delta = std::min(0.15, std::max(-0.15, delta));
            }
            else if (isPressureIdx_(pvIdx)) {
                // dampen to at most 15% change in pressure per
                // iteration
                delta = std::min(0.15*fluidState.pressure(0), std::max(-0.15*fluidState.pressure(0), delta));
            };

            setQuantityRaw_(fluidState, pvIdx, tmp - delta);
        }

        // make sure all saturations, pressures and mole fractions are non-negative
        Scalar sumSat = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar value = fluidState.saturation(phaseIdx);
            if (value < -0.05) {
                value = -0.05;
                fluidState.setSaturation(phaseIdx, value);
            }
            sumSat += value;

            value = fluidState.pressure(phaseIdx);
            if (value < 0)
                fluidState.setPressure(phaseIdx, 0.0);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                value = fluidState.moleFraction(phaseIdx, compIdx);
                if (value < 0)
                    fluidState.setMoleFraction(phaseIdx, compIdx, 0.0);
            }
        };

        // last saturation
        if (sumSat > 1.05) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar value = fluidState.saturation(phaseIdx)/(0.95*sumSat);
                fluidState.setSaturation(phaseIdx, value);
            };
        };

        completeFluidState_<MaterialLaw>(fluidState, paramCache, matParams);

        return relError;
    }

    template <class MaterialLaw, class FluidState>
    static void completeFluidState_(FluidState &fluidState,
                                    ParameterCache &paramCache,
                                    const typename MaterialLaw::Params &matParams)
    {
        // calculate the saturation of the last phase as a function of
        // the other saturations
        Scalar sumSat = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            sumSat += fluidState.saturation(phaseIdx);
        fluidState.setSaturation(/*phaseIdx=*/numPhases - 1, 1.0 - sumSat);

        // update the pressures using the material law (saturations
        // and first pressure are already set because it is implicitly
        // solved for.)
        ComponentVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fluidState);
        for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
            fluidState.setPressure(phaseIdx,
                                   fluidState.pressure(0)
                                   + (pC[phaseIdx] - pC[0]));

        // update the parameter cache
        paramCache.updateAll(fluidState, /*except=*/ParameterCache::Temperature);

        // update all densities and fugacity coefficients
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                Scalar phi = FluidSystem::fugacityCoefficient( fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }
    }

    static bool isPressureIdx_(int pvIdx)
    { return pvIdx == 0; }

    static bool isSaturationIdx_(int pvIdx)
    { return 1 <= pvIdx && pvIdx < numPhases; }

    static bool isMoleFracIdx_(int pvIdx)
    { return numPhases <= pvIdx && pvIdx < numPhases + numPhases*numComponents; }

    // retrieves a quantity from the fluid state
    template <class FluidState>
    static Scalar getQuantity_(const FluidState &fs, int pvIdx)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            int phaseIdx = 0;
            return fs.pressure(phaseIdx);
        }
        // first M - 1 saturations
        else if (pvIdx < numPhases) {
            int phaseIdx = pvIdx - 1;
            return fs.saturation(phaseIdx);
        }
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
        {
            int phaseIdx = (pvIdx - numPhases)/numComponents;
            int compIdx = (pvIdx - numPhases)%numComponents;
            return fs.moleFraction(phaseIdx, compIdx);
        }
    }

    // set a quantity in the fluid state
    template <class MaterialLaw, class FluidState>
    static void setQuantity_(FluidState &fs,
                             ParameterCache &paramCache,
                             const typename MaterialLaw::Params &matParams,
                             int pvIdx,
                             Scalar value)
    {
        assert(0 <= pvIdx && pvIdx < numEq);

        if (pvIdx < 1) { // <- first pressure
            Scalar delta = value - fs.pressure(0);

            // set all pressures. here we assume that the capillary
            // pressure does not depend on absolute pressure.
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx, fs.pressure(phaseIdx) + delta);
            paramCache.updateAllPressures(fs);

            // update all densities and fugacity coefficients
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
                fs.setDensity(phaseIdx, rho);

                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    Scalar phi = FluidSystem::fugacityCoefficient(fs, paramCache, phaseIdx, compIdx);
                    fs.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }
        else if (pvIdx < numPhases) { // <- first M - 1 saturations
            Scalar delta = value - fs.saturation(/*phaseIdx=*/pvIdx - 1);
            fs.setSaturation(/*phaseIdx=*/pvIdx - 1, value);

            // set last saturation (-> minus the change of the
            // satuation of the other phase)
            fs.setSaturation(/*phaseIdx=*/numPhases - 1,
                             fs.saturation(numPhases - 1) - delta);

            // update all fluid pressures using the capillary pressure
            // law
            ComponentVector pC;
            MaterialLaw::capillaryPressures(pC, matParams, fs);
            for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx,
                               fs.pressure(0)
                               + (pC[phaseIdx] - pC[0]));
            paramCache.updateAllPressures(fs);

            // update all densities and fugacity coefficients
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
                fs.setDensity(phaseIdx, rho);

                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    Scalar phi = FluidSystem::fugacityCoefficient(fs, paramCache, phaseIdx, compIdx);
                    fs.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }
        else if (pvIdx < numPhases + numPhases*numComponents) // <- mole fractions
        {
            int phaseIdx = (pvIdx - numPhases)/numComponents;
            int compIdx = (pvIdx - numPhases)%numComponents;

            fs.setMoleFraction(phaseIdx, compIdx, value);
            paramCache.updateSingleMoleFraction(fs, phaseIdx, compIdx);

            // update the density of the phase
            Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
            fs.setDensity(phaseIdx, rho);

            // if the phase's fugacity coefficients are composition
            // dependent, update them as well.
            if (!FluidSystem::isIdealMixture(phaseIdx)) {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    Scalar phi = FluidSystem::fugacityCoefficient(fs, paramCache, phaseIdx, compIdx);
                    fs.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }
        else {
            assert(false);
        }
    }

    // set a quantity in the fluid state
    template <class FluidState>
    static void setQuantityRaw_(FluidState &fs, int pvIdx, Scalar value)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            int phaseIdx = 0;
            fs.setPressure(phaseIdx, value);
        }
        // first M - 1 saturations
        else if (pvIdx < numPhases) {
            int phaseIdx = pvIdx - 1;
            fs.setSaturation(phaseIdx, value);
        }
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
        {
            int phaseIdx = (pvIdx - numPhases)/numComponents;
            int compIdx = (pvIdx - numPhases)%numComponents;
            fs.setMoleFraction(phaseIdx, compIdx, value);
        }
    }

    template <class FluidState>
    static Scalar quantityWeight_(const FluidState &fs, int pvIdx)
    {
        // first pressure
        if (pvIdx < 1)
            return 1.0/1e5;
        // first M - 1 saturations
        else if (pvIdx < numPhases)
            return 1.0;
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
            return 1.0;
    }
};

} // end namespace Dumux

#endif
