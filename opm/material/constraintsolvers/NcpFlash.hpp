// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2013 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \copydoc Opm::NcpFlash
 */
#ifndef OPM_NCP_FLASH_HPP
#define OPM_NCP_FLASH_HPP

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/ErrorMacros.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Means.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <limits>
#include <iostream>

namespace Opm {

/*!
 * \brief Determines the phase compositions, pressures and saturations
 *        given the total mass of all components.
 *
 * In a M-phase, N-component context, we have the following
 * unknowns:
 *
 * - M pressures
 * - M saturations
 * - M*N mole fractions
 *
 * This sums up to M*(N + 2). On the equations side of things,
 * we have:
 *
 * - (M - 1)*N equation stemming from the fact that the
 *   fugacity of any component is the same in all phases
 * - 1 equation from the closure condition of all saturations
 *   (they sum up to 1)
 * - M - 1 constraints from the capillary pressures
 *   \f$(-> p_\beta = p_\alpha + p_c\alpha,\beta)\f$
 * - N constraints from the fact that the total mass of each
 *   component is given \f$(-> sum_\alpha rhoMolar_\alpha *
 *   x_\alpha^\kappa = const)\f$
 * - M model constraints. Here we use the NCP constraints
 *   (-> 0 = min \f$ {S_\alpha, 1 - \sum_\kappa x_\alpha^\kappa}\f$)
 *
 * this also sums up to M*(N + 2).
 *
 * We use the following catches: Capillary pressures are taken into
 * account explicitly, so that only the pressure of the first phase is
 * solved implicitly, also the closure condition for the saturations
 * is taken into account explicitly, which means that we don't need to
 * implicitly solve for the last saturation. These two measures reduce
 * the number of unknowns to M*(N + 1), namely:
 *
 * - 1 pressure
 * - M - 1 saturations
 * - M*N mole fractions
 */
template <class Scalar, class FluidSystem>
class NcpFlash
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    typedef typename FluidSystem::ParameterCache ParameterCache;

    static const int numEq = numPhases*(numComponents + 1);

public:
    /*!
     * \brief Guess initial values for all quantities.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static void guessInitial(FluidState &fluidState,
                             ParameterCache &paramCache,
                             const Dune::FieldVector<Evaluation, numComponents>& globalMolarities)
    {
        // the sum of all molarities
        Evaluation sumMoles = 0;
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
            sumMoles += globalMolarities[compIdx];

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // composition
            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
                fluidState.setMoleFraction(phaseIdx,
                                           compIdx,
                                           globalMolarities[compIdx]/sumMoles);

            // pressure. use atmospheric pressure as initial guess
            fluidState.setPressure(phaseIdx, 1.0135e5);

            // saturation. assume all fluids to be equally distributed
            fluidState.setSaturation(phaseIdx, 1.0/numPhases);
        }

        // set the fugacity coefficients of all components in all phases
        paramCache.updateAll(fluidState);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
                const typename FluidState::Scalar phi =
                    FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
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
                      const Dune::FieldVector<typename FluidState::Scalar, numComponents>& globalMolarities,
                      Scalar tolerance = 0.0)
    {
        typedef typename FluidState::Scalar Evaluation;
        typedef Dune::FieldMatrix<Evaluation, numEq, numEq> Matrix;
        typedef Dune::FieldVector<Evaluation, numEq> Vector;

        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-35);

        if (tolerance <= 0.0) {
            tolerance = std::min<Scalar>(1e-10,
                                         Opm::geometricMean(Scalar(1.0),
                                                            std::numeric_limits<Scalar>::epsilon()));
        }

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
        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
            std::cout << globalMolarities[compIdx] << " ";
        std::cout << "\n";
        */
        const int nMax = 50; // <- maximum number of newton iterations
        for (int nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate Jacobian matrix and right hand side
            linearize_<MaterialLaw>(J,
                                    b,
                                    fluidState,
                                    paramCache,
                                    matParams,
                                    globalMolarities);
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

                throw Opm::NumericalIssue(e.what());
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
            }

            std::cout << "deltaX: " << deltaX << "\n";
            std::cout << "---------------\n";
            */

            // update the fluid quantities.
            //update_<MaterialLaw>(fluidState, paramCache, matParams, deltaX);
            Scalar relError = update_<MaterialLaw>(fluidState, paramCache, matParams, deltaX);

            if (relError < 1e-9)
                return;
        }

        /*
        printFluidState_(fluidState);
        std::cout << "globalMolarities: ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
            std::cout << globalMolarities[compIdx] << " ";
        std::cout << "\n";
        */

        OPM_THROW(NumericalIssue,
                  "Flash calculation failed."
                  " {c_alpha^kappa} = {" << globalMolarities << "}, T = "
                  << fluidState.temperature(/*phaseIdx=*/0));
    }

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     *
     * This is a convenience method which assumes that the capillary pressure is
     * zero...
     */
    template <class FluidState, class ComponentVector>
    static void solve(FluidState &fluidState,
                      const ComponentVector &globalMolarities,
                      Scalar tolerance = 0.0)
    {
        ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        typedef NullMaterialTraits<Scalar, numPhases> MaterialTraits;
        typedef NullMaterial<MaterialTraits> MaterialLaw;
        typedef typename MaterialLaw::Params MaterialLawParams;

        MaterialLawParams matParams;
        solve<MaterialLaw>(fluidState, paramCache, matParams, globalMolarities, tolerance);
    }


protected:
    template <class FluidState>
    static void printFluidState_(const FluidState &fluidState)
    {
        std::cout << "saturations: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.saturation(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "pressures: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.pressure(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "densities: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fluidState.density(phaseIdx) << " ";
        std::cout << "\n";

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            std::cout << "composition " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                std::cout << fluidState.moleFraction(phaseIdx, compIdx) << " ";
            }
            std::cout << "\n";
        }

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            std::cout << "fugacities " << FluidSystem::phaseName(phaseIdx) << "Phase: ";
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                std::cout << fluidState.fugacity(phaseIdx, compIdx) << " ";
            }
            std::cout << "\n";
        }

        std::cout << "global component molarities: ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            Scalar sum = 0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                sum += fluidState.saturation(phaseIdx)*fluidState.molarity(phaseIdx, compIdx);
            }
            std::cout << sum << " ";
        }
        std::cout << "\n";
    }

    template <class MaterialLaw,
              class FluidState,
              class Matrix,
              class Vector,
              class ComponentVector>
    static void linearize_(Matrix &J,
                           Vector &b,
                           FluidState &fluidState,
                           ParameterCache &paramCache,
                           const typename MaterialLaw::Params &matParams,
                           const ComponentVector &globalMolarities)
    {
        typedef typename FluidState::Scalar Evaluation;

        FluidState origFluidState(fluidState);
        ParameterCache origParamCache(paramCache);

        Vector tmp;

        // reset jacobian
        J = 0;

        Valgrind::SetUndefined(b);
        calculateDefect_(b, fluidState, fluidState, globalMolarities);
        Valgrind::CheckDefined(b);

        ///////
        // assemble jacobian matrix
        ///////
        for (unsigned pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            ////////
            // approximately calculate partial derivatives of the
            // fugacity defect of all components in regard to the mole
            // fraction of the i-th component. This is done via
            // forward differences

            // deviate the mole fraction of the i-th component
            const Evaluation& x_i = getQuantity_(fluidState, pvIdx);
            const Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e7/(quantityWeight_(fluidState, pvIdx));

            setQuantity_<MaterialLaw>(fluidState, paramCache, matParams, pvIdx, x_i + eps);

            // compute derivative of the defect
            calculateDefect_(tmp, origFluidState, fluidState, globalMolarities);
            tmp -= b;
            for (unsigned i = 0; i < numEq; ++ i)
                tmp[i] /= eps;

            // store derivative in jacobian matrix
            for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx)
                J[eqIdx][pvIdx] = tmp[eqIdx];

            // fluid state and parameter cache to their original values
            fluidState = origFluidState;
            paramCache = origParamCache;

            // end forward differences
            ////////
        }
    }

    template <class FluidState, class Vector, class ComponentVector>
    static void calculateDefect_(Vector &b,
                                 const FluidState &fluidStateEval,
                                 const FluidState &fluidState,
                                 const ComponentVector &globalMolarities)
    {
        typedef typename FluidState::Scalar Evaluation;

        unsigned eqIdx = 0;

        // fugacity of any component must be equal in all phases
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
                b[eqIdx] =
                    fluidState.fugacity(/*phaseIdx=*/0, compIdx)
                    - fluidState.fugacity(phaseIdx, compIdx);
                ++eqIdx;
            }
        }

        assert(eqIdx == numComponents*(numPhases - 1));

        // the fact saturations must sum up to 1 is included implicitly and also,
        // capillary pressures are treated implicitly!

        // global molarities are given
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            b[eqIdx] = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                b[eqIdx] +=
                    fluidState.saturation(phaseIdx)
                    * fluidState.molarity(phaseIdx, compIdx);
            }

            b[eqIdx] -= globalMolarities[compIdx];
            ++eqIdx;
        }

        // model assumptions (-> non-linear complementarity functions) must be adhered
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Evaluation sumMoleFracEval = 0.0;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                sumMoleFracEval += fluidStateEval.moleFraction(phaseIdx, compIdx);

            if (1.0 - sumMoleFracEval > fluidStateEval.saturation(phaseIdx)) {
                b[eqIdx] = fluidState.saturation(phaseIdx);
            }
            else {
                Evaluation sumMoleFrac = 0.0;
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                    sumMoleFrac += fluidState.moleFraction(phaseIdx, compIdx);
                b[eqIdx] = 1.0 - sumMoleFrac;
            }

            ++eqIdx;
        }
    }

    template <class MaterialLaw, class FluidState, class Vector>
    static Scalar update_(FluidState &fluidState,
                          ParameterCache &paramCache,
                          const typename MaterialLaw::Params &matParams,
                          const Vector &deltaX)
    {
        typedef typename FluidState::Scalar Evaluation;
        typedef Opm::MathToolbox<Evaluation> Toolbox;

        // make sure we don't swallow non-finite update vectors
#ifndef NDEBUG
        assert(deltaX.dimension == numEq);
        for (unsigned i = 0; i < numEq; ++i)
            assert(std::isfinite(Toolbox::value(deltaX[i])));
#endif

        Scalar relError = 0;
        for (unsigned pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            const Evaluation& tmp = getQuantity_(fluidState, pvIdx);
            Evaluation delta = deltaX[pvIdx];

            relError = std::max<Scalar>(relError,
                                        std::abs(Toolbox::value(delta))
                                        * quantityWeight_(fluidState, pvIdx));

            if (isSaturationIdx_(pvIdx)) {
                // dampen to at most 25% change in saturation per iteration
                delta = Toolbox::min(0.25, Toolbox::max(-0.25, delta));
            }
            else if (isMoleFracIdx_(pvIdx)) {
                // dampen to at most 20% change in mole fraction per iteration
                delta = Toolbox::min(0.20, Toolbox::max(-0.20, delta));
            }
            else if (isPressureIdx_(pvIdx)) {
                // dampen to at most 50% change in pressure per iteration
                delta = Toolbox::min(0.5*fluidState.pressure(0),
                                     Toolbox::max(-0.5*fluidState.pressure(0), delta));
            }

            setQuantityRaw_(fluidState, pvIdx, tmp - delta);
        }

        /*
        // make sure all saturations, pressures and mole fractions are non-negative
        Scalar sumSat = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar value = fluidState.saturation(phaseIdx);
            if (value < -0.05) {
                value = -0.05;
                fluidState.setSaturation(phaseIdx, value);
            }
            sumSat += value;

            value = fluidState.pressure(phaseIdx);
            if (value < 0)
                fluidState.setPressure(phaseIdx, 0.0);

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                value = fluidState.moleFraction(phaseIdx, compIdx);
                if (value < 0)
                    fluidState.setMoleFraction(phaseIdx, compIdx, 0.0);
            }
        }

        // last saturation
        if (sumSat > 1.05) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar value = fluidState.saturation(phaseIdx)/(0.95*sumSat);
                fluidState.setSaturation(phaseIdx, value);
            }
        }
        */

        completeFluidState_<MaterialLaw>(fluidState, paramCache, matParams);

        return relError;
    }

    template <class MaterialLaw, class FluidState>
    static void completeFluidState_(FluidState &fluidState,
                                    ParameterCache &paramCache,
                                    const typename MaterialLaw::Params &matParams)
    {
        typedef typename FluidState::Scalar Evaluation;

        // calculate the saturation of the last phase as a function of
        // the other saturations
        Evaluation sumSat = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            sumSat += fluidState.saturation(phaseIdx);
        fluidState.setSaturation(/*phaseIdx=*/numPhases - 1, 1.0 - sumSat);

        // update the pressures using the material law (saturations
        // and first pressure are already set because it is implicitly
        // solved for.)
        Dune::FieldVector<Scalar, numPhases> pC;
        MaterialLaw::capillaryPressures(pC, matParams, fluidState);
        for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
            fluidState.setPressure(phaseIdx,
                                   fluidState.pressure(0)
                                   + (pC[phaseIdx] - pC[0]));

        // update the parameter cache
        paramCache.updateAll(fluidState, /*except=*/ParameterCache::Temperature);

        // update all densities and fugacity coefficients
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            const Evaluation& rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
                const Evaluation& phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }
    }

    static bool isPressureIdx_(unsigned pvIdx)
    { return pvIdx == 0; }

    static bool isSaturationIdx_(unsigned pvIdx)
    { return 1 <= pvIdx && pvIdx < numPhases; }

    static bool isMoleFracIdx_(unsigned pvIdx)
    { return numPhases <= pvIdx && pvIdx < numPhases + numPhases*numComponents; }

    // retrieves a quantity from the fluid state
    template <class FluidState>
    static const typename FluidState::Scalar& getQuantity_(const FluidState &fluidState, unsigned pvIdx)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            unsigned phaseIdx = 0;
            return fluidState.pressure(phaseIdx);
        }
        // first M - 1 saturations
        else if (pvIdx < numPhases) {
            unsigned phaseIdx = pvIdx - 1;
            return fluidState.saturation(phaseIdx);
        }
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
        {
            unsigned phaseIdx = (pvIdx - numPhases)/numComponents;
            unsigned compIdx = (pvIdx - numPhases)%numComponents;
            return fluidState.moleFraction(phaseIdx, compIdx);
        }
    }

    // set a quantity in the fluid state
    template <class MaterialLaw, class FluidState>
    static void setQuantity_(FluidState &fluidState,
                             ParameterCache &paramCache,
                             const typename MaterialLaw::Params &matParams,
                             unsigned pvIdx,
                             const typename FluidState::Scalar& value)
    {
        typedef typename FluidState::Scalar Evaluation;

        assert(0 <= pvIdx && pvIdx < numEq);

        if (pvIdx < 1) { // <- first pressure
            Evaluation delta = value - fluidState.pressure(0);

            // set all pressures. here we assume that the capillary
            // pressure does not depend on absolute pressure.
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState.setPressure(phaseIdx, fluidState.pressure(phaseIdx) + delta);
            paramCache.updateAllPressures(fluidState);

            // update all densities and fugacity coefficients
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                const Evaluation& rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                Valgrind::CheckDefined(rho);
                fluidState.setDensity(phaseIdx, rho);

                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Evaluation& phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                    Valgrind::CheckDefined(phi);
                    fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }
        else if (pvIdx < numPhases) { // <- first M - 1 saturations
            const Evaluation& delta = value - fluidState.saturation(/*phaseIdx=*/pvIdx - 1);
            Valgrind::CheckDefined(delta);
            fluidState.setSaturation(/*phaseIdx=*/pvIdx - 1, value);

            // set last saturation (-> minus the change of the saturation of the other
            // phase)
            fluidState.setSaturation(/*phaseIdx=*/numPhases - 1,
                                     fluidState.saturation(numPhases - 1) - delta);

            // update all fluid pressures using the capillary pressure law
            Dune::FieldVector<Evaluation, numPhases> pC;
            MaterialLaw::capillaryPressures(pC, matParams, fluidState);
            Valgrind::CheckDefined(pC);
            for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
                fluidState.setPressure(phaseIdx,
                                       fluidState.pressure(0)
                                       + (pC[phaseIdx] - pC[0]));
            paramCache.updateAllPressures(fluidState);

            // update all densities and fugacity coefficients
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                const Evaluation& rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                Valgrind::CheckDefined(rho);
                fluidState.setDensity(phaseIdx, rho);

                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    const Evaluation& phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                    Valgrind::CheckDefined(phi);
                    fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
                }
            }
        }
        else if (pvIdx < numPhases + numPhases*numComponents) // <- mole fractions
        {
            unsigned phaseIdx = (pvIdx - numPhases)/numComponents;
            unsigned compIdx = (pvIdx - numPhases)%numComponents;

            Valgrind::CheckDefined(value);
            fluidState.setMoleFraction(phaseIdx, compIdx, value);
            paramCache.updateSingleMoleFraction(fluidState, phaseIdx, compIdx);

            // update the density of the phase
            const Evaluation& rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            Valgrind::CheckDefined(rho);
            fluidState.setDensity(phaseIdx, rho);

            // if the phase's fugacity coefficients are composition
            // dependent, update them as well.
            if (!FluidSystem::isIdealMixture(phaseIdx)) {
                for (unsigned fugCompIdx = 0; fugCompIdx < numComponents; ++fugCompIdx) {
                    const Evaluation& phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, fugCompIdx);
                    Valgrind::CheckDefined(phi);
                    fluidState.setFugacityCoefficient(phaseIdx, fugCompIdx, phi);
                }
            }
        }
        else {
            assert(false);
        }
    }

    // set a quantity in the fluid state
    template <class FluidState>
    static void setQuantityRaw_(FluidState &fluidState,
                                unsigned pvIdx,
                                const typename FluidState::Scalar& value)
    {
        assert(pvIdx < numEq);

        Valgrind::CheckDefined(value);
        // first pressure
        if (pvIdx < 1) {
            unsigned phaseIdx = 0;
            fluidState.setPressure(phaseIdx, value);
        }
        // first M - 1 saturations
        else if (pvIdx < numPhases) {
            unsigned phaseIdx = pvIdx - 1;
            fluidState.setSaturation(phaseIdx, value);
        }
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
        {
            unsigned phaseIdx = (pvIdx - numPhases)/numComponents;
            unsigned compIdx = (pvIdx - numPhases)%numComponents;
            fluidState.setMoleFraction(phaseIdx, compIdx, value);
        }
    }

    template <class FluidState>
    static Scalar quantityWeight_(const FluidState &/*fluidState*/, unsigned pvIdx)
    {
        // first pressure
        if (pvIdx < 1)
            return 1/1e5;
        // first M - 1 saturations
        else if (pvIdx < numPhases)
            return 1.0;
        // mole fractions
        else // if (pvIdx < numPhases + numPhases*numComponents)
            return 1.0;
    }
};

} // namespace Opm

#endif
