// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::ImmiscibleFlash
 */
#ifndef OPM_IMMISCIBLE_FLASH_HPP
#define OPM_IMMISCIBLE_FLASH_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <limits>
#include <iostream>

namespace Opm {

/*!
 * \brief Determines the pressures and saturations of all fluid phases
 *        given the total mass of all components.
 *
 * In a N-phase, N-component context, we have the following
 * unknowns if assuming immiscibility:
 *
 * - N pressures
 * - N saturations
 *
 * This sums up to 2*N unknowns. On the equations side of things, we
 * have:
 *
 * - N total component molarities
 * - 1 The sum of all saturations is 1
 * - N-1 Relations from capillary pressure
 *
 * this also sums up to 2*N. We include the capillary pressures and
 * the sum of the saturations explicitly. This means that we only
 * solve for the first pressure and N-1 saturations.
 *
 * If a fluid phase is incompressible, the pressure cannot determined
 * by this, though. In this case the original pressure is kept, and
 * the saturation of the phase is calculated by dividing the global
 * molarity of the component by the phase density.
 */
template <class Scalar, class FluidSystem>
class ImmiscibleFlash
{
    static const int numPhases = FluidSystem::numPhases;
    static const int numComponents = FluidSystem::numComponents;
    static_assert(numPhases == numComponents,
                  "Immiscibility assumes that the number of phases"
                  " is equal to the number of components");


    static const int numEq = numPhases;

    typedef Dune::FieldMatrix<Scalar, numEq, numEq> Matrix;
    typedef Dune::FieldVector<Scalar, numEq> Vector;

public:
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    /*!
     * \brief Guess initial values for all quantities.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static void guessInitial(FluidState &fluidState,
                             const Dune::FieldVector<Evaluation, numComponents>& /*globalMolarities*/)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // pressure. use atmospheric pressure as initial guess
            fluidState.setPressure(phaseIdx, 2e5);

            // saturation. assume all fluids to be equally distributed
            fluidState.setSaturation(phaseIdx, 1.0/numPhases);
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
                      const typename MaterialLaw::Params &matParams,
                      typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                      const Dune::FieldVector<typename FluidState::Scalar, numComponents>& globalMolarities)
    {
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-25);

        paramCache.updateAll(fluidState);

        /////////////////////////
        // Check if all fluid phases are incompressible
        /////////////////////////
        bool allIncompressible = true;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::isCompressible(phaseIdx)) {
                allIncompressible = false;
                break;
            }
        };

        if (allIncompressible) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);

                Scalar saturation =
                    globalMolarities[/*compIdx=*/phaseIdx]
                    / fluidState.molarDensity(phaseIdx);
                fluidState.setSaturation(phaseIdx, saturation);
            }
        };

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

        completeFluidState_<MaterialLaw>(fluidState, matParams, paramCache);

        const int nMax = 50; // <- maximum number of newton iterations
        for (int nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate Jacobian matrix and right hand side
            linearize_<MaterialLaw>(J, b, fluidState, matParams, paramCache, globalMolarities);
            Valgrind::CheckDefined(J);
            Valgrind::CheckDefined(b);

            // Solve J*x = b
            deltaX = 0;

            try { J.solve(deltaX, b); }
            catch (Dune::FMatrixError e)
            {
                /*
                std::cout << "error: " << e << "\n";
                std::cout << "b: " << b << "\n";
                */

                throw Opm::NumericalProblem(e.what());
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
            std::cout << "residual: " << b << "\n";
            std::cout << "deltaX: " << deltaX << "\n";
            std::cout << "---------------\n";
            */

            // update the fluid quantities.
            Scalar relError = update_<MaterialLaw>(fluidState, matParams, paramCache, deltaX);

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

        OPM_THROW(Opm::NumericalProblem,
                  "Flash calculation failed."
                  " {c_alpha^kappa} = {" << globalMolarities << "}, T = "
                  << fluidState.temperature(/*phaseIdx=*/0));
    }


protected:
    template <class FluidState>
    static void printFluidState_(const FluidState &fs)
    {
        std::cout << "saturations: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fs.saturation(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "pressures: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fs.pressure(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "densities: ";
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            std::cout << fs.density(phaseIdx) << " ";
        std::cout << "\n";

        std::cout << "global component molarities: ";
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            Scalar sum = 0;
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
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
                           const typename MaterialLaw::Params &matParams,
                           typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                           const ComponentVector &globalMolarities)
    {
        // copy the undisturbed fluid state and parameter cache
        FluidState origFluidState(fluidState);
        auto origParamCache(paramCache);

        Vector tmp;

        // reset jacobian
        J = 0;

        Valgrind::SetUndefined(b);
        calculateDefect_(b, fluidState, fluidState, globalMolarities);
        Valgrind::CheckDefined(b);

        // assemble jacobian matrix
        for (unsigned pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            ////////
            // approximately calculate partial derivatives of the
            // fugacity defect of all components in regard to the mole
            // fraction of the i-th component. This is done via
            // forward differences

            // deviate the mole fraction of the i-th component
            Scalar xI = getQuantity_(fluidState, pvIdx);
            const Scalar eps = 1e-10/quantityWeight_(fluidState, pvIdx);
            setQuantity_<MaterialLaw>(fluidState, matParams, paramCache, pvIdx, xI + eps);
            assert(std::abs(getQuantity_(fluidState, pvIdx) - (xI + eps))
                   <= std::max<Scalar>(1.0, std::abs(xI))*std::numeric_limits<Scalar>::epsilon()*100);

            // compute derivative of the defect
            calculateDefect_(tmp, origFluidState, fluidState, globalMolarities);
            tmp -= b;
            tmp /= eps;
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

    template <class FluidState>
    static void calculateDefect_(Vector &b,
                                 const FluidState &/*fluidStateEval*/,
                                 const FluidState &fluidState,
                                 const ComponentVector &globalMolarities)
    {
        // global molarities are given
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            b[phaseIdx] =
                fluidState.saturation(phaseIdx)
                * fluidState.molarity(phaseIdx, /*compIdx=*/phaseIdx);

            b[phaseIdx] -= globalMolarities[/*compIdx=*/phaseIdx];
        }
    }

    template <class MaterialLaw, class FlashFluidState, class EvalVector>
    static Scalar update_(FlashFluidState& fluidState,
                          const typename MaterialLaw::Params& matParams,
                          typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>& paramCache,
                          const EvalVector& deltaX)
    {
        Scalar relError = 0;
        for (unsigned pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            Scalar tmp = getQuantity_(fluidState, pvIdx);
            Scalar delta = deltaX[pvIdx];

            relError = std::max(relError, std::abs(delta)*quantityWeight_(fluidState, pvIdx));

            if (isSaturationIdx_(pvIdx)) {
                // dampen to at most 20% change in saturation per
                // iteration
                delta = std::min<Scalar>(0.2, std::max<Scalar>(-0.2, delta));
            }
            else if (isPressureIdx_(pvIdx)) {
                // dampen to at most 30% change in pressure per
                // iteration
                delta = std::min<Scalar>(0.30*fluidState.pressure(0),
                                         std::max<Scalar>(-0.30*fluidState.pressure(0), delta));
            };

            setQuantityRaw_(fluidState, pvIdx, tmp - delta);
        }

        completeFluidState_<MaterialLaw>(fluidState, matParams, paramCache);

        return relError;
    }

    template <class MaterialLaw, class FluidState>
    static void completeFluidState_(FluidState &fluidState,
                                    const typename MaterialLaw::Params &matParams,
                                    typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache)
    {
        typedef typename FluidSystem::template ParameterCache<typename FluidState::Scalar> ParamCache;

        // calculate the saturation of the last phase as a function of
        // the other saturations
        Scalar sumSat = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            sumSat += fluidState.saturation(phaseIdx);

        if (sumSat > 1.0) {
            // make sure that the last saturation does not become
            // negative
            for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            {
                Scalar S = fluidState.saturation(phaseIdx);
                fluidState.setSaturation(phaseIdx, S/sumSat);
            }
            sumSat = 1;
        }

        // set the last saturation
        fluidState.setSaturation(/*phaseIdx=*/numPhases - 1, 1.0 - sumSat);

        // update the pressures using the material law (saturations
        // and first pressure are already set because it is implicitly
        // solved for.)
        ComponentVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fluidState);
        for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
            fluidState.setPressure(phaseIdx,
                                   fluidState.pressure(0)
                                   + (pC[phaseIdx] - pC[0]));

        // update the parameter cache
        paramCache.updateAll(fluidState, /*except=*/ParamCache::Temperature|ParamCache::Composition);

        // update all densities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);
        }
    }

    static bool isPressureIdx_(unsigned pvIdx)
    { return pvIdx == 0; }

    static bool isSaturationIdx_(unsigned pvIdx)
    { return 1 <= pvIdx && pvIdx < numPhases; }

    // retrieves a quantity from the fluid state
    template <class FluidState>
    static Scalar getQuantity_(const FluidState &fs, unsigned pvIdx)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            unsigned phaseIdx = 0;
            return fs.pressure(phaseIdx);
        }
        // saturations
        else {
            assert(pvIdx < numPhases);
            unsigned phaseIdx = pvIdx - 1;
            return fs.saturation(phaseIdx);
        }
    }

    // set a quantity in the fluid state
    template <class MaterialLaw, class FluidState>
    static void setQuantity_(FluidState &fs,
                             const typename MaterialLaw::Params &matParams,
                             typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                             unsigned pvIdx,
                             Scalar value)
    {
        assert(0 <= pvIdx && pvIdx < numEq);

        if (pvIdx < 1) {
            // -> first pressure
            Scalar delta = value - fs.pressure(0);

            // set all pressures. here we assume that the capillary
            // pressure does not depend on absolute pressure.
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx, fs.pressure(phaseIdx) + delta);

            // update the parameter cache
            paramCache.updateAllPressures(fs);

            // update all densities
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
                fs.setDensity(phaseIdx, rho);
            }
        }
        else {
            assert(pvIdx < numPhases);

            // -> first M-1 saturations
            Scalar delta = value - fs.saturation(/*phaseIdx=*/pvIdx - 1);
            fs.setSaturation(/*phaseIdx=*/pvIdx - 1, value);

            // set last saturation (-> minus the change of the
            // saturation of the other phases)
            fs.setSaturation(/*phaseIdx=*/numPhases - 1,
                             fs.saturation(numPhases - 1) - delta);

            // update all fluid pressures using the capillary pressure
            // law
            ComponentVector pC;
            MaterialLaw::capillaryPressures(pC, matParams, fs);
            for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx,
                               fs.pressure(0)
                               + (pC[phaseIdx] - pC[0]));

            // update the parameter cache
            paramCache.updateAllPressures(fs);

            // update all densities
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
                fs.setDensity(phaseIdx, rho);
            }
        }
    }

    // set a quantity in the fluid state
    template <class FluidState>
    static void setQuantityRaw_(FluidState &fs, unsigned pvIdx, Scalar value)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            unsigned phaseIdx = 0;
            fs.setPressure(phaseIdx, value);
        }
        // saturations
        else {
            assert(pvIdx < numPhases);
            unsigned phaseIdx = pvIdx - 1;

            // make sure that the first M-1 saturations does not get
            // negative
            value = std::max<Scalar>(0.0, value);
            fs.setSaturation(phaseIdx, value);
        }
    }

    template <class FluidState>
    static Scalar quantityWeight_(const FluidState &/*fs*/, unsigned pvIdx)
    {
        // first pressure
        if (pvIdx < 1)
            return 1.0/1e5;
        // first M - 1 saturations
        else {
            assert(pvIdx < numPhases);
            return 1.0;
        }
    }
};

} // namespace Opm

#endif
