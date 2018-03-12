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

#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

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

    enum {
        p0PvIdx = 0,
        S0PvIdx = 1
    };

    static const int numEq = numPhases;

public:
    /*!
     * \brief Guess initial values for all quantities.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static void guessInitial(FluidState& fluidState,
                             const Dune::FieldVector<Evaluation, numComponents>& /*globalMolarities*/)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            // pressure. use 1 bar as initial guess
            fluidState.setPressure(phaseIdx, 1e5);

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
    static void solve(FluidState& fluidState,
                      const typename MaterialLaw::Params& matParams,
                      typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                      const Dune::FieldVector<typename FluidState::Scalar, numComponents>& globalMolarities,
                      Scalar tolerance = -1)
    {
        typedef typename FluidState::Scalar InputEval;

        /////////////////////////
        // Check if all fluid phases are incompressible
        /////////////////////////
        bool allIncompressible = true;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::isCompressible(phaseIdx)) {
                allIncompressible = false;
                break;
            }
        }

        if (allIncompressible) {
            // yes, all fluid phases are incompressible. in this case the flash solver
            // can only determine the saturations, not the pressures. (but this
            // determination is much simpler than a full flash calculation.)
            paramCache.updateAll(fluidState);
            solveAllIncompressible_(fluidState, paramCache, globalMolarities);
            return;
        }

        typedef Dune::FieldMatrix<InputEval, numEq, numEq> Matrix;
        typedef Dune::FieldVector<InputEval, numEq> Vector;

        typedef Opm::DenseAd::Evaluation<InputEval, numEq> FlashEval;
        typedef Dune::FieldVector<FlashEval, numEq> FlashDefectVector;
        typedef Opm::ImmiscibleFluidState<FlashEval, FluidSystem> FlashFluidState;

#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2,7)
        Dune::FMatrixPrecision<InputEval>::set_singular_limit(1e-35);
#endif

        if (tolerance <= 0)
            tolerance = std::min<Scalar>(1e-5,
                                         1e8*std::numeric_limits<Scalar>::epsilon());

        typename FluidSystem::template ParameterCache<FlashEval> flashParamCache;
        flashParamCache.assignPersistentData(paramCache);

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

        FlashFluidState flashFluidState;
        assignFlashFluidState_<MaterialLaw>(fluidState, flashFluidState, matParams, flashParamCache);

        // copy the global molarities to a vector of evaluations. Remember that the
        // global molarities are constants. (but we need to copy them to a vector of
        // FlashEvals anyway in order to avoid getting into hell's kitchen.)
        Dune::FieldVector<FlashEval, numComponents> flashGlobalMolarities;
        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
            flashGlobalMolarities[compIdx] = globalMolarities[compIdx];

        FlashDefectVector defect;
        const unsigned nMax = 50; // <- maximum number of newton iterations
        for (unsigned nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate Jacobian matrix and right hand side
            evalDefect_(defect, flashFluidState, flashGlobalMolarities);
            Valgrind::CheckDefined(defect);

            // create field matrices and vectors out of the evaluation vector to solve
            // the linear system of equations.
            for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
                for (unsigned pvIdx = 0; pvIdx < numEq; ++ pvIdx)
                    J[eqIdx][pvIdx] = defect[eqIdx].derivative(pvIdx);

                b[eqIdx] = defect[eqIdx].value();
            }
            Valgrind::CheckDefined(J);
            Valgrind::CheckDefined(b);

            // Solve J*x = b
            deltaX = 0;

            try { J.solve(deltaX, b); }
            catch (Dune::FMatrixError e) {
                throw Opm::NumericalIssue(e.what());
            }
            Valgrind::CheckDefined(deltaX);

            // update the fluid quantities.
            Scalar relError = update_<MaterialLaw>(flashFluidState, flashParamCache, matParams, deltaX);

            if (relError < tolerance) {
                assignOutputFluidState_(flashFluidState, fluidState);
                return;
            }
        }

        std::ostringstream oss;
        oss << "ImmiscibleFlash solver failed:"
            << " {c_alpha^kappa} = {" << globalMolarities << "},"
            << " T = " << fluidState.temperature(/*phaseIdx=*/0);
    throw Opm::NumericalIssue(oss.str());
    }


protected:
    template <class FluidState>
    static void printFluidState_(const FluidState& fs)
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

    template <class MaterialLaw, class InputFluidState, class FlashFluidState>
    static void assignFlashFluidState_(const InputFluidState& inputFluidState,
                                       FlashFluidState& flashFluidState,
                                       const typename MaterialLaw::Params& matParams,
                                       typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>& flashParamCache)
    {
        typedef typename FlashFluidState::Scalar FlashEval;

        // copy the temperature: even though the model which uses the flash solver might
        // be non-isothermal, the flash solver does not consider energy. (it could be
        // modified to do so relatively easily, but it would come at increased
        // computational cost and normally temperature instead of "total internal energy
        // of the fluids" is specified.)
        flashFluidState.setTemperature(inputFluidState.temperature(/*phaseIdx=*/0));

        // copy the saturations: the first N-1 phases are primary variables, the last one
        // is one minus the sum of the former.
        FlashEval Slast = 1.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx) {
            FlashEval S = inputFluidState.saturation(phaseIdx);
            S.setDerivative(S0PvIdx + phaseIdx, 1.0);

            Slast -= S;

            flashFluidState.setSaturation(phaseIdx, S);
        }
        flashFluidState.setSaturation(numPhases - 1, Slast);

        // copy the pressures: the first pressure is the first primary variable, the
        // remaining ones are given as p_beta = p_alpha + p_calpha,beta
        FlashEval p0 = inputFluidState.pressure(0);
        p0.setDerivative(p0PvIdx, 1.0);

        std::array<FlashEval, numPhases> pc;
        MaterialLaw::capillaryPressures(pc, matParams, flashFluidState);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            flashFluidState.setPressure(phaseIdx, p0 + (pc[phaseIdx] - pc[0]));

        flashParamCache.updateAll(flashFluidState);

        // compute the density of each phase.
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const FlashEval& rho = FluidSystem::density(flashFluidState, flashParamCache, phaseIdx);
            flashFluidState.setDensity(phaseIdx, rho);
        }
    }

    template <class FlashFluidState, class OutputFluidState>
    static void assignOutputFluidState_(const FlashFluidState& flashFluidState,
                                        OutputFluidState& outputFluidState)
    {
        outputFluidState.setTemperature(flashFluidState.temperature(/*phaseIdx=*/0).value());

        // copy the saturations, pressures and densities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const auto& S = flashFluidState.saturation(phaseIdx).value();
            outputFluidState.setSaturation(phaseIdx, S);

            const auto& p = flashFluidState.pressure(phaseIdx).value();
            outputFluidState.setPressure(phaseIdx, p);

            const auto& rho = flashFluidState.density(phaseIdx).value();
            outputFluidState.setDensity(phaseIdx, rho);
        }
    }

    template <class FluidState>
    static void solveAllIncompressible_(FluidState& fluidState,
                                        typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                                        const Dune::FieldVector<typename FluidState::Scalar, numComponents>& globalMolarities)
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            Scalar saturation =
                globalMolarities[/*compIdx=*/phaseIdx]
                / fluidState.molarDensity(phaseIdx);
            fluidState.setSaturation(phaseIdx, saturation);
        }
    }

    template <class FluidState, class FlashDefectVector, class FlashComponentVector>
    static void evalDefect_(FlashDefectVector& b,
                            const FluidState& fluidState,
                            const FlashComponentVector& globalMolarities)
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
                          typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>& paramCache,
                          const typename MaterialLaw::Params& matParams,
                          const EvalVector& deltaX)
    {
        // note that it is possible that FlashEval::Scalar is an Evaluation itself
        typedef typename FlashFluidState::Scalar FlashEval;
        typedef Opm::MathToolbox<FlashEval> FlashEvalToolbox;

        typedef typename FlashEvalToolbox::ValueType InnerEval;

#ifndef NDEBUG
        // make sure we don't swallow non-finite update vectors
        assert(deltaX.dimension == numEq);
        for (unsigned i = 0; i < numEq; ++i)
            assert(std::isfinite(Opm::scalarValue(deltaX[i])));
#endif

        Scalar relError = 0;
        for (unsigned pvIdx = 0; pvIdx < numEq; ++ pvIdx) {
            FlashEval tmp = getQuantity_(fluidState, pvIdx);
            InnerEval delta = deltaX[pvIdx];

            relError = std::max(relError,
                                std::abs(Opm::scalarValue(delta))
                                * quantityWeight_(fluidState, pvIdx));

            if (isSaturationIdx_(pvIdx)) {
                // dampen to at most 20% change in saturation per
                // iteration
                delta = Opm::min(0.25, Opm::max(-0.25, delta));
            }
            else if (isPressureIdx_(pvIdx)) {
                // dampen to at most 30% change in pressure per
                // iteration
                delta = Opm::min(0.5*fluidState.pressure(0).value(),
                                 Opm::max(-0.50*fluidState.pressure(0).value(), delta));
            }

            tmp -= delta;
            setQuantity_(fluidState, pvIdx, tmp);
        }

        completeFluidState_<MaterialLaw>(fluidState, paramCache, matParams);

        return relError;
    }

    template <class MaterialLaw, class FlashFluidState>
    static void completeFluidState_(FlashFluidState& flashFluidState,
                                    typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar>& paramCache,
                                    const typename MaterialLaw::Params& matParams)
    {
        typedef typename FluidSystem::template ParameterCache<typename FlashFluidState::Scalar> ParamCache;

        typedef typename FlashFluidState::Scalar FlashEval;

        // calculate the saturation of the last phase as a function of
        // the other saturations
        FlashEval sumSat = 0.0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            sumSat += flashFluidState.saturation(phaseIdx);
        flashFluidState.setSaturation(/*phaseIdx=*/numPhases - 1, 1.0 - sumSat);

        // update the pressures using the material law (saturations
        // and first pressure are already set because it is implicitly
        // solved for.)
        Dune::FieldVector<FlashEval, numPhases> pC;
        MaterialLaw::capillaryPressures(pC, matParams, flashFluidState);
        for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
            flashFluidState.setPressure(phaseIdx,
                                        flashFluidState.pressure(0)
                                        + (pC[phaseIdx] - pC[0]));

        // update the parameter cache
        paramCache.updateAll(flashFluidState, /*except=*/ParamCache::Temperature|ParamCache::Composition);

        // update all densities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            const FlashEval& rho = FluidSystem::density(flashFluidState, paramCache, phaseIdx);
            flashFluidState.setDensity(phaseIdx, rho);
        }
    }

    static bool isPressureIdx_(unsigned pvIdx)
    { return pvIdx == 0; }

    static bool isSaturationIdx_(unsigned pvIdx)
    { return 1 <= pvIdx && pvIdx < numPhases; }

    // retrieves a quantity from the fluid state
    template <class FluidState>
    static const typename FluidState::Scalar& getQuantity_(const FluidState& fluidState, unsigned pvIdx)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            unsigned phaseIdx = 0;
            return fluidState.pressure(phaseIdx);
        }
        // saturations
        else {
            assert(pvIdx < numPhases);
            unsigned phaseIdx = pvIdx - 1;
            return fluidState.saturation(phaseIdx);
        }
    }

    // set a quantity in the fluid state
    template <class FluidState>
    static void setQuantity_(FluidState& fluidState,
                             unsigned pvIdx,
                             const typename FluidState::Scalar& value)
    {
        assert(pvIdx < numEq);

        // first pressure
        if (pvIdx < 1) {
            unsigned phaseIdx = 0;
            fluidState.setPressure(phaseIdx, value);
        }
        // saturations
        else {
            assert(pvIdx < numPhases);
            unsigned phaseIdx = pvIdx - 1;
            fluidState.setSaturation(phaseIdx, value);
        }
    }

    template <class FluidState>
    static Scalar quantityWeight_(const FluidState& /*fluidState*/, unsigned pvIdx)
    {
        // first pressure
        if (pvIdx < 1)
            return 1e-8;
        // first M - 1 saturations
        else {
            assert(pvIdx < numPhases);
            return 1.0;
        }
    }
};

} // namespace Opm

#endif
