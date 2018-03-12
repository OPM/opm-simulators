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
 * \copydoc Opm::CompositionFromFugacities
 */
#ifndef OPM_COMPOSITION_FROM_FUGACITIES_HPP
#define OPM_COMPOSITION_FROM_FUGACITIES_HPP

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <limits>

namespace Opm {

/*!
 * \brief Calculates the chemical equilibrium from the component
 *        fugacities in a phase.
 */
template <class Scalar, class FluidSystem, class Evaluation = Scalar>
class CompositionFromFugacities
{
    enum { numComponents = FluidSystem::numComponents };

public:
    typedef Dune::FieldVector<Evaluation, numComponents> ComponentVector;

    /*!
     * \brief Guess an initial value for the composition of the phase.
     */
    template <class FluidState>
    static void guessInitial(FluidState& fluidState,
                             unsigned phaseIdx,
                             const ComponentVector& /*fugVec*/)
    {
        if (FluidSystem::isIdealMixture(phaseIdx))
            return;

        // Pure component fugacities
        for (unsigned i = 0; i < numComponents; ++ i) {
            //std::cout << f << " -> " << mutParams.fugacity(phaseIdx, i)/f << "\n";
            fluidState.setMoleFraction(phaseIdx,
                                       i,
                                       1.0/numComponents);
        }
    }

    /*!
     * \brief Calculates the chemical equilibrium from the component
     *        fugacities in a phase.
     *
     * The phase's fugacities must already be set.
     */
    template <class FluidState>
    static void solve(FluidState& fluidState,
                      typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                      unsigned phaseIdx,
                      const ComponentVector& targetFug)
    {
        // use a much more efficient method in case the phase is an
        // ideal mixture
        if (FluidSystem::isIdealMixture(phaseIdx)) {
            solveIdealMix_(fluidState, paramCache, phaseIdx, targetFug);
            return;
        }

        // save initial composition in case something goes wrong
        Dune::FieldVector<Evaluation, numComponents> xInit;
        for (unsigned i = 0; i < numComponents; ++i) {
            xInit[i] = fluidState.moleFraction(phaseIdx, i);
        }

        /////////////////////////
        // Newton method
        /////////////////////////

        // Jacobian matrix
        Dune::FieldMatrix<Evaluation, numComponents, numComponents> J;
        // solution, i.e. phase composition
        Dune::FieldVector<Evaluation, numComponents> x;
        // right hand side
        Dune::FieldVector<Evaluation, numComponents> b;

        paramCache.updatePhase(fluidState, phaseIdx);

        // maximum number of iterations
        const int nMax = 25;
        for (int nIdx = 0; nIdx < nMax; ++nIdx) {
            // calculate Jacobian matrix and right hand side
            linearize_(J, b, fluidState, paramCache, phaseIdx, targetFug);
            Valgrind::CheckDefined(J);
            Valgrind::CheckDefined(b);

            /*
            std::cout << FluidSystem::phaseName(phaseIdx) << "Phase composition: ";
            for (unsigned i = 0; i < FluidSystem::numComponents; ++i)
                std::cout << fluidState.moleFraction(phaseIdx, i) << " ";
            std::cout << "\n";
            std::cout << FluidSystem::phaseName(phaseIdx) << "Phase phi: ";
            for (unsigned i = 0; i < FluidSystem::numComponents; ++i)
                std::cout << fluidState.fugacityCoefficient(phaseIdx, i) << " ";
            std::cout << "\n";
            */

            // Solve J*x = b
            x = 0.0;
            try { J.solve(x, b); }
            catch (Dune::FMatrixError e)
            { throw Opm::NumericalIssue(e.what()); }

            //std::cout << "original delta: " << x << "\n";

            Valgrind::CheckDefined(x);

            /*
            std::cout << FluidSystem::phaseName(phaseIdx) << "Phase composition: ";
            for (unsigned i = 0; i < FluidSystem::numComponents; ++i)
                std::cout << fluidState.moleFraction(phaseIdx, i) << " ";
            std::cout << "\n";
            std::cout << "J: " << J << "\n";
            std::cout << "rho: " << fluidState.density(phaseIdx) << "\n";
            std::cout << "delta: " << x << "\n";
            std::cout << "defect: " << b << "\n";

            std::cout << "J: " << J << "\n";

            std::cout << "---------------------------\n";
            */

            // update the fluid composition. b is also used to store
            // the defect for the next iteration.
            Scalar relError = update_(fluidState, paramCache, x, b, phaseIdx, targetFug);

            if (relError < 1e-9) {
                const Evaluation& rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
                fluidState.setDensity(phaseIdx, rho);

                //std::cout << "num iterations: " << nIdx << "\n";
                return;
            }
        }

        std::ostringstream oss;
        oss << "Calculating the " << FluidSystem::phaseName(phaseIdx)
            << "Phase composition failed. Initial {x} = {"
            << xInit
            << "}, {fug_t} = {" << targetFug << "}, p = " << fluidState.pressure(phaseIdx)
            << ", T = " << fluidState.temperature(phaseIdx);
        throw Opm::NumericalIssue(oss.str());
    }


protected:
    // update the phase composition in case the phase is an ideal
    // mixture, i.e. the component's fugacity coefficients are
    // independent of the phase's composition.
    template <class FluidState>
    static void solveIdealMix_(FluidState& fluidState,
                               typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                               unsigned phaseIdx,
                               const ComponentVector& fugacities)
    {
        for (unsigned i = 0; i < numComponents; ++ i) {
            const Evaluation& phi = FluidSystem::fugacityCoefficient(fluidState,
                                                                     paramCache,
                                                                     phaseIdx,
                                                                     i);
            const Evaluation& gamma = phi * fluidState.pressure(phaseIdx);
            Valgrind::CheckDefined(phi);
            Valgrind::CheckDefined(gamma);
            Valgrind::CheckDefined(fugacities[i]);
            fluidState.setFugacityCoefficient(phaseIdx, i, phi);
            fluidState.setMoleFraction(phaseIdx, i, fugacities[i]/gamma);
        };

        paramCache.updatePhase(fluidState, phaseIdx);

        const Evaluation& rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, rho);
        return;
    }

    template <class FluidState>
    static Scalar linearize_(Dune::FieldMatrix<Evaluation, numComponents, numComponents>& J,
                             Dune::FieldVector<Evaluation, numComponents>& defect,
                             FluidState& fluidState,
                             typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                             unsigned phaseIdx,
                             const ComponentVector& targetFug)
    {
        // reset jacobian
        J = 0;

        Scalar absError = 0;
        // calculate the defect (deviation of the current fugacities
        // from the target fugacities)
        for (unsigned i = 0; i < numComponents; ++ i) {
            const Evaluation& phi = FluidSystem::fugacityCoefficient(fluidState,
                                                          paramCache,
                                                          phaseIdx,
                                                          i);
            const Evaluation& f = phi*fluidState.pressure(phaseIdx)*fluidState.moleFraction(phaseIdx, i);
            fluidState.setFugacityCoefficient(phaseIdx, i, phi);

            defect[i] = targetFug[i] - f;
            absError = std::max(absError, std::abs(Opm::scalarValue(defect[i])));
        }

        // assemble jacobian matrix of the constraints for the composition
        static const Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e6;
        for (unsigned i = 0; i < numComponents; ++ i) {
            ////////
            // approximately calculate partial derivatives of the
            // fugacity defect of all components in regard to the mole
            // fraction of the i-th component. This is done via
            // forward differences

            // deviate the mole fraction of the i-th component
            Evaluation xI = fluidState.moleFraction(phaseIdx, i);
            fluidState.setMoleFraction(phaseIdx, i, xI + eps);
            paramCache.updateSingleMoleFraction(fluidState, phaseIdx, i);

            // compute new defect and derivative for all component
            // fugacities
            for (unsigned j = 0; j < numComponents; ++j) {
                // compute the j-th component's fugacity coefficient ...
                const Evaluation& phi = FluidSystem::fugacityCoefficient(fluidState,
                                                                         paramCache,
                                                                         phaseIdx,
                                                                         j);
                // ... and its fugacity ...
                const Evaluation& f =
                    phi *
                    fluidState.pressure(phaseIdx) *
                    fluidState.moleFraction(phaseIdx, j);
                // as well as the defect for this fugacity
                const Evaluation& defJPlusEps = targetFug[j] - f;

                // use forward differences to calculate the defect's
                // derivative
                J[j][i] = (defJPlusEps - defect[j])/eps;
            }

            // reset composition to original value
            fluidState.setMoleFraction(phaseIdx, i, xI);
            paramCache.updateSingleMoleFraction(fluidState, phaseIdx, i);

            // end forward differences
            ////////
        }

        return absError;
    }

    template <class FluidState>
    static Scalar update_(FluidState& fluidState,
                          typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache,
                          Dune::FieldVector<Evaluation, numComponents>& x,
                          Dune::FieldVector<Evaluation, numComponents>& /*b*/,
                          unsigned phaseIdx,
                          const Dune::FieldVector<Evaluation, numComponents>& targetFug)
    {
        // store original composition and calculate relative error
        Dune::FieldVector<Evaluation, numComponents> origComp;
        Scalar relError = 0;
        Evaluation sumDelta = 0.0;
        Evaluation sumx = 0.0;
        for (unsigned i = 0; i < numComponents; ++i) {
            origComp[i] = fluidState.moleFraction(phaseIdx, i);
            relError = std::max(relError, std::abs(Opm::scalarValue(x[i])));

            sumx += Opm::abs(fluidState.moleFraction(phaseIdx, i));
            sumDelta += Opm::abs(x[i]);
        }

        // chop update to at most 20% change in composition
        const Scalar maxDelta = 0.2;
        if (sumDelta > maxDelta)
            x /= (sumDelta/maxDelta);

        // change composition
        for (unsigned i = 0; i < numComponents; ++i) {
            Evaluation newx = origComp[i] - x[i];
            // only allow negative mole fractions if the target fugacity is negative
            if (targetFug[i] > 0)
                newx = Opm::max(0.0, newx);
            // only allow positive mole fractions if the target fugacity is positive
            else if (targetFug[i] < 0)
                newx = Opm::min(0.0, newx);
            // if the target fugacity is zero, the mole fraction must also be zero
            else
                newx = 0;

            fluidState.setMoleFraction(phaseIdx, i, newx);
        }

        paramCache.updateComposition(fluidState, phaseIdx);

        return relError;
    }

    template <class FluidState>
    static Scalar calculateDefect_(const FluidState& params,
                                   unsigned phaseIdx,
                                   const ComponentVector& targetFug)
    {
        Scalar result = 0.0;
        for (unsigned i = 0; i < numComponents; ++i) {
            // sum of the fugacity defect weighted by the inverse
            // fugacity coefficient
            result += std::abs(
                (targetFug[i] - params.fugacity(phaseIdx, i))
                /
                params.fugacityCoefficient(phaseIdx, i) );
        };
        return result;
    }
}; // namespace Opm

} // end namespace Opm

#endif
