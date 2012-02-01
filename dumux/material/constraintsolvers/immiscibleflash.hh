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
 * \brief Determines the pressures and saturations of all fluid phases
 *        given the total mass of all components.
 */
#ifndef DUMUX_IMMISCIBLE_FLASH_HH
#define DUMUX_IMMISCIBLE_FLASH_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {

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
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;
    static_assert(numPhases == numComponents, 
                  "Immiscibility assumes that the number of phases"
                  " is equal to the number of components");

    typedef typename FluidSystem::ParameterCache ParameterCache;

    static constexpr int numEq = numPhases;

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
                      ParameterCache &paramCache,
                      const typename MaterialLaw::Params &matParams,
                      const ComponentVector &globalMolarities)
    {
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-25);

        /////////////////////////
        // Check if all fluid phases are incompressible
        /////////////////////////
        bool allIncompressible = true;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::isCompressible(phaseIdx)) {
                allIncompressible = false;
                break;
            }
        };
        
        if (allIncompressible) {
            paramCache.updateAll(fluidState);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
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

        completeFluidState_<MaterialLaw>(fluidState, paramCache, matParams);

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
                std::cout << "error: " << e << "\n";
                std::cout << "b: " << b << "\n";
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
            std::cout << "residual: " << b << "\n";
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
            const Scalar eps = 1e-10/quantityWeight_(fluidState, pvIdx);
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
        // global molarities are given
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            b[phaseIdx] = 
                fluidState.saturation(phaseIdx)
                * fluidState.molarity(phaseIdx, /*compIdx=*/phaseIdx);
            
            b[phaseIdx] -= globalMolarities[/*compIdx=*/phaseIdx];
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
            else if (isPressureIdx_(pvIdx)) {
                // dampen to at most 30% change in pressure per
                // iteration
                delta = std::min(0.30*fluidState.pressure(0), std::max(-0.30*fluidState.pressure(0), delta));
            };
            
            setQuantityRaw_(fluidState, pvIdx, tmp - delta);
        }
       
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

        if (sumSat > 1.0) {
            // make sure that the last saturation does not become
            // negative
            for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
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
        for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
            fluidState.setPressure(phaseIdx, 
                                   fluidState.pressure(0)
                                   + (pC[phaseIdx] - pC[0]));

        // update the parameter cache
        paramCache.updateAll(fluidState, /*except=*/ParameterCache::Temperature|ParameterCache::Composition);

        // update all densities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);
        }
    }

    static bool isPressureIdx_(int pvIdx)
    { return pvIdx == 0; }

    static bool isSaturationIdx_(int pvIdx)
    { return 1 <= pvIdx && pvIdx < numPhases; }

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
        // saturations
        else { 
            assert(pvIdx < numPhases);
            int phaseIdx = pvIdx - 1;
            return fs.saturation(phaseIdx);
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

        if (pvIdx < 1) {
            // -> first pressure
            Scalar delta = value - fs.pressure(0);

            // set all pressures. here we assume that the capillary
            // pressure does not depend on absolute pressure.
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx, fs.pressure(phaseIdx) + delta);
            paramCache.updateAllPressures(fs);

            // update all densities
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
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
            for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx, 
                               fs.pressure(0)
                               + (pC[phaseIdx] - pC[0]));
            paramCache.updateAllPressures(fs);

            // update all densities
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
                fs.setDensity(phaseIdx, rho);
            }
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
        // saturations
        else {
            assert(pvIdx < numPhases);
            int phaseIdx = pvIdx - 1;

            // make sure that the first M-1 saturations does not get
            // negative
            value = std::max(0.0, value);
            fs.setSaturation(phaseIdx, value);
        }
    }

    template <class FluidState>
    static Scalar quantityWeight_(const FluidState &fs, int pvIdx)
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

} // end namespace Dumux

#endif
