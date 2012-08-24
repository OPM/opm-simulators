// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Klaus Mosthaf                                     *
 *   Copyright (C) 2012 by Bernd Flemisch                                    *
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
 * \brief Computes the composition of all phases of a N-phase,
 *        N-component fluid system assuming that all N phases are
 *        present
 */
#ifndef DUMUX_MISCIBLE_MULTIPHASE_COMPOSITION_HH
#define DUMUX_MISCIBLE_MULTIPHASE_COMPOSITION_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {

/*!
 * \brief Specifies an auxiliary constraint for the
 *        MiscibleMultiPhaseComposition constraint solver.
 *
 * For this constraint solver, an auxiliary constraint is defined as a
 * fixed mole fraction of a component in a fluid phase.
 */
template <class Scalar>
class MMPCAuxConstraint
{
public:
    MMPCAuxConstraint()
    {}

    MMPCAuxConstraint(int phaseIdx, int compIdx, Scalar value)
        : phaseIdx_(phaseIdx)
        , compIdx_(compIdx)
        , value_(value)
    {}

    /*!
     * \brief Specify the auxiliary constraint.
     */
    void set(int phaseIdx, int compIdx, Scalar value)
    {
        phaseIdx_ = phaseIdx;
        compIdx_ = compIdx;
        value_ = value;
    }
    
    /*!
     * \brief Returns the index of the fluid phase for which the
     *        auxiliary constraint is specified.
     */
    int phaseIdx() const
    { return phaseIdx_; }

    /*!
     * \brief Returns the index of the component for which the
     *        auxiliary constraint is specified.
     */
    int compIdx() const
    { return compIdx_; }

    /*!
     * \brief Returns value of the mole fraction of the auxiliary
     *        constraint.
     */
    Scalar value() const
    { return value_; }

private:
    int phaseIdx_;
    int compIdx_;
    Scalar value_;
};

/*!
 * \brief Computes the composition of all phases of a N-phase,
 *        N-component fluid system assuming that all N phases are
 *        present
 *
 * The constraint solver assumes the following quantities to be set:
 *
 * - temperatures of *all* phases
 * - saturations of *all* phases
 * - pressures of *all* phases
 *
 * It also assumes that the mole/mass fractions of all phases sum up
 * to 1. After calling the solve() method the following quantities
 * are calculated in addition:
 *
 * - temperature of *all* phases
 * - density, molar density, molar volume of *all* phases
 * - composition in mole and mass fractions and molarities of *all* phases
 * - mean molar masses of *all* phases
 * - fugacity coefficients of *all* components in *all* phases
 * - if the setViscosity parameter is true, also dynamic viscosities of *all* phases
 * - if the setInternalEnergy parameter is true, also specific enthalpies and internal energies of *all* phases
 */
template <class Scalar, class FluidSystem>
class MiscibleMultiPhaseComposition
{
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;

    static_assert(numPhases <= numComponents,
                  "This solver requires that the number fluid phases is smaller or equal "
                  "to the number of components");


public:
    /*!
     * \brief Computes the composition of all phases of a N-phase,
     *        N-component fluid system assuming that all N phases are
     *        present
     *
     * The constraint solver assumes the following quantities to be set:
     *
     * - temperatures of *all* phases
     * - saturations of *all* phases
     * - pressures of *all* phases
     *
     * It also assumes that the mole/mass fractions of all phases sum up
     * to 1. After calling the solve() method the following quantities
     * are calculated in addition:
     *
     * - temperature of *all* phases
     * - density, molar density, molar volume of *all* phases
     * - composition in mole and mass fractions and molarities of *all* phases
     * - mean molar masses of *all* phases
     * - fugacity coefficients of *all* components in *all* phases
     * - if the setViscosity parameter is true, also dynamic viscosities of *all* phases
     * - if the setInternalEnergy parameter is true, also specific enthalpies and internal energies of *all* phases
     */
    template <class FluidState, class ParameterCache>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      int phaseState,
                      const MMPCAuxConstraint<Scalar> *auxConstraints,
                      int numAuxConstraints,
                      bool setViscosity,
                      bool setInternalEnergy)
    {
#ifndef NDEBUG
        // currently this solver can only handle fluid systems which
        // assume ideal mixtures of all fluids. TODO: relax this
        // (requires solving a non-linear system of equations, i.e. using
        // newton method.)
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(FluidSystem::isIdealMixture(phaseIdx));
        }
#endif

        // compute all fugacity coefficients
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            paramCache.updatePhase(fluidState, phaseIdx);

            // since we assume ideal mixtures, the fugacity
            // coefficients of the components cannot depend on
            // composition, i.e. the parameters in the cache are valid
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar fugCoeff = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
            }
        }

        // create the linear system of equations which defines the
        // mole fractions
        static constexpr int numEq = numComponents*numPhases;
        Dune::FieldMatrix<Scalar, numEq, numEq> M(0.0);
        Dune::FieldVector<Scalar, numEq> x(0.0);
        Dune::FieldVector<Scalar, numEq> b(0.0);

        // assemble the equations expressing the fact that the
        // fugacities divided by the pressure of each component is
        // equal in all phases (assuming thermal equilibrium, that is)
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            Scalar entryCol1 =
                fluidState.fugacityCoefficient(/*phaseIdx=*/0, compIdx);
            int col1Idx = compIdx;

            for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
                int rowIdx = (phaseIdx - 1)*numComponents + compIdx;
                int col2Idx = phaseIdx*numComponents + compIdx;

                Scalar entryCol2 =
                    fluidState.fugacityCoefficient(phaseIdx, compIdx);

                M[rowIdx][col1Idx] = entryCol1;
                M[rowIdx][col2Idx] = -entryCol2;
            }
        }

        // assemble the equations expressing the assumption that the
        // sum of all mole fractions in each phase must be 1 for the
        // phases present.
        int presentPhases = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!(phaseState & (1 << phaseIdx)))
                continue;

            int rowIdx = numComponents*(numPhases - 1) + presentPhases;
            presentPhases += 1;

            b[rowIdx] = 1.0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                int colIdx = phaseIdx*numComponents + compIdx;

                M[rowIdx][colIdx] = 1.0;
            }
        }

        assert(presentPhases + numAuxConstraints == numComponents);
        
        // assemble the equations expressing the assumption that the
        // sum of all mole fractions in each phase must be 1 for the
        // phases present.
        for (int auxEqIdx = 0; auxEqIdx < numAuxConstraints; ++auxEqIdx) {
            int rowIdx = numComponents*(numPhases - 1) + presentPhases + auxEqIdx;
            b[rowIdx] = auxConstraints[auxEqIdx].value();

            int colIdx = auxConstraints[auxEqIdx].phaseIdx()*numComponents + auxConstraints[auxEqIdx].compIdx();
            M[rowIdx][colIdx] = 1.0;
        }

        // solve for all mole fractions
        try {
            Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-50);
            M.solve(x, b);
        }
        catch (const Dune::FMatrixError &e) {
            DUNE_THROW(NumericalProblem,
                       "Numerical problem in MiscibleMultiPhaseComposition::solve(): " << NumericalProblem(e.what()) << "; M="<<M);
        }
        catch (...) {
            throw;
        }


        // set all mole fractions and the additional quantities in
        // the fluid state
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                int rowIdx = phaseIdx*numComponents + compIdx;
                fluidState.setMoleFraction(phaseIdx, compIdx, x[rowIdx]);
            }
            paramCache.updateComposition(fluidState, phaseIdx);

            Scalar value = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, value);

            if (setViscosity) {
                value = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
                fluidState.setViscosity(phaseIdx, value);
            }

            if (setInternalEnergy) {
                value = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
                fluidState.setEnthalpy(phaseIdx, value);
            }
        }
    }

    /*!
     * \brief Computes the composition of all phases of a N-phase,
     *        N-component fluid system assuming that all N phases are
     *        present
     *
     * This is a convenience method where no auxiliary constraints are used.
     */
    template <class FluidState, class ParameterCache>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      bool setViscosity,
                      bool setInternalEnergy)
    {
        solve(fluidState,
              paramCache,
              /*phaseState=*/0xffffff,
              /*numAuxConstraints=*/0,
              /*auxConstraints=*/0,
              setViscosity,
              setInternalEnergy);
    }
};

} // end namespace Dumux

#endif
