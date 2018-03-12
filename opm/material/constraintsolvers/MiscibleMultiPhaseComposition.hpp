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
 * \copydoc Opm::MiscibleMultiPhaseComposition
 */
#ifndef OPM_MISCIBLE_MULTIPHASE_COMPOSITION_HPP
#define OPM_MISCIBLE_MULTIPHASE_COMPOSITION_HPP

#include <opm/material/common/MathToolbox.hpp>

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

namespace Opm {

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

    MMPCAuxConstraint(unsigned phaseIndex, unsigned compIndex, Scalar val)
        : phaseIdx_(phaseIndex)
        , compIdx_(compIndex)
        , value_(val)
    {}

    /*!
     * \brief Specify the auxiliary constraint.
     */
    void set(unsigned phaseIndex, unsigned compIndex, Scalar val)
    {
        phaseIdx_ = phaseIndex;
        compIdx_ = compIndex;
        value_ = val;
    }

    /*!
     * \brief Returns the index of the fluid phase for which the
     *        auxiliary constraint is specified.
     */
    unsigned phaseIdx() const
    { return phaseIdx_; }

    /*!
     * \brief Returns the index of the component for which the
     *        auxiliary constraint is specified.
     */
    unsigned compIdx() const
    { return compIdx_; }

    /*!
     * \brief Returns value of the mole fraction of the auxiliary
     *        constraint.
     */
    Scalar value() const
    { return value_; }

private:
    unsigned phaseIdx_;
    unsigned compIdx_;
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
template <class Scalar, class FluidSystem, class Evaluation = Scalar>
class MiscibleMultiPhaseComposition
{
    static const int numPhases = FluidSystem::numPhases;
    static const int numComponents = FluidSystem::numComponents;

    typedef MathToolbox<Evaluation> Toolbox;

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
    static void solve(FluidState& fluidState,
                      ParameterCache& paramCache,
                      int phasePresence,
                      const MMPCAuxConstraint<Evaluation>* auxConstraints,
                      unsigned numAuxConstraints,
                      bool setViscosity,
                      bool setInternalEnergy)
    {
        static_assert(std::is_same<typename FluidState::Scalar, Evaluation>::value,
                      "The scalar type of the fluid state must be 'Evaluation'");

#ifndef NDEBUG
        // currently this solver can only handle fluid systems which
        // assume ideal mixtures of all fluids. TODO: relax this
        // (requires solving a non-linear system of equations, i.e. using
        // newton method.)
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(FluidSystem::isIdealMixture(phaseIdx));
        }
#endif

        // compute all fugacity coefficients
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            paramCache.updatePhase(fluidState, phaseIdx);

            // since we assume ideal mixtures, the fugacity
            // coefficients of the components cannot depend on
            // composition, i.e. the parameters in the cache are valid
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                Evaluation fugCoeff = Opm::decay<Evaluation>(
                    FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx));
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
            }
        }

        // create the linear system of equations which defines the
        // mole fractions
        static const int numEq = numComponents*numPhases;
        Dune::FieldMatrix<Evaluation, numEq, numEq> M(0.0);
        Dune::FieldVector<Evaluation, numEq> x(0.0);
        Dune::FieldVector<Evaluation, numEq> b(0.0);

        // assemble the equations expressing the fact that the
        // fugacities of each component are equal in all phases
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            const Evaluation& entryCol1 =
                fluidState.fugacityCoefficient(/*phaseIdx=*/0, compIdx)
                *fluidState.pressure(/*phaseIdx=*/0);
            unsigned col1Idx = compIdx;

            for (unsigned phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
                unsigned rowIdx = (phaseIdx - 1)*numComponents + compIdx;
                unsigned col2Idx = phaseIdx*numComponents + compIdx;

                const Evaluation& entryCol2 =
                    fluidState.fugacityCoefficient(phaseIdx, compIdx)
                    *fluidState.pressure(phaseIdx);

                M[rowIdx][col1Idx] = entryCol1;
                M[rowIdx][col2Idx] = -entryCol2;
            }
        }

        // assemble the equations expressing the assumption that the
        // sum of all mole fractions in each phase must be 1 for the
        // phases present.
        unsigned presentPhases = 0;
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!(phasePresence&  (1 << phaseIdx)))
                continue;

            unsigned rowIdx = numComponents*(numPhases - 1) + presentPhases;
            presentPhases += 1;

            b[rowIdx] = 1.0;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                unsigned colIdx = phaseIdx*numComponents + compIdx;

                M[rowIdx][colIdx] = 1.0;
            }
        }

        assert(presentPhases + numAuxConstraints == numComponents);

        // incorperate the auxiliary equations, i.e., the explicitly given mole fractions
        for (unsigned auxEqIdx = 0; auxEqIdx < numAuxConstraints; ++auxEqIdx) {
            unsigned rowIdx = numComponents*(numPhases - 1) + presentPhases + auxEqIdx;
            b[rowIdx] = auxConstraints[auxEqIdx].value();

            unsigned colIdx = auxConstraints[auxEqIdx].phaseIdx()*numComponents + auxConstraints[auxEqIdx].compIdx();
            M[rowIdx][colIdx] = 1.0;
        }

        // solve for all mole fractions
        try {
#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2,7)
            static constexpr Scalar eps = std::numeric_limits<Scalar>::min()*1000.0;
            Dune::FMatrixPrecision<Scalar>::set_singular_limit(eps);
#endif
            M.solve(x, b);
        }
        catch (const Dune::FMatrixError& e) {
            std::ostringstream oss;
            oss << "Numerical problem in MiscibleMultiPhaseComposition::solve(): " << e.what() << "; M="<<M;
            throw NumericalIssue(oss.str());
        }
        catch (...) {
            throw;
        }


        // set all mole fractions and the additional quantities in
        // the fluid state
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                unsigned rowIdx = phaseIdx*numComponents + compIdx;
                fluidState.setMoleFraction(phaseIdx, compIdx, x[rowIdx]);
            }
            paramCache.updateComposition(fluidState, phaseIdx);

            const Evaluation& rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, rho);

            if (setViscosity) {
                const Evaluation& mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
                fluidState.setViscosity(phaseIdx, mu);
            }

            if (setInternalEnergy) {
                const Evaluation& h =  FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
                fluidState.setEnthalpy(phaseIdx, h);
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
    static void solve(FluidState& fluidState,
                      ParameterCache& paramCache,
                      bool setViscosity,
                      bool setInternalEnergy)
    {
        solve(fluidState,
              paramCache,
              /*phasePresence=*/0xffffff,
              /*numAuxConstraints=*/0,
              /*auxConstraints=*/0,
              setViscosity,
              setInternalEnergy);
    }
};

} // namespace Opm

#endif
