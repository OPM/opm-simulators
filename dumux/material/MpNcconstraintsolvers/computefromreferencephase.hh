/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \brief Computes all quantities of a generic fluid state if a
 *        reference phase has been specified.
 *
 * This makes it is possible to specify just one phase and let the
 * remaining ones be calculated by the constraint solver. This
 * constraint solver assumes thermodynamic equilibrium
 */
#ifndef DUMUX_COMPUTE_FROM_REFERENCE_PHASE_HH
#define DUMUX_COMPUTE_FROM_REFERENCE_PHASE_HH

#include "../MpNcfluidstates/genericfluidstate.hh"
#include <dumux/material/MpNcconstraintsolvers/compositionfromfugacities.hh>

namespace Dumux {

/*!
 * \brief Computes all quantities of a generic fluid state if a
 *        reference phase has been specified.
 *
 * This makes it is possible to specify just one phase and let the
 * remaining ones be calculated by the constraint solver. This
 * constraint solver assumes thermodynamic equilibrium. It assumes the
 * following quantities to be set:
 *
 * - composition (mole+mass fractions) of the *reference* phase
 * - temperature of the *reference* phase
 * - saturations of *all* phases
 * - pressures of *all* phases
 *
 * after calling the solve() method the following quantities are
 * calculated in addition:
 *
 * - temperature of *all* phases
 * - density, molar density, molar volume of *all* phases
 * - composition in mole and mass fractions and molaries of *all* phases
 * - mean molar masses of *all* phases
 * - fugacity coefficients of *all* components in *all* phases
 * - if the setViscosity parameter is true, also dynamic viscosities of *all* phases
 * - if the setEnthalpy parameter is true, also specific enthalpies and internal energies of *all* phases
 */
template <class Scalar, class FluidSystem>
class ComputeFromReferencePhase
{
    typedef typename FluidSystem::MutableParameters MutableParameters;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    typedef Dumux::CompositionFromFugacities<Scalar, FluidSystem> CompositionFromFugacities;
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

public:
    /*!
     * \brief Computes all quantities of a generic fluid state if a
     *        reference phase has been specified.
     *
     * This makes it is possible to specify just one phase and let the
     * remaining ones be calculated by the constraint solver. This
     * constraint solver assumes thermodynamic equilibrium. It assumes the
     * following quantities to be set:
     *
     * - composition (mole+mass fractions) of the *reference* phase
     * - temperature of the *reference* phase
     * - saturations of *all* phases
     * - pressures of *all* phases
     *
     * after calling the solve() method the following quantities are
     * calculated in addition:
     *
     * - temperature of *all* phases
     * - density, molar density, molar volume of *all* phases
     * - composition in mole and mass fractions and molaries of *all* phases
     * - mean molar masses of *all* phases
     * - fugacity coefficients of *all* components in *all* phases
     * - if the setViscosity parameter is true, also dynamic viscosities of *all* phases
     * - if the setEnthalpy parameter is true, also specific enthalpies and internal energies of *all* phases
     *
     * \param mutParams The mutable parameters object which ought to be set
     * \param refPhaseIdx The phase index of the reference phase
     * \param setViscosity Specify whether the dynamic viscosity of
     *                     each phase should also be set.
     * \param setEnthalpy Specify whether the specific
     *                    enthalpy/internal energy of each phase
     *                    should also be set.
     */
    static void solve(MutableParameters &mutParams,
                      int refPhaseIdx,
                      bool setViscosity,
                      bool setEnthalpy)
    {
        ComponentVector fugVec;

        // compute the density and enthalpy of the
        // reference phase
        mutParams.updateMeanMolarMass(refPhaseIdx);
        mutParams.setMolarVolume(refPhaseIdx,
                                 FluidSystem::computeMolarVolume(mutParams, refPhaseIdx));

        if (setEnthalpy)
            mutParams.setEnthalpy(refPhaseIdx,
                                  FluidSystem::computeEnthalpy(mutParams, refPhaseIdx));

        if (setViscosity)
            mutParams.setViscosity(refPhaseIdx,
                                   FluidSystem::computeViscosity(mutParams, refPhaseIdx));

        // compute the fugacities of all components in the reference phase
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            mutParams.setFugacityCoeff(refPhaseIdx, compIdx,
                                       FluidSystem::computeFugacityCoeff(mutParams, refPhaseIdx, compIdx));
            fugVec[compIdx] = mutParams.fugacity(refPhaseIdx, compIdx);
        }

        // compute all quantities for the non-reference phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (phaseIdx == refPhaseIdx)
                continue; // reference phase is already calculated

            mutParams.setTemperature(phaseIdx, mutParams.temperature(refPhaseIdx));

            CompositionFromFugacities::guessInitial(mutParams, phaseIdx, fugVec);
            CompositionFromFugacities::solve(mutParams, phaseIdx, fugVec);

            if (setViscosity)
                mutParams.setViscosity(phaseIdx,
                                       FluidSystem::computeViscosity(mutParams, phaseIdx));

            if (setEnthalpy)
                mutParams.setEnthalpy(phaseIdx,
                                      FluidSystem::computeEnthalpy(mutParams, phaseIdx));
        }
    };
};

} // end namespace Dumux

#endif
