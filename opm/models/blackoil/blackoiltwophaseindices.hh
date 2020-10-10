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
 *
 * \copydoc Opm::BlackOilTwoPhaseIndices
 */
#ifndef EWOMS_BLACK_OIL_TWO_PHASE_INDICES_HH
#define EWOMS_BLACK_OIL_TWO_PHASE_INDICES_HH

#include <cassert>

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief The primary variable and equation indices for the black-oil model.
 */
template <unsigned numSolventsV, unsigned numExtbosV, unsigned numPolymersV, unsigned numEnergyV, bool enableFoam, bool enableBrine, unsigned PVOffset, unsigned disabledCanonicalCompIdx>
struct BlackOilTwoPhaseIndices
{
    //! Is phase enabled or not
    static const bool oilEnabled = (disabledCanonicalCompIdx != 0);
    static const bool waterEnabled = (disabledCanonicalCompIdx != 1);
    static const bool gasEnabled = (disabledCanonicalCompIdx != 2);

    //! Are solvents involved?
    static const bool enableSolvent = numSolventsV > 0;

    //! Is extbo invoked?
    static const bool enableExtbo = numExtbosV > 0;

    //! Are polymers involved?
    static const bool enablePolymer = numPolymersV > 0;

    //! Shall energy be conserved?
    static const bool enableEnergy = numEnergyV > 0;

    //! Number of solvent components to be considered
    static const int numSolvents = enableSolvent ? numSolventsV : 0;

    //! Number of components to be considered for extbo
    static const int numExtbos = enableExtbo ? numExtbosV : 0;

    //! Number of polymer components to be considered
    static const int numPolymers = enablePolymer ? numPolymersV : 0;

    //! Number of energy equations to be considered
    static const int numEnergy = enableEnergy ? numEnergyV : 0;

    //! Number of foam equations to be considered
    static const int numFoam = enableFoam? 1 : 0;

    //! Number of salt equations to be considered
    static const int numBrine = enableBrine? 1 : 0;

    //! The number of fluid phases
    static const int numPhases = 2;

    //! The number of equations
    static const int numEq = numPhases + numSolvents + numExtbos + numPolymers + numEnergy + numFoam + numBrine;

    //////////////////////////////
    // Primary variable indices
    //////////////////////////////

    //! The index of the water saturation. For two-phase oil gas models this is disabled.
    static const int waterSaturationIdx  = waterEnabled ? PVOffset + 0 : -10000;

    //! Index of the oil pressure in a vector of primary variables
    static const int pressureSwitchIdx  = waterEnabled ? PVOffset + 1 : PVOffset + 0;

    /*!
     * \brief Index of the switching variable which determines the composition of the
     *        hydrocarbon phases.
     *
     * \note For two-phase water oil models this is disabled.
     */
    static const int compositionSwitchIdx = gasEnabled ? PVOffset + 1 : -10000;

    //! Index of the primary variable for the first solvent
    static const int solventSaturationIdx =
        enableSolvent ? PVOffset + numPhases : -1000;

    //! Index of the primary variable for the first extbo component
    static const int zFractionIdx =
        enableExtbo ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the primary variable for the first polymer
    static const int polymerConcentrationIdx =
        enablePolymer ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the primary variable for the second polymer primary variable (molecular weight)
    static const int polymerMoleWeightIdx =
        numPolymers > 1 ? polymerConcentrationIdx + 1 : -1000;

    //! Index of the primary variable for the foam
    static const int foamConcentrationIdx =
        enableFoam ? PVOffset + numPhases + numSolvents + numPolymers : -1000;

    //! Index of the primary variable for the salt
    static const int saltConcentrationIdx =
        enableBrine ? PVOffset + numPhases + numSolvents + numPolymers + numFoam : -1000;

    //! Index of the primary variable for temperature
    static const int temperatureIdx  =
        enableEnergy ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers + numFoam + numBrine : - 1000;

    //////////////////////
    // Equation indices
    //////////////////////

    //! \brief returns the index of "active" component
    static unsigned canonicalToActiveComponentIndex(unsigned compIdx)
    {
        // assumes canonical oil = 0, water = 1, gas = 2;
        if(!gasEnabled) {
            assert(compIdx != 2);
            // oil = 0, water = 1
            return compIdx;
        } else if (!waterEnabled) {
            assert(compIdx != 1);
            // oil = 0, gas = 1
            return compIdx / 2;
        } else {
            assert(!oilEnabled);
            assert(compIdx != 0);
        }
        // water = 0, gas = 1;
        return compIdx-1;
    }

    static unsigned activeToCanonicalComponentIndex(unsigned compIdx)
    {
        // assumes canonical oil = 0, water = 1, gas = 2;
        assert(compIdx < 2);
        if(!gasEnabled) {
            // oil = 0, water = 1
            return compIdx;
        } else if (!waterEnabled) {
            // oil = 0, gas = 1
            return compIdx * 2;
        } else {
            assert(!oilEnabled);
        }
        // water = 0, gas = 1;
        return compIdx+1;
    }

    //! Index of the continuity equation of the first phase
    static const int conti0EqIdx = PVOffset + 0;
    // one continuity equation follows

    //! Index of the continuity equation for the first solvent component
    static const int contiSolventEqIdx =
        enableSolvent ? PVOffset + numPhases : -1000;

    //! Index of the continuity equation for the first extbo component
    static const int contiZfracEqIdx =
        enableExtbo ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the continuity equation for the first polymer component
    static const int contiPolymerEqIdx =
        enablePolymer ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the continuity equation for the second polymer component (molecular weight)
    static const int contiPolymerMWEqIdx =
        numPolymers > 1 ? contiPolymerEqIdx + 1 : -1000;

    //! Index of the continuity equation for the foam component
    static const int contiFoamEqIdx =
        enableFoam ? PVOffset + numPhases + numSolvents + numPolymers : -1000;

    //! Index of the continuity equation for the salt component
    static const int contiBrineEqIdx =
        enableBrine ? PVOffset + numPhases + numSolvents + numPolymers + numFoam : -1000;

    //! Index of the continuity equation for energy
    static const int contiEnergyEqIdx =
        enableEnergy ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers + numFoam + numBrine : -1000;
};

} // namespace Opm

#endif
