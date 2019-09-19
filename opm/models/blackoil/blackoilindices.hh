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
 * \copydoc Opm::BlackOilIndices
 */
#ifndef EWOMS_BLACK_OIL_INDICES_HH
#define EWOMS_BLACK_OIL_INDICES_HH

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief The primary variable and equation indices for the black-oil model.
 */
template <unsigned numSolventsV, unsigned numPolymersV, unsigned numEnergyV, bool enableFoam, unsigned PVOffset>
struct BlackOilIndices
{
    //! Number of phases active at all times
    static const int numPhases = 3;

    //! All phases are enabled
    static const bool oilEnabled = true;
    static const bool waterEnabled = true;
    static const bool gasEnabled = true;

    //! Are solvents involved?
    static const bool enableSolvent = numSolventsV > 0;

    //! Are polymers involved?
    static const bool enablePolymer = numPolymersV > 0;

    //! Shall energy be conserved?
    static const bool enableEnergy = numEnergyV > 0;

    //! Number of solvent components to be considered
    static const int numSolvents = enableSolvent ? numSolventsV : 0;

    //! Number of polymer components to be considered
    static const int numPolymers = enablePolymer ? numPolymersV : 0;

    //! Number of energy equations to be considered
    static const int numEnergy = enableEnergy ? numEnergyV : 0;

    //! Number of foam equations to be considered
    static const int numFoam = enableFoam? 1 : 0;

    //! The number of equations
    static const int numEq = numPhases + numSolvents + numPolymers + numEnergy + numFoam;

    //! \brief returns the index of "active" component
    static constexpr unsigned canonicalToActiveComponentIndex(unsigned compIdx)
    { return compIdx; }

    static constexpr unsigned activeToCanonicalComponentIndex(unsigned compIdx)
    { return compIdx; }

    ////////
    // Primary variable indices
    ////////

    //! The index of the water saturation
    static const int waterSaturationIdx = PVOffset + 0;

    //! Index of the oil pressure in a vector of primary variables
    static const int pressureSwitchIdx = PVOffset + 1;

    /*!
     * \brief Index of the switching variable which determines the composition of the
     *        hydrocarbon phases.
     *
     * Depending on the phases present, this variable is either interpreted as the
     * saturation of the gas phase, as the mole fraction of the gas component in the oil
     * phase or as the mole fraction of the oil component in the gas phase.
     */
    static const int compositionSwitchIdx = PVOffset + 2;

    //! Index of the primary variable for the first solvent
    static const int solventSaturationIdx =
        enableSolvent ? PVOffset + numPhases : -1000;

    //! Index of the primary variable for the first polymer
    static const int polymerConcentrationIdx =
        enablePolymer ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the primary variable for the second polymer primary variable (molecular weight)
    static const int polymerMoleWeightIdx =
        numPolymers > 1 ? polymerConcentrationIdx + 1 : -1000;

    //! Index of the primary variable for the foam
    static const int foamConcentrationIdx =
        enableFoam ? PVOffset + numPhases + numSolvents + numPolymers : -1000;

    //! Index of the primary variable for temperature
    static const int temperatureIdx  =
        enableEnergy ? PVOffset + numPhases + numSolvents + numPolymers + numFoam : - 1000;


    ////////
    // Equation indices
    ////////

    //! Index of the continuity equation of the first phase
    static const int conti0EqIdx = PVOffset + 0;
    // two continuity equations follow

    //! Index of the continuity equation for the first solvent component
    static const int contiSolventEqIdx =
        enableSolvent ? PVOffset + numPhases : -1000;

    //! Index of the continuity equation for the first polymer component
    static const int contiPolymerEqIdx =
        enablePolymer ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the continuity equation for the second polymer component (molecular weight)
    static const int contiPolymerMWEqIdx =
        numPolymers > 1 ? contiPolymerEqIdx + 1 : -1000;

    //! Index of the continuity equation for the foam component
    static const int contiFoamEqIdx =
        enableFoam ? PVOffset + numPhases + numSolvents + numPolymers : -1000;


    //! Index of the continuity equation for energy
    static const int contiEnergyEqIdx =
        enableEnergy ? PVOffset + numPhases + numSolvents + numPolymers + numFoam : -1000;
};

} // namespace Opm

#endif
