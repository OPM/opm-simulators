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
 * \copydoc Opm::BlackOilVariableAndEquationIndices
 */
#ifndef OPM_BLACK_OIL_VARIABLE_AND_EQUATION_INDICES_HH
#define OPM_BLACK_OIL_VARIABLE_AND_EQUATION_INDICES_HH

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief The primary variable and equation indices for the three-phase black-oil model.
 */
template<unsigned numSolventsV,
         unsigned numExtbosV,
         unsigned numPolymersV,
         unsigned numEnergyV,
         bool enableTemperature,
         bool enableFoam,
         bool enableBrine,
         unsigned PVOffset,
         unsigned numBioCompV>
struct BlackOilVariableAndEquationIndices
{
    //! Number of phases active at all times
    static constexpr int numPhases = 3;

    //! All phases are enabled
    static constexpr bool oilEnabled = true;
    static constexpr bool waterEnabled = true;
    static constexpr bool gasEnabled = true;

    //! Are solvents involved?
    static constexpr bool enableSolvent = numSolventsV > 0;

    //! Is extbo invoked?
    static constexpr bool enableExtbo = numExtbosV > 0;

    //! Are polymers involved?
    static constexpr bool enablePolymer = numPolymersV > 0;

    //! Shall energy be conserved?
    static constexpr bool enableEnergy = numEnergyV > 0;

    //! No bioeffects for three phase indices
    static constexpr bool enableMICP = false;
    static constexpr bool enableBiofilm = false;
    static constexpr int numBioInWat = 0;

    //! Number of solvent components to be considered
    static constexpr int numSolvents = enableSolvent ? numSolventsV : 0;

    //! Number of components to be considered for extbo
    static constexpr int numExtbos = enableExtbo ? numExtbosV : 0;

    //! Number of polymer components to be considered
    static constexpr int numPolymers = enablePolymer ? numPolymersV : 0;

    //! Number of energy equations to be considered
    static constexpr int numEnergy = enableEnergy ? numEnergyV : 0;

    //! Number of foam equations to be considered
    static constexpr int numFoam = enableFoam? 1 : 0;

    //! Number of salt equations to be considered
    static constexpr int numBrine = enableBrine? 1 : 0;

    //! The number of equations
    static constexpr int numEq = numPhases + numSolvents + numExtbos + numPolymers +
                                 numEnergy + numFoam + numBrine;

    static constexpr int numDerivatives = numEq + enableTemperature;


    ////////
    // Primary variable indices
    ////////

    /*!
     * \brief Index of the switching variable which determines the composistion of the water phase
     *
     * Depending on the phases present, this variable is either interpreted as
     * water saturation or vapporized water in gas phase
     */
    static constexpr int waterSwitchIdx = PVOffset + 0;

    /*!
     * \brief Index of the switching variable which determines the pressure
     *
     * Depending on the phases present, this variable is either interpreted as the
     * pressure of the oil phase, gas phase (if no oil) or water phase (if only water)
     */
    static constexpr int pressureSwitchIdx = PVOffset + 1;

    /*!
     * \brief Index of the switching variable which determines the composition of the
     *        hydrocarbon phases.
     *
     * Depending on the phases present, this variable is either interpreted as the
     * saturation of the gas phase, as the mole fraction of the gas component in the oil
     * phase or as the mole fraction of the oil component in the gas phase.
     */
    static constexpr int compositionSwitchIdx = PVOffset + 2;

    //! Index of the primary variable for the first solvent
    static constexpr int solventSaturationIdx =
        enableSolvent ? PVOffset + numPhases : -1000;

    //! Index of the primary variable for the first extbo component
    static constexpr int zFractionIdx =
        enableExtbo ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the primary variable for the first polymer
    static constexpr int polymerConcentrationIdx =
        enablePolymer ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the primary variable for the second polymer primary variable (molecular weight)
    static constexpr int polymerMoleWeightIdx =
        numPolymers > 1 ? polymerConcentrationIdx + 1 : -1000;

    //! No bioeffects for three phase indices
    static constexpr int microbialConcentrationIdx = -1000;
    static constexpr int oxygenConcentrationIdx = -1000;
    static constexpr int ureaConcentrationIdx = -1000;
    static constexpr int biofilmVolumeFractionIdx = -1000;
    static constexpr int calciteVolumeFractionIdx = -1000;

    //! Index of the primary variable for the foam
    static constexpr int foamConcentrationIdx =
        enableFoam ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers : -1000;

    //! Index of the primary variable for the brine
    static constexpr int saltConcentrationIdx =
        enableBrine ? PVOffset + numPhases + numSolvents + numExtbos + numExtbos + numPolymers + numFoam : -1000;

    //! Index of the primary variable for temperature
    static constexpr int temperatureIdx  =
        (enableTemperature || enableEnergy) ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers + numFoam + numBrine : - 1000;


    ////////
    // Equation indices
    ////////

    //! Index of the continuity equation of the first phase
    static constexpr int conti0EqIdx = PVOffset + 0;
    // two continuity equations follow

    //! Index of the continuity equation for the first solvent component
    static constexpr int contiSolventEqIdx =
        enableSolvent ? PVOffset + numPhases : -1000;

    //! Index of the continuity equation for the first extbo component
    static constexpr int contiZfracEqIdx =
        enableExtbo ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the continuity equation for the first polymer component
    static constexpr int contiPolymerEqIdx =
        enablePolymer ? PVOffset + numPhases + numSolvents + numExtbos : -1000;

    //! Index of the continuity equation for the second polymer component (molecular weight)
    static constexpr int contiPolymerMWEqIdx =
        numPolymers > 1 ? contiPolymerEqIdx + 1 : -1000;
    
    //! No bioeffects for three phase indices
    static constexpr int contiMicrobialEqIdx = -1000;
    static constexpr int contiOxygenEqIdx = -1000;
    static constexpr int contiUreaEqIdx = -1000;
    static constexpr int contiBiofilmEqIdx = -1000;
    static constexpr int contiCalciteEqIdx = -1000;

    //! Index of the continuity equation for the foam component
    static constexpr int contiFoamEqIdx =
        enableFoam ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers : -1000;

    //! Index of the continuity equation for the salt water component
    static constexpr int contiBrineEqIdx =
        enableBrine ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers + numFoam : -1000;

    //! Index of the continuity equation for energy
    static constexpr int contiEnergyEqIdx =
        enableEnergy ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers + numFoam + numBrine: -1000;
};

} // namespace Opm

#endif
