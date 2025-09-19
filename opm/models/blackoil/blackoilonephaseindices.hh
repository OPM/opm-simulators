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
 * \copydoc Ewoms::BlackOilOnePhaseIndices
 */
#ifndef EWOMS_BLACK_OIL_ONE_PHASE_INDICES_HH
#define EWOMS_BLACK_OIL_ONE_PHASE_INDICES_HH

#include <cassert>

#include <opm/common/utility/ConstexprAssert.hpp>

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief The primary variable and equation indices for the black-oil model.
 */
template<unsigned numSolventsV,
         unsigned numExtbosV,
         unsigned numPolymersV,
         unsigned numEnergyV,
         bool enableFoam,
         bool enableBrine,
         unsigned PVOffset,
         unsigned canonicalCompIdx,
         unsigned numBioCompV>
struct BlackOilOnePhaseIndices
{
    //! Is phase enabled or not
    static constexpr bool oilEnabled = canonicalCompIdx == 0;
    static constexpr bool waterEnabled = canonicalCompIdx == 1;
    static constexpr bool gasEnabled = canonicalCompIdx == 2;

    //! Are solvents involved?
    static constexpr bool enableSolvent = numSolventsV > 0;

    //! Is extbo invoked?
    static constexpr bool enableExtbo = numExtbosV > 0;

    //! Are polymers involved?
    static constexpr bool enablePolymer = numPolymersV > 0;

    //! Shall energy be conserved?
    static constexpr bool enableEnergy = numEnergyV > 0;

    //! Is MICP involved? (microbes, oxygen, urea, biofilm, and calcite)
    static constexpr bool enableMICP = numBioCompV == 5;

    //! Biofilm effects on co2/h2store only for two phase indices  
    static constexpr bool enableBiofilm = false;

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

    //! The number of fluid phases
    static constexpr int numPhases = 1;

    //! Number of MICP components to be considered
    static constexpr int numMICPs = enableMICP ? numBioCompV : 0;

    //! Number of biocomponents in the water phase
    static constexpr int numBioInWat = enableMICP ? 3 : 0;

    //! The number of equations
    static constexpr int numEq = numPhases + numSolvents + numExtbos + numPolymers +
                                 numEnergy + numFoam + numBrine + numMICPs;

    //////////////////////////////
    // Primary variable indices
    //////////////////////////////

    /*!
     * \brief Index of the switching variable which determines the composistion of the water phase
     *
     * Depending on the phases present, this variable is either interpreted as
     * water saturation or vapporized water in gas phase.
     *
     * \note For one-phase models this is disabled.
     */
    static constexpr int waterSwitchIdx  = -10000;

    /*!
     * \brief Index of the switching variable which determines the pressure
     *
     * Depending on the phases present, this variable is either interpreted as the
     * pressure of the oil phase, gas phase (if no oil) or water phase (if only water)
     */
     static constexpr int pressureSwitchIdx  = PVOffset + 0;

    /*!
     * \brief Index of the switching variable which determines the composition of the
     *        hydrocarbon phases.
     *
     * \note For one-phase models this is disabled.
     */
    static constexpr int compositionSwitchIdx = -10000;

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

    //! Index of the primary variable for the first MICP component
    static constexpr int microbialConcentrationIdx =
        enableMICP ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the primary variable for the second MICP component
    static constexpr int oxygenConcentrationIdx =
        numMICPs > 1 ? microbialConcentrationIdx + 1 : -1000;

    //! Index of the primary variable for the third MICP component
    static constexpr int ureaConcentrationIdx =
        numMICPs > 2 ? oxygenConcentrationIdx + 1 : -1000;

    //! Index of the primary variable for the fourth MICP component
    static constexpr int biofilmVolumeFractionIdx =
        numMICPs > 3 ? ureaConcentrationIdx + 1 : -1000;

    //! Index of the primary variable for the fifth MICP component
    static constexpr int calciteVolumeFractionIdx =
        numMICPs > 4 ? biofilmVolumeFractionIdx + 1 : -1000;

    //! Index of the primary variable for the foam
    static constexpr int foamConcentrationIdx =
        enableFoam ? PVOffset + numPhases + numSolvents + numPolymers + numMICPs : -1000;

    //! Index of the primary variable for the salt
    static constexpr int saltConcentrationIdx =
        enableBrine ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers + numMICPs + numFoam : -1000;

    //! Index of the primary variable for temperature
    static constexpr int temperatureIdx  =
        enableEnergy ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers + numMICPs + numFoam + numBrine: - 1000;

    //////////////////////
    // Equation indices
    //////////////////////

    //! Index of the continuity equation of the first (and only) phase
    static constexpr int conti0EqIdx = PVOffset + 0;

    //! Index of the continuity equation for the first solvent component
    static constexpr int contiSolventEqIdx =
        enableSolvent ? PVOffset + numPhases : -1000;

    //! Index of the continuity equation for the first extbo component
    static constexpr int contiZfracEqIdx =
        enableExtbo ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the continuity equation for the first polymer component
    static constexpr int contiPolymerEqIdx =
        enablePolymer ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the continuity equation for the second polymer component (molecular weight)
    static constexpr int contiPolymerMWEqIdx =
        numPolymers > 1 ? contiPolymerEqIdx + 1 : -1000;

    //! Index of the continuity equation for the first MICP component
    static constexpr int contiMicrobialEqIdx =
        enableMICP ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the continuity equation for the second MICP component
    static constexpr int contiOxygenEqIdx =
        numMICPs > 1 ? contiMicrobialEqIdx + 1 : -1000;

    //! Index of the continuity equation for the third MICP component
    static constexpr int contiUreaEqIdx =
        numMICPs > 2 ? contiOxygenEqIdx + 1 : -1000;

    //! Index of the continuity equation for the fourth MICP component
    static constexpr int contiBiofilmEqIdx =
        numMICPs > 3 ? contiUreaEqIdx + 1 : -1000;

    //! Index of the continuity equation for the fifth MICP component
    static constexpr int contiCalciteEqIdx =
        numMICPs > 4 ? contiBiofilmEqIdx + 1 : -1000;

    //! Index of the continuity equation for the foam component
    static constexpr int contiFoamEqIdx =
        enableFoam ? PVOffset + numPhases + numSolvents + numPolymers + numMICPs : -1000;

    //! Index of the continuity equation for the salt component
    static constexpr int contiBrineEqIdx =
        enableBrine ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers + numMICPs + numFoam : -1000;

    //! Index of the continuity equation for energy
    static constexpr int contiEnergyEqIdx =
        enableEnergy ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers + numMICPs + numFoam + numBrine: -1000;
};

} // namespace Opm

#endif
