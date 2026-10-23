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
 * \copydoc Opm::BlackOilDynTwoPhaseIndices
 */
#ifndef OPM_BLACK_OIL_DYN_TWO_PHASE_INDICES_HH
#define OPM_BLACK_OIL_DYN_TWO_PHASE_INDICES_HH

#include <opm/common/utility/ConstexprAssert.hpp>

#include <opm/models/blackoil/blackoildynamicphases.hh>

#include <cassert>
#include <limits>

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief The primary variable and equation indices for the two-phase black-oil model
 *        where the pair of active phases is chosen at run time.
 *
 * Same layout rules as BlackOilTwoPhaseIndices, but oil-water, gas-water and gas-oil
 * share one type. Which phases are active is read from BlackOilDynamicPhases, so
 * the indices that depend on it are functions and not constants. The numbers of
 * equations and variables do not depend on the phases and stay compile-time.
 */
template<unsigned numSolventsV,
         unsigned numExtbosV,
         unsigned numPolymersV,
         unsigned numEnergyV,
         bool enableFoam,
         bool enableBrine,
         int PVOffset,
         unsigned numBioCompV>
struct BlackOilDynTwoPhaseIndices
{
    //! Is phase enabled or not (decided at run time)
    static bool oilIsActive() { return BlackOilDynamicPhases::oilActive(); }
    static bool waterIsActive() { return BlackOilDynamicPhases::waterActive(); }
    static bool gasIsActive() { return BlackOilDynamicPhases::gasActive(); }

    //! Are solvents involved?
    static constexpr bool enableSolvent = numSolventsV > 0;

    //! Is extbo invoked?
    static constexpr bool enableExtbo = numExtbosV > 0;

    //! Are polymers involved?
    static constexpr bool enablePolymer = numPolymersV > 0;

    //! Shall energy be conserved?
    static constexpr bool enableFullyImplicitThermal = numEnergyV > 0;

    //! MICP only available for one phase indices
    static constexpr bool enableMICP = false;

    //! Are biofilms involved?
    static constexpr bool enableBiofilm = numBioCompV > 0;

    //! Number of solvent components to be considered
    static constexpr int numSolvents = enableSolvent ? numSolventsV : 0;

    //! Number of components to be considered for extbo
    static constexpr int numExtbos = enableExtbo ? numExtbosV : 0;

    //! Number of polymer components to be considered
    static constexpr int numPolymers = enablePolymer ? numPolymersV : 0;

    //! Number of energy equations to be considered
    static constexpr int numEnergy = enableFullyImplicitThermal ? numEnergyV : 0;

    //! Number of foam equations to be considered
    static constexpr int numFoam = enableFoam? 1 : 0;

    //! Number of salt equations to be considered
    static constexpr int numBrine = enableBrine? 1 : 0;

    //! Number of biofilm equations to be considered
    static constexpr int numBioComp = enableBiofilm ? numBioCompV : 0;

    //! Number of biocomponents in the water phase
    static constexpr int numBioInWat = enableBiofilm ? 1 : 0;

    //! The number of fluid phases
    static constexpr int numPhases = 2;

    //! The number of equations
    static constexpr int numEq = numPhases + numSolvents + numExtbos + numPolymers +
                                 numEnergy + numFoam + numBrine + numBioComp;

    //! We may want to evaluate with more derivatives that equations, but not this time
    static constexpr int numDerivatives = numEq;

    //////////////////////////////
    // Primary variable indices
    //////////////////////////////

    //! Index of the primary variable for the first solvent
    static constexpr unsigned solventSaturationIdx =
        enableSolvent ? PVOffset + numPhases
                      : std::numeric_limits<unsigned>::max();

    //! Index of the primary variable for the first extbo component
    static constexpr unsigned zFractionIdx =
        enableExtbo ? PVOffset + numPhases + numSolvents
                    : std::numeric_limits<unsigned>::max();

    //! Index of the primary variable for the first polymer
    static constexpr unsigned polymerConcentrationIdx =
        enablePolymer ? PVOffset + numPhases + numSolvents
                      : std::numeric_limits<unsigned>::max();

    //! Index of the primary variable for the second polymer primary variable (molecular weight)
    static constexpr unsigned polymerMoleWeightIdx =
        numPolymers > 1 ? polymerConcentrationIdx + 1
                        : std::numeric_limits<unsigned>::max();

    //! Index of the primary variable for the first microbial component
    static constexpr unsigned microbialConcentrationIdx =
        enableBiofilm ? PVOffset + numPhases + numSolvents
                      : std::numeric_limits<unsigned>::max();

    //! Index of the primary variable for the biofilm component
    static constexpr unsigned biofilmVolumeFractionIdx =
        enableBiofilm ? PVOffset + numPhases + numSolvents + 1
                      : std::numeric_limits<unsigned>::max();

    //! MICP only available for one phase indices
    static constexpr unsigned oxygenConcentrationIdx = std::numeric_limits<unsigned>::max();
    static constexpr unsigned ureaConcentrationIdx = std::numeric_limits<unsigned>::max();
    static constexpr unsigned calciteVolumeFractionIdx = std::numeric_limits<unsigned>::max();

    //! Index of the primary variable for the foam
    static constexpr unsigned foamConcentrationIdx =
        enableFoam ? PVOffset + numPhases + numSolvents + numPolymers + numBioComp
                   : std::numeric_limits<unsigned>::max();

    //! Index of the primary variable for the salt
    static constexpr unsigned saltConcentrationIdx =
        enableBrine ? PVOffset + numPhases + numSolvents + numPolymers + numBioComp + numFoam
                    : std::numeric_limits<unsigned>::max();

    //! Index of the primary variable for temperature
    static constexpr unsigned temperatureIdx  =
        enableFullyImplicitThermal ? PVOffset + numPhases + numSolvents + numExtbos +
                                     numPolymers + numBioComp + numFoam + numBrine
                                   : std::numeric_limits<unsigned>::max();

    //////////////////////
    // Equation indices
    //////////////////////

    //! Index of the continuity equation of the first phase
    static constexpr int conti0EqIdx = PVOffset + 0;
    // one continuity equation follows

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

    //! Index of the continuity equation for the first microbial component
    static constexpr int contiMicrobialEqIdx =
        enableBiofilm ? PVOffset + numPhases + numSolvents : -1000;

    //! Index of the continuity equation for the biofilm component
    static constexpr int contiBiofilmEqIdx =
        enableBiofilm ? PVOffset + numPhases + numSolvents + 1 : -1000;

    //! MICP only available for one phase indices
    static constexpr int contiOxygenEqIdx = -1000;
    static constexpr int contiUreaEqIdx = -1000;
    static constexpr int contiCalciteEqIdx = -1000;

    //! Index of the continuity equation for the foam component
    static constexpr int contiFoamEqIdx =
        enableFoam ? PVOffset + numPhases + numSolvents + numPolymers + numBioComp : -1000;

    //! Index of the continuity equation for the salt component
    static constexpr int contiBrineEqIdx =
        enableBrine ? PVOffset + numPhases + numSolvents + numPolymers + numBioComp + numFoam : -1000;

    //! Index of the continuity equation for energy
    static constexpr int contiEnergyEqIdx =
        enableFullyImplicitThermal ? PVOffset + numPhases + numSolvents + numExtbos + numPolymers + numBioComp + numFoam + numBrine: -1000;

    //! Primary variable indices. The position of the water saturation, the
    //! pressure and the hydrocarbon composition switch depend on the active phases.
    //! Indices that are not in use are std::numeric_limits<unsigned>::max().
    static unsigned activeWaterSwitchIdx()
    {
        return waterIsActive() ? PVOffset + 0 : std::numeric_limits<unsigned>::max();
    }

    static unsigned activePressureSwitchIdx()
    {
        return waterIsActive() ? PVOffset + 1 : PVOffset + 0;
    }

    static unsigned activeCompositionSwitchIdx()
    {
        return gasIsActive() && oilIsActive() ? PVOffset + 1
                                              : std::numeric_limits<unsigned>::max();
    }

    static constexpr unsigned activeSolventSaturationIdx() { return solventSaturationIdx; }

    //! A composition switch can exist (gas-oil), so the fluid state must store Rs/Rv.
    static constexpr bool hasCompositionSwitch = true;

    //! The gas phase can be active.
    static constexpr bool mayHaveGas = true;
};

} // namespace Opm

#endif
