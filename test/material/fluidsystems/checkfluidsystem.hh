// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \brief This file provides the actual code for the fluid systems
 *        test.
 *
 * It is not directly in test_fluidsystems.cc so that external modules
 * like dumux-devel can use it easily
 */
#ifndef DUMUX_CHECK_FLUIDSYSTEM_HH
#define DUMUX_CHECK_FLUIDSYSTEM_HH

// include all fluid systems in dumux-stable
#include <dumux/material/fluidsystems/1pfluidsystem.hh>
#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/material/fluidsystems/h2oairmesitylenefluidsystem.hh>
#include <dumux/material/fluidsystems/h2oairxylenefluidsystem.hh>
#include <dumux/material/fluidsystems/spe5fluidsystem.hh>

// include all fluid states
#include <dumux/material/fluidstates/pressureoverlayfluidstate.hh>
#include <dumux/material/fluidstates/saturationoverlayfluidstate.hh>
#include <dumux/material/fluidstates/temperatureoverlayfluidstate.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/material/fluidstates/nonequilibriumfluidstate.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

#include <dune/common/classname.hh>

#include <iostream>
#include <string>

// this is a fluid state which makes sure that only the quantities
// allowed are accessed
template <class Scalar,
          class FluidSystem,
          class BaseFluidState = Dumux::CompositionalFluidState<Scalar, FluidSystem> >
class HairSplittingFluidState
    : protected BaseFluidState
{
public:
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    HairSplittingFluidState()
    {
        // initially, do not allow anything
        allowTemperature(false);
        allowPressure(false);
        allowComposition(false);
        allowDensity(false);

        // do not allow accessing any phase
        restrictToPhase(1000);
    }

    void allowTemperature(bool yesno)
    { allowTemperature_ = yesno; }

    void allowPressure(bool yesno)
    { allowPressure_ = yesno; }

    void allowComposition(bool yesno)
    { allowComposition_ = yesno; }

    void allowDensity(bool yesno)
    { allowDensity_ = yesno; }

    void restrictToPhase(int phaseIdx)
    { restrictPhaseIdx_ = phaseIdx; }

    Scalar temperature(int phaseIdx) const
    {
        assert(allowTemperature_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        DUNE_UNUSED Scalar tmp = BaseFluidState::temperature(phaseIdx);
        return 1e100;
    }

    Scalar pressure(int phaseIdx) const
    {
        assert(allowPressure_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        DUNE_UNUSED Scalar tmp = BaseFluidState::pressure(phaseIdx);
        return 1e100;
    }

    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        DUNE_UNUSED Scalar tmp = BaseFluidState::moleFraction(phaseIdx, compIdx);
        return 1e100;
    }

    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        DUNE_UNUSED Scalar tmp = BaseFluidState::massFraction(phaseIdx, compIdx);
        return 1e100;
    }

    Scalar averageMolarMass(int phaseIdx) const
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        DUNE_UNUSED Scalar tmp = BaseFluidState::averageMolarMass(phaseIdx);
        return 1e100;
    }

    Scalar molarity(int phaseIdx, int compIdx) const
    {
        assert(allowDensity_ && allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        DUNE_UNUSED Scalar tmp = BaseFluidState::molarity(phaseIdx, compIdx);
        return 1e100;
    }

    Scalar molarDensity(int phaseIdx) const
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        DUNE_UNUSED Scalar tmp = BaseFluidState::molarDensity(phaseIdx);
        return 1e100;
    }

    Scalar molarVolume(int phaseIdx) const
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        DUNE_UNUSED Scalar tmp = BaseFluidState::molarVolume(phaseIdx);
        return 1e100;
    }

    Scalar density(int phaseIdx) const
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == phaseIdx);
        DUNE_UNUSED Scalar tmp = BaseFluidState::density(phaseIdx);
        return 1e100;
    }

    Scalar saturation(int phaseIdx) const
    {
        assert(false);
        DUNE_UNUSED Scalar tmp =  BaseFluidState::saturation(phaseIdx);
        return 1e100;
    }

    Scalar fugacity(int phaseIdx, int compIdx) const
    {
        assert(false);
        DUNE_UNUSED Scalar tmp = BaseFluidState::fugacity(phaseIdx, compIdx);
        return 1e100;
    }

    Scalar fugacityCoefficient(int phaseIdx, int compIdx) const
    {
        assert(false);
        DUNE_UNUSED Scalar tmp = BaseFluidState::fugacityCoefficient(phaseIdx, compIdx);
        return 1e100;
    }

    Scalar enthalpy(int phaseIdx) const
    {
        assert(false);
        DUNE_UNUSED Scalar tmp = BaseFluidState::enthalpy(phaseIdx);
        return 1e100;
    }

    Scalar internalEnergy(int phaseIdx) const
    {
        assert(false);
        DUNE_UNUSED Scalar tmp = BaseFluidState::internalEnergy(phaseIdx);
        return 1e100;
    }

    Scalar viscosity(int phaseIdx) const
    {
        assert(false);
        DUNE_UNUSED Scalar tmp = BaseFluidState::viscosity(phaseIdx);
        return 1e100;
    }

private:
    bool allowSaturation_;
    bool allowTemperature_;
    bool allowPressure_;
    bool allowComposition_;
    bool allowDensity_;
    int restrictPhaseIdx_;
};

template <class Scalar, class BaseFluidState>
void checkFluidState(const BaseFluidState &fs)
{
    // fluid states must be copy-able
    BaseFluidState tmpFs(fs);
    tmpFs = fs;

    // a fluid state must provide a checkDefined() method
    fs.checkDefined();

    // make sure the fluid state provides all mandatory methods
    while (false) {
        Scalar DUNE_UNUSED val;

        val = fs.temperature(/*phaseIdx=*/0);
        val = fs.pressure(/*phaseIdx=*/0);
        val = fs.moleFraction(/*phaseIdx=*/0, /*compIdx=*/0);
        val = fs.massFraction(/*phaseIdx=*/0, /*compIdx=*/0);
        val = fs.averageMolarMass(/*phaseIdx=*/0);
        val = fs.molarity(/*phaseIdx=*/0, /*compIdx=*/0);
        val = fs.molarDensity(/*phaseIdx=*/0);
        val = fs.molarVolume(/*phaseIdx=*/0);
        val = fs.density(/*phaseIdx=*/0);
        val = fs.saturation(/*phaseIdx=*/0);
        val = fs.fugacity(/*phaseIdx=*/0, /*compIdx=*/0);
        val = fs.fugacityCoefficient(/*phaseIdx=*/0, /*compIdx=*/0);
        val = fs.enthalpy(/*phaseIdx=*/0);
        val = fs.internalEnergy(/*phaseIdx=*/0);
        val = fs.viscosity(/*phaseIdx=*/0);
    };
}

template <class Scalar, class FluidSystem>
void checkFluidSystem()
{
    std::cout << "Testing fluid system '" << Dune::className<FluidSystem>() << "'\n";

    // make sure the fluid system provides the number of phases and
    // the number of components
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    HairSplittingFluidState<Scalar, FluidSystem> fs;
    fs.allowTemperature(true);
    fs.allowPressure(true);
    fs.allowComposition(true);
    fs.restrictToPhase(-1);

    // check whether the parameter cache adheres to the API
    typedef typename FluidSystem::ParameterCache PC;
    PC paramCache;
    try { paramCache.updateAll(fs); } catch (...) {};
    try { paramCache.updateAll(fs, /*except=*/PC::None); } catch (...) {};
    try { paramCache.updateAll(fs, /*except=*/PC::Temperature | PC::Pressure | PC::Composition); } catch (...) {};
    try { paramCache.updateAllPressures(fs); } catch (...) {};

    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        fs.restrictToPhase(phaseIdx);
        try { paramCache.updatePhase(fs, phaseIdx); } catch (...) {};
        try { paramCache.updatePhase(fs, phaseIdx, /*except=*/PC::None); } catch (...) {};
        try { paramCache.updatePhase(fs, phaseIdx, /*except=*/PC::Temperature | PC::Pressure | PC::Composition); } catch (...) {};
        try { paramCache.updateTemperature(fs, phaseIdx); } catch (...) {};
        try { paramCache.updatePressure(fs, phaseIdx); } catch (...) {};
        try { paramCache.updateComposition(fs, phaseIdx); } catch (...) {};
        try { paramCache.updateSingleMoleFraction(fs, phaseIdx, /*compIdx=*/0); } catch (...) {};
    }

    // some value to make sure the return values of the fluid system
    // are convertible to scalars
    Scalar DUNE_UNUSED val;

    // actually check the fluid system API
    try { FluidSystem::init(); } catch (...) {};
    for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        fs.restrictToPhase(phaseIdx);
        fs.allowPressure(FluidSystem::isCompressible(phaseIdx));
        fs.allowComposition(true);
        fs.allowDensity(false);
        try { val = FluidSystem::density(fs, paramCache, phaseIdx); } catch (...) {};
            
        fs.allowPressure(true);
        fs.allowDensity(true);
        try { val = FluidSystem::viscosity(fs, paramCache, phaseIdx); } catch (...) {};
        try { val = FluidSystem::enthalpy(fs, paramCache, phaseIdx); } catch (...) {};
        try { val = FluidSystem::heatCapacity(fs, paramCache, phaseIdx); } catch (...) {};
        try { val = FluidSystem::thermalConductivity(fs, paramCache, phaseIdx); } catch (...) {};
            
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            fs.allowComposition(!FluidSystem::isIdealMixture(phaseIdx));
            try { val = FluidSystem::fugacityCoefficient(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
            fs.allowComposition(true);
            try { val = FluidSystem::diffusionCoefficient(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
        }
            
    }

    // test for phaseName(), isLiquid() and isIdealGas()
    for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        std::string DUNE_UNUSED name = FluidSystem::phaseName(phaseIdx);
        bool DUNE_UNUSED bVal = FluidSystem::isLiquid(phaseIdx);
        bVal = FluidSystem::isIdealGas(phaseIdx);
    }

    // test for componentName()
    for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
        val = FluidSystem::molarMass(compIdx);
        std::string DUNE_UNUSED name = FluidSystem::componentName(compIdx);
    }

    std::cout << "----------------------------------\n";
}

#endif
