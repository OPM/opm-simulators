// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
*/
/*!
 * \file
 * \copydoc checkFluidSystem
 */
#ifndef OPM_CHECK_FLUIDSYSTEM_HPP
#define OPM_CHECK_FLUIDSYSTEM_HPP

// include all fluid systems in opm-material
#include <opm/material/fluidsystems/SinglePhaseFluidSystem.hpp>
#include <opm/material/fluidsystems/TwoPhaseImmiscibleFluidSystem.hpp>
#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidsystems/H2ON2LiquidPhaseFluidSystem.hpp>
#include <opm/material/fluidsystems/H2OAirFluidSystem.hpp>
#include <opm/material/fluidsystems/H2OAirMesityleneFluidSystem.hpp>
#include <opm/material/fluidsystems/H2OAirXyleneFluidSystem.hpp>
#include <opm/material/fluidsystems/Spe5FluidSystem.hpp>

// include all fluid states
#include <opm/material/fluidstates/PressureOverlayFluidState.hpp>
#include <opm/material/fluidstates/SaturationOverlayFluidState.hpp>
#include <opm/material/fluidstates/TemperatureOverlayFluidState.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/NonEquilibriumFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>

#include <opm/material/common/Unused.hpp>
#include <opm/material/common/ClassName.hpp>

#include <iostream>
#include <string>

/*!
 * \brief This is a fluid state which makes sure that only the quantities
 *        allowed are accessed.
 */
template <class ScalarT,
          class FluidSystem,
          class BaseFluidState = Opm::CompositionalFluidState<ScalarT, FluidSystem> >
class HairSplittingFluidState
    : protected BaseFluidState
{
public:
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    typedef ScalarT Scalar;

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

    Scalar temperature(unsigned phaseIdx) const
    {
        assert(allowTemperature_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        OPM_UNUSED Scalar tmp = BaseFluidState::temperature(phaseIdx);
        return 1e100;
    }

    Scalar pressure(unsigned phaseIdx) const
    {
        assert(allowPressure_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        OPM_UNUSED Scalar tmp = BaseFluidState::pressure(phaseIdx);
        return 1e100;
    }

    Scalar moleFraction(unsigned phaseIdx, unsigned compIdx) const
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        OPM_UNUSED Scalar tmp = BaseFluidState::moleFraction(phaseIdx, compIdx);
        return 1e100;
    }

    Scalar massFraction(unsigned phaseIdx, unsigned compIdx) const
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        OPM_UNUSED Scalar tmp = BaseFluidState::massFraction(phaseIdx, compIdx);
        return 1e100;
    }

    Scalar averageMolarMass(unsigned phaseIdx) const
    {
        assert(allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        OPM_UNUSED Scalar tmp = BaseFluidState::averageMolarMass(phaseIdx);
        return 1e100;
    }

    Scalar molarity(unsigned phaseIdx, unsigned compIdx) const
    {
        assert(allowDensity_ && allowComposition_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        OPM_UNUSED Scalar tmp = BaseFluidState::molarity(phaseIdx, compIdx);
        return 1e100;
    }

    Scalar molarDensity(unsigned phaseIdx) const
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        OPM_UNUSED Scalar tmp = BaseFluidState::molarDensity(phaseIdx);
        return 1e100;
    }

    Scalar molarVolume(unsigned phaseIdx) const
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        OPM_UNUSED Scalar tmp = BaseFluidState::molarVolume(phaseIdx);
        return 1e100;
    }

    Scalar density(unsigned phaseIdx) const
    {
        assert(allowDensity_);
        assert(restrictPhaseIdx_ < 0 || restrictPhaseIdx_ == static_cast<int>(phaseIdx));
        OPM_UNUSED Scalar tmp = BaseFluidState::density(phaseIdx);
        return 1e100;
    }

    Scalar saturation(unsigned phaseIdx) const
    {
        assert(false);
        OPM_UNUSED Scalar tmp =  BaseFluidState::saturation(phaseIdx);
        return 1e100;
    }

    Scalar fugacity(unsigned phaseIdx, unsigned compIdx) const
    {
        assert(false);
        OPM_UNUSED Scalar tmp = BaseFluidState::fugacity(phaseIdx, compIdx);
        return 1e100;
    }

    Scalar fugacityCoefficient(unsigned phaseIdx, unsigned compIdx) const
    {
        assert(false);
        OPM_UNUSED Scalar tmp = BaseFluidState::fugacityCoefficient(phaseIdx, compIdx);
        return 1e100;
    }

    Scalar enthalpy(unsigned phaseIdx) const
    {
        assert(false);
        OPM_UNUSED Scalar tmp = BaseFluidState::enthalpy(phaseIdx);
        return 1e100;
    }

    Scalar internalEnergy(unsigned phaseIdx) const
    {
        assert(false);
        OPM_UNUSED Scalar tmp = BaseFluidState::internalEnergy(phaseIdx);
        return 1e100;
    }

    Scalar viscosity(unsigned phaseIdx) const
    {
        assert(false);
        OPM_UNUSED Scalar tmp = BaseFluidState::viscosity(phaseIdx);
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

    // fluid states must export the types which they use as Scalars
    typedef typename BaseFluidState::Scalar FsScalar;
    static_assert(std::is_same<FsScalar, Scalar>::value,
                  "Fluid states must export the type they are given as scalar in an unmodified way");

    // make sure the fluid state provides all mandatory methods
    while (false) {
        Scalar val = 1.0;

        val = 2*val; // get rid of GCC warning (only occurs with paranoid warning flags)

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

/*!
 * \brief Checks whether a fluid system adheres to the specification.
 */
template <class Scalar, class FluidSystem, class RhsEval, class LhsEval>
void checkFluidSystem()
{
    std::cout << "Testing fluid system '" << Opm::className<FluidSystem>() << "'\n";

    // make sure the fluid system provides the number of phases and
    // the number of components
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    typedef HairSplittingFluidState<RhsEval, FluidSystem> FluidState;
    FluidState fs;
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

    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        fs.restrictToPhase(static_cast<int>(phaseIdx));
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
    LhsEval val = 0.0;
    Scalar scalarVal = 0.0;

    scalarVal = 2*scalarVal; // get rid of GCC warning (only occurs with paranoid warning flags)
    val = 2*val; // get rid of GCC warning (only occurs with paranoid warning flags)

    // actually check the fluid system API
    try { FluidSystem::init(); } catch (...) {};
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        fs.restrictToPhase(static_cast<int>(phaseIdx));
        fs.allowPressure(FluidSystem::isCompressible(phaseIdx));
        fs.allowComposition(true);
        fs.allowDensity(false);
        try { auto tmpVal = FluidSystem::density(fs, paramCache, phaseIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
        try { val = FluidSystem::template density<FluidState, LhsEval>(fs, paramCache, phaseIdx); } catch (...) {};
        try { scalarVal = FluidSystem::template density<FluidState, Scalar>(fs, paramCache, phaseIdx); } catch (...) {};

        fs.allowPressure(true);
        fs.allowDensity(true);
        try { auto tmpVal = FluidSystem::viscosity(fs, paramCache, phaseIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
        try { auto tmpVal = FluidSystem::enthalpy(fs, paramCache, phaseIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
        try { auto tmpVal = FluidSystem::heatCapacity(fs, paramCache, phaseIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
        try { auto tmpVal = FluidSystem::thermalConductivity(fs, paramCache, phaseIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
        try { val = FluidSystem::template viscosity<FluidState, LhsEval>(fs, paramCache, phaseIdx); } catch (...) {};
        try { val = FluidSystem::template enthalpy<FluidState, LhsEval>(fs, paramCache, phaseIdx); } catch (...) {};
        try { val = FluidSystem::template heatCapacity<FluidState, LhsEval>(fs, paramCache, phaseIdx); } catch (...) {};
        try { val = FluidSystem::template thermalConductivity<FluidState, LhsEval>(fs, paramCache, phaseIdx); } catch (...) {};
        try { scalarVal = FluidSystem::template viscosity<FluidState, Scalar>(fs, paramCache, phaseIdx); } catch (...) {};
        try { scalarVal = FluidSystem::template enthalpy<FluidState, Scalar>(fs, paramCache, phaseIdx); } catch (...) {};
        try { scalarVal = FluidSystem::template heatCapacity<FluidState, Scalar>(fs, paramCache, phaseIdx); } catch (...) {};
        try { scalarVal = FluidSystem::template thermalConductivity<FluidState, Scalar>(fs, paramCache, phaseIdx); } catch (...) {};

        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
            fs.allowComposition(!FluidSystem::isIdealMixture(phaseIdx));
            try { auto tmpVal = FluidSystem::fugacityCoefficient(fs, paramCache, phaseIdx, compIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
            try { val = FluidSystem::template fugacityCoefficient<FluidState, LhsEval>(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
            try { scalarVal = FluidSystem::template fugacityCoefficient<FluidState, Scalar>(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
            fs.allowComposition(true);
            try { auto tmpVal = FluidSystem::diffusionCoefficient(fs, paramCache, phaseIdx, compIdx); static_assert(std::is_same<decltype(tmpVal), RhsEval>::value, "The default return value must be the scalar used by the fluid state!"); } catch (...) {};
            try { val = FluidSystem::template diffusionCoefficient<FluidState, LhsEval>(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
            try { scalarVal = FluidSystem::template fugacityCoefficient<FluidState, Scalar>(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
        }
    }

    // test for phaseName(), isLiquid() and isIdealGas()
    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        std::string OPM_UNUSED name = FluidSystem::phaseName(phaseIdx);
        bool bVal = FluidSystem::isLiquid(phaseIdx);
        bVal = FluidSystem::isIdealGas(phaseIdx);
        bVal = !bVal; // get rid of GCC warning (only occurs with paranoid warning flags)
    }

    // test for molarMass() and componentName()
    for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx) {
        val = FluidSystem::molarMass(compIdx);
        std::string name = FluidSystem::componentName(compIdx);
    }

    std::cout << "----------------------------------\n";
}

#endif
