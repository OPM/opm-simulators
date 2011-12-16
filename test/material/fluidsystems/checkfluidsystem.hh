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
 * \brief This file provides the actual code for the fluid systems
 *        test. 
 *
 * It is not directly in test_fluidsystems.cc so that external modules
 * like dumux-devel can use it easily
 */
#ifndef DUMUX_CHECK_FLUIDSYSTEM_HH
#define DUMUX_CHECK_FLUIDSYSTEM_HH

// include all fluid systems in dumux-stable
#include <test/boxmodels/1p2c/interstitialfluidtrailfluidsystem.hh>
#include <dumux/material/fluidsystems/1pfluidsystem.hh>
#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <appl/lecture/msm/1p2cvs2p/watercontaminantfluidsystem.hh>

// include the compositional fluid state
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include <dune/common/classname.hh>

// this is a fluid state which makes sure that only the quantities
// allowed are accessed
template <class Scalar, class FluidSystem>
class HairSplittingFluidState 
    : protected Dumux::CompositionalFluidState<Scalar, FluidSystem>
{
    typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> ParentType;

public:
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    HairSplittingFluidState()
    {
        // set some fake values
        ParentType::setTemperature(293.15);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            ParentType::setSaturation(phaseIdx, 1.0 / numPhases);
            ParentType::setDensity(phaseIdx, 1.0);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                ParentType::setMoleFraction(phaseIdx, compIdx, 1.0 / numComponents);
                
            }
        }
        
        // initially, do not allow anything
        allowSaturation(false);
        allowTemperature(false);
        allowPressure(false);
        allowComposition(false);
        allowDensity(false);
    }

    void allowSaturation(bool yesno)
    { allowSaturation_ = yesno; }

    void allowPressure(bool yesno)
    { allowPressure_ = yesno; }

    void allowComposition(bool yesno)
    { allowComposition_ = yesno; }
    
    void allowDensity(bool yesno) 
    { allowDensity_ = yesno; }

    Scalar saturation(int phaseIdx) const
    { assert(allowSaturation_); return ParentType::saturation(phaseIdx); }

    Scalar temperature(int phaseIdx) const
    { assert(allowTemperature_); return ParentType::temperature(phaseIdx); }

    void allowTemperature(bool yesno)
    { allowTemperature_ = yesno; }
    
    Scalar pressure(int phaseIdx) const
    { assert(allowPressure_); return ParentType::pressure(phaseIdx); }

    Scalar moleFraction(int phaseIdx, int compIdx) const
    { assert(allowComposition_); return ParentType::moleFraction(phaseIdx, compIdx); }

    Scalar massFraction(int phaseIdx, int compIdx) const
    { assert(allowComposition_); return ParentType::massFraction(phaseIdx, compIdx); }

    Scalar averageMolarMass(int phaseIdx) const
    { assert(allowComposition_); return ParentType::averageMolarMass(phaseIdx); }

    Scalar molarity(int phaseIdx, int compIdx) const
    { assert(allowDensity_ && allowComposition_); return ParentType::molarity(phaseIdx, compIdx); }

    Scalar molarDensity(int phaseIdx) const
    { assert(allowDensity_); return ParentType::molarDensity(phaseIdx); }

    Scalar molarVolume(int phaseIdx) const
    { assert(allowDensity_); return ParentType::molarVolume(phaseIdx); }

    Scalar density(int phaseIdx) const
    { assert(allowDensity_); return ParentType::density(phaseIdx); }

private:
    bool allowSaturation_;
    bool allowTemperature_;
    bool allowPressure_;
    bool allowComposition_;
    bool allowDensity_;
};

template <class Scalar, class FluidSystem>
void checkFluidSystem()
{
    std::cout << "Testing fluid system '" << Dune::className<FluidSystem>() << "'\n";

    HairSplittingFluidState<Scalar, FluidSystem> fs;

    fs.allowTemperature(true);
    fs.allowPressure(true);
    fs.allowComposition(true);

    // check whether the parameter cache adheres to the API
    typedef typename FluidSystem::ParameterCache PC;
    PC paramCache;
    try { paramCache.updateAll(fs); } catch (...) {};
    try { paramCache.updateAll(fs, /*except=*/PC::None); } catch (...) {};
    try { paramCache.updateAll(fs, /*except=*/PC::Temperature | PC::Pressure | PC::Composition); } catch (...) {};
    try { paramCache.updatePhase(fs, /*phaseIdx=*/0); } catch (...) {};
    try { paramCache.updatePhase(fs, /*phaseIdx=*/0, /*except=*/PC::None); } catch (...) {};
    try { paramCache.updatePhase(fs, /*phaseIdx=*/0, /*except=*/PC::Temperature | PC::Pressure | PC::Composition); } catch (...) {};
    try { paramCache.updateTemperature(fs, /*phaseIdx=*/0); } catch (...) {};
    try { paramCache.updateAllPressures(fs); } catch (...) {};
    try { paramCache.updateSinglePressure(fs, /*phaseIdx=*/0); } catch (...) {};
    try { paramCache.updateComposition(fs, /*phaseIdx=*/0); } catch (...) {};
    try { paramCache.updateSingleMoleFraction(fs, /*phaseIdx=*/0, /*compIdx=*/0); } catch (...) {};

    // some value to make sure the return values of the fluid system
    // are convertible to scalars
    Scalar __attribute__((unused)) val;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    // actually check the fluid system API
    FluidSystem::init();
    for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
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
                
        fs.allowComposition(!FluidSystem::isIdealMixture(phaseIdx));
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            try { val = FluidSystem::fugacityCoefficient(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
            try { val = FluidSystem::diffusionCoefficient(fs, paramCache, phaseIdx, compIdx); } catch (...) {};
            for (int comp2Idx = 0; comp2Idx < numComponents; ++ comp2Idx) {
                try { val = FluidSystem::binaryDiffusionCoefficient(fs, paramCache, phaseIdx, compIdx, comp2Idx); } catch (...) {};
            }
        }
    }

    for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
        std::string __attribute__((unused)) name = FluidSystem::phaseName(phaseIdx);
        bool __attribute__((unused)) bVal = FluidSystem::isLiquid(phaseIdx);
    }
    
    for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
        val = FluidSystem::molarMass(compIdx);
        std::string __attribute__((unused)) name = FluidSystem::componentName(compIdx);
    }

    std::cout << "----------------------------------\n";
};

#endif
