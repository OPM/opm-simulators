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
 * \copydoc Opm::FlashIntensiveQuantities
 */
#ifndef EWOMS_FLASH_INTENSIVE_QUANTITIES_HH
#define EWOMS_FLASH_INTENSIVE_QUANTITIES_HH

#include "flashproperties.hh"
#include "flashindices.hh"

#include <opm/models/common/energymodule.hh>
#include <opm/models/common/diffusionmodule.hh>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Opm {

/*!
 * \ingroup FlashModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the intensive quantities of the flash-based compositional multi-phase model
 */
template <class TypeTag>
class FlashIntensiveQuantities
    : public GetPropType<TypeTag, Properties::DiscIntensiveQuantities>
    , public DiffusionIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDiffusion>() >
    , public EnergyIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableEnergy>() >
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxIntensiveQuantities
{
    using ParentType = GetPropType<TypeTag, Properties::DiscIntensiveQuantities>;

    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;

    // primary variable indices
    enum { z0Idx = Indices::z0Idx };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { dimWorld = GridView::dimensionworld };
    enum { pressure0Idx = Indices::pressure0Idx };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FlashSolver = GetPropType<TypeTag, Properties::FlashSolver>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using ComponentVector = Dune::FieldVector<Evaluation, numComponents>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;
    using DiffusionIntensiveQuantities = Opm::DiffusionIntensiveQuantities<TypeTag, enableDiffusion>;
    using EnergyIntensiveQuantities = Opm::EnergyIntensiveQuantities<TypeTag, enableEnergy>;

public:
    //! The type of the object returned by the fluidState() method
    using FluidState = Opm::CompositionalFluidState<Evaluation, FluidSystem, enableEnergy>;

    FlashIntensiveQuantities()
    { }

    FlashIntensiveQuantities(const FlashIntensiveQuantities& other) = default;

    FlashIntensiveQuantities& operator=(const FlashIntensiveQuantities& other) = default;

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        EnergyIntensiveQuantities::updateTemperatures_(fluidState_, elemCtx, dofIdx, timeIdx);

        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto& problem = elemCtx.problem();
        
        Scalar flashTolerance = 1.e-12;//EWOMS_GET_PARAM(TypeTag, Scalar, FlashTolerance);
        int flashVerbosity = 0;//EWOMS_GET_PARAM(TypeTag, int, FlashVerbosity);
        std::string flashTwoPhaseMethod = EWOMS_GET_PARAM(TypeTag, std::string, FlashTwoPhaseMethod);
        
        // extract the total molar densities of the components
        ComponentVector z(0.);
        {
            Evaluation lastZ = 1.0;
            for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
                z[compIdx] = priVars.makeEvaluation(z0Idx + compIdx, timeIdx);
                lastZ -= z[compIdx];
            }
            z[numComponents - 1] = lastZ;

            Evaluation sumz = 0.0;
            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                z[compIdx] = Opm::max(z[compIdx], 1e-8);
                sumz +=z[compIdx];
            }
            z /= sumz;

        }

        Evaluation p = priVars.makeEvaluation(pressure0Idx, timeIdx);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fluidState_.setPressure(phaseIdx, p);

        // Get initial K and L from storage initially (if enabled)
        const auto *hint = elemCtx.thermodynamicHint(dofIdx, timeIdx);
        const auto *hint2 = elemCtx.thermodynamicHint(dofIdx, 1);
        if (hint) {
             for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                 const Evaluation& Ktmp = hint->fluidState().K(compIdx);
                 fluidState_.setKvalue(compIdx, Ktmp);
             }
             const Evaluation& Ltmp = hint->fluidState().L();
             fluidState_.setLvalue(Ltmp);
        }
        else if (hint2) {
             for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                 const Evaluation& Ktmp = hint2->fluidState().K(compIdx);
                 fluidState_.setKvalue(compIdx, Ktmp);
             }
             const Evaluation& Ltmp = hint2->fluidState().L();
             fluidState_.setLvalue(Ltmp);
        }
        else {
             for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                 const Evaluation Ktmp = fluidState_.wilsonK_(compIdx);
                 fluidState_.setKvalue(compIdx, Ktmp);
             }
             const Evaluation& Ltmp = -1.0;
             fluidState_.setLvalue(Ltmp);
         }

        /////////////
        // Compute the phase compositions and densities 
        /////////////
        int spatialIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        //FlashSolver::solve(fluidState_, z, spatialIdx, flashVerbosity, flashTwoPhaseMethod, flashTolerance);
        //Flash::solve(fluidState_, z, spatialIdx, flashVerbosity, flashTwoPhaseMethod, flashTolerance);

        
        using Flash = Opm::PTFlash<double, FluidSystem>;
        FlashSolver::solve(fluidState_, z, spatialIdx, flashTwoPhaseMethod, flashTolerance, flashVerbosity);

//printing of flashresult after solve
    // std::cout << " \n After flashsolve for cell " << spatialIdx << std::endl;
    // ComponentVector x, y;
    // Evaluation L0 = fluidState_.L();
    // for (unsigned comp_idx = 0; comp_idx < numComponents; ++comp_idx) {
    //     x[comp_idx] = fluidState_.moleFraction(FluidSystem::oilPhaseIdx, comp_idx);
    //     y[comp_idx] = fluidState_.moleFraction(FluidSystem::gasPhaseIdx, comp_idx);
    // }
    //         for (unsigned comp_idx = 0; comp_idx < numComponents; ++comp_idx) {
    //     std::cout << " x for component: " << comp_idx << " is " << x[comp_idx] << std::endl;
    //      for (int i = 0; i < 3; ++i) {
    //          std::cout << " x deriv " << i << " is: " << x[comp_idx].derivative(i) << std::endl;
    //      }

    //     std::cout << " y for component: " << comp_idx << "is " << y[comp_idx] << std::endl;
    //      for (int i = 0; i < 3; ++i) {
    //          std::cout << " y deriv " << i << " is: " << y[comp_idx].derivative(i) << std::endl;
    //      }
    // }
    // std::cout << " L is " << L0 << std::endl;
    //  for (int i = 0; i < L0.size(); ++i) {
    //          std::cout << " L deriv " << i << " is: " << L0.derivative(i) << std::endl;
    //  }
    //  //end printinting 1

        // Update phases        
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updatePhase(fluidState_, FluidSystem::oilPhaseIdx);
        
        const Scalar R = Opm::Constants<Scalar>::R;
        Evaluation Z_L = (paramCache.molarVolume(FluidSystem::oilPhaseIdx) * fluidState_.pressure(FluidSystem::oilPhaseIdx) )/
        (R * fluidState_.temperature(FluidSystem::oilPhaseIdx));
        paramCache.updatePhase(fluidState_, FluidSystem::gasPhaseIdx);
        Evaluation Z_V = (paramCache.molarVolume(FluidSystem::gasPhaseIdx) * fluidState_.pressure(FluidSystem::gasPhaseIdx) )/
        (R * fluidState_.temperature(FluidSystem::gasPhaseIdx));


        // Update saturation
        Evaluation L = fluidState_.L();
        Evaluation So = Opm::max((L*Z_L/(L*Z_L+(1-L)*Z_V)), 0.0);
        Evaluation Sg = Opm::max(1-So, 0.0);
        Scalar sumS = Opm::getValue(So) + Opm::getValue(Sg);
        So /= sumS;
        Sg /= sumS;
        
        fluidState_.setSaturation(0, So);
        fluidState_.setSaturation(1, Sg);

        fluidState_.setCompressFactor(0, Z_L);
        fluidState_.setCompressFactor(1, Z_V);

        // Print saturation
         if (flashVerbosity == 5) {
             std::cout << "So = " << So <<std::endl;
             std::cout << "Sg = " << Sg <<std::endl;
           }

        // Print saturation
         if (flashVerbosity == 5) {
             std::cout << "So = " << So <<std::endl;
             std::cout << "Sg = " << Sg <<std::endl;
             std::cout << "Z_L = " << Z_L <<std::endl;
             std::cout << "Z_V = " << Z_V <<std::endl;
         }
   
        /////////////
        // Compute rel. perm and viscosities and densities
        /////////////
        const MaterialLawParams& materialParams = problem.materialLawParams(elemCtx, dofIdx, timeIdx);

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_,
                                            materialParams, fluidState_);
        Opm::Valgrind::CheckDefined(relativePermeability_);

        // set the phase viscosities and density
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            paramCache.updatePhase(fluidState_, phaseIdx);

            const Evaluation& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            
            fluidState_.setViscosity(phaseIdx, mu);
            //fluidState_.setViscosity(phaseIdx, newmu);

            mobility_[phaseIdx] = relativePermeability_[phaseIdx] / mu;
            Opm::Valgrind::CheckDefined(mobility_[phaseIdx]);

            const Evaluation& rho = FluidSystem::density(fluidState_, paramCache, phaseIdx);
            fluidState_.setDensity(phaseIdx, rho);
        }

        /////////////
        // Compute the remaining quantities
        /////////////

        // porosity
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        Opm::Valgrind::CheckDefined(porosity_);

        // intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities specific for the velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // energy related quantities
        EnergyIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        DiffusionIntensiveQuantities::update_(fluidState_, paramCache, elemCtx, dofIdx, timeIdx);
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::intrinsicPermeability
     */
    const DimMatrix& intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    const Evaluation& relativePermeability(unsigned phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    {
        return mobility_[phaseIdx];
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

private:
    DimMatrix intrinsicPerm_;
    FluidState fluidState_;
    Evaluation porosity_;
    Evaluation relativePermeability_[numPhases];
    Evaluation mobility_[numPhases];
};

} // namespace Opm

#endif
