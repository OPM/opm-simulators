// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026, SINTEF Digital

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
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
 * \copydoc Opm::FlashLocalResidualTPFA
 */
#ifndef OPM_PTFLASH_LOCAL_RESIDUAL_TPFA_HH
#define OPM_PTFLASH_LOCAL_RESIDUAL_TPFA_HH

#include <opm/models/ptflash/flashlocalresidual.hh>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/Schedule/BCState.hpp>

#include <opm/material/common/ConditionalStorage.hpp>
#include <opm/material/common/MathToolbox.hpp>

#include <stdexcept>

namespace Opm {

/*!
 * \ingroup FlashModel
 *
 * \brief Local residual of the PT-flash model for the TPFA linearizer.
 *
 * The TPFA linearizer evaluates the storage, face flux and source terms one
 * cell or face at a time, without an element context. The face flux here is
 * the same transmissibility-based flux that the element-context assembly
 * computes. The element-context versions are still inherited from
 * FlashLocalResidual.
 */
template <class TypeTag>
class FlashLocalResidualTPFA : public FlashLocalResidual<TypeTag>
{
    using ParentType = FlashLocalResidual<TypeTag>;

    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr unsigned numPhases = getPropValue<TypeTag, Properties::NumPhases>();
    static constexpr unsigned numComponents = getPropValue<TypeTag, Properties::NumComponents>();
    static constexpr unsigned conti0EqIdx = Indices::conti0EqIdx;
    static constexpr bool waterEnabled = Indices::waterEnabled;

    static_assert(!getPropValue<TypeTag, Properties::EnableEnergy>(),
                  "The TPFA flash residual does not support the energy equation");
    static_assert(!getPropValue<TypeTag, Properties::EnableDiffusion>(),
                  "The TPFA flash residual does not support molecular diffusion");

    using Toolbox = MathToolbox<Evaluation>;

public:
    //! Face data. The TPFA linearizer fills it by aggregate initialization,
    //! so the members must keep this order. The last five, for heat
    //! conduction, diffusion, dispersion and polymer, are unused.
    struct ResidualNBInfo
    {
        double trans;
        double faceArea;
        double thpresInToEx;
        double thpresExToIn;
        double dZg;
        FaceDir::DirEnum faceDir;
        double Vin;
        double Vex;
        ConditionalStorage<false, double> inAlpha;
        ConditionalStorage<false, double> outAlpha;
        ConditionalStorage<false, double> diffusivity;
        ConditionalStorage<false, double> dispersivity;
        ConditionalStorage<false, double> dist;
    };

    using ParentType::computeFlux;
    using ParentType::computeSource;

    /*!
     * \brief Mass flux of each component from the inside to the outside cell,
     *        per unit face area.
     *
     * Derivatives are taken with respect to the inside cell only. \p darcy
     * receives the volume rate of each phase, indexed by phase.
     */
    template <class ModuleParamsT, class RateVectorT, class ResidualNBInfoT>
    static void computeFlux(RateVectorT& flux,
                            RateVectorT& darcy,
                            const unsigned globalIndexIn,
                            const unsigned globalIndexEx,
                            const IntensiveQuantities& intQuantsIn,
                            const IntensiveQuantities& intQuantsEx,
                            const ResidualNBInfoT& nbInfo,
                            const ModuleParamsT& moduleParams)
    {
        flux = 0.0;
        darcy = 0.0;

        constexpr unsigned interiorIdx = 0;
        constexpr unsigned exteriorIdx = 1;

        // Average the rock compaction multipliers of the two cells, as the
        // element-context flux does.
        const Evaluation transMult = (intQuantsIn.rockCompTransMultiplier() +
                                      Toolbox::value(intQuantsEx.rockCompTransMultiplier())) / 2;
        const Scalar transFactor = -nbInfo.trans / nbInfo.faceArea;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            short upIdx = interiorIdx;
            short dnIdx = exteriorIdx;
            Evaluation pressureDifference;
            ExtensiveQuantities::calculatePhasePressureDiff_(upIdx,
                                                             dnIdx,
                                                             pressureDifference,
                                                             intQuantsIn,
                                                             intQuantsEx,
                                                             phaseIdx,
                                                             interiorIdx,
                                                             exteriorIdx,
                                                             nbInfo.Vin,
                                                             nbInfo.Vex,
                                                             globalIndexIn,
                                                             globalIndexEx,
                                                             nbInfo.dZg,
                                                             nbInfo.thpresInToEx,
                                                             nbInfo.thpresExToIn,
                                                             moduleParams);

            const bool upIsIn = static_cast<unsigned>(upIdx) == interiorIdx;
            const IntensiveQuantities& up = upIsIn ? intQuantsIn : intQuantsEx;
            const auto& fs = up.fluidState();

            const Evaluation volumeFlux = upIsIn
                ? pressureDifference * up.mobility(phaseIdx, nbInfo.faceDir) * transMult * transFactor
                : pressureDifference * (Toolbox::value(up.mobility(phaseIdx, nbInfo.faceDir)) *
                                        transMult * transFactor);
            darcy[phaseIdx] = volumeFlux.value() * nbInfo.faceArea;

            const Evaluation massFlux = upIsIn
                ? fs.density(phaseIdx) * volumeFlux
                : Toolbox::value(fs.density(phaseIdx)) * volumeFlux;

            if (waterEnabled && phaseIdx == static_cast<unsigned>(FluidSystem::waterPhaseIdx)) {
                flux[conti0EqIdx + numComponents] = massFlux;
                continue;
            }

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                flux[conti0EqIdx + compIdx] += upIsIn
                    ? massFlux * fs.massFraction(phaseIdx, compIdx)
                    : massFlux * Toolbox::value(fs.massFraction(phaseIdx, compIdx));
            }
        }
    }

    //! Source term of a single cell, per unit bulk volume.
    static void computeSource(RateVector& source,
                              const Problem& problem,
                              const IntensiveQuantities& /*intQuants*/,
                              const unsigned globalSpaceIdx,
                              const unsigned timeIdx)
    {
        problem.source(source, globalSpaceIdx, timeIdx);
    }

    //! Source term of a single cell, per unit bulk volume, without the wells.
    static void computeSourceDense(RateVector& source,
                                   const Problem& problem,
                                   const IntensiveQuantities& /*intQuants*/,
                                   const unsigned globalSpaceIdx,
                                   const unsigned timeIdx)
    {
        source = 0.0;
        problem.addToSourceDense(source, globalSpaceIdx, timeIdx);
    }

    //! Only no-flow boundaries are supported.
    template <class RateVectorT, class BoundaryConditionData>
    static void computeBoundaryFlux(RateVectorT& bdyFlux,
                                    const Problem& /*problem*/,
                                    const BoundaryConditionData& bdyInfo,
                                    const IntensiveQuantities& /*insideIntQuants*/,
                                    const unsigned /*globalSpaceIdx*/)
    {
        if (bdyInfo.type != BCType::NONE) {
            throw std::logic_error("Boundary conditions are not supported "
                                   "by the compositional model yet");
        }
        bdyFlux = 0.0;
    }
};

} // namespace Opm

#endif // OPM_PTFLASH_LOCAL_RESIDUAL_TPFA_HH
