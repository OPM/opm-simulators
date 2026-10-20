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
 * \file Shared per-cell black-oil Newton update physics.
 */
#ifndef OPM_BLACK_OIL_NEWTON_UPDATE_HPP
#define OPM_BLACK_OIL_NEWTON_UPDATE_HPP

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/models/blackoil/blackoilnewtonmethodparams.hpp>
#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/material/common/Valgrind.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace Opm {

// Keep all Newton chopping, bounds and phase switches in one implementation.
// The owning nonlinear method supplies the problem, fluid system and history;
// this operation only writes one candidate cell and returns its switch history.
template<class TypeTag>
struct BlackOilNewtonUpdate
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr bool enableBioeffects = getPropValue<TypeTag, Properties::EnableBioeffects>();
    static constexpr bool enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();

    template<class Problem, class FluidSystem, class PrimaryVariables, class Update>
    OPM_HOST_DEVICE static bool update(const Problem& problem,
                                       const FluidSystem& fluidSystem,
                                       unsigned globalDofIdx,
                                       PrimaryVariables& nextValue,
                                       const PrimaryVariables& currentValue,
                                       const Update& update,
                                       const BlackoilNewtonParams<Scalar>& params,
                                       bool wasSwitched)
    {
        static constexpr bool enableSolvent =
            Indices::solventSaturationIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enableExtbo =
            Indices::zFractionIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enablePolymer =
            Indices::polymerConcentrationIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enablePolymerWeight =
            Indices::polymerMoleWeightIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enableFullyImplicitThermal =
            Indices::temperatureIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enableFoam =
            Indices::foamConcentrationIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enableBrine =
            Indices::saltConcentrationIdx != std::numeric_limits<unsigned>::max();
        static constexpr bool enableMICP = Indices::enableMICP;

        currentValue.checkDefined();
        Valgrind::CheckDefined(update);

        // saturation delta for each phase
        Scalar deltaSw = 0.0;
        Scalar deltaSo = 0.0;
        Scalar deltaSg = 0.0;
        Scalar deltaSs = 0.0;

        if (currentValue.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw)
        {
            if constexpr (Indices::waterSwitchIdx != std::numeric_limits<unsigned>::max()) {
                deltaSw = update[Indices::waterSwitchIdx];
                deltaSo -= deltaSw;
            }
        }
        if (currentValue.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg)
        {
            if constexpr (Indices::compositionSwitchIdx != std::numeric_limits<unsigned>::max()) {
                deltaSg = update[Indices::compositionSwitchIdx];
                deltaSo -= deltaSg;
            }
        }
        if (currentValue.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
            if constexpr (Indices::solventSaturationIdx != std::numeric_limits<unsigned>::max()) {
                deltaSs = update[Indices::solventSaturationIdx];
                deltaSo -= deltaSs;
            }
        }

        // maximum saturation delta
        Scalar maxSatDelta = std::max(std::abs(deltaSg), std::abs(deltaSo));
        maxSatDelta = std::max(maxSatDelta, std::abs(deltaSw));
        maxSatDelta = std::max(maxSatDelta, std::abs(deltaSs));

        // scaling factor for saturation deltas to make sure that none of them exceeds
        // the specified threshold value.
        Scalar satAlpha = 1.0;
        if (maxSatDelta > params.dsMax_) {
            satAlpha = params.dsMax_ / maxSatDelta;
        }

        for (unsigned pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            // calculate the update of the current primary variable. For the black-oil
            // model we limit the pressure delta relative to the pressure's current
            // absolute value (Default: 30%) and saturation deltas to an absolute change
            // (Default: 20%). Further, we ensure that the R factors, solvent
            // "saturation" and polymer concentration do not become negative after the
            // update.
            Scalar delta = update[pvIdx];

            // limit pressure delta
            if (pvIdx == Indices::pressureSwitchIdx) {
                if (std::abs(delta) > params.dpMaxRel_ * currentValue[pvIdx]) {
                    delta = ((Scalar{0} < delta) - (delta < Scalar{0})) * params.dpMaxRel_ * currentValue[pvIdx];
                }
            }
            // water saturation delta
            else if (pvIdx == Indices::waterSwitchIdx)
                if (currentValue.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                    delta *= satAlpha;
                }
                else {
                    //Ensure Rvw and Rsw factor does not become negative
                    if (delta > currentValue[ Indices::waterSwitchIdx]) {
                        delta = currentValue[ Indices::waterSwitchIdx];
                    }
                }
            else if (pvIdx == Indices::compositionSwitchIdx) {
                // the switching primary variable for composition is tricky because the
                // "reasonable" value ranges it exhibits vary widely depending on its
                // interpretation since it can represent Sg, Rs or Rv. For now, we only
                // limit saturation deltas and ensure that the R factors do not become
                // negative.
                if (currentValue.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg) {
                    delta *= satAlpha;
                }
                else {
                    // Ensure Rv and Rs factor does not become negative
                    if (delta > currentValue[Indices::compositionSwitchIdx]) {
                        delta = currentValue[Indices::compositionSwitchIdx];
                    }
                }
            }
            else if (enableSolvent && pvIdx == Indices::solventSaturationIdx) {
                // solvent saturation updates are also subject to the Appleyard chop
                if (currentValue.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
                    delta *= satAlpha;
                }
                else {
                    // Ensure Rssolw factor does not become negative
                    if (delta > currentValue[Indices::solventSaturationIdx]) {
                        delta = currentValue[Indices::solventSaturationIdx];
                    }
                }
            }
            else if (enableExtbo && pvIdx == Indices::zFractionIdx) {
                // z fraction updates are also subject to the Appleyard chop
                const auto& curr = currentValue[Indices::zFractionIdx]; // or currentValue[pvIdx] given the block condition
                delta = std::clamp(delta, curr - Scalar{1.0}, curr);
            }
            else if (enablePolymerWeight && pvIdx == Indices::polymerMoleWeightIdx) {
                const double sign = delta >= 0. ? 1. : -1.;
                // maximum change of polymer molecular weight, the unit is MDa.
                // applying this limit to stabilize the simulation. The value itself is still experimental.
                const Scalar maxMolarWeightChange = 100.0;
                delta = sign * std::min(std::abs(delta), maxMolarWeightChange);
                delta *= satAlpha;
            }
            else if (enableFullyImplicitThermal && pvIdx == Indices::temperatureIdx) {
                const double sign = delta >= 0. ? 1. : -1.;
                delta = sign * std::min(std::abs(delta), params.maxTempChange_);
            }
            else if (enableBrine && pvIdx == Indices::saltConcentrationIdx &&
                     enableSaltPrecipitation &&
                     currentValue.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp)
            {
                const Scalar maxSaltSaturationChange = 0.1;
                const Scalar sign = delta >= 0. ? 1. : -1.;
                delta = sign * std::min(std::abs(delta), maxSaltSaturationChange);
            }

            // do the actual update
            nextValue[pvIdx] = currentValue[pvIdx] - delta;

            // keep the solvent saturation between 0 and 1
            if (enableSolvent && pvIdx == Indices::solventSaturationIdx) {
                if (currentValue.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
                    nextValue[pvIdx] = std::min(std::max(nextValue[pvIdx], Scalar{0.0}), Scalar{1.0});
                }
            }

            // keep the z fraction between 0 and 1
            if (enableExtbo && pvIdx == Indices::zFractionIdx) {
                nextValue[pvIdx] = std::min(std::max(nextValue[pvIdx], Scalar{0.0}), Scalar{1.0});
            }

            // keep the polymer concentration above 0
            if (enablePolymer && pvIdx == Indices::polymerConcentrationIdx) {
                nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
            }

            if (enablePolymerWeight && pvIdx == Indices::polymerMoleWeightIdx) {
                nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                const double polymerConcentration = nextValue[Indices::polymerConcentrationIdx];
                if (polymerConcentration < 1.e-10) {
                    nextValue[pvIdx] = 0.0;
                }
            }

            // keep the foam concentration above 0
            if (enableFoam && pvIdx == Indices::foamConcentrationIdx) {
                nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
            }

            if (enableBrine && pvIdx == Indices::saltConcentrationIdx) {
               // keep the salt concentration above 0
                if (!enableSaltPrecipitation ||
                    currentValue.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Cs)
               {
                   nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                }
               // keep the salt saturation below upperlimit
                if (enableSaltPrecipitation &&
                    currentValue.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp)
                {
                   nextValue[pvIdx] = std::min(nextValue[pvIdx], Scalar{1.0-1.e-8});
                }
            }

            // keep the temperature within given values
            if (enableFullyImplicitThermal && pvIdx == Indices::temperatureIdx) {
                nextValue[pvIdx] = std::clamp(nextValue[pvIdx], params.tempMin_, params.tempMax_);
            }

            if (pvIdx == Indices::pressureSwitchIdx) {
                nextValue[pvIdx] = std::clamp(nextValue[pvIdx], params.pressMin_, params.pressMax_);
            }

            // keep the values above 0
            // for the biofilm and calcite, we set an upper limit equal to the initial porosity
            // minus 1e-8. This prevents singularities (e.g., one of the calcite source term is
            // evaluated at 1/(iniPoro - calcite)). The value 1e-8 is taken from the salt precipitation
            // clapping above.
            if constexpr (enableBioeffects) {
                if (pvIdx == Indices::microbialConcentrationIdx) {
                    nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                }
                if (pvIdx == Indices::biofilmVolumeFractionIdx) {
                    nextValue[pvIdx] = std::clamp(nextValue[pvIdx],
                                                  Scalar{0.0},
                                                  problem.referencePorosity(globalDofIdx, 0) - 1e-8);
                }
                if constexpr (enableMICP) {
                    if (pvIdx == Indices::oxygenConcentrationIdx) {
                        nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                    }
                    if (pvIdx == Indices::ureaConcentrationIdx) {
                        nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                    }
                    if (pvIdx == Indices::calciteVolumeFractionIdx) {
                        nextValue[pvIdx] = std::clamp(nextValue[pvIdx], Scalar{0.0},
                                                                        problem.referencePorosity(globalDofIdx, 0) - 1e-8);
                    }
                }
            }
        }

        // switch the new primary variables to something which is physically meaningful.
        // use a threshold value after a switch to make it harder to switch back
        // immediately.
        const bool switched = nextValue.adaptPrimaryVariables(
            problem, fluidSystem, globalDofIdx, params.waterSaturationMax_,
            params.waterOnlyThreshold_, wasSwitched ? params.priVarOscilationThreshold_ : Scalar{0});

        if (params.projectSaturations_) {
            nextValue.chopAndNormalizeSaturations();
        }

        nextValue.checkDefined();
        return switched;
    }
};

} // namespace Opm
#endif
