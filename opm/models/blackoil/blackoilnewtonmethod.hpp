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
 * \copydoc Opm::BlackOilNewtonMethod
 */
#ifndef OPM_BLACK_OIL_NEWTON_METHOD_HPP
#define OPM_BLACK_OIL_NEWTON_METHOD_HPP

#include <opm/common/Exceptions.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/blackoil/blackoilmicpmodules.hh>
#include <opm/models/blackoil/blackoilnewtonmethodparams.hpp>

#include <opm/models/nonlinear/newtonmethod.hh>

#include <opm/models/utils/signum.hh>

#include <algorithm>
#include <cmath>
#include <vector>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct DiscNewtonMethod;

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief A newton solver which is specific to the black oil model.
 */
template <class TypeTag>
class BlackOilNewtonMethod : public GetPropType<TypeTag, Properties::DiscNewtonMethod>
{
    using ParentType = GetPropType<TypeTag, Properties::DiscNewtonMethod>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Linearizer = GetPropType<TypeTag, Properties::Linearizer>;
    using MICPModule = BlackOilMICPModule<TypeTag>;

    static const unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr bool enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();

public:
    explicit BlackOilNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {
        bparams_.read();
    }

    /*!
     * \copydoc NewtonMethod::finishInit()
     */
    void finishInit()
    {
        ParentType::finishInit();

        wasSwitched_.resize(this->model().numTotalDof(), false);
    }

    /*!
     * \brief Register all run-time parameters for the blackoil newton method.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
        BlackoilNewtonParams<Scalar>::registerParameters();
    }

    /*!
     * \brief Returns the number of degrees of freedom for which the
     *        interpretation has changed for the most recent iteration.
     */
    unsigned numPriVarsSwitched() const
    { return numPriVarsSwitched_; }

protected:
    friend NewtonMethod<TypeTag>;
    friend ParentType;

    /*!
     * \copydoc FvBaseNewtonMethod::beginIteration_
     */
    void beginIteration_()
    {
        numPriVarsSwitched_ = 0;
        ParentType::beginIteration_();
    }

    /*!
     * \copydoc FvBaseNewtonMethod::endIteration_
     *
     * \param uCurrentIter Current solution iterator
     * \param uLastIter Last solution iterator
     */
    void endIteration_(SolutionVector& uCurrentIter,
                       const SolutionVector& uLastIter)
    {
#if HAVE_MPI
        // in the MPI enabled case we need to add up the number of DOF
        // for which the interpretation changed over all processes.
        const int localSwitched = numPriVarsSwitched_;
        MPI_Allreduce(&localSwitched,
                      &numPriVarsSwitched_,
                      /*num=*/1,
                      MPI_INT,
                      MPI_SUM,
                      MPI_COMM_WORLD);
#endif // HAVE_MPI

        this->simulator_.model().newtonMethod().endIterMsg()
            << ", num switched=" << numPriVarsSwitched_;

        ParentType::endIteration_(uCurrentIter, uLastIter);
    }

public:
    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual)
    {
        const auto& comm = this->simulator_.gridView().comm();

        int succeeded;
        try {
            ParentType::update_(nextSolution,
                                currentSolution,
                                solutionUpdate,
                                currentResidual);
            succeeded = 1;
        }
        catch (...) {
            succeeded = 0;
        }
        succeeded = comm.min(succeeded);

        if (!succeeded) {
            throw NumericalProblem("A process did not succeed in adapting the primary variables");
        }

        numPriVarsSwitched_ = comm.sum(numPriVarsSwitched_);
    }

    template <class DofIndices>
    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual,
                 const DofIndices& dofIndices)
    {
        const auto zero = 0.0 * solutionUpdate[0];
        for (auto dofIdx : dofIndices) {
            if (solutionUpdate[dofIdx] == zero) {
                continue;
            }
            updatePrimaryVariables_(dofIdx,
                                    nextSolution[dofIdx],
                                    currentSolution[dofIdx],
                                    solutionUpdate[dofIdx],
                                    currentResidual[dofIdx]);
        }
    }

protected:
    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(unsigned globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual)
    {
        static constexpr bool enableSolvent = Indices::solventSaturationIdx >= 0;
        static constexpr bool enableExtbo = Indices::zFractionIdx >= 0;
        static constexpr bool enablePolymer = Indices::polymerConcentrationIdx >= 0;
        static constexpr bool enablePolymerWeight = Indices::polymerMoleWeightIdx >= 0;
        static constexpr bool enableEnergy = Indices::temperatureIdx >= 0;
        static constexpr bool enableFoam = Indices::foamConcentrationIdx >= 0;
        static constexpr bool enableBrine = Indices::saltConcentrationIdx >= 0;
        static constexpr bool enableMICP = Indices::microbialConcentrationIdx >= 0;

        currentValue.checkDefined();
        Valgrind::CheckDefined(update);
        Valgrind::CheckDefined(currentResidual);

        // saturation delta for each phase
        Scalar deltaSw = 0.0;
        Scalar deltaSo = 0.0;
        Scalar deltaSg = 0.0;
        Scalar deltaSs = 0.0;

        if (currentValue.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw)
        {
            deltaSw = update[Indices::waterSwitchIdx];
            deltaSo -= deltaSw;
        }
        if (currentValue.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg)
        {
            deltaSg = update[Indices::compositionSwitchIdx];
            deltaSo -= deltaSg;
        }
        if (currentValue.primaryVarsMeaningSolvent() == PrimaryVariables::SolventMeaning::Ss) {
            deltaSs = update[Indices::solventSaturationIdx];
            deltaSo -= deltaSs;
        }

        // maximum saturation delta
        Scalar maxSatDelta = std::max(std::abs(deltaSg), std::abs(deltaSo));
        maxSatDelta = std::max(maxSatDelta, std::abs(deltaSw));
        maxSatDelta = std::max(maxSatDelta, std::abs(deltaSs));

        // scaling factor for saturation deltas to make sure that none of them exceeds
        // the specified threshold value.
        Scalar satAlpha = 1.0;
        if (maxSatDelta > bparams_.dsMax_) {
            satAlpha = bparams_.dsMax_ / maxSatDelta;
        }

        for (int pvIdx = 0; pvIdx < int(numEq); ++pvIdx) {
            // calculate the update of the current primary variable. For the black-oil
            // model we limit the pressure delta relative to the pressure's current
            // absolute value (Default: 30%) and saturation deltas to an absolute change
            // (Default: 20%). Further, we ensure that the R factors, solvent
            // "saturation" and polymer concentration do not become negative after the
            // update.
            Scalar delta = update[pvIdx];

            // limit pressure delta
            if (pvIdx == Indices::pressureSwitchIdx) {
                if (std::abs(delta) > bparams_.dpMaxRel_ * currentValue[pvIdx]) {
                    delta = signum(delta) * bparams_.dpMaxRel_ * currentValue[pvIdx];
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
            else if (enableEnergy && pvIdx == Indices::temperatureIdx) {
                const double sign = delta >= 0. ? 1. : -1.;
                delta = sign * std::min(std::abs(delta), bparams_.maxTempChange_);
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
            if (enableEnergy && pvIdx == Indices::temperatureIdx) {
                nextValue[pvIdx] = std::clamp(nextValue[pvIdx], bparams_.tempMin_, bparams_.tempMax_);
            }

            if (pvIdx == Indices::pressureSwitchIdx) {
                nextValue[pvIdx] = std::clamp(nextValue[pvIdx], bparams_.pressMin_, bparams_.pressMax_);
            }

            // keep the values above 0
            // for the biofilm and calcite, we set an upper limit equal to the initial porosity
            // minus 1e-8. This prevents singularities (e.g., one of the calcite source term is 
            // evaluated at 1/(iniPoro - calcite)). The value 1e-8 is taken from the salt precipitation
            // clapping above. 
            if constexpr (enableMICP) {
                if (pvIdx == Indices::microbialConcentrationIdx) {
                    nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                }
                if (pvIdx == Indices::oxygenConcentrationIdx) {
                    nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                }
                if (pvIdx == Indices::ureaConcentrationIdx) {
                    nextValue[pvIdx] = std::max(nextValue[pvIdx], Scalar{0.0});
                }
                if (pvIdx == Indices::biofilmConcentrationIdx) {
                    nextValue[pvIdx] = std::clamp(nextValue[pvIdx],
                                                  Scalar{0.0},
                                                  this->problem().referencePorosity(globalDofIdx, 0) - 1e-8);
                }
                if (pvIdx == Indices::calciteConcentrationIdx) {
                    nextValue[pvIdx] = std::clamp(nextValue[pvIdx],
                                                  Scalar{0.0},
                                                  this->problem().referencePorosity(globalDofIdx, 0) - 1e-8);
                }
            }
        }

        // switch the new primary variables to something which is physically meaningful.
        // use a threshold value after a switch to make it harder to switch back
        // immediately.
        if (wasSwitched_[globalDofIdx]) {
            wasSwitched_[globalDofIdx] = nextValue.adaptPrimaryVariables(this->problem(),
                                                                         globalDofIdx,
                                                                         bparams_.waterSaturationMax_,
                                                                         bparams_.waterOnlyThreshold_,
                                                                         bparams_.priVarOscilationThreshold_);
        }
        else {
            wasSwitched_[globalDofIdx] = nextValue.adaptPrimaryVariables(this->problem(),
                                                                         globalDofIdx,
                                                                         bparams_.waterSaturationMax_,
                                                                         bparams_.waterOnlyThreshold_);
        }

        if (wasSwitched_[globalDofIdx]) {
            ++numPriVarsSwitched_;
        }
        if (bparams_.projectSaturations_) {
            nextValue.chopAndNormalizeSaturations();
        }

        nextValue.checkDefined();
    }

private:
    int numPriVarsSwitched_{};

    BlackoilNewtonParams<Scalar> bparams_{};

    // keep track of cells where the primary variable meaning has changed
    // to detect and hinder oscillations
    std::vector<bool> wasSwitched_{};
};

} // namespace Opm

#endif // OPM_BLACK_OIL_NEWTHON_METHOD_HPP
