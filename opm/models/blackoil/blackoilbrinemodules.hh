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
 * \brief Contains the classes required to extend the black-oil model by brine.
 */
#ifndef EWOMS_BLACK_OIL_BRINE_MODULE_HH
#define EWOMS_BLACK_OIL_BRINE_MODULE_HH

#include <dune/common/fvector.hh>

#include <opm/models/blackoil/blackoilbrineparams.hpp>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/discretization/common/linearizationtype.hh>

#include <opm/models/utils/basicproperties.hh>

#include <cassert>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace Opm {

/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by brine.
 */
template <class TypeTag, bool enableBrineV = getPropValue<TypeTag, Properties::EnableBrine>()>
class BlackOilBrineModule
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using Toolbox = MathToolbox<Evaluation>;

    using TabulatedFunction = typename BlackOilBrineParams<Scalar>::TabulatedFunction;

    static constexpr unsigned saltConcentrationIdx = Indices::saltConcentrationIdx;
    static constexpr unsigned contiBrineEqIdx = Indices::contiBrineEqIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr bool gasEnabled = Indices::gasEnabled;
    static constexpr bool oilEnabled = Indices::oilEnabled;
    static constexpr bool enableBrine = enableBrineV;
    static constexpr bool enableSaltPrecipitation =
        getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
    //! \brief Set parameters.
    static void setParams(BlackOilBrineParams<Scalar>&& params)
    {
        params_ = params;
    }

    /*!
     * \brief Register all run-time parameters for the black-oil brine module.
     */
    static void registerParameters()
    {}

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if constexpr (enableBrine) {
            return pvIdx == saltConcentrationIdx;
        }
        else {
            return false;
        }
    }

    /*!
     * \brief Assign the brine specific primary variables to a PrimaryVariables object
     */
    template <class FluidState>
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  const FluidState& fluidState)
    {
        if constexpr (enableBrine) {
            priVars[saltConcentrationIdx] = fluidState.saltConcentration();
        }
    }

    static std::string primaryVarName([[maybe_unused]] unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        return "saltConcentration";
    }

    static Scalar primaryVarWeight([[maybe_unused]] unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to choose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if constexpr (enableBrine) {
            return eqIdx == contiBrineEqIdx;
        }
        else {
            return false;
        }
    }

    static std::string eqName([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        return "conti^brine";
    }

    static Scalar eqWeight([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to choose this differently.
        return static_cast<Scalar>(1.0);
    }

    // must be called after water storage is computed
    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if constexpr (enableBrine) {
            const auto& fs = intQuants.fluidState();

            // avoid singular matrix if no water is present.
            const LhsEval surfaceVolumeWater =
                    max(Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx)) *
                        Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx)) *
                        Toolbox::template decay<LhsEval>(intQuants.porosity()),
                        1e-10);

            // Brine in water phase
            const LhsEval massBrine = surfaceVolumeWater *
                                      Toolbox::template decay<LhsEval>(fs.saltConcentration());

            if constexpr (enableSaltPrecipitation) {
                const double saltDensity = intQuants.saltDensity(); // Solid salt density kg/m3
                const LhsEval solidSalt =
                              Toolbox::template decay<LhsEval>(intQuants.porosity()) /
                              (1.0 - Toolbox::template decay<LhsEval>(intQuants.saltSaturation()) + 1.e-8) *
                              saltDensity *
                              Toolbox::template decay<LhsEval>(intQuants.saltSaturation());

                storage[contiBrineEqIdx] += massBrine + solidSalt;
            }
            else {
                storage[contiBrineEqIdx] += massBrine;
            }
        }
    }

    static void computeFlux([[maybe_unused]] RateVector& flux,
                            [[maybe_unused]] const ElementContext& elemCtx,
                            [[maybe_unused]] unsigned scvfIdx,
                            [[maybe_unused]] unsigned timeIdx)
    {
        if constexpr (enableBrine) {
            const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
            unsigned focusIdx = elemCtx.focusDofIndex();
            unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
            flux[contiBrineEqIdx] = 0.0;
            if (upIdx == focusIdx)
                addBrineFluxes_<Evaluation>(flux, elemCtx, scvfIdx, timeIdx);
            else
                addBrineFluxes_<Scalar>(flux, elemCtx, scvfIdx, timeIdx);
        }
    }

    template <class UpstreamEval>
    static void addBrineFluxes_(RateVector& flux,
                                const ElementContext& elemCtx,
                                unsigned scvfIdx,
                                unsigned timeIdx)
    {
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& upFs = up.fluidState();
        const auto& volFlux = extQuants.volumeFlux(waterPhaseIdx);
        addBrineFluxes_<UpstreamEval>(flux, waterPhaseIdx, volFlux, upFs);
    }

    template <class UpEval, class FluidState>
    static void addBrineFluxes_(RateVector& flux,
                                unsigned phaseIdx,
                                const Evaluation& volFlux,
                                const FluidState& upFs)
    {
        if constexpr (enableBrine) { 
            if (phaseIdx == waterPhaseIdx) {
                flux[contiBrineEqIdx] =
                    decay<UpEval>(upFs.saltConcentration())
                    * decay<UpEval>(upFs.invB(waterPhaseIdx))
                    * volFlux;
            }
        }
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables&,
                                     const EqVector&)
    {
        // do not consider the change of Brine primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if constexpr (enableBrine) {
            const unsigned dofIdx = model.dofMapper().index(dof);
            const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
            outstream << priVars[saltConcentrationIdx];
        }
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if constexpr (enableBrine) {
            const unsigned dofIdx = model.dofMapper().index(dof);
            PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
            PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

            instream >> priVars0[saltConcentrationIdx];

            // set the primary variables for the beginning of the current time step.
            priVars1[saltConcentrationIdx] = priVars0[saltConcentrationIdx];
        }
    }

    static Scalar referencePressure(const ElementContext& elemCtx,
                                    unsigned scvIdx,
                                    unsigned timeIdx)
    {
        const unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.referencePressure_[pvtnumRegionIdx];
    }

    static const TabulatedFunction& bdensityTable(const ElementContext& elemCtx,
                                                  unsigned scvIdx,
                                                  unsigned timeIdx)
    {
        const unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.bdensityTable_[pvtnumRegionIdx];
    }

    static const TabulatedFunction& pcfactTable(unsigned satnumRegionIdx)
    { return params_.pcfactTable_[satnumRegionIdx]; }

    static const TabulatedFunction& permfactTable(const ElementContext& elemCtx,
                                                  unsigned scvIdx,
                                                  unsigned timeIdx)
    {
        const unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.permfactTable_[satnumRegionIdx];
    }

    static const TabulatedFunction& permfactTable(unsigned satnumRegionIdx)
    { return params_.permfactTable_[satnumRegionIdx]; }

    static Scalar saltsolTable(const ElementContext& elemCtx,
                               unsigned scvIdx,
                               unsigned timeIdx)
    {
        const unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.saltsolTable_[pvtnumRegionIdx];
    }

    static Scalar saltsolTable(const unsigned pvtnumRegionIdx)
    {
        return params_.saltsolTable_[pvtnumRegionIdx];
    }

    static Scalar saltdenTable(const ElementContext& elemCtx,
                               unsigned scvIdx,
                               unsigned timeIdx)
    {
        const unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.saltdenTable_[pvtnumRegionIdx];
    }

    static Scalar saltdenTable(const unsigned pvtnumRegionIdx)
    {
        return params_.saltdenTable_[pvtnumRegionIdx];
    }

    static bool hasBDensityTables()
    { return !params_.bdensityTable_.empty(); }

    static bool hasSaltsolTables()
    { return !params_.saltsolTable_.empty(); }

    static bool hasPcfactTables()
    {
        if constexpr (enableSaltPrecipitation) {
            return !params_.pcfactTable_.empty();
        }
        else {
            return false;
        }
    }

    static Scalar saltSol(unsigned regionIdx)
    { return params_.saltsolTable_[regionIdx]; }

private:
    static BlackOilBrineParams<Scalar> params_;
};

template <class TypeTag, bool enableBrineV>
BlackOilBrineParams<typename BlackOilBrineModule<TypeTag, enableBrineV>::Scalar>
BlackOilBrineModule<TypeTag, enableBrineV>::params_;

template <class TypeTag, bool enableBrineV>
class BlackOilBrineIntensiveQuantities;

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilBrineIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        brine extension of the black-oil model.
 */
template <class TypeTag>
class BlackOilBrineIntensiveQuantities<TypeTag, /*enableBrineV=*/true>
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using BrineModule = BlackOilBrineModule<TypeTag>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    static constexpr int saltConcentrationIdx = Indices::saltConcentrationIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr bool enableBrine = true;
    static constexpr bool enableSaltPrecipitation =
        getPropValue<TypeTag, Properties::EnableSaltPrecipitation>();
    static constexpr int contiBrineEqIdx = Indices::contiBrineEqIdx;

public:
    /*!
     * \brief Update the intensive properties needed to handle brine from the
     *        primary variables
     *
     */
    void updateSaltConcentration_(const ElementContext& elemCtx,
                                  unsigned dofIdx,
                                  unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const LinearizationType lintype = elemCtx.linearizationType();
        updateSaltConcentration_(priVars, timeIdx, lintype);
    }

    void updateSaltConcentration_(const PrimaryVariables& priVars,
                                  const unsigned timeIdx,
                                  const LinearizationType lintype)
    {
        const unsigned pvtnumRegionIdx = priVars.pvtRegionIndex();
        auto& fs = asImp_().fluidState_;

        if constexpr (enableSaltPrecipitation) {
            saltSolubility_ = BrineModule::saltsolTable(pvtnumRegionIdx);
            saltDensity_ = BrineModule::saltdenTable(pvtnumRegionIdx);

            if (priVars.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp) {
                saltSaturation_ = priVars.makeEvaluation(saltConcentrationIdx, timeIdx, lintype);
                fs.setSaltConcentration(saltSolubility_);
            }
            else {
                saltConcentration_ = priVars.makeEvaluation(saltConcentrationIdx, timeIdx, lintype);
                fs.setSaltConcentration(saltConcentration_);
                saltSaturation_ = 0.0;
            }
            fs.setSaltSaturation(saltSaturation_);
        }
        else {
            saltConcentration_ = priVars.makeEvaluation(saltConcentrationIdx, timeIdx, lintype);
            fs.setSaltConcentration(saltConcentration_);
        }
    }

    void saltPropertiesUpdate_([[maybe_unused]] const ElementContext& elemCtx,
                               [[maybe_unused]] unsigned dofIdx,
                               [[maybe_unused]] unsigned timeIdx)
    {
        if constexpr (enableSaltPrecipitation) {
            const Evaluation porosityFactor  = min(1.0 - saltSaturation(), 1.0); //phi/phi_0

            const auto& permfactTable = BrineModule::permfactTable(elemCtx, dofIdx, timeIdx);

            permFactor_ = permfactTable.eval(porosityFactor);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                asImp_().mobility_[phaseIdx] *= permFactor_;
            }
        }
    }

    const Evaluation& saltConcentration() const
    { return saltConcentration_; }

    const Evaluation& brineRefDensity() const
    { return refDensity_; }

    const Evaluation& saltSaturation() const
    { return saltSaturation_; }

    Scalar saltSolubility() const
    { return saltSolubility_; }

    Scalar saltDensity() const
    { return saltDensity_; }

    const Evaluation& permFactor() const
    { return permFactor_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation saltConcentration_;
    Evaluation refDensity_;
    Evaluation saltSaturation_;
    Evaluation permFactor_;
    Scalar saltSolubility_;
    Scalar saltDensity_;
};

template <class TypeTag>
class BlackOilBrineIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    void updateSaltConcentration_(const ElementContext&,
                                  unsigned,
                                  unsigned)
    {}

    void saltPropertiesUpdate_(const ElementContext&,
                                  unsigned,
                                  unsigned)
    {}

    const Evaluation& saltConcentration() const
    { throw std::runtime_error("saltConcentration() called but brine are disabled"); }

    const Evaluation& brineRefDensity() const
    { throw std::runtime_error("brineRefDensity() called but brine are disabled"); }

    const Evaluation& saltSaturation() const
    { throw std::logic_error("saltSaturation() called but salt precipitation is disabled"); }

    const Scalar saltSolubility() const
    { throw std::logic_error("saltSolubility() called but salt precipitation is disabled"); }

    const Scalar saltDensity() const
    { throw std::logic_error("saltDensity() called but salt precipitation is disabled"); }

    const Evaluation& permFactor() const
    { throw std::logic_error("permFactor() called but salt precipitation is disabled"); }
};

} // namespace Opm

#endif
