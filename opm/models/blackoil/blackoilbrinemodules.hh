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

#include "blackoilproperties.hh"
#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/IntervalTabulated2DFunction.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PvtwsaltTable.hpp>
#endif

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/common/fvector.hh>

#include <string>
#include <math.h>

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

    using TabulatedFunction = Tabulated1DFunction<Scalar>;
    using TabulatedTwoDFunction = IntervalTabulated2DFunction<Scalar>;

    static constexpr unsigned saltConcentrationIdx = Indices::saltConcentrationIdx;
    static constexpr unsigned contiBrineEqIdx = Indices::contiBrineEqIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr unsigned enableBrine = enableBrineV;

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the brine module
     */
    static void initFromState(const EclipseState& eclState)
    {
        // some sanity checks: if brine are enabled, the BRINE keyword must be
        // present, if brine are disabled the keyword must not be present.
        if (enableBrine && !eclState.runspec().phases().active(Phase::BRINE)) {
            throw std::runtime_error("Non-trivial brine treatment requested at compile time, but "
                                     "the deck does not contain the BRINE keyword");
        }
        else if (!enableBrine && eclState.runspec().phases().active(Phase::BRINE)) {
            throw std::runtime_error("Brine treatment disabled at compile time, but the deck "
                                     "contains the BRINE keyword");
        }


        if (!eclState.runspec().phases().active(Phase::BRINE))
            return; // brine treatment is supposed to be disabled

        const auto& tableManager = eclState.getTableManager();

        unsigned numPvtRegions = tableManager.getTabdims().getNumPVTTables();
        referencePressure_.resize(numPvtRegions);

        const auto& pvtwsaltTables = tableManager.getPvtwSaltTables();

        // initialize the objects which deal with the BDENSITY keyword
        const auto& bdensityTables = tableManager.getBrineDensityTables();
        if (!bdensityTables.empty()) {
            bdensityTable_.resize(numPvtRegions);
            assert(numPvtRegions == bdensityTables.size());
            for (unsigned pvtRegionIdx = 0; pvtRegionIdx < numPvtRegions; ++ pvtRegionIdx) {
                const auto& bdensityTable = bdensityTables[pvtRegionIdx];
                const auto& pvtwsaltTable = pvtwsaltTables[pvtRegionIdx];
                const auto& c = pvtwsaltTable.getSaltConcentrationColumn();
                bdensityTable_[pvtRegionIdx].setXYContainers(c, bdensityTable);
            }
        }
    }
#endif

    /*!
     * \brief Register all run-time parameters for the black-oil brine module.
     */
    static void registerParameters()
    {
        if (!enableBrine)
            // brine have been disabled at compile time
            return;
    }


    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enableBrine)
            // brine have been disabled at compile time
            return false;

        return pvIdx == saltConcentrationIdx;
    }

    /*!
     * \brief Assign the brine specific primary variables to a PrimaryVariables object
     */
    template <class FluidState>
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  const FluidState& fluidState)
    {
        if (!enableBrine)
            return;

        priVars[saltConcentrationIdx] = fluidState.saltConcentration();
    }

    static std::string primaryVarName([[maybe_unused]] unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        return "saltConcentration";
    }

    static Scalar primaryVarWeight([[maybe_unused]] unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableBrine)
            return false;

        return eqIdx == contiBrineEqIdx;
    }

    static std::string eqName([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        return "conti^brine";
    }

    static Scalar eqWeight([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    // must be called after water storage is computed
    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enableBrine)
            return;

        const auto& fs = intQuants.fluidState();

        LhsEval surfaceVolumeWater =
                Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx))
                * Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.porosity());

        // avoid singular matrix if no water is present.
        surfaceVolumeWater = max(surfaceVolumeWater, 1e-10);

        // Brine in water phase
        const LhsEval massBrine = surfaceVolumeWater
                * Toolbox::template decay<LhsEval>(fs.saltConcentration());

        storage[contiBrineEqIdx] += massBrine;
    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)

    {
        if (!enableBrine)
            return;

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        const unsigned upIdx = extQuants.upstreamIndex(FluidSystem::waterPhaseIdx);
        const unsigned inIdx = extQuants.interiorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

        if (upIdx == inIdx) {
            flux[contiBrineEqIdx] =
                    extQuants.volumeFlux(waterPhaseIdx)
                    *up.fluidState().invB(waterPhaseIdx)
                    *up.fluidState().saltConcentration();
        }
        else {
            flux[contiBrineEqIdx] =
                    extQuants.volumeFlux(waterPhaseIdx)
                    *decay<Scalar>(up.fluidState().invB(waterPhaseIdx))
                    *decay<Scalar>(up.fluidState().saltConcentration());
        }
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
    {
        // do not consider consider the change of Brine primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enableBrine)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        outstream << priVars[saltConcentrationIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enableBrine)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        instream >> priVars0[saltConcentrationIdx];

        // set the primary variables for the beginning of the current time step.
        priVars1[saltConcentrationIdx] = priVars0[saltConcentrationIdx];
    }

    static const Scalar& referencePressure(const ElementContext& elemCtx,
                                           unsigned scvIdx,
                                           unsigned timeIdx)
    {
        unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return referencePressure_[pvtnumRegionIdx];
    }


    static const TabulatedFunction& bdensityTable(const ElementContext& elemCtx,
                                                  unsigned scvIdx,
                                                  unsigned timeIdx)
    {
        unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
        return bdensityTable_[pvtnumRegionIdx];
    }

    static bool hasBDensityTables()
    {
        return !bdensityTable_.empty();
    }

private:
    static std::vector<TabulatedFunction> bdensityTable_;
    static std::vector<Scalar> referencePressure_;
};


template <class TypeTag, bool enableBrineV>
std::vector<typename BlackOilBrineModule<TypeTag, enableBrineV>::TabulatedFunction>
BlackOilBrineModule<TypeTag, enableBrineV>::bdensityTable_;

template <class TypeTag, bool enableBrineV>
std::vector<typename BlackOilBrineModule<TypeTag, enableBrineV>::Scalar>
BlackOilBrineModule<TypeTag, enableBrineV>::referencePressure_;

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilBrineIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        brine extension of the black-oil model.
 */
template <class TypeTag, bool enableBrineV = getPropValue<TypeTag, Properties::EnableBrine>()>
class BlackOilBrineIntensiveQuantities
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
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr unsigned enableBrine = enableBrineV;
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

        auto& fs = asImp_().fluidState_;
        // set saltconcentration
        fs.setSaltConcentration(priVars.makeEvaluation(saltConcentrationIdx, timeIdx));

    }

    const Evaluation& saltConcentration() const
    { return saltConcentration_; }

    const Evaluation& brineRefDensity() const
    { return refDensity_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation saltConcentration_;
    Evaluation refDensity_;

};

template <class TypeTag>
class BlackOilBrineIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    void updateSaltConcentration_(const ElementContext& elemCtx OPM_UNUSED,
                                  unsigned dofIdx OPM_UNUSED,
                                  unsigned timeIdx OPM_UNUSED)
    { }

    const Evaluation& saltConcentration() const
    { throw std::runtime_error("saltConcentration() called but brine are disabled"); }

    const Evaluation& brineRefDensity() const
    { throw std::runtime_error("brineRefDensity() called but brine are disabled"); }

};





} // namespace Ewoms

#endif
