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
 * \copydoc Opm::PvsModel
 */
#ifndef EWOMS_PVS_MODEL_HH
#define EWOMS_PVS_MODEL_HH

#include <opm/material/densead/Math.hpp>

#include "pvsproperties.hh"
#include "pvslocalresidual.hh"
#include "pvsnewtonmethod.hh"
#include "pvsprimaryvariables.hh"
#include "pvsratevector.hh"
#include "pvsboundaryratevector.hh"
#include "pvsintensivequantities.hh"
#include "pvsextensivequantities.hh"
#include "pvsindices.hh"

#include <opm/models/common/multiphasebasemodel.hh>
#include <opm/models/common/diffusionmodule.hh>
#include <opm/models/common/energymodule.hh>
#include <opm/models/io/vtkcompositionmodule.hh>
#include <opm/models/io/vtkenergymodule.hh>
#include <opm/models/io/vtkdiffusionmodule.hh>

#include <opm/material/fluidmatrixinteractions/NullMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace Opm {
template <class TypeTag>
class PvsModel;
}

BEGIN_PROPERTIES

//! The type tag for the isothermal single phase problems
NEW_TYPE_TAG(PvsModel, INHERITS_FROM(MultiPhaseBaseModel,
                                     VtkPhasePresence,
                                     VtkComposition,
                                     VtkEnergy,
                                     VtkDiffusion));

//! Use the PVS local jacobian operator for the PVS model
SET_TYPE_PROP(PvsModel,
              LocalResidual,
              Opm::PvsLocalResidual<TypeTag>);

//! Use the PVS specific newton method for the PVS model
SET_TYPE_PROP(PvsModel, NewtonMethod, Opm::PvsNewtonMethod<TypeTag>);

//! the Model property
SET_TYPE_PROP(PvsModel, Model, Opm::PvsModel<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(PvsModel, PrimaryVariables, Opm::PvsPrimaryVariables<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(PvsModel, RateVector, Opm::PvsRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(PvsModel, BoundaryRateVector, Opm::PvsBoundaryRateVector<TypeTag>);

//! the IntensiveQuantities property
SET_TYPE_PROP(PvsModel, IntensiveQuantities, Opm::PvsIntensiveQuantities<TypeTag>);

//! the ExtensiveQuantities property
SET_TYPE_PROP(PvsModel, ExtensiveQuantities, Opm::PvsExtensiveQuantities<TypeTag>);

//! The indices required by the isothermal PVS model
SET_TYPE_PROP(PvsModel, Indices, Opm::PvsIndices<TypeTag, /*PVIdx=*/0>);

// set the model to a medium verbosity
SET_INT_PROP(PvsModel, PvsVerbosity, 1);

//! Disable the energy equation by default
SET_BOOL_PROP(PvsModel, EnableEnergy, false);

// disable molecular diffusion by default
SET_BOOL_PROP(PvsModel, EnableDiffusion, false);

//! The basis value for the weight of the pressure primary variable
SET_SCALAR_PROP(PvsModel, PvsPressureBaseWeight, 1.0);

//! The basis value for the weight of the saturation primary variables
SET_SCALAR_PROP(PvsModel, PvsSaturationsBaseWeight, 1.0);

//! The basis value for the weight of the mole fraction primary variables
SET_SCALAR_PROP(PvsModel, PvsMoleFractionsBaseWeight, 1.0);

END_PROPERTIES

namespace Opm {

/*!
 * \ingroup PvsModel
 *
 * \brief A generic compositional multi-phase model using primary-variable
 *        switching.
 *
 * This model assumes a flow of \f$M \geq 1\f$ fluid phases
 * \f$\alpha\f$, each of which is assumed to be a mixture \f$N \geq
 * M\f$ chemical species \f$\kappa\f$.
 *
 * By default, the standard multi-phase Darcy approach is used to determine
 * the velocity, i.e.
 * \f[
 *  \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 *  \left(\mathbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;,
 * \f]
 * although the actual approach which is used can be specified via the
 * \c FluxModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * SET_TYPE_PROP(MyProblemTypeTag, FluxModule, Opm::ForchheimerFluxModule<TypeTag>);
 * \endcode
 *
 * The core of the model is the conservation mass of each component by
 * means of the equation
 * \f[
 * \sum_\alpha \frac{\partial\;\phi c_\alpha^\kappa S_\alpha }{\partial t}
 * - \sum_\alpha \mathrm{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha \right\}
 * - q^\kappa = 0 \;.
 * \f]
 *
 * To close the system mathematically, \f$M\f$ model equations are
 * also required. This model uses the primary variable switching
 * assumptions, which are given by:
 * \f[
 * 0 \stackrel{!}{=}
 *  f_\alpha = \left\{
 *  \begin{array}{cl}
 *    S_\alpha&  \quad \text{if phase }\alpha\text{ is not present} \    \
 *    1 - \sum_\kappa x_\alpha^\kappa&  \quad \text{else}
 *  \end{array}
 *  \right.
 * \f]
 *
 * To make this approach applicable, a pseudo primary variable
 * <em>phase presence</em> has to be introduced. Its purpose is to
 * specify for each phase whether it is present or not. It is a
 * <em>pseudo</em> primary variable because it is not directly considered when
 * linearizing the system in the Newton method, but after each Newton
 * iteration, it gets updated like the "real" primary variables.  The
 * following rules are used for this update procedure:
 *
 * <ul>
 * <li>If phase \f$\alpha\f$ is present according to the pseudo
 *     primary variable, but \f$S_\alpha < 0\f$ after the Newton
 *     update, consider the phase \f$\alpha\f$ disappeared for the
 *     next iteration and use the set of primary variables which
 *     correspond to the new phase presence.</li>
 * <li>If phase \f$\alpha\f$ is not present according to the pseudo
 *     primary variable, but the sum of the component mole fractions
 *     in the phase is larger than 1, i.e. \f$\sum_\kappa
 *     x_\alpha^\kappa > 1\f$, consider the phase \f$\alpha\f$ present
 *     in the the next iteration and update the set of primary
 *     variables to make it consistent with the new phase
 *     presence.</li>
 * <li>In all other cases don't modify the phase presence for phase
 *     \f$\alpha\f$.</li>
 *
 * </ul>
 *
 * The model always requires \f$N\f$ primary variables, but their
 * interpretation is dependent on the phase presence:
 *
 * <ul>
 *
 * <li>The first primary variable is always interpreted as the
 *      pressure of the phase with the lowest index \f$PV_0 =
 *      p_0\f$.</li>
 *
 * <li>Then, \f$M - 1\f$ "switching primary variables" follow, which
 *     are interpreted depending in the presence of the first
 *     \f$M-1\f$ phases: If phase \f$\alpha\f$ is present, its
 *     saturation \f$S_\alpha = PV_i\f$ is used as primary variable;
 *     if it is not present, the mole fraction \f$PV_i =
 *     x_{\alpha^\star}^\alpha\f$ of the component with index
 *     \f$\alpha\f$ in the phase with the lowest index that is present
 *     \f$\alpha^\star\f$ is used instead.</li>
 *
 * <li>Finally, the mole fractions of the \f$N-M\f$ components with
 *     the largest index in the phase with the lowest index that is
 *     present \f$x_{\alpha^\star}^\kappa\f$ are used as primary
 *     variables.</li>
 *
 * </ul>
 */
template <class TypeTag>
class PvsModel
    : public MultiPhaseBaseModel<TypeTag>
{
    typedef MultiPhaseBaseModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Opm::EnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    PvsModel(Simulator& simulator)
        : ParentType(simulator)
    {
        verbosity_ = EWOMS_GET_PARAM(TypeTag, int, PvsVerbosity);
        numSwitched_ = 0;
    }

    /*!
     * \brief Register all run-time parameters for the PVS compositional model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        Opm::VtkPhasePresenceModule<TypeTag>::registerParameters();
        Opm::VtkCompositionModule<TypeTag>::registerParameters();

        if (enableDiffusion)
            Opm::VtkDiffusionModule<TypeTag>::registerParameters();

        if (enableEnergy)
            Opm::VtkEnergyModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, int, PvsVerbosity,
                             "The verbosity level of the primary variable "
                             "switching model");
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "pvs"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(unsigned pvIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::primaryVarName(pvIdx)).empty())
            return s;

        std::ostringstream oss;
        if (pvIdx == Indices::pressure0Idx)
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        else if (Indices::switch0Idx <= pvIdx
                 && pvIdx < Indices::switch0Idx + numPhases - 1)
            oss << "switch_" << pvIdx - Indices::switch0Idx;
        else if (Indices::switch0Idx + numPhases - 1 <= pvIdx
                 && pvIdx < Indices::switch0Idx + numComponents - 1)
            oss << "auxMoleFrac^" << FluidSystem::componentName(pvIdx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::eqName
     */
    std::string eqName(unsigned eqIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::eqName(eqIdx)).empty())
            return s;

        std::ostringstream oss;
        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx
                                                     + numComponents) {
            unsigned compIdx = eqIdx - Indices::conti0EqIdx;
            oss << "continuity^" << FluidSystem::componentName(compIdx);
        }
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::updateFailed
     */
    void updateFailed()
    {
        ParentType::updateFailed();
        numSwitched_ = 0;
    }

    /*!
     * \copydoc FvBaseDiscretization::updateBegin
     */
    void updateBegin()
    {
        ParentType::updateBegin();

        // find the a reference pressure. The first degree of freedom
        // might correspond to non-interior entities which would lead
        // to an undefined value, so we have to iterate...
        size_t nDof = this->numTotalDof();
        for (unsigned dofIdx = 0; dofIdx < nDof; ++ dofIdx) {
            if (this->dofTotalVolume(dofIdx) > 0.0) {
                referencePressure_ =
                    this->solution(/*timeIdx=*/0)[dofIdx][/*pvIdx=*/Indices::pressure0Idx];
                if (referencePressure_ > 0.0)
                    break;
            }
        }
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(unsigned globalDofIdx, unsigned pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(*this, globalDofIdx, pvIdx);
        if (tmp > 0)
            // energy related quantity
            return tmp;

        if (Indices::pressure0Idx == pvIdx) {
            return 10 / referencePressure_;
        }

        if (Indices::switch0Idx <= pvIdx && pvIdx < Indices::switch0Idx
                                                    + numPhases - 1) {
            unsigned phaseIdx = pvIdx - Indices::switch0Idx;

            if (!this->solution(/*timeIdx=*/0)[globalDofIdx].phaseIsPresent(phaseIdx))
                // for saturations, the weight is always 1
                return 1;

            // for saturations, the PvsMoleSaturationsBaseWeight
            // property determines the weight
            return GET_PROP_VALUE(TypeTag, PvsSaturationsBaseWeight);
        }

        // for mole fractions, the PvsMoleFractionsBaseWeight
        // property determines the weight
        return GET_PROP_VALUE(TypeTag, PvsMoleFractionsBaseWeight);
    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(unsigned globalDofIdx, unsigned eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(*this, globalDofIdx, eqIdx);
        if (tmp > 0)
            // energy related equation
            return tmp;

        unsigned compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numComponents);

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
    }

    /*!
     * \copydoc FvBaseDiscretization::advanceTimeLevel
     */
    void advanceTimeLevel()
    {
        ParentType::advanceTimeLevel();
        numSwitched_ = 0;
    }

    /*!
     * \brief Return true if the primary variables were switched for
     *        at least one vertex after the last timestep.
     */
    bool switched() const
    { return numSwitched_ > 0; }

    /*!
     * \copydoc FvBaseDiscretization::serializeEntity
     */
    template <class DofEntity>
    void serializeEntity(std::ostream& outstream, const DofEntity& dofEntity)
    {
        // write primary variables
        ParentType::serializeEntity(outstream, dofEntity);

        unsigned dofIdx = static_cast<unsigned>(this->dofMapper().index(dofEntity));
        if (!outstream.good())
            throw std::runtime_error("Could not serialize DOF "+std::to_string(dofIdx));

        outstream << this->solution(/*timeIdx=*/0)[dofIdx].phasePresence() << " ";
    }

    /*!
     * \copydoc FvBaseDiscretization::deserializeEntity
     */
    template <class DofEntity>
    void deserializeEntity(std::istream& instream, const DofEntity& dofEntity)
    {
        // read primary variables
        ParentType::deserializeEntity(instream, dofEntity);

        // read phase presence
        unsigned dofIdx = static_cast<unsigned>(this->dofMapper().index(dofEntity));
        if (!instream.good())
            throw std::runtime_error("Could not deserialize DOF "+std::to_string(dofIdx));

        short tmp;
        instream >> tmp;
        this->solution(/*timeIdx=*/0)[dofIdx].setPhasePresence(tmp);
        this->solution(/*timeIdx=*/1)[dofIdx].setPhasePresence(tmp);
    }

    /*!
     * \internal
     * \brief Do the primary variable switching after a Newton iteration.
     *
     * This is an internal method that needs to be public because it
     * gets called by the Newton method after an update.
     */
    void switchPrimaryVars_()
    {
        numSwitched_ = 0;

        int succeeded;
        try {
            std::vector<bool> visited(this->numGridDof(), false);
            ElementContext elemCtx(this->simulator_);

            ElementIterator elemIt = this->gridView_.template begin<0>();
            ElementIterator elemEndIt = this->gridView_.template end<0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const Element& elem = *elemIt;
                if (elem.partitionType() != Dune::InteriorEntity)
                    continue;
                elemCtx.updateStencil(elem);

                size_t numLocalDof = elemCtx.stencil(/*timeIdx=*/0).numPrimaryDof();
                for (unsigned dofIdx = 0; dofIdx < numLocalDof; ++dofIdx) {
                    unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                    if (visited[globalIdx])
                        continue;
                    visited[globalIdx] = true;

                    // compute the intensive quantities of the current degree of freedom
                    auto& priVars = this->solution(/*timeIdx=*/0)[globalIdx];
                    elemCtx.updateIntensiveQuantities(priVars, dofIdx, /*timeIdx=*/0);
                    const IntensiveQuantities& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);

                    // evaluate primary variable switch
                    short oldPhasePresence = priVars.phasePresence();

                    // set the primary variables and the new phase state
                    // from the current fluid state
                    priVars.assignNaive(intQuants.fluidState());

                    if (oldPhasePresence != priVars.phasePresence()) {
                        if (verbosity_ > 1)
                            printSwitchedPhases_(elemCtx,
                                                 dofIdx,
                                                 intQuants.fluidState(),
                                                 oldPhasePresence,
                                                 priVars);
                        ++numSwitched_;
                    }
                }
            }

            succeeded = 1;
        }
        catch (...)
        {
            std::cout << "rank " << this->simulator_.gridView().comm().rank()
                      << " caught an exception during primary variable switching"
                      << "\n"  << std::flush;
            succeeded = 0;
        }
        succeeded = this->simulator_.gridView().comm().min(succeeded);

        if (!succeeded)
            throw Opm::NumericalIssue("A process did not succeed in adapting the primary variables");

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        numSwitched_ = this->gridView_.comm().sum(numSwitched_);

        if (verbosity_ > 0)
            this->simulator_.model().newtonMethod().endIterMsg()
                << ", num switched=" << numSwitched_;
    }

    template <class FluidState>
    void printSwitchedPhases_(const ElementContext& elemCtx,
                              unsigned dofIdx,
                              const FluidState& fs,
                              short oldPhasePresence,
                              const PrimaryVariables& newPv) const
    {
        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            bool oldPhasePresent = (oldPhasePresence&  (1 << phaseIdx)) > 0;
            bool newPhasePresent = newPv.phaseIsPresent(phaseIdx);
            if (oldPhasePresent == newPhasePresent)
                continue;

            const auto& pos = elemCtx.pos(dofIdx, /*timeIdx=*/0);
            if (oldPhasePresent && !newPhasePresent) {
                std::cout << "'" << FluidSystem::phaseName(phaseIdx)
                          << "' phase disappears at position " << pos
                          << ". saturation=" << fs.saturation(phaseIdx)
                          << std::flush;
            }
            else {
                Scalar sumx = 0;
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                    sumx += FsToolbox::value(fs.moleFraction(phaseIdx, compIdx));

                std::cout << "'" << FluidSystem::phaseName(phaseIdx)
                          << "' phase appears at position " << pos
                          << " sum x = " << sumx  << std::flush;
            }
        }

        std::cout << ", new primary variables: ";
        newPv.print();
        std::cout << "\n"  << std::flush;
    }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        // add the VTK output modules which are meaningful for the model
        this->addOutputModule(new Opm::VtkPhasePresenceModule<TypeTag>(this->simulator_));
        this->addOutputModule(new Opm::VtkCompositionModule<TypeTag>(this->simulator_));
        if (enableDiffusion)
            this->addOutputModule(new Opm::VtkDiffusionModule<TypeTag>(this->simulator_));
        if (enableEnergy)
            this->addOutputModule(new Opm::VtkEnergyModule<TypeTag>(this->simulator_));
    }

    mutable Scalar referencePressure_;

    // number of switches of the phase state in the last Newton
    // iteration
    unsigned numSwitched_;

    // verbosity of the model
    int verbosity_;
};
}

#endif
