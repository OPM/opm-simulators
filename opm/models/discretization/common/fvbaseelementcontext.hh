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
 * \copydoc Opm::FvBaseElementContext
 */
#ifndef EWOMS_FV_BASE_ELEMENT_CONTEXT_HH
#define EWOMS_FV_BASE_ELEMENT_CONTEXT_HH

#include "fvbaseproperties.hh"

#include <opm/models/discretization/common/linearizationtype.hh>
#include <opm/models/utils/alignedallocator.hh>

#include <opm/material/common/Unused.hpp>

#include <dune/common/fvector.hh>

#include <vector>

namespace Opm {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief This class stores an array of IntensiveQuantities objects, one
 *        intensive quantities object for each of the element's vertices
 */
template<class TypeTag>
class FvBaseElementContext
{
    using Implementation = GetPropType<TypeTag, Properties::ElementContext>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

    // the history size of the time discretization in number of steps
    enum { timeDiscHistorySize = getPropValue<TypeTag, Properties::TimeDiscHistorySize>() };

    struct DofStore_ {
        IntensiveQuantities intensiveQuantities[timeDiscHistorySize];
        PrimaryVariables priVars[timeDiscHistorySize];
        const IntensiveQuantities *thermodynamicHint[timeDiscHistorySize];
    };
    using DofVarsVector = std::vector<DofStore_>;
    using ExtensiveQuantitiesVector = std::vector<ExtensiveQuantities>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Stencil = GetPropType<TypeTag, Properties::Stencil>;
    using GradientCalculator = GetPropType<TypeTag, Properties::GradientCalculator>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;

    static const unsigned dimWorld = GridView::dimensionworld;
    static const unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();

    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    // we don't allow copies of element contexts!
    FvBaseElementContext(const FvBaseElementContext& ) = delete;

public:
    /*!
     * \brief The constructor.
     */
    explicit FvBaseElementContext(const Simulator& simulator)
        : gridView_(simulator.gridView())
        , stencil_(gridView_, simulator.model().dofMapper() )
    {
        // remember the simulator object
        simulatorPtr_ = &simulator;
        enableStorageCache_ = EWOMS_GET_PARAM(TypeTag, bool, EnableStorageCache);
        stashedDofIdx_ = -1;
        focusDofIdx_ = -1;
    }

    static void *operator new(size_t size)
    { return aligned_alloc(alignof(FvBaseElementContext), size); }

    static void operator delete(void *ptr)
    { aligned_free(ptr); }

    /*!
     * \brief Construct all volume and extensive quantities of an element
     *        from scratch.
     *
     * \param elem The DUNE Codim<0> entity for which the volume
     *             variables ought to be calculated
     */
    void updateAll(const Element& elem)
    {
        asImp_().updateStencil(elem);
        asImp_().updateAllIntensiveQuantities();
        asImp_().updateAllExtensiveQuantities();
    }

    /*!
     * \brief Compute the finite volume geometry for an element.
     *
     * \param elem The grid element for which the finite volume geometry ought to be
     *             computed.
     */
    void updateStencil(const Element& elem)
    {
        // remember the current element
        elemPtr_ = &elem;

        // update the stencil. the center gradients are quite expensive to calculate and
        // most models don't need them, so that we only do this if the model explicitly
        // enables them
        stencil_.update(elem);

        // resize the arrays containing the flux and the volume variables
        dofVars_.resize(stencil_.numDof());
        extensiveQuantities_.resize(stencil_.numInteriorFaces());
    }

    /*!
     * \brief Update the primary topological part of the stencil, but nothing else.
     *
     * \param elem The grid element for which the finite volume geometry ought to be
     *             computed.
     */
    void updatePrimaryStencil(const Element& elem)
    {
        // remember the current element
        elemPtr_ = &elem;

        // update the finite element geometry
        stencil_.updatePrimaryTopology(elem);

        dofVars_.resize(stencil_.numPrimaryDof());
    }

    /*!
     * \brief Update the topological part of the stencil, but nothing else.
     *
     * \param elem The grid element for which the finite volume geometry ought to be
     *             computed.
     */
    void updateStencilTopology(const Element& elem)
    {
        // remember the current element
        elemPtr_ = &elem;

        // update the finite element geometry
        stencil_.updateTopology(elem);
    }

    /*!
     * \brief Compute the intensive quantities of all sub-control volumes of the current
     *        element for all time indices.
     */
    void updateAllIntensiveQuantities()
    {
        if (!enableStorageCache_) {
            // if the storage cache is disabled, we need to calculate the storage term
            // from scratch, i.e. we need the intensive quantities of all of the history.
            for (unsigned timeIdx = 0; timeIdx < timeDiscHistorySize; ++ timeIdx)
                asImp_().updateIntensiveQuantities(timeIdx);
        }
        else
            // if the storage cache is enabled, we only need to recalculate the storage
            // term for the most recent point of history (i.e., for the current iterative
            // solution)
            asImp_().updateIntensiveQuantities(/*timeIdx=*/0);
    }

    /*!
     * \brief Compute the intensive quantities of all sub-control volumes of the current
     *        element for a single time index.
     *
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    void updateIntensiveQuantities(unsigned timeIdx)
    { updateIntensiveQuantities_(timeIdx, numDof(timeIdx)); }

    /*!
     * \brief Compute the intensive quantities of all sub-control volumes of the current
     *        element for a single time index.
     *
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    void updatePrimaryIntensiveQuantities(unsigned timeIdx)
    { updateIntensiveQuantities_(timeIdx, numPrimaryDof(timeIdx)); }

    /*!
     * \brief Compute the intensive quantities of a single sub-control volume of the
     *        current element for a single time index.
     *
     * \param priVars The PrimaryVariables which should be used to calculate the
     *                intensive quantities.
     * \param dofIdx The local index in the current element of the sub-control volume
     *               which should be updated.
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    void updateIntensiveQuantities(const PrimaryVariables& priVars, unsigned dofIdx, unsigned timeIdx)
    { asImp_().updateSingleIntQuants_(priVars, dofIdx, timeIdx); }

    /*!
     * \brief Compute the extensive quantities of all sub-control volume
     *        faces of the current element for all time indices.
     */
    void updateAllExtensiveQuantities()
    { asImp_().updateExtensiveQuantities(/*timeIdx=*/0); }

    /*!
     * \brief Compute the extensive quantities of all sub-control volume
     *        faces of the current element for a single time index.
     *
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    void updateExtensiveQuantities(unsigned timeIdx)
    {
        gradientCalculator_.prepare(/*context=*/asImp_(), timeIdx);

        for (unsigned fluxIdx = 0; fluxIdx < numInteriorFaces(timeIdx); fluxIdx++) {
            extensiveQuantities_[fluxIdx].update(/*context=*/asImp_(),
                                                 /*localIndex=*/fluxIdx,
                                                 timeIdx);
        }
    }

    /*!
     * \brief Sets the degree of freedom on which the simulator is currently "focused" on
     *
     * I.e., in the case of automatic differentiation, all derivatives are with regard to
     * the primary variables of that degree of freedom. Only "primary" DOFs can be
     * focused on.
     */
    void setFocusDofIndex(unsigned dofIdx)
    { focusDofIdx_ = dofIdx; }

    /*!
     * \brief Returns the degree of freedom on which the simulator is currently "focused" on
     *
     * \copydetails setFocusDof()
     */
    unsigned focusDofIndex() const
    { return focusDofIdx_; }

    /*!
     * \brief Returns the linearization type.
     *
     * \copydetails setLinearizationType()
     */
    LinearizationType linearizationType() const
    { return this->model().linearizer().getLinearizationType(); }

    /*!
     * \brief Return a reference to the simulator.
     */
    const Simulator& simulator() const
    { return *simulatorPtr_; }

    /*!
     * \brief Return a reference to the problem.
     */
    const Problem& problem() const
    { return simulatorPtr_->problem(); }

    /*!
     * \brief Return a reference to the model.
     */
    const Model& model() const
    { return simulatorPtr_->model(); }

    /*!
     * \brief Return a reference to the grid view.
     */
    const GridView& gridView() const
    { return gridView_; }

    /*!
     * \brief Return the current element.
     */
    const Element& element() const
    { return *elemPtr_; }

    /*!
     * \brief Return the number of sub-control volumes of the current element.
     */
    size_t numDof(unsigned timeIdx) const
    { return stencil(timeIdx).numDof(); }

    /*!
     * \brief Return the number of primary degrees of freedom of the current element.
     */
    size_t numPrimaryDof(unsigned timeIdx) const
    { return stencil(timeIdx).numPrimaryDof(); }

    /*!
     * \brief Return the number of non-boundary faces which need to be
     *        considered for the flux apporixmation.
     */
    size_t numInteriorFaces(unsigned timeIdx) const
    { return stencil(timeIdx).numInteriorFaces(); }

    /*!
     * \brief Return the number of boundary faces which need to be
     *        considered for the flux apporixmation.
     */
    size_t numBoundaryFaces(unsigned timeIdx) const
    { return stencil(timeIdx).numBoundaryFaces(); }

    /*!
     * \brief Return the current finite element geometry.
     *
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    const Stencil& stencil(unsigned timeIdx OPM_UNUSED) const
    { return stencil_; }

    /*!
     * \brief Return the position of a local entities in global coordinates
     *
     * \param dofIdx The local index of the degree of freedom
     *               in the current element.
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    const GlobalPosition& pos(unsigned dofIdx, unsigned timeIdx OPM_UNUSED) const
    { return stencil_.subControlVolume(dofIdx).globalPos(); }

    /*!
     * \brief Return the global spatial index for a sub-control volume
     *
     * \param dofIdx The local index of the degree of freedom
     *               in the current element.
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    unsigned globalSpaceIndex(unsigned dofIdx, unsigned timeIdx) const
    { return stencil(timeIdx).globalSpaceIndex(dofIdx); }


    /*!
     * \brief Return the element-local volume associated with a degree of freedom
     *
     * In the case of the vertex-centered finite volume method, this is different from
     * the total volume because a finite volume usually spans multiple elements...
     *
     * \param dofIdx The local index of the degree of freedom in the current element.
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    Scalar dofVolume(unsigned dofIdx, unsigned timeIdx) const
    { return stencil(timeIdx).subControlVolume(dofIdx).volume();  }

    /*!
     * \brief Return the total volume associated with a degree of freedom
     *
     * "Total" means the volume controlled by a degree of freedom disregarding the
     * element. (For example in the vertex-centered finite volume method, a control
     * volume typically encompasses parts of multiple elements.)
     *
     * \param dofIdx The local index of the degree of freedom in the current element.
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    Scalar dofTotalVolume(unsigned dofIdx, unsigned timeIdx) const
    { return model().dofTotalVolume(globalSpaceIndex(dofIdx, timeIdx)); }

    /*!
     * \brief Returns whether the current element is on the domain's
     *        boundary.
     */
    bool onBoundary() const
    { return element().hasBoundaryIntersections(); }

    /*!
     * \brief Return a reference to the intensive quantities of a
     *        sub-control volume at a given time.
     *
     * If the time step index is not given, return the volume
     * variables for the current time.
     *
     * \param dofIdx The local index of the degree of freedom in the current element.
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    const IntensiveQuantities& intensiveQuantities(unsigned dofIdx, unsigned timeIdx) const
    {
#ifndef NDEBUG
        assert(dofIdx < numDof(timeIdx));

        if (enableStorageCache_ && timeIdx != 0 && problem().recycleFirstIterationStorage())
            throw std::logic_error("If caching of the storage term is enabled, only the intensive quantities "
                                   "for the most-recent substep (i.e. time index 0) are available!");
#endif

        return dofVars_[dofIdx].intensiveQuantities[timeIdx];
    }

    /*!
     * \brief Return the thermodynamic hint for a given local index.
     *
     * \sa Discretization::thermodynamicHint(int, int)
     *
     * \param dofIdx The local index of the degree of freedom in the current element.
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    const IntensiveQuantities *thermodynamicHint(unsigned dofIdx, unsigned timeIdx) const
    {
        assert(dofIdx < numDof(timeIdx));
        return dofVars_[dofIdx].thermodynamicHint[timeIdx];
    }
    /*!
     * \copydoc intensiveQuantities()
     */
    IntensiveQuantities& intensiveQuantities(unsigned dofIdx, unsigned timeIdx)
    {
        assert(dofIdx < numDof(timeIdx));
        return dofVars_[dofIdx].intensiveQuantities[timeIdx];
    }

    /*!
     * \brief Return the primary variables for a given local index.
     *
     * \param dofIdx The local index of the degree of freedom
     *               in the current element.
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    PrimaryVariables& primaryVars(unsigned dofIdx, unsigned timeIdx)
    {
        assert(dofIdx < numDof(timeIdx));
        return dofVars_[dofIdx].priVars[timeIdx];
    }
    /*!
     * \copydoc primaryVars()
     */
    const PrimaryVariables& primaryVars(unsigned dofIdx, unsigned timeIdx) const
    {
        assert(dofIdx < numDof(timeIdx));
        return dofVars_[dofIdx].priVars[timeIdx];
    }

    /*!
     * \brief Returns true if no intensive quanties are stashed
     *
     * In most cases quantities are stashed only if a partial derivative is to be
     * calculated via finite difference methods.
     */
    bool haveStashedIntensiveQuantities() const
    { return stashedDofIdx_ != -1; }

    /*!
     * \brief Return the (local) index of the DOF for which the primary variables were
     *        stashed
     *
     * If none, then this returns -1.
     */
    int stashedDofIdx() const
    { return stashedDofIdx_; }

    /*!
     * \brief Stash the intensive quantities for a degree of freedom on internal memory.
     *
     * \param dofIdx The local index of the degree of freedom in the current element.
     */
    void stashIntensiveQuantities(unsigned dofIdx)
    {
        assert(dofIdx < numDof(/*timeIdx=*/0));

        intensiveQuantitiesStashed_ = dofVars_[dofIdx].intensiveQuantities[/*timeIdx=*/0];
        priVarsStashed_ = dofVars_[dofIdx].priVars[/*timeIdx=*/0];
        stashedDofIdx_ = static_cast<int>(dofIdx);
    }

    /*!
     * \brief Restores the intensive quantities for a degree of freedom from internal memory.
     *
     * \param dofIdx The local index of the degree of freedom in the current element.
     */
    void restoreIntensiveQuantities(unsigned dofIdx)
    {
        dofVars_[dofIdx].priVars[/*timeIdx=*/0] = priVarsStashed_;
        dofVars_[dofIdx].intensiveQuantities[/*timeIdx=*/0] = intensiveQuantitiesStashed_;
        stashedDofIdx_ = -1;
    }

    /*!
     * \brief Return a reference to the gradient calculation class of
     *        the chosen spatial discretization.
     */
    const GradientCalculator& gradientCalculator() const
    { return gradientCalculator_; }

    /*!
     * \brief Return a reference to the extensive quantities of a
     *        sub-control volume face.
     *
     * \param fluxIdx The local index of the sub-control volume face for which the
     *               extensive quantities are requested
     * \param timeIdx The index of the solution vector used by the time discretization.
     */
    const ExtensiveQuantities& extensiveQuantities(unsigned fluxIdx, unsigned timeIdx OPM_UNUSED) const
    { return extensiveQuantities_[fluxIdx]; }

    /*!
     * \brief Returns true iff the cache for the storage term ought to be used for this context.
     *
     * If it is used, intensive quantities can only be accessed for the most recent time
     * index. (time index 0.)
     */
    bool enableStorageCache() const
    { return enableStorageCache_; }

    /*!
     * \brief Specifies if the cache for the storage term ought to be used for this context.
     */
    void setEnableStorageCache(bool yesno)
    { enableStorageCache_ = yesno; }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

protected:
    /*!
     * \brief Update the first 'n' intensive quantities objects from the primary variables.
     *
     * This method considers the intensive quantities cache.
     */
    void updateIntensiveQuantities_(unsigned timeIdx, size_t numDof)
    {
        // update the intensive quantities for the whole history
        const SolutionVector& globalSol = model().solution(timeIdx);

        // update the non-gradient quantities
        for (unsigned dofIdx = 0; dofIdx < numDof; dofIdx++) {
            unsigned globalIdx = globalSpaceIndex(dofIdx, timeIdx);
            const PrimaryVariables& dofSol = globalSol[globalIdx];
            dofVars_[dofIdx].priVars[timeIdx] = dofSol;

            dofVars_[dofIdx].thermodynamicHint[timeIdx] =
                model().thermodynamicHint(globalIdx, timeIdx);

            const auto *cachedIntQuants = model().cachedIntensiveQuantities(globalIdx, timeIdx);
            if (cachedIntQuants) {
                dofVars_[dofIdx].intensiveQuantities[timeIdx] = *cachedIntQuants;
            }
            else {
                updateSingleIntQuants_(dofSol, dofIdx, timeIdx);
                model().updateCachedIntensiveQuantities(dofVars_[dofIdx].intensiveQuantities[timeIdx],
                                                        globalIdx,
                                                        timeIdx);
            }
        }
    }

    void updateSingleIntQuants_(const PrimaryVariables& priVars, unsigned dofIdx, unsigned timeIdx)
    {
#ifndef NDEBUG
        if (enableStorageCache_ && timeIdx != 0 && problem().recycleFirstIterationStorage())
            throw std::logic_error("If caching of the storage term is enabled, only the intensive quantities "
                                   "for the most-recent substep (i.e. time index 0) are available!");
#endif

        dofVars_[dofIdx].priVars[timeIdx] = priVars;
        dofVars_[dofIdx].intensiveQuantities[timeIdx].update(/*context=*/asImp_(), dofIdx, timeIdx);
    }

    IntensiveQuantities intensiveQuantitiesStashed_;
    PrimaryVariables priVarsStashed_;

    GradientCalculator gradientCalculator_;

    std::vector<DofStore_, aligned_allocator<DofStore_, alignof(DofStore_)> > dofVars_;
    std::vector<ExtensiveQuantities, aligned_allocator<ExtensiveQuantities, alignof(ExtensiveQuantities)> > extensiveQuantities_;

    const Simulator *simulatorPtr_;
    const Element *elemPtr_;
    const GridView gridView_;
    Stencil stencil_;

    int stashedDofIdx_;
    int focusDofIdx_;
    bool enableStorageCache_;
};

} // namespace Opm

#endif
