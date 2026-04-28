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
 * \brief GPU-compatible stub of FvBaseElementContext.
 *
 * All methods are annotated with OPM_HOST_DEVICE and return dummy/default
 * values. This class is intended as a placeholder so that GPU kernels that
 * receive an element context can compile; it does not implement any real
 * behaviour.
 *
 * The heavy Dune/grid infrastructure present in FvBaseElementContext is
 * intentionally absent here because those types are not usable inside a
 * CUDA/HIP kernel.
 */
#ifndef OPM_FV_BASE_ELEMENT_CONTEXT_GPU_HH
#define OPM_FV_BASE_ELEMENT_CONTEXT_GPU_HH

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/discretization/common/linearizationtype.hh>

#include <cstddef>

namespace Opm {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief GPU-compatible stub mirroring the public interface of FvBaseElementContext.
 *
 * Uses the same TypeTag-based template parameter as FvBaseElementContext; all
 * concrete types are extracted via GetPropType at compile time.
 *
 * \tparam TypeTag  The type tag from which Scalar, PrimaryVariables, etc. are derived.
 */
template<class TypeTag>
class FvBaseElementContextGpu
{
public:
    using Scalar              = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables    = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using Problem             = GetPropType<TypeTag, Properties::Problem>;
    using Model               = GetPropType<TypeTag, Properties::Model>;
    using GradientCalculator  = GetPropType<TypeTag, Properties::GradientCalculator>;

    // -----------------------------------------------------------------------
    // Construction / assignment
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE FvBaseElementContextGpu() = default;
    OPM_HOST_DEVICE FvBaseElementContextGpu(const FvBaseElementContextGpu&) = default;
    OPM_HOST_DEVICE FvBaseElementContextGpu& operator=(const FvBaseElementContextGpu&) = default;

    // -----------------------------------------------------------------------
    // Stencil / element update stubs
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE void updateAll() {}
    OPM_HOST_DEVICE void updateStencil() {}
    OPM_HOST_DEVICE void updatePrimaryStencil() {}
    OPM_HOST_DEVICE void updateStencilTopology() {}
    OPM_HOST_DEVICE void updateAllIntensiveQuantities() {}
    OPM_HOST_DEVICE void updateIntensiveQuantities(unsigned /*timeIdx*/) {}
    OPM_HOST_DEVICE void updatePrimaryIntensiveQuantities(unsigned /*timeIdx*/) {}
    OPM_HOST_DEVICE void updateIntensiveQuantities(const PrimaryVariables& /*priVars*/,
                                                   unsigned /*dofIdx*/,
                                                   unsigned /*timeIdx*/) {}
    OPM_HOST_DEVICE void updateAllExtensiveQuantities() {}
    OPM_HOST_DEVICE void updateExtensiveQuantities(unsigned /*timeIdx*/) {}

    // -----------------------------------------------------------------------
    // Focus DOF
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE void setFocusDofIndex(unsigned dofIdx) { focusDofIdx_ = dofIdx; }
    OPM_HOST_DEVICE unsigned focusDofIndex() const { return focusDofIdx_; }

    // -----------------------------------------------------------------------
    // Linearization type
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE LinearizationType linearizationType() const
    { return LinearizationType{}; }

    // -----------------------------------------------------------------------
    // Problem / model accessors — return references to dummy stored objects
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE const Problem& problem() const { return problem_; }
    OPM_HOST_DEVICE const Model&   model()   const { return model_;   }

    // -----------------------------------------------------------------------
    // DOF / face counts
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE std::size_t numDof(unsigned /*timeIdx*/) const        { return 0; }
    OPM_HOST_DEVICE std::size_t numPrimaryDof(unsigned /*timeIdx*/) const  { return 0; }
    OPM_HOST_DEVICE std::size_t numInteriorFaces(unsigned /*timeIdx*/) const { return 0; }
    OPM_HOST_DEVICE std::size_t numBoundaryFaces(unsigned /*timeIdx*/) const { return 0; }

    // -----------------------------------------------------------------------
    // Global space index / volume
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE unsigned globalSpaceIndex(unsigned /*dofIdx*/, unsigned /*timeIdx*/) const
    { return 0; }

    OPM_HOST_DEVICE Scalar dofVolume(unsigned /*dofIdx*/, unsigned /*timeIdx*/) const
    { return Scalar{0}; }

    OPM_HOST_DEVICE Scalar dofTotalVolume(unsigned /*dofIdx*/, unsigned /*timeIdx*/) const
    { return Scalar{0}; }

    // -----------------------------------------------------------------------
    // Boundary
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE bool onBoundary() const { return false; }

    // -----------------------------------------------------------------------
    // Intensive quantities
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE const IntensiveQuantities& intensiveQuantities(unsigned /*dofIdx*/,
                                                                   unsigned /*timeIdx*/) const
    { return intensiveQuantitiesStashed_; }

    OPM_HOST_DEVICE IntensiveQuantities& intensiveQuantities(unsigned /*dofIdx*/,
                                                             unsigned /*timeIdx*/)
    { return intensiveQuantitiesStashed_; }

    OPM_HOST_DEVICE const IntensiveQuantities* thermodynamicHint(unsigned /*dofIdx*/,
                                                                  unsigned /*timeIdx*/) const
    { return nullptr; }

    // -----------------------------------------------------------------------
    // Primary variables
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE const PrimaryVariables& primaryVars(unsigned /*dofIdx*/,
                                                        unsigned /*timeIdx*/) const
    { return priVarsStashed_; }

    // -----------------------------------------------------------------------
    // Stash / restore
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE bool haveStashedIntensiveQuantities() const
    { return stashedDofIdx_ != -1; }

    OPM_HOST_DEVICE int stashedDofIdx() const { return stashedDofIdx_; }

    OPM_HOST_DEVICE void stashIntensiveQuantities(unsigned dofIdx)
    { stashedDofIdx_ = static_cast<int>(dofIdx); }

    OPM_HOST_DEVICE void restoreIntensiveQuantities(unsigned /*dofIdx*/)
    { stashedDofIdx_ = -1; }

    // -----------------------------------------------------------------------
    // Gradient calculator
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE const GradientCalculator& gradientCalculator() const
    { return gradientCalculator_; }

    // -----------------------------------------------------------------------
    // Extensive quantities
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE const ExtensiveQuantities& extensiveQuantities(unsigned /*fluxIdx*/,
                                                                   unsigned /*timeIdx*/) const
    { return extensiveQuantitiesStashed_; }

    // -----------------------------------------------------------------------
    // Storage cache flag
    // -----------------------------------------------------------------------

    OPM_HOST_DEVICE bool enableStorageCache() const        { return enableStorageCache_; }
    OPM_HOST_DEVICE void setEnableStorageCache(bool yesno) { enableStorageCache_ = yesno; }

private:
    // Dummy stored objects returned by reference accessors.
    // On device these are default-constructed and carry no meaningful data.
    Problem            problem_{};
    Model              model_{};
    GradientCalculator gradientCalculator_{};
    IntensiveQuantities intensiveQuantitiesStashed_{};
    ExtensiveQuantities extensiveQuantitiesStashed_{};
    PrimaryVariables   priVarsStashed_{};

    int  stashedDofIdx_{-1};
    int  focusDofIdx_{-1};
    bool enableStorageCache_{false};
};

} // namespace Opm

#endif // OPM_FV_BASE_ELEMENT_CONTEXT_GPU_HH
