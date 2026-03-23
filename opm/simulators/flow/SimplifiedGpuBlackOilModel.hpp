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
 * \copydoc Opm::FIBlackOilModel
 */
#ifndef SIMPLIFIED_BLACK_OIL_MODEL_HPP
#define SIMPLIFIED_BLACK_OIL_MODEL_HPP

#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/utility/ElementChunks.hpp>

#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoilmoduleparams.hh>
#include <opm/models/parallel/threadmanager.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterialParams.hpp>

#include <cstddef>
#include <stdexcept>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
    This file implements a simplified version of the FIBlackOilModel
    which can be used in a GPU environment.
    This is developed while parallelizing the tpfa matrix assembly, so
    only functionality needed there is implemented.
    To start with this means just being able to access intensiveQuantities.
*/

namespace Opm
{

template <typename TypeTag, template <class> class Storage = Opm::VectorWithDefaultAllocator>
class SimplifiedGpuFIBlackOilModel
{
public:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using TypeTagPublic = TypeTag;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    SimplifiedGpuFIBlackOilModel(unsigned int nCells)
        : cachedIntensiveQuantities0_(nCells)
    {
    }

    // TODO: copy for now, but should be move!
    SimplifiedGpuFIBlackOilModel(
        const Storage<BlackOilIntensiveQuantities<TypeTag>>& cachedIntensiveQuantities0,
        const Storage<BlackOilIntensiveQuantities<TypeTag>>& cachedIntensiveQuantities1)
        : cachedIntensiveQuantities0_(cachedIntensiveQuantities0)
        , cachedIntensiveQuantities1_(cachedIntensiveQuantities1)
    {
    }

    SimplifiedGpuFIBlackOilModel() = default;

    void invalidateAndUpdateIntensiveQuantities(unsigned /* timeIdx */) const
    {
        OPM_THROW(std::logic_error, "invalidateAndUpdateIntensiveQuantities should not be called in SimplifiedGpuFIBlackOilModel");
    }

    void invalidateAndUpdateIntensiveQuantitiesOverlap(unsigned /* timeIdx */) const
    {
        OPM_THROW(std::logic_error, "invalidateAndUpdateIntensiveQuantitiesOverlap should not be called in SimplifiedGpuFIBlackOilModel");
    }

    template <class GridSubDomain>
    void
    invalidateAndUpdateIntensiveQuantities(unsigned /* timeIdx */,
                                           const GridSubDomain& /* gridSubDomain */) const
    {
        OPM_THROW(std::logic_error, "invalidateAndUpdateIntensiveQuantities with GridSubDomain should not be called in SimplifiedGpuFIBlackOilModel");
    }

    void updateFailed()
    {
        OPM_THROW(std::logic_error, "updateFailed should not be called in SimplifiedGpuFIBlackOilModel");
    }

    OPM_HOST_DEVICE const auto& intensiveQuantities(unsigned int globalIdx,
                                                    unsigned int timeIdx) const
    {
        assert(timeIdx == 0
               || timeIdx
                   == 1); // TODO: do not use assert for this because it can be run in GPU kernel...
        if (timeIdx == 0) {
            return cachedIntensiveQuantities0_[globalIdx];
        } else {
            return cachedIntensiveQuantities1_[globalIdx];
        }
    }

protected:
    void updateCachedIntQuants(const unsigned /* timeIdx */) const
    {
        OPM_THROW(std::logic_error, "updateCachedIntQuants should not be called in SimplifiedGpuFIBlackOilModel");
    }

    template <class EMDArg>
    void updateCachedIntQuants1(const unsigned /* timeIdx */) const
    {
        OPM_THROW(std::logic_error, "updateCachedIntQuants1 should not be called in SimplifiedGpuFIBlackOilModel");
    }

    template <class... Args>
    void updateCachedIntQuantsLoop(const unsigned /* timeIdx */) const
    {
        OPM_THROW(std::logic_error, "updateCachedIntQuantsLoop should not be called in SimplifiedGpuFIBlackOilModel");
    }

    template <class... Args>
    void updateSingleCachedIntQuantUnchecked(const unsigned /* globalIdx */,
                                             const unsigned /* timeIdx */) const
    {
        OPM_THROW(std::logic_error, "updateSingleCachedIntQuantUnchecked should not be called in SimplifiedGpuFIBlackOilModel");
    }

public:
    // assumes --enable-storage-cache=false so that we need not worry about updating cpu-caches from
    // gpu kernels
    Storage<BlackOilIntensiveQuantities<TypeTag>> cachedIntensiveQuantities0_;
    Storage<BlackOilIntensiveQuantities<TypeTag>> cachedIntensiveQuantities1_;
};

namespace gpuistl
{
    // For now we make a copy of the simplified CPU model because we need to set the fluid system
    // pointer without breaking copy semantics
    template <class TypeTag, typename GpuFluidSystem>
    auto copy_to_gpu(const SimplifiedGpuFIBlackOilModel<TypeTag>& cpuModel, const GpuFluidSystem& fsys)
    {
        using Scalar = typename SimplifiedGpuFIBlackOilModel<TypeTag>::Scalar;
        // In this copy_to_gpu we need to take care of two things:
        // 1) Copy the cachedIntensiveQuantities_ to GPU (go from vector to GpuBuffer to allocate
        // things on GPU) 2) Ensure that the FluidSystem pointer inside each IntensiveQuantities
        // points to a GPU FluidSystem
        //    The pointer should be set before we move the entire thing to the GPU.

        // set pointers
        using CorrectTypeTagView =
            typename ::Opm::Properties::TTag::to_gpu_type_t<TypeTag, GpuView>;
        using CorrectBOIQ = BlackOilIntensiveQuantities<CorrectTypeTagView>;
        using SimplifiedGpuBufferModel
            = SimplifiedGpuFIBlackOilModel<CorrectTypeTagView, GpuBuffer>;

        std::size_t nCells0 = cpuModel.cachedIntensiveQuantities0_.size();
        std::size_t nCells1 = cpuModel.cachedIntensiveQuantities1_.size();

        // CorrectBOIQ is not default-constructible (non-static FluidSystem requires a pointer at
        // construction time). To enable parallel fill, we create a prototype from the first element
        // and use it to pre-size the vector with valid objects, then overwrite in parallel.
        // Assumption: all cells share the same fsys pointer, so the prototype object is a valid
        // stand-in until it is overwritten.
        std::vector<CorrectBOIQ> cpuIntQuantsWithGpuPtr0;
        if (nCells0 > 0) {
            CorrectBOIQ proto0 =
                cpuModel.cachedIntensiveQuantities0_[0].template withOtherFluidSystem<CorrectTypeTagView>(fsys);
            cpuIntQuantsWithGpuPtr0.resize(nCells0, proto0);


#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (std::size_t i = 1; i < nCells0; ++i) {
                cpuIntQuantsWithGpuPtr0[i] =
                    cpuModel.cachedIntensiveQuantities0_[i].template withOtherFluidSystem<CorrectTypeTagView>(fsys);
            }
        }

        std::vector<CorrectBOIQ> cpuIntQuantsWithGpuPtr1;
        if (nCells1 > 0) {
            CorrectBOIQ proto1 =
                cpuModel.cachedIntensiveQuantities1_[0].template withOtherFluidSystem<CorrectTypeTagView>(fsys);
            cpuIntQuantsWithGpuPtr1.resize(nCells1, proto1);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (std::size_t i = 1; i < nCells1; ++i) {
                cpuIntQuantsWithGpuPtr1[i] =
                    cpuModel.cachedIntensiveQuantities1_[i].template withOtherFluidSystem<CorrectTypeTagView>(fsys);
            }
        }

        GpuBuffer<CorrectBOIQ> gpuBuf0(cpuIntQuantsWithGpuPtr0);
        GpuBuffer<CorrectBOIQ> gpuBuf1(cpuIntQuantsWithGpuPtr1);

        // create new blackoil model providing a gpubuffer allocated version of the intQuants
        auto result = SimplifiedGpuBufferModel(std::move(gpuBuf0), std::move(gpuBuf1));

        return result;
    }

    template <typename LocalTypeTag>
    auto make_view(SimplifiedGpuFIBlackOilModel<LocalTypeTag, GpuBuffer>& gpuSimplifiedGpuFIBlackOilModel)
    {
        return SimplifiedGpuFIBlackOilModel<LocalTypeTag, gpuistl::GpuView>(
            make_view(gpuSimplifiedGpuFIBlackOilModel.cachedIntensiveQuantities0_),
            make_view(gpuSimplifiedGpuFIBlackOilModel.cachedIntensiveQuantities1_));
    }
} // namespace gpuistl

} // namespace Opm

#endif // SIMPLIFIED_BLACK_OIL_MODEL_HPP
