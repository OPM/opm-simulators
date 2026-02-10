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

#include <opm/models/parallel/threadmanager.hpp>
#include <opm/models/blackoil/blackoilconvectivemixingmodule.hh>
#include <opm/models/blackoil/blackoilmoduleparams.hh>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterialParams.hpp>

#include <cstddef>
#include <stdexcept>
#include <type_traits>

/*
    This file implements a simplified version of the FIBlackOilModel
    which can easily be used in a GPU environment.
    This is developed wihle parallelizing the tpfa matrix assembly, so
    only functionality needed there is implemented.
    To start with this means just being able to access intensiveQuantities
*/

namespace Opm {

template <typename TypeTag, template <class> class Storage = Opm::VectorWithDefaultAllocator>
class SimplifiedGpuFIBlackOilModel
{
public:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using TypeTagPublic = TypeTag;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModuleParams = BlackoilModuleParams<ConvectiveMixingModuleParam<Scalar, Storage>>;

    OPM_HOST_DEVICE SimplifiedGpuFIBlackOilModel(unsigned int nCells)
        : cachedIntensiveQuantities0_(nCells)
    {}

    // TODO: copy for now, but should be move!
    OPM_HOST_DEVICE SimplifiedGpuFIBlackOilModel(
        Storage<BlackOilIntensiveQuantities<TypeTag>> cachedIntensiveQuantities0, Storage<BlackOilIntensiveQuantities<TypeTag>> cachedIntensiveQuantities1, ModuleParams moduleParams)
        : cachedIntensiveQuantities0_(cachedIntensiveQuantities0)
        , cachedIntensiveQuantities1_(cachedIntensiveQuantities1)
        , moduleParams_(moduleParams)
    {
    }

    OPM_HOST_DEVICE ModuleParams& moduleParams() {
        return moduleParams_;
    }

    OPM_HOST_DEVICE SimplifiedGpuFIBlackOilModel() = default;

    OPM_HOST_DEVICE void invalidateAndUpdateIntensiveQuantities(unsigned /* timeIdx */) const {}
    OPM_HOST_DEVICE void invalidateAndUpdateIntensiveQuantitiesOverlap(unsigned /* timeIdx */) const {}
    template <class GridSubDomain>
    OPM_HOST_DEVICE void invalidateAndUpdateIntensiveQuantities(unsigned /* timeIdx */, const GridSubDomain& /* gridSubDomain */) const {}
    OPM_HOST_DEVICE void updateFailed() {}
    OPM_HOST_DEVICE const auto& intensiveQuantities(unsigned int globalIdx, unsigned int timeIdx) const {
        assert(timeIdx == 0 || timeIdx == 1); // TODO: do not use assert for this because it can be run in GPU kernel...
        if (timeIdx == 0) {
            return cachedIntensiveQuantities0_[globalIdx];
        } else {
            return cachedIntensiveQuantities1_[globalIdx];
        }
    }

protected:
    OPM_HOST_DEVICE void updateCachedIntQuants(const unsigned /* timeIdx */) const {}
    template <class EMDArg>
    OPM_HOST_DEVICE void updateCachedIntQuants1(const unsigned /* timeIdx */) const {}
    template <class ...Args>
    OPM_HOST_DEVICE void updateCachedIntQuantsLoop(const unsigned /* timeIdx */) const {}
    template <class ...Args>
    OPM_HOST_DEVICE void updateSingleCachedIntQuantUnchecked(const unsigned /* globalIdx */, const unsigned /* timeIdx */) const {}

public:
    // assumes --enable-storage-cache=false so that we need not worry about updating cpu-caches from gpu kernels
    Storage<BlackOilIntensiveQuantities<TypeTag>> cachedIntensiveQuantities0_;
    Storage<BlackOilIntensiveQuantities<TypeTag>> cachedIntensiveQuantities1_;
    ModuleParams moduleParams_;
};

namespace gpuistl {
    // For now we make a copy of the simplified CPU model because we need to set the fluid system pointer
    // without breaking copy semantics
    template <class TypeTag, typename CpuSimplifiedGpuFIBlackOilModel, typename GpuFluidSystemPtr>
    auto copy_to_gpu_just_find_me(CpuSimplifiedGpuFIBlackOilModel cpuModel, GpuFluidSystemPtr* ptrFluidSystem)
    {
        using Scalar = typename CpuSimplifiedGpuFIBlackOilModel::Scalar;
        // In this copy_to_gpu we need to take care of two things:
        // 1) Copy the cachedIntensiveQuantities_ to GPU (go from vector to GpuBuffer to allocate things on GPU)
        // 2) Ensure that the FluidSystem pointer inside each IntensiveQuantities points to a GPU FluidSystem
        //    The pointer should be set before we move the entire thing to the GPU.

        using IntQuantsType = decltype(cpuModel.cachedIntensiveQuantities0_[0]);

        // Here I have to declare the new type of BOIQ that I want?
        // Because I cannot use the existing boiq and just set another type of poiner
        // in the fluid state.
        // the BOIQ is currently just defined from the TypeTag, so I should
        // probably add another template argument such that I can adjust
        // what type of fluidstate I am using (template it on the gpufluidsystem!)

        // set pointers
        using CorrectTypeTag = typename ::Opm::Properties::TTag::to_gpu_type_t<TypeTag, GpuBuffer>;
        using CorrectTypeTagView = typename ::Opm::Properties::TTag::to_gpu_type_t<TypeTag, GpuView>;
        using CorrectBOIQ = BlackOilIntensiveQuantities<CorrectTypeTagView>;
        using SimplifiedGpuBufferModel = SimplifiedGpuFIBlackOilModel<CorrectTypeTagView, GpuBuffer>;

        std::vector<CorrectBOIQ> cpuIntQuantsWithGpuPtr0;
        for (auto& iq : cpuModel.cachedIntensiveQuantities0_) {
            cpuIntQuantsWithGpuPtr0.push_back(iq.template withOtherFluidSystem<CorrectTypeTagView>(ptrFluidSystem));
        }

        std::vector<CorrectBOIQ> cpuIntQuantsWithGpuPtr1;
        for (auto& iq : cpuModel.cachedIntensiveQuantities1_) {
            cpuIntQuantsWithGpuPtr1.push_back(iq.template withOtherFluidSystem<CorrectTypeTagView>(ptrFluidSystem));
        }

        using ModuleParams = BlackoilModuleParams<ConvectiveMixingModuleParam<Scalar, GpuBuffer>>;
        ModuleParams moduleParams {
            gpuistl::copy_to_gpu(cpuModel.moduleParams_.convectiveMixingModuleParam)
        };

        // create new blackoil model providing a gpubuffer allocated version of the intQuants
        return SimplifiedGpuBufferModel(GpuBuffer<CorrectBOIQ>(cpuIntQuantsWithGpuPtr0), GpuBuffer<CorrectBOIQ>(cpuIntQuantsWithGpuPtr1), moduleParams);

    }

    template <typename GpuSimplifiedGpuFIBlackOilModel>
    auto make_view_just_find_me(GpuSimplifiedGpuFIBlackOilModel& gpuSimplifiedGpuFIBlackOilModel)
    {
        using Scalar = typename GpuSimplifiedGpuFIBlackOilModel::Scalar;
        using ModuleParams = BlackoilModuleParams<ConvectiveMixingModuleParam<Scalar, GpuView>>;

        return SimplifiedGpuFIBlackOilModel<
            typename GpuSimplifiedGpuFIBlackOilModel::TypeTagPublic,
            gpuistl::GpuView>
            (make_view(gpuSimplifiedGpuFIBlackOilModel.cachedIntensiveQuantities0_),
             make_view(gpuSimplifiedGpuFIBlackOilModel.cachedIntensiveQuantities1_),
             ModuleParams{
                 make_view(gpuSimplifiedGpuFIBlackOilModel.moduleParams_.convectiveMixingModuleParam)
             });
    }
} // namespace gpuistl

} // namespace Opm

#endif // SIMPLIFIED_BLACK_OIL_MODEL_HPP
