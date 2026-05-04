/*
  Copyright 2026 Equinor ASA
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
*/

#ifndef OPM_SIMPLIFIED_FLOW_PROBLEM_GPU_HPP
#define OPM_SIMPLIFIED_FLOW_PROBLEM_GPU_HPP

#include <opm/common/utility/gpuistl_if_available.hpp>
#include <opm/common/utility/VectorWithDefaultAllocator.hpp>

namespace Opm {

template<class Scalar, template <class> class Storage = Opm::VectorWithDefaultAllocator>
class SimplifiedFlowProblemGPU
{
public:
    using ModuleParams = BlackoilModuleParams<ConvectiveMixingModuleParam<Scalar, Storage>>;

    SimplifiedFlowProblemGPU() = default;

    SimplifiedFlowProblemGPU(Storage<Scalar> alpha0,
                             Storage<Scalar> alpha1,
                             Storage<Scalar> alpha2,
                             ModuleParams moduleParams)
        : alpha0_(alpha0)
        , alpha1_(alpha1)
        , alpha2_(alpha2)
        , moduleParams_(moduleParams)
    {}

    OPM_HOST_DEVICE Scalar getAlpha(unsigned globalIndex, unsigned boundaryFaceIndex) const
    {
        assert(boundaryFaceIndex < 3 && "SOMETHING IS WRONG WITH BOUNDARYFACEINDEX");
        if (boundaryFaceIndex == 0) {
            return alpha0_[globalIndex];
        } else if (boundaryFaceIndex == 1) {
            return alpha1_[globalIndex];
        } else if (boundaryFaceIndex == 2) {
            return alpha2_[globalIndex];
        } else {
            OPM_THROW(std::logic_error, "Invalid boundary face index: " + std::to_string(boundaryFaceIndex));
        }
    }

    OPM_HOST_DEVICE const ModuleParams& moduleParams() const { return moduleParams_; }
    OPM_HOST_DEVICE ModuleParams& moduleParams() { return moduleParams_; }

    Storage<Scalar>& alpha0() {
        return alpha0_;
    }

    Storage<Scalar>& alpha1() {
        return alpha1_;
    }

    Storage<Scalar>& alpha2() {
        return alpha2_;
    }

private:
    Storage<Scalar> alpha0_;
    Storage<Scalar> alpha1_;
    Storage<Scalar> alpha2_;
    ModuleParams moduleParams_;
};

namespace gpuistl {

    template<class Scalar>
    SimplifiedFlowProblemGPU<Scalar, GpuBuffer>
    copy_to_gpu(SimplifiedFlowProblemGPU<Scalar, Opm::VectorWithDefaultAllocator>& cpuProblem)
    {
        using ModuleParams = BlackoilModuleParams<ConvectiveMixingModuleParam<Scalar, GpuBuffer>>;
        return SimplifiedFlowProblemGPU<Scalar, GpuBuffer>(
            GpuBuffer<Scalar>(cpuProblem.alpha0()),
            GpuBuffer<Scalar>(cpuProblem.alpha1()),
            GpuBuffer<Scalar>(cpuProblem.alpha2()),
            ModuleParams { copy_to_gpu(cpuProblem.moduleParams().convectiveMixingModuleParam) }
        );
    }

    template<class Scalar>
    SimplifiedFlowProblemGPU<Scalar, GpuView>
    make_view(SimplifiedFlowProblemGPU<Scalar, GpuBuffer>& buffer)
    {
        using ModuleParams = BlackoilModuleParams<ConvectiveMixingModuleParam<Scalar, GpuView>>;
        return SimplifiedFlowProblemGPU<Scalar, GpuView>(
            make_view(buffer.alpha0()),
            make_view(buffer.alpha1()),
            make_view(buffer.alpha2()),
            ModuleParams { make_view(buffer.moduleParams().convectiveMixingModuleParam) }
        );
    }

} // namespace gpuistl

} // namespace Opm

#endif // OPM_SIMPLIFIED_FLOW_PROBLEM_GPU_HPP
