#include <config.h>
#include <opm/material/components/CO2.hpp>
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>

#include <opm/material/components/co2tables.inc>

namespace Opm {

/*
    TODO: instantiate GPUBuf and GPUView versions of the CO2 static object
    This file does not really need to be a .cu file, but currently it is to avoid
    possible collisions with the existing CO2.cpp in opm/common
*/

using GpuBufD = const gpuistl::GpuBuffer<double>;
using GpuViewD = gpuistl::GpuView<const double>;

// instantiate double-precision GpuBuffer CO2
template<>
const UniformTabulated2DFunction<double, GpuBufD>&
CO2<double, GpuBufD>::tabulatedEnthalpy = gpuistl::move_to_gpu<double, GpuBufD>(CO2Tables::tabulatedEnthalpy);

template<>
const UniformTabulated2DFunction<double, GpuBufD>&
CO2<double, GpuBufD>::tabulatedDensity = gpuistl::move_to_gpu<double, GpuBufD>(CO2Tables::tabulatedDensity);

template<>
const double CO2<double, GpuBufD>::brineSalinity = CO2Tables::brineSalinity;

// instantiate double-precision GpuView CO2
template<>
const UniformTabulated2DFunction<double, GpuViewD>&
CO2<double, GpuViewD>::tabulatedEnthalpy = gpuistl::make_view<double, GpuBufD, GpuViewD>(CO2<double, GpuBufD>::getEnthalpy());

template<>
const UniformTabulated2DFunction<double, GpuViewD>&
CO2<double, GpuViewD>::tabulatedDensity = gpuistl::make_view<double, GpuBufD, GpuViewD>(CO2<double, GpuBufD>::getDensity());

template<>
const double CO2<double, GpuViewD>::brineSalinity = CO2Tables::brineSalinity;


// instantiate single precision CO2
template<>
const UniformTabulated2DFunction<double, GpuBufD>&
CO2<float, GpuBufD>::tabulatedEnthalpy = gpuistl::move_to_gpu<double, GpuBufD>(CO2Tables::tabulatedEnthalpy);

template<>
const UniformTabulated2DFunction<double, GpuBufD>&
CO2<float, GpuBufD>::tabulatedDensity = gpuistl::move_to_gpu<double, GpuBufD>(CO2Tables::tabulatedDensity);

template<>
const float CO2<float, GpuBufD>::brineSalinity = CO2Tables::brineSalinity;

// instantiate single-precision GpuView CO2
template<>
const UniformTabulated2DFunction<double, GpuViewD>&
CO2<float, GpuViewD>::tabulatedEnthalpy = gpuistl::make_view<double, GpuBufD, GpuViewD>(CO2<float, GpuBufD>::getEnthalpy());

template<>
const UniformTabulated2DFunction<double, GpuViewD>&
CO2<float, GpuViewD>::tabulatedDensity = gpuistl::make_view<double, GpuBufD, GpuViewD>(CO2<float, GpuBufD>::getDensity());

template<>
const float CO2<float, GpuViewD>::brineSalinity = CO2Tables::brineSalinity;

} // namespace Opm
