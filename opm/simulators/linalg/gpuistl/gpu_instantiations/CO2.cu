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

using GpuBufD = gpuistl::GpuBuffer<double>;
using GpuViewD = gpuistl::GpuView<double>;

// instantiate double precision CO2
template<>
const UniformTabulated2DFunction<double, GpuBufD>&
CO2<double, GpuBufD>::tabulatedEnthalpy = gpuistl::move_to_gpu<double, GpuBufD>(CO2Tables::tabulatedEnthalpy);

template<>
const UniformTabulated2DFunction<double, GpuBufD>&
CO2<double, GpuBufD>::tabulatedDensity = gpuistl::move_to_gpu<double, GpuBufD>(CO2Tables::tabulatedDensity);

template<>
const double CO2<double, GpuBufD>::brineSalinity = CO2Tables::brineSalinity;


// instantiate single precision CO2
template<>
const UniformTabulated2DFunction<double, GpuBufD>&
CO2<float, GpuBufD>::tabulatedEnthalpy = gpuistl::move_to_gpu<double, GpuBufD>(CO2Tables::tabulatedEnthalpy);

template<>
const UniformTabulated2DFunction<double, GpuBufD>&
CO2<float, GpuBufD>::tabulatedDensity = gpuistl::move_to_gpu<double, GpuBufD>(CO2Tables::tabulatedDensity);

template<>
const float CO2<float, GpuBufD>::brineSalinity = CO2Tables::brineSalinity;

// template<>
// const UniformTabulated2DFunction<double>&
// CO2<float>::tabulatedEnthalpy = CO2Tables::tabulatedEnthalpy;
// template<>
// const UniformTabulated2DFunction<double>&
// CO2<float>::tabulatedDensity = CO2Tables::tabulatedDensity;
// template<>
// const float CO2<float>::brineSalinity = CO2Tables::brineSalinity;

} // namespace Opm
