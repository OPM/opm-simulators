// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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

/*!
 * \file
 *
 * \brief Unit tests for the black-oil fluid system on the GPU with the live oil PVT class
 *        and the dry or wet gas PVT class.
 *
 * Two configurations are tested, each with two PVT regions:
 * - an SPE1-like deck (PVDG, PVTO, PVTW and DISGAS), copied to the GPU with
 *   copy_to_gpu<LiveOilDryGasPvtApproaches>(), and
 * - a deck with vaporized oil as in Norne (PVTG, PVTO, PVTW, DISGAS and VAPOIL), copied
 *   with copy_to_gpu<LiveOilWetGasPvtApproaches>().
 * The GPU fluid systems are evaluated in kernels, and all GPU results are compared to the
 * same evaluations of the CPU fluid system.
 *
 * The test also checks the compile-time mapping from PVT approaches to PVT classes, and
 * that copy_to_gpu() rejects a CPU fluid system whose PVT approaches differ from the
 * requested ones.
 */
#include <config.h>

#define BOOST_TEST_MODULE TestGpuBlackOilFluidSystemSpe1
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystemNonStatic.hpp>
#include <opm/material/fluidsystems/blackoilpvt/BlackOilPvtApproaches.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>

#include <cuda_runtime.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using Opm::gpuistl::GpuBuffer;
using Opm::gpuistl::GpuView;
using Opm::VectorWithDefaultAllocator;
using Indices = Opm::BlackOilDefaultFluidSystemIndices;

//! The AD variables of the fluid state, see makeFluidState().
enum VarIdx : int { pressureVarIdx, rsVarIdx, rvVarIdx, sgVarIdx, numVars };

template <class Scalar>
using Evaluation = Opm::DenseAd::Evaluation<Scalar, numVars>;

template <class Scalar>
using CpuFluidSystem = Opm::BlackOilFluidSystem<Scalar>;
template <class Scalar, class PvtApproaches = Opm::LiveOilDryGasPvtApproaches>
using GpuBufFluidSystem = Opm::BlackOilFluidSystemNonStatic<Scalar, Indices, GpuBuffer, PvtApproaches>;
template <class Scalar, class PvtApproaches = Opm::LiveOilDryGasPvtApproaches>
using GpuViewFluidSystem = Opm::BlackOilFluidSystemNonStatic<Scalar, Indices, GpuView, PvtApproaches>;

// BlackOilFluidSystem<float> is explicitly instantiated in opm-common independently of
// FLOW_INSTANTIATE_FLOAT, so both types are always tested.
using Scalars = boost::mpl::list<double, float>;

// ---------------------------------------------------------------------------------------
// Compile-time checks of the mapping from PVT approaches to PVT classes.
// ---------------------------------------------------------------------------------------

template <class Scalar, template <typename> typename Storage, class PvtApproaches>
using PvtTypes = Opm::BlackOilPvtTypes<Scalar, Storage, PvtApproaches>;

template <class Scalar>
struct CheckPvtTypes {
    using LiveOilDryGas = Opm::LiveOilDryGasPvtApproaches;
    using LiveOilWetGas = Opm::LiveOilWetGasPvtApproaches;
    using Co2 = Opm::Co2StorePvtApproaches;

    // Live oil, dry gas (SPE1)
    static_assert(std::is_same_v<typename PvtTypes<Scalar, VectorWithDefaultAllocator, LiveOilDryGas>::GasPvt,
                                 Opm::DryGasPvt<Scalar>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuBuffer, LiveOilDryGas>::GasPvt,
                                 Opm::DryGasPvt<Scalar, GpuBuffer, GpuView>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuView, LiveOilDryGas>::GasPvt,
                                 Opm::DryGasPvt<Scalar, GpuView, GpuView>>);

    static_assert(std::is_same_v<typename PvtTypes<Scalar, VectorWithDefaultAllocator, LiveOilDryGas>::OilPvt,
                                 Opm::LiveOilPvt<Scalar>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuBuffer, LiveOilDryGas>::OilPvt,
                                 Opm::LiveOilPvt<Scalar, GpuBuffer, GpuView>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuView, LiveOilDryGas>::OilPvt,
                                 Opm::LiveOilPvt<Scalar, GpuView, GpuView>>);

    static_assert(std::is_same_v<typename PvtTypes<Scalar, VectorWithDefaultAllocator, LiveOilDryGas>::WaterPvt,
                                 Opm::ConstantCompressibilityWaterPvt<Scalar>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuBuffer, LiveOilDryGas>::WaterPvt,
                                 Opm::ConstantCompressibilityWaterPvt<Scalar, GpuBuffer>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuView, LiveOilDryGas>::WaterPvt,
                                 Opm::ConstantCompressibilityWaterPvt<Scalar, GpuView>>);

    // Live oil, wet gas (Norne). Oil and water are as for live oil and dry gas.
    static_assert(std::is_same_v<typename PvtTypes<Scalar, VectorWithDefaultAllocator, LiveOilWetGas>::GasPvt,
                                 Opm::WetGasPvt<Scalar>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuBuffer, LiveOilWetGas>::GasPvt,
                                 Opm::WetGasPvt<Scalar, GpuBuffer, GpuView>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuView, LiveOilWetGas>::GasPvt,
                                 Opm::WetGasPvt<Scalar, GpuView, GpuView>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuView, LiveOilWetGas>::OilPvt,
                                 typename PvtTypes<Scalar, GpuView, LiveOilDryGas>::OilPvt>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuView, LiveOilWetGas>::WaterPvt,
                                 typename PvtTypes<Scalar, GpuView, LiveOilDryGas>::WaterPvt>);

    // CO2STORE
    static_assert(std::is_same_v<typename PvtTypes<Scalar, VectorWithDefaultAllocator, Co2>::GasPvt,
                                 Opm::Co2GasPvt<Scalar>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuBuffer, Co2>::GasPvt,
                                 Opm::Co2GasPvt<Scalar, GpuBuffer>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuView, Co2>::GasPvt,
                                 Opm::Co2GasPvt<Scalar, GpuView>>);

    static_assert(std::is_same_v<typename PvtTypes<Scalar, VectorWithDefaultAllocator, Co2>::OilPvt,
                                 Opm::NullOilPvt<Scalar>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuBuffer, Co2>::OilPvt,
                                 Opm::NullOilPvt<Scalar>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuView, Co2>::OilPvt,
                                 Opm::NullOilPvt<Scalar>>);

    static_assert(std::is_same_v<typename PvtTypes<Scalar, VectorWithDefaultAllocator, Co2>::WaterPvt,
                                 Opm::BrineCo2Pvt<Scalar>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuBuffer, Co2>::WaterPvt,
                                 Opm::BrineCo2Pvt<Scalar, GpuBuffer>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, GpuView, Co2>::WaterPvt,
                                 Opm::BrineCo2Pvt<Scalar, GpuView>>);

    // Multiplexers (CPU only)
    using Multiplexer = Opm::MultiplexerPvtApproaches;
    static_assert(std::is_same_v<typename PvtTypes<Scalar, VectorWithDefaultAllocator, Multiplexer>::GasPvt,
                                 Opm::GasPvtMultiplexer<Scalar>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, VectorWithDefaultAllocator, Multiplexer>::OilPvt,
                                 Opm::OilPvtMultiplexer<Scalar>>);
    static_assert(std::is_same_v<typename PvtTypes<Scalar, VectorWithDefaultAllocator, Multiplexer>::WaterPvt,
                                 Opm::WaterPvtMultiplexer<Scalar>>);

    // The fluid system holds the classes given by its PVT approaches, for any storage. The
    // default is the multiplexers.
    static_assert(std::is_same_v<typename GpuViewFluidSystem<Scalar>::GasPvt,
                                 Opm::DryGasPvt<Scalar, GpuView, GpuView>>);
    static_assert(std::is_same_v<typename GpuViewFluidSystem<Scalar, Co2>::GasPvt,
                                 Opm::Co2GasPvt<Scalar, GpuView>>);
    static_assert(std::is_same_v<typename Opm::BlackOilFluidSystemNonStatic<Scalar, Indices,
                                     VectorWithDefaultAllocator, LiveOilDryGas>::GasPvt,
                                 Opm::DryGasPvt<Scalar>>);
    static_assert(std::is_same_v<typename CpuFluidSystem<Scalar>::GasPvt,
                                 Opm::GasPvtMultiplexer<Scalar>>);
    static_assert(std::is_same_v<Opm::BlackOilFluidSystemNonStatic<Scalar, Indices>,
                                 Opm::BlackOilFluidSystemNonStatic<Scalar, Indices,
                                     VectorWithDefaultAllocator, Multiplexer>>);

    // copy_to_gpu() returns the GPU fluid system for the requested approaches, and
    // make_view() keeps them.
    using CpuNonStatic = Opm::BlackOilFluidSystemNonStatic<Scalar, Indices>;
    static_assert(std::is_same_v<decltype(Opm::gpuistl::copy_to_gpu<LiveOilDryGas>(std::declval<const CpuNonStatic&>())),
                                 GpuBufFluidSystem<Scalar>>);
    static_assert(std::is_same_v<decltype(Opm::gpuistl::copy_to_gpu<LiveOilWetGas>(
                                     std::declval<const CpuNonStatic&>())),
                                 GpuBufFluidSystem<Scalar, LiveOilWetGas>>);
    static_assert(std::is_same_v<decltype(Opm::gpuistl::copy_to_gpu<Co2>(std::declval<const CpuNonStatic&>())),
                                 GpuBufFluidSystem<Scalar, Co2>>);
    static_assert(std::is_same_v<decltype(Opm::gpuistl::make_view(std::declval<GpuBufFluidSystem<Scalar>&>())),
                                 GpuViewFluidSystem<Scalar>>);
    using WetGasBuf = GpuBufFluidSystem<Scalar, LiveOilWetGas>;
    static_assert(std::is_same_v<decltype(Opm::gpuistl::make_view(std::declval<WetGasBuf&>())),
                                 GpuViewFluidSystem<Scalar, LiveOilWetGas>>);

    // The owning GPU PVT objects of DryGasPvt, WetGasPvt and LiveOilPvt are move-only, so
    // the owning GPU fluid systems are as well.
    static_assert(!std::is_copy_constructible_v<GpuBufFluidSystem<Scalar>>);
    static_assert(std::is_move_constructible_v<GpuBufFluidSystem<Scalar>>);
    static_assert(!std::is_copy_constructible_v<GpuBufFluidSystem<Scalar, LiveOilWetGas>>);
    static_assert(std::is_move_constructible_v<GpuBufFluidSystem<Scalar, LiveOilWetGas>>);
};
template struct CheckPvtTypes<double>;
template struct CheckPvtTypes<float>;

// ---------------------------------------------------------------------------------------
// Decks
// ---------------------------------------------------------------------------------------

// Two PVT regions with different table sizes. The PVTO tables have rows with only the
// saturated point (extended by LiveOilPvt) and rows with undersaturated points.
constexpr const char* spe1DeckString = R"(
RUNSPEC
DIMENS
 1 1 1 /
TABDIMS
 * 2 /
OIL
GAS
WATER
DISGAS
METRIC
GRID
DX
 1*100 /
DY
 1*100 /
DZ
 1*10 /
TOPS
 1*1000 /
PORO
 1*0.2 /
PROPS
DENSITY
 850.0 1030.0 0.90 /
 870.0 1020.0 0.80 /
PVTW
 250.0 1.03 4.5E-5 0.40 0.0 /
 200.0 1.02 4.0E-5 0.50 1.0E-5 /
PVDG
  20.0 0.0600 0.012
 100.0 0.0120 0.016
 300.0 0.0045 0.025 /
  10.0 0.1100 0.011
  60.0 0.0200 0.013
 150.0 0.0080 0.018
 350.0 0.0040 0.027 /
PVTO
-- RS     P      BO     MU
   20.0   50.0  1.10   1.20 /
   50.0  100.0  1.20   1.00 /
   90.0  150.0  1.30   0.80
         200.0  1.28   0.85
         250.0  1.26   0.90 /
/
   10.0   40.0  1.05   1.50
          80.0  1.04   1.55 /
   30.0   70.0  1.12   1.30 /
   60.0  110.0  1.21   1.10
         160.0  1.20   1.12
         210.0  1.19   1.15 /
  120.0  170.0  1.40   0.85
         230.0  1.38   0.88 /
/
)";

// Vaporized oil as in Norne: PVTO, PVTW and DENSITY as in the SPE1 deck, and the PVTG
// tables of test_gpuWetGasPvt.cu, which are checked there to be unchanged by the parser's
// low-pressure padding. One region has many pressure rows, the other a row with many Rv
// points. At low pressures region 1 has Rv in Norne's range (about 1e-5).
constexpr const char* vapoilDeckString = R"(
RUNSPEC
DIMENS
 1 1 1 /
TABDIMS
 * 2 /
OIL
GAS
WATER
DISGAS
VAPOIL
METRIC
GRID
DX
 1*100 /
DY
 1*100 /
DZ
 1*10 /
TOPS
 1*1000 /
PORO
 1*0.2 /
PROPS
DENSITY
 850.0 1030.0 0.90 /
 870.0 1020.0 0.80 /
PVTW
 250.0 1.03 4.5E-5 0.40 0.0 /
 200.0 1.02 4.0E-5 0.50 1.0E-5 /
PVTG
-- PG     RV        BG       MUG
   50.0   0.00010   0.0200   0.0140
          0.00000   0.0198   0.0136 /
  100.0   0.00025   0.0120   0.0160 /
  150.0   0.00040   0.0085   0.0185
          0.00030   0.0084   0.0182
          0.00020   0.0083   0.0179 /
/
   40.0   0.00002   0.0250   0.0130
          0.00000   0.0249   0.0129 /
   60.0   0.00012   0.0180   0.0138
          0.00000   0.0178   0.0135 /
   90.0   0.00025   0.0125   0.0150
          0.00010   0.0123   0.0147 /
  120.0   0.00038   0.0098   0.0163
          0.00030   0.0097   0.0161
          0.00022   0.0096   0.0159
          0.00016   0.0095   0.0157
          0.00010   0.0094   0.0155
          0.00006   0.0093   0.0153
          0.00003   0.0092   0.0151
          0.00000   0.0091   0.0149 /
  150.0   0.00052   0.0082   0.0178 /
  180.0   0.00065   0.0072   0.0192
          0.00030   0.0070   0.0185 /
  210.0   0.00075   0.0065   0.0205 /
  240.0   0.00082   0.0060   0.0218
          0.00050   0.0058   0.0210
          0.00020   0.0056   0.0202 /
/
PVTO
-- RS     P      BO     MU
   20.0   50.0  1.10   1.20 /
   50.0  100.0  1.20   1.00 /
   90.0  150.0  1.30   0.80
         200.0  1.28   0.85
         250.0  1.26   0.90 /
/
   10.0   40.0  1.05   1.50
          80.0  1.04   1.55 /
   30.0   70.0  1.12   1.30 /
   60.0  110.0  1.21   1.10
         160.0  1.20   1.12
         210.0  1.19   1.15 /
  120.0  170.0  1.40   0.85
         230.0  1.38   0.88 /
/
)";

// The deck of test_gpuBlackOilFluidSystem.cu, used to check the rejection of a CO2STORE
// fluid system by copy_to_gpu<LiveOilDryGasPvtApproaches>().
constexpr const char* co2StoreDeckString = R"(
RUNSPEC
DIMENS
3 3 3 /
EQLDIMS
/
TABDIMS
/
WATER
GAS
CO2STORE
METRIC
GRID
GRIDFILE
0 0 /
DX
27*1 /
DY
27*1 /
DZ
27*1 /
BOX
1 3 1 3 1 1 /
TOPS
9*0 /
ENDBOX
PERMX
27*1013.25 /
PORO
27*0.25 /
COPY
PERMX PERMY /
PERMX PERMZ /
/
PROPS
SGWFN
0.000000E+00 0.000000E+00 1.000000E+00 3.060000E-02
1.000000E+00 1.000000E+00 0.000000E+00 3.060000E-01 /
SOLUTION
EQUIL
0 300 100 0 0 0 1 1 0 /
SCHEDULE
TSTEP
1 /
)";

constexpr unsigned numRegions = 2;

template <class Scalar>
void initCpuFluidSystem(const char* deckString)
{
    const auto deck = Opm::Parser{}.parseString(deckString);
    const Opm::EclipseState eclState(deck);
    const Opm::Schedule schedule(deck, eclState, std::make_shared<Opm::Python>());
    CpuFluidSystem<Scalar>::initFromState(eclState, schedule);
}

// ---------------------------------------------------------------------------------------
// Evaluation points, fluid state and functors
// ---------------------------------------------------------------------------------------

// Relative tolerance for comparing CPU and GPU results. The GPU may contract
// multiplications and additions into FMAs, so results can differ in the last bits.
template <class Scalar>
constexpr Scalar relTol = std::is_same_v<Scalar, float> ? 1e-4 : 1e-10;

//! One evaluation point. Pressure in [Pa], temperature in [K], Rs and Rv in [sm3/sm3].
//! The oil saturation is given explicitly, so that So = 0 is exact.
template <class Scalar>
struct TestPoint {
    unsigned regionIdx;
    Scalar p;
    Scalar T;
    Scalar rs;
    Scalar rv;
    Scalar sw;
    Scalar sg;
    Scalar so;
    //! The size of Rv used to scale its derivatives in isClose(). Rv is far below 1
    //! (about 1e-5 in Norne), so the saturated Rv at the point is used.
    Scalar rvScale;
};

//! The parts of a fluid state used by the tested fluid system methods. Providing Rs(),
//! Rv(), Rvw(), Rsw() and saltConcentration() means that the fluid system does not fall
//! back to computing them from mass fractions.
template <class Scalar>
struct TestFluidState {
    using ValueType = Evaluation<Scalar>;

    OPM_HOST_DEVICE const ValueType& pressure(unsigned) const
    { return p; }

    OPM_HOST_DEVICE const ValueType& temperature(unsigned) const
    { return T; }

    OPM_HOST_DEVICE const ValueType& saturation(unsigned phaseIdx) const
    { return s[phaseIdx]; }

    OPM_HOST_DEVICE const ValueType& Rs() const
    { return rs; }

    OPM_HOST_DEVICE const ValueType& Rv() const
    { return rv; }

    OPM_HOST_DEVICE ValueType Rvw() const
    { return ValueType(0.0); }

    OPM_HOST_DEVICE ValueType Rsw() const
    { return ValueType(0.0); }

    OPM_HOST_DEVICE ValueType saltConcentration() const
    { return ValueType(0.0); }

    ValueType p;
    ValueType T;
    ValueType rs;
    ValueType rv;
    ValueType s[3]; // indexed by the canonical phase index
};

//! The fluid state of \p pt, with pressure, Rs, Rv and Sg as AD variables (see VarIdx).
//! So = 1 - Sw - Sg, but its value is taken from \p pt.
template <class Scalar>
OPM_HOST_DEVICE TestFluidState<Scalar> makeFluidState(const TestPoint<Scalar>& pt)
{
    using Eval = Evaluation<Scalar>;
    TestFluidState<Scalar> fs;
    fs.p = Eval::createVariable(pt.p, pressureVarIdx);
    fs.T = Eval(pt.T);
    fs.rs = Eval::createVariable(pt.rs, rsVarIdx);
    fs.rv = Eval::createVariable(pt.rv, rvVarIdx);
    fs.s[Indices::waterPhaseIdx] = Eval(pt.sw);
    fs.s[Indices::gasPhaseIdx] = Eval::createVariable(pt.sg, sgVarIdx);
    Eval so(pt.so);
    so.setDerivative(sgVarIdx, Scalar(-1.0));
    fs.s[Indices::oilPhaseIdx] = so;
    return fs;
}

enum class Quantity {
    Density,
    InverseFormationVolumeFactor,
    Viscosity,
    SaturatedDissolutionFactor,
    SaturatedDissolutionFactorMaxSo,
    // Rvw: zero for DryGasPvt and WetGasPvt. Kept to check that the call compiles on the
    // device. The saturated Rv is checked by SaturatedDissolutionFactor of the gas phase.
    SaturatedVaporizationFactor,
};

std::ostream& operator<<(std::ostream& os, Quantity q)
{
    switch (q) {
    case Quantity::Density: return os << "density";
    case Quantity::InverseFormationVolumeFactor: return os << "inverseFormationVolumeFactor";
    case Quantity::Viscosity: return os << "viscosity";
    case Quantity::SaturatedDissolutionFactor: return os << "saturatedDissolutionFactor";
    case Quantity::SaturatedDissolutionFactorMaxSo: return os << "saturatedDissolutionFactor(maxSo)";
    case Quantity::SaturatedVaporizationFactor: return os << "saturatedVaporizationFactor";
    }
    return os;
}

//! Evaluates one fluid system method for one phase. Templated on the fluid system, so
//! the same functor is used for the CPU reference and in the kernel.
struct FluidSystemQuantity {
    Quantity quantity;
    unsigned phaseIdx;

    template <class FluidSystem, class Scalar>
    OPM_HOST_DEVICE Evaluation<Scalar> operator()(const FluidSystem& fluidSystem,
                                                  const TestPoint<Scalar>& pt) const
    {
        const auto fs = makeFluidState(pt);
        switch (quantity) {
        case Quantity::Density:
            return fluidSystem.density(fs, phaseIdx, pt.regionIdx);
        case Quantity::InverseFormationVolumeFactor:
            return fluidSystem.inverseFormationVolumeFactor(fs, phaseIdx, pt.regionIdx);
        case Quantity::Viscosity:
            return fluidSystem.viscosity(fs, phaseIdx, pt.regionIdx);
        case Quantity::SaturatedDissolutionFactor:
            return fluidSystem.saturatedDissolutionFactor(fs, phaseIdx, pt.regionIdx);
        case Quantity::SaturatedDissolutionFactorMaxSo:
            return fluidSystem.saturatedDissolutionFactor(fs, phaseIdx, pt.regionIdx,
                                                          Evaluation<Scalar>(1.0));
        case Quantity::SaturatedVaporizationFactor:
            return fluidSystem.saturatedVaporizationFactor(fs, phaseIdx, pt.regionIdx);
        }
        return Evaluation<Scalar>(0.0);
    }
};

struct ReferenceDensity {
    unsigned phaseIdx;

    template <class FluidSystem, class Scalar>
    OPM_HOST_DEVICE Scalar operator()(const FluidSystem& fluidSystem,
                                      const TestPoint<Scalar>& pt) const
    { return fluidSystem.referenceDensity(phaseIdx, pt.regionIdx); }
};

template <class Functor, class FluidSystem, class Point, class Result>
__global__ void evaluateKernel(Functor f, FluidSystem fluidSystem, const Point* points,
                               Result* results, std::size_t n)
{
    const std::size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        results[i] = f(fluidSystem, points[i]);
    }
}

template <class Functor, class FluidSystem, class Point>
auto evaluateOnCpu(Functor f, const FluidSystem& fluidSystem, const std::vector<Point>& points)
{
    using Result = decltype(f(fluidSystem, points.front()));
    std::vector<Result> results;
    results.reserve(points.size());
    for (const auto& pt : points) {
        results.push_back(f(fluidSystem, pt));
    }
    return results;
}

template <class Functor, class FluidSystem, class Point>
auto evaluateOnGpu(Functor f, FluidSystem fluidSystem, const std::vector<Point>& points)
{
    using Result = decltype(f(fluidSystem, points.front()));
    const GpuBuffer<Point> gpuPoints(points);
    GpuBuffer<Result> gpuResults(points.size());

    constexpr unsigned threadsPerBlock = 128;
    const auto numBlocks
        = static_cast<unsigned>((points.size() + threadsPerBlock - 1) / threadsPerBlock);
    evaluateKernel<<<numBlocks, threadsPerBlock>>>(f, fluidSystem, gpuPoints.data(),
                                                   gpuResults.data(), points.size());
    OPM_GPU_SAFE_CALL(cudaPeekAtLastError());
    OPM_GPU_SAFE_CALL(cudaDeviceSynchronize());

    return gpuResults.asStdVector();
}

template <std::floating_point Scalar>
bool isClose(Scalar a, Scalar b)
{
    if (a == b) {
        return true;
    }
    return std::abs(a - b) <= relTol<Scalar> * std::max(std::abs(a), std::abs(b));
}

//! A derivative is also close if the difference, times the size of its variable at \p pt,
//! is small compared to the value. Derivatives that are zero up to round-off (for example
//! the water viscosity with zero viscosibility) would otherwise fail the relative check.
template <class Scalar>
bool isClose(const Evaluation<Scalar>& a, const Evaluation<Scalar>& b, const TestPoint<Scalar>& pt)
{
    if (!isClose(a.value(), b.value())) {
        return false;
    }

    std::array<Scalar, numVars> varScale{};
    varScale[pressureVarIdx] = pt.p;
    varScale[rsVarIdx] = std::max(std::abs(pt.rs), Scalar(1));
    varScale[rvVarIdx] = pt.rvScale;
    varScale[sgVarIdx] = Scalar(1);
    const Scalar valueScale = std::max(std::abs(a.value()), std::abs(b.value()));
    for (int varIdx = 0; varIdx < Evaluation<Scalar>::numVars; ++varIdx) {
        const Scalar da = a.derivative(varIdx);
        const Scalar db = b.derivative(varIdx);
        if (!isClose(da, db) && std::abs(da - db) * varScale[varIdx] > relTol<Scalar> * valueScale) {
            return false;
        }
    }
    return true;
}

template <class Result, class Scalar>
void checkResults(const std::vector<Result>& cpu,
                  const std::vector<Result>& gpu,
                  const std::vector<TestPoint<Scalar>>& points,
                  const std::string& name)
{
    BOOST_REQUIRE_EQUAL(cpu.size(), gpu.size());
    for (std::size_t i = 0; i < cpu.size(); ++i) {
        bool close = false;
        if constexpr (std::is_floating_point_v<Result>) {
            close = isClose(cpu[i], gpu[i]);
        } else {
            close = isClose(cpu[i], gpu[i], points[i]);
        }
        BOOST_CHECK_MESSAGE(close,
                            name << " at point " << i << " (region " << points[i].regionIdx
                                 << ", p = " << points[i].p << ", Rs = " << points[i].rs
                                 << ", Rv = " << points[i].rv << ", Sw = " << points[i].sw
                                 << ", Sg = " << points[i].sg << ", So = " << points[i].so
                                 << "): CPU = " << cpu[i] << ", GPU = " << gpu[i]);
    }
}

//! Evaluates \p f with the CPU and the GPU fluid system, and compares the results.
template <class Functor, class CpuFs, class GpuFs, class Scalar>
void checkOnGpu(Functor f,
                const CpuFs& cpuFluidSystem,
                const GpuFs& gpuFluidSystem,
                const std::vector<TestPoint<Scalar>>& points,
                const std::string& name)
{
    checkResults(evaluateOnCpu(f, cpuFluidSystem, points),
                 evaluateOnGpu(f, gpuFluidSystem, points),
                 points, name);
}

template <class CpuFs>
std::string quantityName(const CpuFs& cpuFluidSystem, Quantity quantity, unsigned phaseIdx)
{
    std::ostringstream name;
    name << quantity << " (" << cpuFluidSystem.phaseName(phaseIdx) << ")";
    return name.str();
}

constexpr std::array allPhases{Indices::waterPhaseIdx, Indices::oilPhaseIdx, Indices::gasPhaseIdx};

/*!
 * \brief Copies the current CPU fluid system to the GPU with \p PvtApproaches and
 *        compares all tested quantities of all phases at \p points.
 */
template <class PvtApproaches, class Scalar>
void checkFluidSystemOnGpu(const std::vector<TestPoint<Scalar>>& points)
{
    const auto& cpuFluidSystem = CpuFluidSystem<Scalar>::getNonStaticInstance();
    auto gpuFluidSystemBuffer = Opm::gpuistl::copy_to_gpu<PvtApproaches>(cpuFluidSystem);
    const auto gpuFluidSystemView = Opm::gpuistl::make_view(gpuFluidSystemBuffer);

    constexpr std::array quantities{
        Quantity::Density,
        Quantity::InverseFormationVolumeFactor,
        Quantity::Viscosity,
        Quantity::SaturatedDissolutionFactor,
        Quantity::SaturatedDissolutionFactorMaxSo,
        Quantity::SaturatedVaporizationFactor,
    };

    for (const unsigned phaseIdx : allPhases) {
        for (const auto quantity : quantities) {
            checkOnGpu(FluidSystemQuantity{quantity, phaseIdx}, cpuFluidSystem, gpuFluidSystemView,
                       points, quantityName(cpuFluidSystem, quantity, phaseIdx));
        }
        checkOnGpu(ReferenceDensity{phaseIdx}, cpuFluidSystem, gpuFluidSystemView, points,
                   "referenceDensity (" + std::string(cpuFluidSystem.phaseName(phaseIdx)) + ")");
    }
}

/*!
 * \brief Evaluation points covering both branches of the oil phase and, with vaporized
 *        oil, of the gas phase.
 *
 * For each region and pressure, with RsSat and RvSat the saturated Rs and Rv of the CPU
 * PVT objects (RvSat = 0 without vaporized oil):
 *
 * | Sw  | Sg  | So  | Rs          | Rv          | Branches                              |
 * |-----|-----|-----|-------------|-------------|---------------------------------------|
 * | 0.2 | 0.3 | 0.5 | 0.5 RsSat   | 0.5 RvSat   | oil and gas undersaturated            |
 * | 0.2 | 0.3 | 0.5 | 1.001 RsSat | 1.001 RvSat | oil and gas saturated                 |
 * | 0.2 | 0   | 0.8 | 0.5 RsSat   | 0.5 RvSat   | oil undersaturated (no free gas)      |
 * | 0.2 | 0.8 | 0   | 0.5 RsSat   | 0.5 RvSat   | gas undersaturated (no oil)           |
 * | 0.2 | 0.8 | 0   | 0.5 RsSat   | 1.001 RvSat | gas undersaturated, since So = 0      |
 *
 * The saturated tables are used when Rs (Rv) is at least the saturated value and the
 * other phase is present. Rs and Rv are not set to exactly the saturated values: the
 * branch is chosen by comparing with the saturated value computed on the device, which
 * may differ in the last bits.
 *
 * \param pressuresBar Pressures in [bar], which should cover extrapolation below and
 *                     above the tables, and give positive RsSat and RvSat.
 */
template <class Scalar>
std::vector<TestPoint<Scalar>> makeTestPoints(const std::vector<Scalar>& pressuresBar)
{
    struct Case { Scalar sg; Scalar so; Scalar rsFactor; Scalar rvFactor; };
    constexpr std::array<Case, 5> cases{{
        {0.3, 0.5, 0.5, 0.5},
        {0.3, 0.5, 1.001, 1.001},
        {0.0, 0.8, 0.5, 0.5},
        {0.8, 0.0, 0.5, 0.5},
        {0.8, 0.0, 0.5, 1.001},
    }};

    const Scalar T = CpuFluidSystem<Scalar>::reservoirTemperature();
    std::vector<TestPoint<Scalar>> points;
    for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
        for (const Scalar pBar : pressuresBar) {
            const Scalar p = pBar * 1e5;
            const Scalar rsSat
                = CpuFluidSystem<Scalar>::oilPvt().saturatedGasDissolutionFactor(regionIdx, T, p);
            const Scalar rvSat
                = CpuFluidSystem<Scalar>::gasPvt().saturatedOilVaporizationFactor(regionIdx, T, p);
            for (const auto& c : cases) {
                points.push_back({.regionIdx = regionIdx, .p = p, .T = T,
                                  .rs = c.rsFactor * rsSat, .rv = c.rvFactor * rvSat,
                                  .sw = 0.2, .sg = c.sg, .so = c.so,
                                  // Without vaporized oil Rv and its derivatives are zero.
                                  .rvScale = rvSat > 0.0 ? rvSat : Scalar(1.0)});
            }
        }
    }
    return points;
}

//! The number of points that take the saturated branch of the oil and of the gas phase.
template <class Scalar>
std::pair<std::size_t, std::size_t> countSaturated(const std::vector<TestPoint<Scalar>>& points)
{
    std::size_t numOil = 0;
    std::size_t numGas = 0;
    for (const auto& pt : points) {
        const Scalar rsSat = CpuFluidSystem<Scalar>::oilPvt()
                                 .saturatedGasDissolutionFactor(pt.regionIdx, pt.T, pt.p);
        const Scalar rvSat = CpuFluidSystem<Scalar>::gasPvt()
                                 .saturatedOilVaporizationFactor(pt.regionIdx, pt.T, pt.p);
        numOil += (pt.sg > 0.0 && pt.rs >= rsSat) ? 1 : 0;
        numGas += (pt.so > 0.0 && pt.rv >= rvSat) ? 1 : 0;
    }
    return {numOil, numGas};
}

} // namespace

// ---------------------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE_TEMPLATE(Spe1FluidSystemOnGpu, Scalar, Scalars)
{
    initCpuFluidSystem<Scalar>(spe1DeckString);
    const auto& cpuFluidSystem = CpuFluidSystem<Scalar>::getNonStaticInstance();

    BOOST_REQUIRE_EQUAL(cpuFluidSystem.numRegions(), std::size_t{numRegions});
    BOOST_REQUIRE_EQUAL(cpuFluidSystem.numActivePhases(), 3u);
    BOOST_REQUIRE(cpuFluidSystem.enableDissolvedGas());
    BOOST_REQUIRE(cpuFluidSystem.useSaturatedTables());

    const auto points = makeTestPoints<Scalar>({30.0, 60.0, 95.0, 140.0, 220.0, 400.0});

    // Make sure the points take both oil branches, so the checks below cover them.
    const std::size_t numSaturatedOil = countSaturated(points).first;
    BOOST_CHECK_GT(numSaturatedOil, 0u);
    BOOST_CHECK_LT(numSaturatedOil, points.size());

    checkFluidSystemOnGpu<Opm::LiveOilDryGasPvtApproaches>(points);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(WetGasFluidSystemOnGpu, Scalar, Scalars)
{
    using Opm::GasPvtApproach;

    initCpuFluidSystem<Scalar>(vapoilDeckString);
    const auto& cpuFluidSystem = CpuFluidSystem<Scalar>::getNonStaticInstance();

    BOOST_REQUIRE_EQUAL(cpuFluidSystem.numRegions(), std::size_t{numRegions});
    BOOST_REQUIRE_EQUAL(cpuFluidSystem.numActivePhases(), 3u);
    BOOST_REQUIRE(cpuFluidSystem.enableDissolvedGas());
    BOOST_REQUIRE(cpuFluidSystem.enableVaporizedOil());
    BOOST_REQUIRE(!cpuFluidSystem.enableVaporizedWater());
    BOOST_REQUIRE(cpuFluidSystem.useSaturatedTables());
    BOOST_REQUIRE(CpuFluidSystem<Scalar>::gasPvt().approach() == GasPvtApproach::WetGas);

    // Below both the PVTO and the PVTG tables in both regions (but high enough for a
    // positive extrapolated RvSat), inside them, between the end of region 0's PVTG table
    // and the end of the PVTO tables, and above all tables. None is a table row.
    const auto points = makeTestPoints<Scalar>({38.0, 45.0, 75.0, 115.0, 165.0, 220.0, 300.0});

    // Make sure the points take both branches of the oil and the gas phase.
    const auto [numSaturatedOil, numSaturatedGas] = countSaturated(points);
    BOOST_CHECK_GT(numSaturatedOil, 0u);
    BOOST_CHECK_LT(numSaturatedOil, points.size());
    BOOST_CHECK_GT(numSaturatedGas, 0u);
    BOOST_CHECK_LT(numSaturatedGas, points.size());

    // saturationPressure(), and so bubblePointPressure() and dewPointPressure(), are not
    // tested: the Newton iteration of LiveOilPvt and WetGasPvt is CPU only, and the GPU
    // intensive quantities do not use them.
    checkFluidSystemOnGpu<Opm::LiveOilWetGasPvtApproaches>(points);
}

BOOST_AUTO_TEST_CASE(RejectsMismatchingPvtApproaches)
{
    using Scalar = double;
    using Opm::GasPvtApproach;
    using Opm::OilPvtApproach;
    using Opm::WaterPvtApproach;

    // An SPE1 fluid system with the CO2STORE approaches.
    initCpuFluidSystem<Scalar>(spe1DeckString);
    {
        const auto& cpuFluidSystem = CpuFluidSystem<Scalar>::getNonStaticInstance();
        BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu<Opm::Co2StorePvtApproaches>(cpuFluidSystem),
                          std::logic_error);

        // Only the water approach differs.
        using WrongWater = Opm::BlackOilPvtApproaches<GasPvtApproach::DryGas,
                                                      OilPvtApproach::LiveOil,
                                                      WaterPvtApproach::BrineCo2>;
        BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu<WrongWater>(cpuFluidSystem), std::logic_error);

        // The oil phase is active, but NoOil is requested.
        using NoOil = Opm::BlackOilPvtApproaches<GasPvtApproach::DryGas,
                                                 OilPvtApproach::NoOil,
                                                 WaterPvtApproach::ConstantCompressibilityWater>;
        BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu<NoOil>(cpuFluidSystem), std::logic_error);

        // Only the gas approach differs.
        BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu<Opm::LiveOilWetGasPvtApproaches>(cpuFluidSystem),
                          std::logic_error);

        // The matching approaches are accepted.
        BOOST_CHECK_NO_THROW(Opm::gpuistl::copy_to_gpu<Opm::LiveOilDryGasPvtApproaches>(cpuFluidSystem));
    }

    // A fluid system with vaporized oil and the dry gas approaches. initFromState() sets
    // the gas approach from the deck (WetGas for PVTG).
    initCpuFluidSystem<Scalar>(vapoilDeckString);
    {
        const auto& cpuFluidSystem = CpuFluidSystem<Scalar>::getNonStaticInstance();
        BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu<Opm::LiveOilDryGasPvtApproaches>(cpuFluidSystem),
                          std::logic_error);
        BOOST_CHECK_NO_THROW(Opm::gpuistl::copy_to_gpu<Opm::LiveOilWetGasPvtApproaches>(cpuFluidSystem));
    }

    // A CO2STORE fluid system with the live oil and dry gas approaches. The oil phase is
    // inactive in this deck, so initFromState() does not touch the oil multiplexer, which
    // still holds the LiveOilPvt from the previous deck. Reset it, as a fluid system initialized
    // only from the CO2STORE deck would have it.
    CpuFluidSystem<Scalar>::setOilPvt(std::make_shared<typename CpuFluidSystem<Scalar>::OilPvt>());
    initCpuFluidSystem<Scalar>(co2StoreDeckString);
    {
        const auto& cpuFluidSystem = CpuFluidSystem<Scalar>::getNonStaticInstance();
        BOOST_CHECK_THROW(Opm::gpuistl::copy_to_gpu<Opm::LiveOilDryGasPvtApproaches>(cpuFluidSystem),
                          std::logic_error);
        BOOST_CHECK_NO_THROW(Opm::gpuistl::copy_to_gpu<Opm::Co2StorePvtApproaches>(cpuFluidSystem));
    }
}
