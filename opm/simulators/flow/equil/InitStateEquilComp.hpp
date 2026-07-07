// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 SINTEF Digital

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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/**
 * \file
 *
 * \brief Hydrostatic equilibration for the compositional simulator (EQUIL + ZMFVD).
 */
#ifndef OPM_INIT_STATE_EQUIL_COMP_HPP
#define OPM_INIT_STATE_EQUIL_COMP_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/constraintsolvers/SaturationPressure.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <opm/input/eclipse/EclipseState/Compositional/CompositionalConfig.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/Equil.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RtempvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/ZmfvdTable.hpp>

#include <opm/simulators/flow/equil/PressureFunction.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace Opm {
namespace EQUIL {
namespace Comp {

namespace Details {

/// Right-hand side of the hydrostatic ODE dp/ddepth = rho(depth, p) * g for a
/// fluid whose density follows from the cubic equation of state at the given
/// temperature and composition.  The EOS root (liquid or vapour) is selected
/// by the phase index.
template <class FluidSystem>
class EosDensityODE
{
public:
    using Scalar = typename FluidSystem::Scalar;
    using CompVec = std::array<Scalar, FluidSystem::numComponents>;
    using CompositionFunction = std::function<CompVec(const Scalar)>;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;

    EosDensityODE(CompositionFunction composition,
                  const TabulatedFunction& tempVdTable,
                  const unsigned phaseIdx,
                  const CompositionalConfig::EOSType eosType,
                  const Scalar normGrav)
        : composition_(std::move(composition))
        , tempVdTable_(tempVdTable)
        , phaseIdx_(phaseIdx)
        , eosType_(eosType)
        , g_(normGrav)
    {}

    Scalar operator()(const Scalar depth,
                      const Scalar press) const
    {
        const CompVec z = composition_(depth);
        const Scalar temp = tempVdTable_.eval(depth, /*extrapolate=*/true);

        CompositionalFluidState<Scalar, FluidSystem> fs;
        fs.setTemperature(temp);
        fs.setPressure(FluidSystem::oilPhaseIdx, press);
        fs.setPressure(FluidSystem::gasPhaseIdx, press);
        for (unsigned compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx) {
            fs.setMoleFraction(phaseIdx_, compIdx, z[compIdx]);
        }

        typename FluidSystem::template ParameterCache<Scalar> paramCache(eosType_);
        paramCache.updatePhase(fs, phaseIdx_);

        return FluidSystem::density(fs, paramCache, phaseIdx_) * g_;
    }

private:
    CompositionFunction composition_;
    const TabulatedFunction& tempVdTable_;
    unsigned phaseIdx_;
    CompositionalConfig::EOSType eosType_;
    Scalar g_;
};

} // namespace Details

/*!
 * \brief Computes the initial state of a compositional model from hydrostatic
 *        equilibrium (the EQUIL and ZMFVD keywords).
 *
 * The composition versus depth is given by ZMFVD and the temperature by RTEMPVD
 * (or the constant RTEMP).  The phase pressures are obtained by integrating the
 * hydrostatic ODE with the equation-of-state density, reusing the ODE machinery
 * of the black-oil equilibration facility.  Only the total composition, pressure
 * and temperature are needed downstream: the phase split and the saturations are
 * recomputed by the flash from these quantities.
 *
 * The supported initialization procedures (EQUIL item 10) are
 *  - type 1 (default): ZMFVD provides the total composition and the fluid is
 *    treated as a single phase throughout the column;
 *  - type 3: ZMFVD provides the liquid composition below the gas-oil contact.
 *    The contact acts as the datum, where the pressure is the saturation
 *    (bubble-point) pressure of the contact liquid unless EQUIL item 11
 *    requests the given datum pressure.  Above the contact the gas has the
 *    constant composition of the equilibrium vapour at the contact.
 */
template <class FluidSystem>
class InitialStateComputer
{
public:
    using Scalar = typename FluidSystem::Scalar;
    using FluidState = CompositionalFluidState<Scalar, FluidSystem>;

    /// \param[in] eclipseState    Input state, provides EQUIL, ZMFVD, RTEMP(VD).
    /// \param[in] eosType         Equation of state used by the fluid system.
    /// \param[in] cellCenterDepth Depth of each cell centre.
    /// \param[in] eqlnum          Zero-based equilibration region of each cell.
    /// \param[in] comm            Communicator for parallel runs.
    /// \param[in] gravity         Norm of the gravity vector.
    /// \param[in] numSamplePoints Sample points in each pressure integration.
    InitialStateComputer(const EclipseState& eclipseState,
                         const CompositionalConfig::EOSType eosType,
                         const std::vector<Scalar>& cellCenterDepth,
                         const std::vector<int>& eqlnum,
                         const Parallel::Communication& comm,
                         const Scalar gravity,
                         const int numSamplePoints)
        : eosType_(eosType)
    {
        const auto& records = eclipseState.getInitConfig().getEquil();
        const auto& tables = eclipseState.getTableManager();

        if (!tables.hasTables("ZMFVD")) {
            // COMPVD is the other accepted way of giving the composition
            // versus depth, but it is not supported here yet; name it so the
            // message does not read as if ZMFVD were the only valid input.
            const std::string msg = tables.hasTables("COMPVD")
                ? "Equilibration of a compositional model with the composition versus "
                  "depth from COMPVD is not supported; use ZMFVD instead."
                : "Equilibration of a compositional model requires the composition "
                  "versus depth from the ZMFVD keyword.";
            OPM_THROW(std::runtime_error, msg);
        }

        std::vector<Region> regions;
        regions.reserve(records.size());
        for (std::size_t r = 0; r < records.size(); ++r) {
            regions.push_back(setupRegion(records.getRecord(r), tables, cellCenterDepth,
                                          eqlnum, comm, gravity, numSamplePoints, r));
        }

        fluidStates_.resize(cellCenterDepth.size());
        for (std::size_t cell = 0; cell < cellCenterDepth.size(); ++cell) {
            const auto region = eqlnum[cell];
            if (region < 0 || std::cmp_greater_equal(region, regions.size())) {
                OPM_THROW(std::runtime_error,
                          fmt::format("Cell {} has EQLNUM {} outside the {} "
                                      "equilibration regions.",
                                      cell, region + 1, regions.size()));
            }
            assignCell(fluidStates_[cell], regions[region], cellCenterDepth[cell]);
        }
    }

    std::vector<FluidState>& fluidStates()
    { return fluidStates_; }

    const std::vector<FluidState>& fluidStates() const
    { return fluidStates_; }

private:
    using CompVec = std::array<Scalar, FluidSystem::numComponents>;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;
    using ODE = Details::EosDensityODE<FluidSystem>;
    using PressFunc = EQUIL::Details::PressureFunction<Scalar, ODE>;

    static constexpr int numComponents = FluidSystem::numComponents;

    /// The equilibrated vertical distributions within one region.
    struct Region {
        int initType{1};                            // EQUIL item 10
        Scalar zgoc{};
        CompVec vaporComposition{};                 // gas above the contact (type 3)
        std::vector<TabulatedFunction> zmfVdTable;  // per-component ZMFVD
        TabulatedFunction tempVdTable;
        std::optional<PressFunc> oilPressure;
        std::optional<PressFunc> gasPressure;       // type 3 only
    };

    static CompVec composition(const Region& reg, const Scalar depth)
    {
        CompVec z{};
        Scalar sum = 0.0;
        for (int c = 0; c < numComponents; ++c) {
            z[c] = std::max(Scalar{0}, reg.zmfVdTable[c].eval(depth, /*extrapolate=*/true));
            sum += z[c];
        }
        if (!(sum > 0.0)) {
            OPM_THROW(std::runtime_error,
                      fmt::format("The ZMFVD composition vanishes at depth {} m.", depth));
        }
        std::ranges::transform(z, z.begin(), [sum](const Scalar zc) { return zc / sum; });
        return z;
    }

    /// Whether the ZMFVD composition differs between the top and the bottom of
    /// \p span, i.e. whether the table carries any variation over the region.
    static bool compositionVariesAcross(const Region& reg,
                                        const std::array<Scalar, 2>& span)
    {
        const CompVec top = composition(reg, span[0]);
        const CompVec bottom = composition(reg, span[1]);
        for (int c = 0; c < numComponents; ++c) {
            if (std::abs(top[c] - bottom[c]) > Scalar{1.0e-10}) {
                return true;
            }
        }
        return false;
    }

    Region setupRegion(const EquilRecord& record,
                       const TableManager& tables,
                       const std::vector<Scalar>& cellCenterDepth,
                       const std::vector<int>& eqlnum,
                       const Parallel::Communication& comm,
                       const Scalar gravity,
                       const int numSamplePoints,
                       const std::size_t regionIdx) const
    {
        Region reg;

        reg.initType = record.compositionalInitType();
        if (reg.initType != 1 && reg.initType != 3) {
            OPM_THROW(std::runtime_error,
                      fmt::format("Compositional initialization type {} (EQUIL item 10) is "
                                  "not supported for region {}; only type 1 (total "
                                  "composition) and type 3 (liquid composition) are.",
                                  reg.initType, regionIdx + 1));
        }

        reg.zgoc = record.gasOilContactDepth();

        const auto& zmfvd = tables.getZmfvdTables().template getTable<ZmfvdTable>(regionIdx);
        reg.zmfVdTable.resize(numComponents);
        std::vector<Scalar> depths(zmfvd.getDepthColumn().begin(),
                                   zmfvd.getDepthColumn().end());
        // A single row means a depth-independent composition; the interpolant
        // needs two sample points, so duplicate it onto an arbitrary interval.
        const bool constantComposition = (depths.size() == 1);
        if (constantComposition) {
            depths.push_back(depths.front() + Scalar{1});
        }
        for (int c = 0; c < numComponents; ++c) {
            const auto& col = zmfvd.getMoleFractionColumn(c);
            std::vector<Scalar> values(col.begin(), col.end());
            if (constantComposition) {
                values.push_back(values.front());
            }
            reg.zmfVdTable[c].setXYContainers(depths, values);
        }

        if (tables.hasTables("RTEMPVD")) {
            const auto& rtempvd = tables.getRtempvdTables().template getTable<RtempvdTable>(regionIdx);
            reg.tempVdTable.setXYContainers(rtempvd.getDepthColumn(),
                                            rtempvd.getTemperatureColumn());
        }
        else {
            const std::vector<Scalar> x{0.0, 1.0};
            const std::vector<Scalar> y(2, tables.rtemp());
            reg.tempVdTable.setXYContainers(x, y);
        }

        // Vertical extent of the region's cells across all processes.
        auto span = std::array{std::numeric_limits<Scalar>::max(),
                               std::numeric_limits<Scalar>::lowest()};
        for (std::size_t cell = 0; cell < cellCenterDepth.size(); ++cell) {
            if (std::cmp_equal(eqlnum[cell], regionIdx)) {
                span[0] = std::min(span[0], cellCenterDepth[cell]);
                span[1] = std::max(span[1], cellCenterDepth[cell]);
            }
        }
        span[0] = comm.min(span[0]);
        span[1] = comm.max(span[1]);
        if (span[0] > span[1]) {
            // No cells anywhere in this region.
            return reg;
        }
        if (span[1] - span[0] < Scalar{1}) {
            // Avoid a degenerate integration interval.
            span = {span[0] - Scalar{1}, span[1] + Scalar{1}};
        }

        // The equilibration covers the hydrocarbon column only.
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            if (record.waterOilContactDepth() < span[1]) {
                OPM_THROW(std::runtime_error,
                          fmt::format("Compositional equilibration does not support a water "
                                      "zone: the water-oil contact at {} m is above the "
                                      "deepest cell centre at {} m of region {}.",
                                      record.waterOilContactDepth(), span[1], regionIdx + 1));
            }
            OpmLog::info(fmt::format("Equilibration region {}: the water phase is "
                                     "initialized with zero saturation.", regionIdx + 1));
        }

        if (reg.initType == 1) {
            setupSinglePhaseRegion(reg, record, span, gravity, numSamplePoints, regionIdx);
        }
        else {
            setupTwoPhaseRegion(reg, record, span, gravity, numSamplePoints, regionIdx);
        }

        return reg;
    }

    /// EQUIL item 10 type 1: ZMFVD is the total composition and the fluid is a
    /// single phase, integrated from the datum with the EOS density.  The EOS
    /// root is the vapour one if the datum lies in the gas zone (above the
    /// gas-oil contact) and the liquid one otherwise.
    ///
    /// The datum-versus-contact test is enough to pick the root because of the
    /// convention on the gas-oil contact (EQUIL item 5): it lies above the top
    /// of the reservoir when there is no initial free gas, and below the bottom
    /// when the region holds only gas.  The defaulted item 5 (0 m, i.e. at the
    /// surface) therefore expresses "no free gas" and correctly yields the
    /// liquid root.
    void setupSinglePhaseRegion(Region& reg,
                                const EquilRecord& record,
                                const std::array<Scalar, 2>& span,
                                const Scalar gravity,
                                const int numSamplePoints,
                                const std::size_t regionIdx) const
    {
        const Scalar datum = record.datumDepth();
        const auto phaseIdx = (datum < reg.zgoc)
            ? FluidSystem::gasPhaseIdx : FluidSystem::oilPhaseIdx;

        // Type 1 describes a continuous hydrocarbon phase, i.e. no gas-oil
        // contact in the region.  A contact placed inside the region is only
        // labelled correctly when ZMFVD varies across it, so say so rather
        // than initializing a two-phase column as a single phase in silence.
        if ((reg.zgoc > span[0]) && (reg.zgoc < span[1]) &&
            !compositionVariesAcross(reg, span))
        {
            OpmLog::warning(fmt::format("Equilibration region {}: EQUIL item 10 is 1 "
                                        "(continuous hydrocarbon phase) but the gas-oil "
                                        "contact at {} m lies inside the region, and ZMFVD "
                                        "gives no compositional variation across it. The "
                                        "phases cannot be labelled reliably; supply a "
                                        "varying ZMFVD or use item 10 = 3.",
                                        regionIdx + 1, reg.zgoc));
        }

        const ODE ode([&reg](const Scalar depth) { return composition(reg, depth); },
                      reg.tempVdTable, phaseIdx, eosType_, gravity);
        reg.oilPressure.emplace(ode,
                                typename PressFunc::InitCond{datum, Scalar(record.datumDepthPressure())},
                                numSamplePoints, span);

        OpmLog::info(fmt::format("Equilibration region {}: single phase, total "
                                 "composition specified (EQUIL item 10 is 1).",
                                 regionIdx + 1));
    }

    /// EQUIL item 10 type 3: ZMFVD is the liquid composition.  The pressure at
    /// the gas-oil contact is the saturation pressure of the contact liquid
    /// whenever the given datum pressure disagrees with it by an atmosphere or
    /// more (or unconditionally the datum pressure when EQUIL item 11 is 1),
    /// and the gas above the contact is the equilibrium vapour of the contact
    /// liquid.
    void setupTwoPhaseRegion(Region& reg,
                             const EquilRecord& record,
                             const std::array<Scalar, 2>& span,
                             const Scalar gravity,
                             const int numSamplePoints,
                             const std::size_t regionIdx) const
    {
        if (std::abs(record.datumDepth() - reg.zgoc) > 0.0) {
            OpmLog::warning(fmt::format("Equilibration region {}: the datum depth {} m "
                                        "must be at the gas-oil contact when EQUIL "
                                        "item 10 is 3; using the contact depth {} m.",
                                        regionIdx + 1, record.datumDepth(), reg.zgoc));
        }

        const CompVec liquid = composition(reg, reg.zgoc);
        const Scalar temp = reg.tempVdTable.eval(reg.zgoc, /*extrapolate=*/true);
        Scalar psat{};
        CompVec vapor{};
        if (!SaturationPressure<Scalar, FluidSystem>::bubblePressure(liquid, temp, eosType_,
                                                                     psat, vapor)) {
            OPM_THROW(std::runtime_error,
                      fmt::format("The saturation pressure calculation at the gas-oil "
                                  "contact of region {} did not converge.", regionIdx + 1));
        }
        reg.vaporComposition = vapor;

        // The datum pressure should already equal the saturation pressure at
        // the contact.  With EQUIL item 11 defaulted, the two are required to
        // agree to within one atmosphere and the datum pressure is reset to the
        // computed saturation pressure when they do not; item 11 = 1 keeps the
        // given datum pressure whatever the outcome of that test.
        constexpr Scalar oneAtmosphere = 101325.0;
        const Scalar pDatum = record.datumDepthPressure();
        const bool resetToPsat = record.setToSaturationPressure()
            && (std::abs(pDatum - psat) >= oneAtmosphere);
        const Scalar pGoc = resetToPsat ? psat : pDatum;

        OpmLog::info(fmt::format("Equilibration region {}: two phases, liquid composition "
                                 "specified (EQUIL item 10 is 3). The saturation pressure "
                                 "at the gas-oil contact ({} m) is {:.6g} bar.",
                                 regionIdx + 1, reg.zgoc, psat / 1e5));

        if (resetToPsat) {
            OpmLog::warning(fmt::format("Equilibration region {}: the datum pressure {:.6g} bar "
                                        "differs from the saturation pressure {:.6g} bar at the "
                                        "gas-oil contact by more than one atmosphere; the "
                                        "saturation pressure is used instead.",
                                        regionIdx + 1, pDatum / 1e5, psat / 1e5));
        }

        const ODE oilOde([&reg](const Scalar depth) { return composition(reg, depth); },
                         reg.tempVdTable, FluidSystem::oilPhaseIdx, eosType_, gravity);
        reg.oilPressure.emplace(oilOde,
                                typename PressFunc::InitCond{reg.zgoc, pGoc},
                                numSamplePoints, span);

        const ODE gasOde([vapor](const Scalar) { return vapor; },
                         reg.tempVdTable, FluidSystem::gasPhaseIdx, eosType_, gravity);
        const Scalar pcgoc = record.gasOilContactCapillaryPressure();
        reg.gasPressure.emplace(gasOde,
                                typename PressFunc::InitCond{reg.zgoc, pGoc + pcgoc},
                                numSamplePoints, span);
    }

    void assignCell(FluidState& fs, const Region& reg, const Scalar depth) const
    {
        const bool inGasZone = (reg.initType == 3) && (depth < reg.zgoc);

        const CompVec z = inGasZone ? reg.vaporComposition : composition(reg, depth);
        const auto& pressFunc = inGasZone ? reg.gasPressure : reg.oilPressure;
        if (!pressFunc.has_value()) {
            OPM_THROW(std::runtime_error,
                      "Evaluating the equilibrated pressure of a region without cells.");
        }
        const Scalar press = pressFunc->value(depth);

        fs.setTemperature(reg.tempVdTable.eval(depth, /*extrapolate=*/true));
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                fs.setPressure(phaseIdx, press);
                fs.setSaturation(phaseIdx, 0.0);
            }
        }
        // Nominal single-phase saturation; the flash recomputes the phase split
        // from the total composition, the pressure and the temperature.
        fs.setSaturation(inGasZone ? FluidSystem::gasPhaseIdx : FluidSystem::oilPhaseIdx, 1.0);

        for (int c = 0; c < numComponents; ++c) {
            fs.setMoleFraction(c, z[c]);
        }
    }

    CompositionalConfig::EOSType eosType_;
    std::vector<FluidState> fluidStates_;
};

} // namespace Comp
} // namespace EQUIL
} // namespace Opm

#endif // OPM_INIT_STATE_EQUIL_COMP_HPP
