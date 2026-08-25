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
#include <opm/input/eclipse/EclipseState/Tables/CompvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SwfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/ZmfvdTable.hpp>

#include <opm/simulators/flow/equil/PressureFunction.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <fmt/format.h>

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


namespace Opm {
namespace EQUIL {
namespace Comp {

namespace Details {

/// Right-hand side of the hydrostatic ODE dp/ddepth = rho(depth, p) * g for a
/// fluid whose density follows from the cubic equation of state at the given
/// The hydrostatic gradient of the water phase, whose density comes from the
/// water PVT rather than from the equation of state.
template <class FluidSystem>
class WaterDensityODE
{
public:
    using Scalar = typename FluidSystem::Scalar;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;

    WaterDensityODE(const TabulatedFunction& tempVdTable,
                    const CompositionalConfig::EOSType eosType,
                    const Scalar normGrav)
        : tempVdTable_(tempVdTable)
        , eosType_(eosType)
        , g_(normGrav)
    {}

    Scalar operator()(const Scalar depth,
                      const Scalar press) const
    {
        CompositionalFluidState<Scalar, FluidSystem> fs;
        fs.setTemperature(tempVdTable_.eval(depth, /*extrapolate=*/true));
        fs.setPressure(FluidSystem::waterPhaseIdx, press);

        typename FluidSystem::template ParameterCache<Scalar> paramCache(eosType_);

        return FluidSystem::density(fs, paramCache, FluidSystem::waterPhaseIdx) * g_;
    }

private:
    const TabulatedFunction& tempVdTable_;
    CompositionalConfig::EOSType eosType_;
    Scalar g_;
};

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

        if (!tables.hasTables("ZMFVD") && !tables.hasTables("COMPVD")) {
            OPM_THROW(std::runtime_error,
                      "Equilibration of a compositional model requires the composition "
                      "versus depth from the ZMFVD or the COMPVD keyword.");
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
    using WaterODE = Details::WaterDensityODE<FluidSystem>;
    using PressFunc = EQUIL::Details::PressureFunction<Scalar, ODE>;
    using WaterPressFunc = EQUIL::Details::PressureFunction<Scalar, WaterODE>;

    static constexpr int numComponents = FluidSystem::numComponents;

    /// The equilibrated vertical distributions within one region.
    struct Region {
        int initType{1};                            // EQUIL item 10
        Scalar zgoc{};
        /// The phase COMPVD states the composition belongs to, when the
        /// composition came from COMPVD.  ZMFVD carries no such column.
        std::optional<unsigned> statedPhaseIdx{};
        /// COMPVD naming both phases describes a gas zone over a liquid one,
        /// each with its own composition and its own hydrostatic column.
        bool twoZone{false};
        std::vector<TabulatedFunction> vaporVdTable;
        CompVec vaporComposition{};                 // gas above the contact (type 3)
        std::vector<TabulatedFunction> zmfVdTable;  // per-component ZMFVD
        TabulatedFunction tempVdTable;
        std::optional<PressFunc> oilPressure;
        std::optional<PressFunc> gasPressure;       // type 3 only

        Scalar zwoc{};                              // water-oil contact
        Scalar connateSw{};                         // water saturation above it
        std::optional<WaterPressFunc> waterPressure;
    };

    /// The vapour-zone composition of a two-zone COMPVD region.
    static CompVec vaporComposition(const Region& reg, const Scalar depth)
    {
        CompVec z{};
        Scalar sum = 0.0;
        for (int c = 0; c < numComponents; ++c) {
            z[c] = std::max(Scalar{0}, reg.vaporVdTable[c].eval(depth, /*extrapolate=*/true));
            sum += z[c];
        }
        if (!(sum > 0.0)) {
            OPM_THROW(std::runtime_error,
                      fmt::format("The COMPVD vapour composition vanishes at depth {} m.", depth));
        }
        std::ranges::transform(z, z.begin(), [sum](const Scalar zc) { return zc / sum; });
        return z;
    }

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

    /// Fills the region's per-component composition interpolants from a
    /// composition-versus-depth table.  ZMFVD and COMPVD differ in their extra
    /// columns, not in the depth and mole-fraction ones this reads.
    template <class Table>
    void setupComposition(Region& reg, const Table& table) const
    {
        reg.zmfVdTable.resize(numComponents);

        std::vector<Scalar> depths(table.getDepthColumn().begin(),
                                   table.getDepthColumn().end());
        // A single row means a depth-independent composition; the interpolant
        // needs two sample points, so duplicate it onto an arbitrary interval.
        const bool constantComposition = (depths.size() == 1);
        if (constantComposition) {
            depths.push_back(depths.front() + Scalar{1});
        }

        for (int c = 0; c < numComponents; ++c) {
            const auto& col = table.getMoleFractionColumn(c);
            std::vector<Scalar> values(col.begin(), col.end());
            if (constantComposition) {
                values.push_back(values.front());
            }
            reg.zmfVdTable[c].setXYContainers(depths, values);
        }
    }

    /// Builds the per-component interpolants from the COMPVD rows carrying
    /// \p phase.
    static void setupCompositionFromRows(std::vector<TabulatedFunction>& table,
                                         const CompvdTable& compvd,
                                         const CompvdTable::Phase phase)
    {
        const auto& flags = compvd.phaseFlags();
        const auto& depthCol = compvd.getDepthColumn();

        std::vector<Scalar> depths;
        for (std::size_t r = 0; r < flags.size(); ++r) {
            if (flags[r] == phase) {
                depths.push_back(depthCol[r]);
            }
        }
        if (depths.empty()) {
            OPM_THROW(std::runtime_error,
                      "COMPVD names both phases but has no row for one of them.");
        }
        const bool single = (depths.size() == 1);
        if (single) {
            depths.push_back(depths.front() + Scalar{1});
        }

        table.resize(numComponents);
        for (int c = 0; c < numComponents; ++c) {
            const auto& col = compvd.getMoleFractionColumn(c);
            std::vector<Scalar> values;
            for (std::size_t r = 0; r < flags.size(); ++r) {
                if (flags[r] == phase) {
                    values.push_back(col[r]);
                }
            }
            if (single) {
                values.push_back(values.front());
            }
            table[c].setXYContainers(depths, values);
        }
    }

    /// The phase COMPVD assigns the composition to.  A table that names both
    /// phases describes a column with a contact in it, which the single-phase
    /// initialization cannot represent, so only a table agreeing on one phase
    /// states one.
    static std::optional<unsigned> statedPhase(const CompvdTable& compvd,
                                               const std::size_t regionIdx)
    {
        const auto& flags = compvd.phaseFlags();
        if (flags.empty()) {
            return std::nullopt;
        }

        const auto first = flags.front();
        if (std::ranges::any_of(flags, [first](const auto f) { return f != first; })) {
            OpmLog::info(fmt::format("Equilibration region {}: COMPVD names both phases, "
                                     "so the composition of each depth decides the phase.",
                                     regionIdx + 1));
            return std::nullopt;
        }

        return (first == CompvdTable::Phase::Vapor)
            ? FluidSystem::gasPhaseIdx : FluidSystem::oilPhaseIdx;
    }

    /// Whether the ZMFVD composition differs between \p depthA and \p depthB.
    static bool compositionVariesBetween(const Region& reg,
                                         const Scalar depthA,
                                         const Scalar depthB)
    {
        const CompVec a = composition(reg, depthA);
        const CompVec b = composition(reg, depthB);
        for (int c = 0; c < numComponents; ++c) {
            if (std::abs(a[c] - b[c]) > Scalar{1.0e-10}) {
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

        if (tables.hasTables("ZMFVD")) {
            const auto& zmfvd = tables.getZmfvdTables().template getTable<ZmfvdTable>(regionIdx);
            setupComposition(reg, zmfvd);
        }
        else {
            const auto& compvd = tables.getCompvdTables().template getTable<CompvdTable>(regionIdx);
            reg.statedPhaseIdx = statedPhase(compvd, regionIdx);
            if (reg.statedPhaseIdx.has_value()) {
                setupComposition(reg, compvd);
            }
            else {
                // Both phases named: the vapour rows describe the gas zone and
                // the liquid rows the one below the contact.
                reg.twoZone = true;
                setupCompositionFromRows(reg.vaporVdTable, compvd, CompvdTable::Phase::Vapor);
                setupCompositionFromRows(reg.zmfVdTable, compvd, CompvdTable::Phase::Liquid);
            }
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

        const bool waterActive = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
        if (waterActive) {
            reg.zwoc = record.waterOilContactDepth();
            reg.connateSw = connateWaterSaturation(tables, regionIdx);
        }

        // A datum below the water-oil contact states the pressure of the water
        // rather than of the hydrocarbon: integrate the water from there and
        // hand the hydrocarbon its pressure at the contact.  Otherwise the
        // hydrocarbon is anchored at the datum and the water follows from the
        // contact.
        Scalar hcDatum = record.datumDepth();
        Scalar hcPressure = record.datumDepthPressure();
        const bool datumInWater = waterActive && (record.datumDepth() > reg.zwoc);

        if (datumInWater) {
            integrateWaterPressure(reg, span, gravity, numSamplePoints,
                                   record.datumDepth(), record.datumDepthPressure());
            hcDatum = reg.zwoc;
            hcPressure = reg.waterPressure->value(reg.zwoc)
                       + record.waterOilContactCapillaryPressure();
            OpmLog::info(fmt::format("Equilibration region {}: the datum at {} m lies below the "
                                     "water-oil contact at {} m, so it gives the water pressure; "
                                     "the hydrocarbon pressure at the contact is {:.5} bar.",
                                     regionIdx + 1, record.datumDepth(), reg.zwoc,
                                     hcPressure / 1.0e5));
        }

        if (reg.initType == 1) {
            setupSinglePhaseRegion(reg, record, span, gravity, numSamplePoints, regionIdx,
                                   hcDatum, hcPressure);
        }
        else {
            setupTwoPhaseRegion(reg, record, span, gravity, numSamplePoints, regionIdx);
        }

        // With the datum in the hydrocarbon the water follows from the contact.
        if (waterActive && !datumInWater) {
            setupWaterZone(reg, record, tables, span, gravity, numSamplePoints, regionIdx);
        }
        else if (waterActive) {
            OpmLog::info(fmt::format("Equilibration region {}: the water-oil contact is at {} m, "
                                     "with a connate water saturation of {} above it.",
                                     regionIdx + 1, reg.zwoc, reg.connateSw));
        }

        return reg;
    }

    /// The water phase of a region: connate above the water-oil contact, fully
    /// water-saturated below it, with its own hydrostatic pressure.
    ///
    /// The water pressure is integrated from the contact rather than derived
    /// from the hydrocarbon pressure: the two only agree there, and away from
    /// it the water gradient is the steeper one.
    void setupWaterZone(Region& reg,
                        const EquilRecord& record,
                        const TableManager& tables,
                        const std::array<Scalar, 2>& span,
                        const Scalar gravity,
                        const int numSamplePoints,
                        const std::size_t regionIdx) const
    {
        // The capillary pressure at the contact (EQUIL item 4) offsets the
        // water pressure from the hydrocarbon pressure there.
        if (!reg.oilPressure.has_value() && !reg.gasPressure.has_value()) {
            return;
        }
        const auto& hcPressure = reg.oilPressure.has_value() ? reg.oilPressure : reg.gasPressure;
        const Scalar pcow = record.waterOilContactCapillaryPressure();
        const Scalar pContact = hcPressure->value(reg.zwoc) - pcow;

        integrateWaterPressure(reg, span, gravity, numSamplePoints, reg.zwoc, pContact);

        OpmLog::info(fmt::format("Equilibration region {}: the water-oil contact is at {} m, "
                                 "with a connate water saturation of {} above it.",
                                 regionIdx + 1, reg.zwoc, reg.connateSw));
    }

    /// Integrates the water pressure over \p span from \p depth, where it is
    /// \p pressure.
    void integrateWaterPressure(Region& reg,
                                const std::array<Scalar, 2>& span,
                                const Scalar gravity,
                                const int numSamplePoints,
                                const Scalar depth,
                                const Scalar pressure) const
    {
        const WaterODE ode(reg.tempVdTable, eosType_, gravity);
        reg.waterPressure.emplace(ode,
                                  typename WaterPressFunc::InitCond{depth, pressure},
                                  numSamplePoints, span);
    }

    /// The smallest water saturation the saturation function defines, i.e. the
    /// connate water left in the hydrocarbon column.
    static Scalar connateWaterSaturation(const TableManager& tables,
                                         const std::size_t regionIdx)
    {
        const auto& swfn = tables.getSwfnTables();
        if (!swfn.empty()) {
            const auto& table = swfn.template getTable<SwfnTable>(
                std::min(regionIdx, swfn.size() - 1));
            const auto& sw = table.getSwColumn();
            if (sw.size() > 0) {
                return sw.front();
            }
        }
        return Scalar{0};
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
                                const std::size_t regionIdx,
                                const Scalar datum,
                                const Scalar datumPressure) const
    {
        // COMPVD states the phase its composition belongs to; without that the
        // datum's side of the gas-oil contact decides the EOS root.
        const auto phaseIdx = reg.statedPhaseIdx.value_or(
            (datum < reg.zgoc) ? FluidSystem::gasPhaseIdx : FluidSystem::oilPhaseIdx);

        // Type 1 describes a continuous hydrocarbon phase, i.e. no gas-oil
        // contact in the region, and the whole column is integrated with the
        // single EOS root chosen from the datum side of the contact.  Warn on
        // a contact inside the region rather than initialize in silence: on
        // the far side of it that root can be the wrong one of a three-root
        // cubic, and a genuinely two-phase column needs item 10 = 3.
        if ((reg.zgoc > span[0]) && (reg.zgoc < span[1])) {
            if (compositionVariesBetween(reg, span[0], reg.zgoc) ||
                compositionVariesBetween(reg, reg.zgoc, span[1])) {
                OpmLog::warning(fmt::format("Equilibration region {}: EQUIL item 10 is 1 "
                                            "(continuous hydrocarbon phase) but the gas-oil "
                                            "contact at {} m lies inside the region. The "
                                            "gas/liquid EOS root is fixed from the datum "
                                            "side, so densities across the contact may use "
                                            "the wrong root; use item 10 = 3 for a "
                                            "two-phase column.",
                                            regionIdx + 1, reg.zgoc));
            }
            else {
                OpmLog::warning(fmt::format("Equilibration region {}: EQUIL item 10 is 1 "
                                            "(continuous hydrocarbon phase) but the gas-oil "
                                            "contact at {} m lies inside the region, and ZMFVD "
                                            "gives no compositional variation across it. The "
                                            "phases cannot be labelled reliably; supply a "
                                            "varying ZMFVD or use item 10 = 3.",
                                            regionIdx + 1, reg.zgoc));
            }
        }

        const ODE ode([&reg](const Scalar depth) { return composition(reg, depth); },
                      reg.tempVdTable, phaseIdx, eosType_, gravity);
        reg.oilPressure.emplace(ode,
                                typename PressFunc::InitCond{datum, datumPressure},
                                numSamplePoints, span);

        if (reg.twoZone) {
            // The gas zone is its own hydrostatic column, continuous with the
            // liquid one at the contact.
            const ODE gasOde([&reg](const Scalar d) { return vaporComposition(reg, d); },
                             reg.tempVdTable, FluidSystem::gasPhaseIdx, eosType_, gravity);
            reg.gasPressure.emplace(gasOde,
                                    typename PressFunc::InitCond{reg.zgoc,
                                                                 reg.oilPressure->value(reg.zgoc)},
                                    numSamplePoints, span);
            OpmLog::info(fmt::format("Equilibration region {}: COMPVD gives a gas zone above the "
                                     "contact at {} m and a liquid one below it.",
                                     regionIdx + 1, reg.zgoc));
        }

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
        const bool inGasZone = ((reg.initType == 3) || reg.twoZone) && (depth < reg.zgoc);

        const CompVec z = !inGasZone      ? composition(reg, depth)
                        : reg.twoZone     ? vaporComposition(reg, depth)
                                          : reg.vaporComposition;
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

        // Below the water-oil contact the pore space holds water alone; above
        // it the hydrocarbon leaves room for the connate water only.
        Scalar sWat = 0.0;
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            sWat = (depth > reg.zwoc) ? Scalar{1} : reg.connateSw;
            fs.setSaturation(FluidSystem::waterPhaseIdx, sWat);
            if (reg.waterPressure.has_value()) {
                fs.setPressure(FluidSystem::waterPhaseIdx, reg.waterPressure->value(depth));
            }
        }

        // Nominal single-phase saturation for the hydrocarbon; the flash
        // recomputes the phase split from the total composition, the pressure
        // and the temperature.
        //
        // The phase is the one the EOS root describes, so the liquid zone of a
        // two-zone region is labelled oil.  The reference simulator reports the
        // hydrocarbon of a single-phase region as gas throughout, whatever its
        // root, which is a reporting convention rather than a statement about
        // the fluid; the saturations here name the phase that is present.
        fs.setSaturation(inGasZone ? FluidSystem::gasPhaseIdx : FluidSystem::oilPhaseIdx,
                         Scalar{1} - sWat);

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
