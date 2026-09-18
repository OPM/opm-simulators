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
 * \brief Hydrostatic equilibration for the compositional simulator
 *        (EQUIL + ZMFVD/COMPVD).
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
#include <opm/input/eclipse/EclipseState/Tables/CompvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RtempvdTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SwfnTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableContainer.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/ZmfvdTable.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/flow/equil/PressureFunction.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>


namespace Opm::EQUIL::Comp {

namespace Details {

/// Evaluate a depth table using constant endpoint values outside its tabulated
/// range. Linear extrapolation would introduce composition or temperature
/// values not specified by the input.
template <class Scalar>
Scalar evalDepthTable(const Tabulated1DFunction<Scalar>& table, const Scalar depth)
{
    return table.eval(std::clamp(depth, table.xMin(), table.xMax()));
}

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
    using CompositionFunction = std::function<CompVec(Scalar)>;
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
        const Scalar temp = evalDepthTable(tempVdTable_, depth);

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
        fs.setTemperature(evalDepthTable(tempVdTable_, depth));
        fs.setPressure(FluidSystem::waterPhaseIdx, press);

        typename FluidSystem::template ParameterCache<Scalar> paramCache(eosType_);

        return FluidSystem::density(fs, paramCache, FluidSystem::waterPhaseIdx) * g_;
    }

private:
    const TabulatedFunction& tempVdTable_;
    CompositionalConfig::EOSType eosType_;
    Scalar g_;
};

} // namespace Details

/*!
 * \brief Computes the initial state of a compositional model from hydrostatic
 *        equilibrium (the EQUIL and ZMFVD or COMPVD keywords).
 *
 * The composition versus depth is given by ZMFVD or COMPVD and the temperature
 * by RTEMPVD (or the constant RTEMP). The pressure is obtained by integrating
 * the hydrostatic ODE with the equation-of-state density. The initializer stores
 * the pressure, temperature, overall composition, and nominal phase saturations;
 * the downstream flash recomputes the equilibrium phase split.
 *
 * The supported initialization procedures (EQUIL item 10) are
 *  - type 1 (default): the table provides the total composition. A single-zone
 *    table uses one EOS root throughout the region. When the gas-oil contact
 *    lies inside such a region, the composition must vary over the region for
 *    the downstream flash to label phases correctly. A COMPVD table naming
 *    both phases instead defines separate gas and liquid zones; the zone
 *    containing the datum is anchored at the datum pressure and the other at
 *    the contact;
 *  - type 3: the table provides the liquid composition below the gas-oil
 *    contact. The contact becomes the reference depth, where the pressure
 *    is the saturation (bubble-point) pressure of the contact liquid unless
 *    EQUIL item 11 retains the input pressure. Above the contact, the gas uses
 *    the COMPVD vapour rows when present, or the equilibrium vapour composition
 *    of the contact liquid otherwise.
 *
 * The COMPVD phase column selects the EOS root. If its rows name both phases,
 * the vapour rows define the gas zone above the gas-oil contact and the liquid
 * rows define the liquid zone below it.
 *
 * Only cell-centre initialization is supported (EQUIL item 9 = 0).
 * Gas-oil contact capillary pressure must be zero: the downstream flash uses
 * a single pressure for all phases.
 */
template <class FluidSystem>
class InitialStateComputer
{
public:
    using Scalar = typename FluidSystem::Scalar;
    using FluidState = CompositionalFluidState<Scalar, FluidSystem>;

    /// \param[in] inputState      Input state, provides EQUIL, ZMFVD or COMPVD,
    ///                            RTEMP(VD).
    /// \param[in] eosType         Equation of state used by the fluid system.
    /// \param[in] cellCenterDepth Depth of each cell centre.
    /// \param[in] eqlnum          Zero-based equilibration region of each cell.
    /// \param[in] comm            Communicator for parallel runs.
    /// \param[in] gravity         Norm of the gravity vector.
    /// \param[in] numSamplePoints Sample points in each pressure integration.
    /// \param[in] connateWater   Scaled connate water saturation of each cell,
    ///                            empty when the water phase is inactive.
    /// \param[in] maxWater       Scaled maximum water saturation of each cell,
    ///                            empty when the water phase is inactive.
    InitialStateComputer(const EclipseState& inputState,
                         const CompositionalConfig::EOSType eosType,
                         const std::vector<Scalar>& cellCenterDepth,
                         const std::vector<int>& eqlnum,
                         const Parallel::Communication& comm,
                         const Scalar gravity,
                         const int numSamplePoints,
                         const std::vector<Scalar>& connateWater = {},
                         const std::vector<Scalar>& maxWater = {})
        : eosType_(eosType)
        , connateWater_(connateWater)
        , maxWater_(maxWater)
    {
        const auto& records = inputState.getInitConfig().getEquil();
        const auto& tables = inputState.getTableManager();

        if (!tables.hasTables("ZMFVD") && !tables.hasTables("COMPVD")) {
            OPM_THROW(std::runtime_error,
                      "Equilibration of a compositional model requires the composition "
                      "versus depth from the ZMFVD or the COMPVD keyword.");
        }

        OPM_BEGIN_PARALLEL_TRY_CATCH();
        if (eqlnum.size() != cellCenterDepth.size()) {
            OPM_THROW(std::runtime_error,
                      fmt::format("EQLNUM contains {} entries for {} cell depths.",
                                  eqlnum.size(), cellCenterDepth.size()));
        }
        for (std::size_t cell = 0; cell < eqlnum.size(); ++cell) {
            const auto region = eqlnum[cell];
            if (region < 0 || std::cmp_greater_equal(region, records.size())) {
                OPM_THROW(std::runtime_error,
                          fmt::format("Cell {} has EQLNUM {} outside the {} "
                                      "equilibration regions.",
                                      cell, region + 1, records.size()));
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("Invalid EQLNUM: ", comm);

        std::vector<Region> regions;
        regions.reserve(records.size());
        for (std::size_t r = 0; r < records.size(); ++r) {
            regions.push_back(setupRegion(records.getRecord(r), tables, cellCenterDepth,
                                          eqlnum, comm, gravity, numSamplePoints, r));
        }

        fluidStates_.resize(cellCenterDepth.size());
        referencePressures_.resize(cellCenterDepth.size());
        for (std::size_t cell = 0; cell < cellCenterDepth.size(); ++cell) {
            referencePressures_[cell] =
                assignCell(fluidStates_[cell], regions[eqlnum[cell]], cellCenterDepth[cell], cell);
        }
    }

    std::vector<FluidState>& fluidStates()
    { return fluidStates_; }

    const std::vector<FluidState>& fluidStates() const
    { return fluidStates_; }

    /// Pressure of the water column below the water-oil contact, or of the
    /// hydrocarbon column above it. The current flash uses this common pressure;
    /// fluidStates() retains the independently integrated phase pressures.
    const std::vector<Scalar>& referencePressures() const
    { return referencePressures_; }

private:
    using CompVec = std::array<Scalar, FluidSystem::numComponents>;
    using TabulatedFunction = Tabulated1DFunction<Scalar>;
    using ODE = Details::EosDensityODE<FluidSystem>;
    using WaterODE = Details::WaterDensityODE<FluidSystem>;
    using PressFunc = EQUIL::Details::PressureFunction<Scalar, ODE>;
    using WaterPressFunc = EQUIL::Details::PressureFunction<Scalar, WaterODE>;

    static constexpr int numComponents = FluidSystem::numComponents;

    /// A depth table with a single row is depth-independent, but the
    /// interpolant still needs two sample points; the duplicate row is placed
    /// this far below the original. The distance does not matter.
    static constexpr Scalar constantTableSpan{1.0};

    /// Regions thinner than this are padded by the same amount on either side,
    /// so the pressure integration never runs on a degenerate interval.
    static constexpr Scalar minimumSpanExtent{1.0};

    /// The equilibrated vertical distributions within one region.
    struct Region {
        int initType{1};                            // EQUIL item 10
        Scalar zgoc{};
        /// Name of the selected composition keyword, used in diagnostics.
        std::string_view compositionKeyword;
        /// EOS phase named by every row of a single-phase COMPVD table.
        std::optional<unsigned> statedPhaseIdx{};
        /// Nominal phase assigned to cells in a single-zone region.
        unsigned nominalPhaseIdx{FluidSystem::oilPhaseIdx};
        /// COMPVD naming both phases describes a gas zone over a liquid one,
        /// each with its own composition and its own hydrostatic column.
        bool twoZone{false};
        /// Gas-zone composition for a COMPVD table that names both phases.
        std::vector<TabulatedFunction> vaporVdTable;
        /// Equilibrium vapour at the contact for type 3 without gas-zone rows.
        CompVec vaporComposition{};
        /// Overall composition for a single-zone type 1 region, liquid composition
        /// for type 3, or liquid-zone composition when vaporVdTable stores the gas zone.
        std::vector<TabulatedFunction> compositionVdTable;
        TabulatedFunction tempVdTable;
        std::optional<PressFunc> oilPressure;
        std::optional<PressFunc> gasPressure;       // two zones, or type 3

        Scalar zwoc{};                              // water-oil contact
        std::optional<WaterPressFunc> waterPressure;
    };

    /// Each pressure column must reach the water-oil contact before another
    /// column can be anchored there. Keep the cell span separate for diagnostics.
    static std::array<Scalar, 2> waterContactSpan(const Region& reg,
                                                const std::array<Scalar, 2>& span)
    {
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            return {std::min(span[0], reg.zwoc), std::max(span[1], reg.zwoc)};
        }
        return span;
    }

    /// Normalized vapour-zone composition of a two-zone COMPVD region at a given depth.
    static CompVec vaporComposition(const Region& reg, const Scalar depth)
    {
        CompVec z{};
        Scalar sum = 0.0;
        for (int c = 0; c < numComponents; ++c) {
            z[c] = std::max(Scalar{0}, Details::evalDepthTable(reg.vaporVdTable[c], depth));
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
        Scalar sum{};
        for (int c = 0; c < numComponents; ++c) {
            z[c] = std::max(Scalar{0}, Details::evalDepthTable(reg.compositionVdTable[c], depth));
            sum += z[c];
        }
        if (!(sum > 0.0)) {
            OPM_THROW(std::runtime_error,
                      fmt::format("The composition vanishes at depth {} m.", depth));
        }
        std::ranges::transform(z, z.begin(), [sum](const Scalar zc) { return zc / sum; });
        return z;
    }

    /// Build per-component composition interpolants from the selected \p rows.
    /// ZMFVD and COMPVD expose the same depth and mole-fraction columns.
    template <class Table>
    static void setupComposition(std::vector<TabulatedFunction>& out,
                                 const Table& table,
                                 const std::vector<std::size_t>& rows)
    {
        const auto& depthCol = table.getDepthColumn();

        std::vector<Scalar> depths;
        depths.reserve(rows.size() + 1);
        for (const auto row : rows) {
            depths.push_back(depthCol[row]);
        }
        // A single row means a depth-independent composition; the interpolant
        // needs two sample points, so duplicate it onto an arbitrary interval.
        const bool constantComposition = (depths.size() == 1);
        if (constantComposition) {
            depths.push_back(depths.front() + constantTableSpan);
        }

        out.resize(numComponents);
        for (int c = 0; c < numComponents; ++c) {
            const auto& col = table.getMoleFractionColumn(c);
            std::vector<Scalar> values;
            values.reserve(depths.size());
            for (const auto row : rows) {
                values.push_back(col[row]);
            }
            if (constantComposition) {
                values.push_back(values.front());
            }
            out[c].setXYContainers(depths, values);
        }
    }

    /// Return the region's table index, or the nearest preceding table index.
    /// Return no value when no table exists at or before \p regionIdx.
    static std::optional<std::size_t> sourceTable(const TableContainer& container,
                                                  const std::size_t regionIdx)
    {
        const auto& byIndex = container.tables();
        const auto after = byIndex.upper_bound(regionIdx);
        if (after == byIndex.begin()) {
            return std::nullopt;
        }
        return std::prev(after)->first;
    }

    /// Every row of a table, in order.
    static std::vector<std::size_t> allRows(const std::size_t count)
    {
        std::vector<std::size_t> rows(count);
        std::iota(rows.begin(), rows.end(), std::size_t{0});
        return rows;
    }

    /// Verify that vapour rows are at or above the gas-oil contact and liquid
    /// rows are at or below it.
    static void checkZonesStraddleContact(const CompvdTable& compvd,
                                          const std::vector<std::size_t>& vaporRows,
                                          const std::vector<std::size_t>& liquidRows,
                                          const Scalar zgoc,
                                          const std::size_t regionIdx)
    {
        if (vaporRows.empty() || liquidRows.empty()) {
            OPM_THROW(std::runtime_error,
                      fmt::format("The COMPVD table of region {} names both phases but has "
                                  "no row for one of them.", regionIdx + 1));
        }

        const auto& depth = compvd.getDepthColumn();
        if ((depth[vaporRows.back()] > zgoc) || (depth[liquidRows.front()] < zgoc)) {
            OPM_THROW(std::runtime_error,
                      fmt::format("The COMPVD table of region {} puts its vapour rows down to "
                                  "{} m and its liquid rows from {} m, which do not meet at "
                                  "the gas-oil contact at {} m.",
                                  regionIdx + 1, depth[vaporRows.back()],
                                  depth[liquidRows.front()], zgoc));
        }
    }

    /// The COMPVD rows carrying \p phase.
    static std::vector<std::size_t> rowsOfPhase(const CompvdTable& compvd,
                                                const CompvdTable::Phase phase)
    {
        const auto& flags = compvd.phaseFlags();
        std::vector<std::size_t> rows;
        for (std::size_t r = 0; r < flags.size(); ++r) {
            if (flags[r] == phase) {
                rows.push_back(r);
            }
        }
        return rows;
    }

    /// Return the EOS phase named by every COMPVD row, or no value for absent or
    /// mixed phase flags.
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

    /// Whether the selected composition differs between \p depthA and \p depthB.
    static bool compositionVariesBetween(const Region& reg,
                                         const Scalar depthA,
                                         const Scalar depthB)
    {
        // Input mole fractions that differ physically across the contact should
        // exceed this round-off tolerance.
        constexpr Scalar sameComposition{1.0e-10};

        const CompVec a = composition(reg, depthA);
        const CompVec b = composition(reg, depthB);
        return !std::ranges::equal(a, b, [](const Scalar x, const Scalar y) {
            return std::abs(x - y) <= sameComposition;
        });
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

        if (record.gasOilContactCapillaryPressure() != 0.0) {
            OPM_THROW(std::runtime_error,
                      fmt::format("Compositional equilibration only supports zero gas-oil "
                                  "contact capillary pressure (EQUIL item 6); region {} "
                                  "specifies {} bar.",
                                  regionIdx + 1,
                                  unit::convert::to(record.gasOilContactCapillaryPressure(),
                                                    unit::barsa)));
        }

        if (const auto accuracy = record.initializationTargetAccuracy(); accuracy != 0) {
            OPM_THROW(std::runtime_error,
                      fmt::format("Compositional equilibration only supports cell-centre "
                                  "initialization (EQUIL item 9 = 0); region {} specifies {}.",
                                  regionIdx + 1, accuracy));
        }

        reg.zgoc = record.gasOilContactDepth();

        // With one composition keyword, a missing regional record inherits the
        // nearest preceding record. If both keywords are present, each region
        // must select one explicitly because inheritance would be ambiguous.
        const bool deckHasZmfvd = tables.hasTables("ZMFVD");
        const bool deckHasCompvd = tables.hasTables("COMPVD");
        const bool statesZmfvd = deckHasZmfvd &&
                                 tables.getZmfvdTables().hasTable(regionIdx);
        const bool statesCompvd = deckHasCompvd &&
                                  tables.getCompvdTables().hasTable(regionIdx);

        if (statesZmfvd && statesCompvd) {
            OPM_THROW(std::runtime_error,
                      fmt::format("Region {} has both a ZMFVD and a COMPVD composition "
                                  "versus depth; give only one of them.", regionIdx + 1));
        }

        std::optional<std::size_t> zmfvdAt;
        std::optional<std::size_t> compvdAt;
        if (deckHasZmfvd && deckHasCompvd) {
            if (!statesZmfvd && !statesCompvd) {
                OPM_THROW(std::runtime_error,
                          fmt::format("Region {} has neither a ZMFVD nor a COMPVD "
                                      "composition versus depth. A deck using both "
                                      "keywords has to give every region a record of "
                                      "its own, as neither can be inherited.",
                                      regionIdx + 1));
            }
            (statesZmfvd ? zmfvdAt : compvdAt) = regionIdx;
        }
        else if (deckHasZmfvd) {
            zmfvdAt = sourceTable(tables.getZmfvdTables(), regionIdx);
        }
        else if (deckHasCompvd) {
            compvdAt = sourceTable(tables.getCompvdTables(), regionIdx);
        }

        if (!zmfvdAt.has_value() && !compvdAt.has_value()) {
            OPM_THROW(std::runtime_error,
                      fmt::format("Region {} has neither a ZMFVD nor a COMPVD composition "
                                  "versus depth.", regionIdx + 1));
        }

        if (zmfvdAt.has_value()) {
            reg.compositionKeyword = "ZMFVD";
            const auto& zmfvd =
                tables.getZmfvdTables().template getTable<ZmfvdTable>(*zmfvdAt);
            setupComposition(reg.compositionVdTable, zmfvd,
                             allRows(zmfvd.getDepthColumn().size()));
        }
        else {
            reg.compositionKeyword = "COMPVD";
            const auto& compvd =
                tables.getCompvdTables().template getTable<CompvdTable>(*compvdAt);
            reg.statedPhaseIdx = statedPhase(compvd, regionIdx);

            if (reg.statedPhaseIdx.has_value()) {
                // One phase named: all rows form one composition profile, and
                // the phase flag selects its EOS root.
                if ((reg.initType == 3) &&
                    (*reg.statedPhaseIdx == FluidSystem::gasPhaseIdx)) {
                    OPM_THROW(std::runtime_error,
                              fmt::format("Region {} states a vapour composition in COMPVD "
                                          "while EQUIL item 10 is 3, which takes the liquid "
                                          "composition at the gas-oil contact.",
                                          regionIdx + 1));
                }
                setupComposition(reg.compositionVdTable, compvd,
                                 allRows(compvd.getDepthColumn().size()));
            }
            else {
                // Both phases named: the vapour rows describe the gas zone and
                // the liquid rows the one below the contact.
                const auto vaporRows = rowsOfPhase(compvd, CompvdTable::Phase::Vapor);
                const auto liquidRows = rowsOfPhase(compvd, CompvdTable::Phase::Liquid);
                checkZonesStraddleContact(compvd, vaporRows, liquidRows,
                                          record.gasOilContactDepth(), regionIdx);
                reg.twoZone = true;
                setupComposition(reg.vaporVdTable, compvd, vaporRows);
                setupComposition(reg.compositionVdTable, compvd, liquidRows);
            }
        }

        if (tables.hasTables("RTEMPVD")) {
            const auto& rtempvd =
                tables.getRtempvdTables().template getTable<RtempvdTable>(regionIdx);
            std::vector<Scalar> tempDepths(rtempvd.getDepthColumn().begin(),
                                           rtempvd.getDepthColumn().end());
            const auto& tempCol = rtempvd.getTemperatureColumn();
            std::vector<Scalar> temps(tempCol.begin(), tempCol.end());
            // As for the composition tables, a single row is a depth-independent
            // temperature and the interpolant needs a second sample point.
            if (tempDepths.size() == 1) {
                tempDepths.push_back(tempDepths.front() + constantTableSpan);
                temps.push_back(temps.front());
            }
            reg.tempVdTable.setXYContainers(tempDepths, temps);
        }
        else {
            const std::vector<Scalar> tempDepths{Scalar{0}, constantTableSpan};
            const std::vector<Scalar> temps(tempDepths.size(), tables.rtemp());
            reg.tempVdTable.setXYContainers(tempDepths, temps);
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
        if (span[1] - span[0] < minimumSpanExtent) {
            span = {span[0] - minimumSpanExtent, span[1] + minimumSpanExtent};
        }

        const bool waterActive = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
        if (waterActive) {
            reg.zwoc = record.waterOilContactDepth();
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
            // Type 3 anchors the hydrocarbon at the gas-oil contact on its own
            // saturation pressure, so it has no use for a datum that states the
            // water pressure: the two columns would not meet at the water-oil
            // contact.
            if (datumInWater) {
                OPM_THROW(std::runtime_error,
                          fmt::format("Compositional equilibration of region {} places the "
                                      "datum at {} m, below the water-oil contact at {} m, "
                                      "while EQUIL item 10 is 3. Put the datum in the "
                                      "hydrocarbon column or use item 10 = 1.",
                                      regionIdx + 1, record.datumDepth(), reg.zwoc));
            }
            setupTwoPhaseRegion(reg, record, span, gravity, numSamplePoints, regionIdx);
        }

        // With the datum in the hydrocarbon the water follows from the contact.
        if (waterActive && !datumInWater) {
            setupWaterZone(reg, record, tables, span, gravity, numSamplePoints, regionIdx);
        }
        else if (waterActive) {
            OpmLog::info(fmt::format("Equilibration region {}: the water-oil contact "
                                 "is at {} m.", regionIdx + 1, reg.zwoc));
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

        OpmLog::info(fmt::format("Equilibration region {}: the water-oil contact "
                                 "is at {} m.", regionIdx + 1, reg.zwoc));
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
                                  numSamplePoints, waterContactSpan(reg, span));
    }

    /// A per-cell saturation endpoint, or \p fallback when the caller supplied
    /// none.
    static Scalar waterLimit(const std::vector<Scalar>& limits,
                             const std::size_t cell,
                             const Scalar fallback)
    {
        return limits.empty() ? fallback : limits[cell];
    }

    /// EQUIL item 10 type 1: the table gives the total composition and the
    /// pressure is integrated from the datum with the EOS density; the
    /// subsequent flash determines the phase split.
    ///
    /// A COMPVD table naming both phases describes two zones instead, each with
    /// its own composition, EOS root and hydrostatic column. The column holding
    /// the datum is anchored there and the other picks its pressure up at the
    /// contact, so the datum pressure is honoured whichever zone it lies in.
    ///
    /// For a single zone the EOS root is the phase COMPVD names, or, when the
    /// table does not name one, the root implied by the datum's side of the
    /// gas-oil contact: EQUIL item 5 lies above the top of the reservoir when
    /// there is no initial free gas and below the bottom when the region holds
    /// only gas, so the defaulted 0 m correctly yields the liquid root.
    void setupSinglePhaseRegion(Region& reg,
                                const EquilRecord& record,
                                const std::array<Scalar, 2>& span,
                                const Scalar gravity,
                                const int numSamplePoints,
                                const std::size_t regionIdx,
                                const Scalar datum,
                                const Scalar datumPressure) const
    {
        if (reg.twoZone) {
            setupTwoZoneRegion(reg, span, gravity, numSamplePoints, regionIdx,
                               datum, datumPressure);
            return;
        }

        // COMPVD states the phase its composition belongs to; without that the
        // datum's side of the gas-oil contact decides the EOS root.
        const auto phaseIdx = reg.statedPhaseIdx.value_or(
            (datum < reg.zgoc) ? FluidSystem::gasPhaseIdx : FluidSystem::oilPhaseIdx);

        // A gas-oil contact inside a single-zone region requires the composition
        // to vary across it, so the flash can label the phases correctly.
        if ((reg.zgoc > span[0]) && (reg.zgoc < span[1]) &&
            !compositionVariesBetween(reg, span[0], reg.zgoc) &&
            !compositionVariesBetween(reg, reg.zgoc, span[1])) {
            OpmLog::warning(fmt::format("Equilibration region {}: the gas-oil contact "
                                        "at {} m lies inside a type-1 region, but the "
                                        "composition does not vary across the contact. "
                                        "Compositional variation is required for proper "
                                        "phase labeling.", regionIdx + 1, reg.zgoc));
        }

        const ODE ode([&reg](const Scalar depth) { return composition(reg, depth); },
                      reg.tempVdTable, phaseIdx, eosType_, gravity);
        reg.oilPressure.emplace(ode,
                                typename PressFunc::InitCond{datum, datumPressure},
                                numSamplePoints, waterContactSpan(reg, span));
        reg.nominalPhaseIdx = phaseIdx;

        OpmLog::info(fmt::format("Equilibration region {}: pressure integrated with one "
                                 "EOS root and the total composition from {} "
                                 "(EQUIL item 10 = 1).",
                                 regionIdx + 1, reg.compositionKeyword));
    }

    /// A COMPVD region naming both phases: a gas zone above the gas-oil contact
    /// and a liquid one below, each integrated with its own composition and EOS
    /// root.
    void setupTwoZoneRegion(Region& reg,
                            const std::array<Scalar, 2>& span,
                            const Scalar gravity,
                            const int numSamplePoints,
                            const std::size_t regionIdx,
                            const Scalar datum,
                            const Scalar datumPressure) const
    {
        const ODE liquidOde([&reg](const Scalar depth) { return composition(reg, depth); },
                            reg.tempVdTable, FluidSystem::oilPhaseIdx, eosType_, gravity);
        const ODE gasOde([&reg](const Scalar depth) { return vaporComposition(reg, depth); },
                         reg.tempVdTable, FluidSystem::gasPhaseIdx, eosType_, gravity);

        // The datum-side pressure function must reach the actual contact even
        // when it lies outside the cell span, because the other pressure function
        // is initialized from its value at the contact.
        const auto pressureSpan = waterContactSpan(reg, span);
        const std::array<Scalar, 2> datumSpan{std::min(pressureSpan[0], reg.zgoc),
                                              std::max(pressureSpan[1], reg.zgoc)};
        if ((reg.zgoc < span[0]) || (reg.zgoc > span[1])) {
            OpmLog::warning(fmt::format("Equilibration region {}: the gas-oil contact at {} m "
                                        "lies outside the cells of the region, so the COMPVD "
                                        "rows of one phase describe no cell.",
                                        regionIdx + 1, reg.zgoc));
        }

        if (datum < reg.zgoc) {
            reg.gasPressure.emplace(gasOde,
                                    typename PressFunc::InitCond{datum, datumPressure},
                                    numSamplePoints, datumSpan);
            reg.oilPressure.emplace(liquidOde,
                                    typename PressFunc::InitCond{
                                        reg.zgoc, reg.gasPressure->value(reg.zgoc)},
                                    numSamplePoints, pressureSpan);
        }
        else {
            reg.oilPressure.emplace(liquidOde,
                                    typename PressFunc::InitCond{datum, datumPressure},
                                    numSamplePoints, datumSpan);
            reg.gasPressure.emplace(gasOde,
                                    typename PressFunc::InitCond{
                                        reg.zgoc, reg.oilPressure->value(reg.zgoc)},
                                    numSamplePoints, pressureSpan);
        }

        OpmLog::info(fmt::format("Equilibration region {}: COMPVD gives a gas zone above the "
                                 "contact at {} m and a liquid one below it "
                                 "(EQUIL item 10 = 1).", regionIdx + 1, reg.zgoc));
    }

    /// EQUIL item 10 type 3: the selected table supplies the liquid composition,
    /// and the gas-oil contact is used as the reference depth. If saturation-
    /// pressure adjustment is enabled and the input pressure differs by one
    /// atmosphere or more, use the saturation pressure at the contact. Item 11 = 1
    /// always preserves the input pressure. Above the contact, a two-zone COMPVD
    /// table supplies the gas composition; otherwise the equilibrium vapour of
    /// the contact liquid is used.
    void setupTwoPhaseRegion(Region& reg,
                             const EquilRecord& record,
                             const std::array<Scalar, 2>& span,
                             const Scalar gravity,
                             const int numSamplePoints,
                             const std::size_t regionIdx) const
    {
        const Scalar inputReferenceDepth = record.datumDepth();
        if (inputReferenceDepth != reg.zgoc) {
            OpmLog::warning(fmt::format("Equilibration region {}: the reference depth {} m "
                                        "does not coincide with the gas-oil contact when "
                                        "EQUIL item 10 is 3; resetting it to the contact "
                                        "depth {} m.",
                                        regionIdx + 1, inputReferenceDepth, reg.zgoc));
        }

        const CompVec liquid = composition(reg, reg.zgoc);
        const Scalar temp = Details::evalDepthTable(reg.tempVdTable, reg.zgoc);
        Scalar psat{};
        CompVec vapor{};
        if (!SaturationPressure<Scalar, FluidSystem>::bubblePressure(liquid, temp, eosType_,
                                                                     psat, vapor)) {
            OPM_THROW(std::runtime_error,
                      fmt::format("The saturation pressure calculation at the gas-oil "
                                  "contact of region {} did not converge.", regionIdx + 1));
        }
        reg.vaporComposition = vapor;

        // For type 3, the contact is the reference depth. With item 11
        // defaulted, the input pressure must agree with the saturation pressure
        // to within one atmosphere and is reset to the computed value otherwise.
        // Item 11 = 1 retains the numeric input pressure at the contact regardless
        // of that test; the result need not be an equilibrium system in that case.
        constexpr Scalar oneAtmosphere = unit::atm;
        const Scalar inputPressure = record.datumDepthPressure();
        const bool resetToPsat = record.setToSaturationPressure()
            && (std::abs(inputPressure - psat) >= oneAtmosphere);
        const Scalar referencePressure = resetToPsat ? psat : inputPressure;

        OpmLog::info(fmt::format("Equilibration region {}: two phases, liquid composition "
                                 "specified (EQUIL item 10 is 3). The saturation pressure "
                                 "at the gas-oil contact ({} m) is {:.6g} bar.",
                                 regionIdx + 1, reg.zgoc,
                                 unit::convert::to(psat, unit::barsa)));

        if (resetToPsat) {
            OpmLog::warning(fmt::format("Equilibration region {}: the datum pressure {:.6g} bar "
                                        "differs from the saturation pressure {:.6g} bar at the "
                                        "gas-oil contact by one atmosphere or more; the "
                                        "saturation pressure is used instead.",
                                        regionIdx + 1,
                                        unit::convert::to(inputPressure, unit::barsa),
                                        unit::convert::to(psat, unit::barsa)));
        }

        const ODE oilOde([&reg](const Scalar depth) { return composition(reg, depth); },
                         reg.tempVdTable, FluidSystem::oilPhaseIdx, eosType_, gravity);
        reg.oilPressure.emplace(oilOde,
                                typename PressFunc::InitCond{reg.zgoc, referencePressure},
                                numSamplePoints, waterContactSpan(reg, span));

        // Integrate a two-zone gas column with the same depth-dependent
        // composition assigned to its cells. Otherwise use the equilibrium
        // vapour composition at the contact throughout the gas column.
        typename ODE::CompositionFunction gasComposition;
        if (reg.twoZone) {
            gasComposition = [&reg](const Scalar depth) { return vaporComposition(reg, depth); };
        }
        else {
            gasComposition = [vapor](const Scalar) { return vapor; };
        }
        const ODE gasOde(gasComposition,
                         reg.tempVdTable, FluidSystem::gasPhaseIdx, eosType_, gravity);
        reg.gasPressure.emplace(gasOde,
                                typename PressFunc::InitCond{reg.zgoc, referencePressure},
                                numSamplePoints, waterContactSpan(reg, span));
    }

    Scalar assignCell(FluidState& fs, const Region& reg, const Scalar depth,
                      const std::size_t cell) const
    {
        const bool inGasZone = ((reg.initType == 3) || reg.twoZone) && (depth < reg.zgoc);

        const CompVec z = [&reg, depth, inGasZone]() {
            if (!inGasZone) {
                return composition(reg, depth);
            }
            // Type 3 holds the contact vapour above the contact; a two-zone
            // COMPVD table gives the gas zone its own composition versus depth.
            return reg.twoZone ? vaporComposition(reg, depth) : reg.vaporComposition;
        }();
        const auto& pressFunc = inGasZone ? reg.gasPressure : reg.oilPressure;
        if (!pressFunc.has_value()) {
            OPM_THROW(std::runtime_error,
                      "Evaluating the equilibrated pressure of a region without cells.");
        }

        // The common pressure used by the current flash follows the water
        // column below the contact, including when the scaled maximum water
        // saturation leaves some residual hydrocarbon. Keep this separate from
        // the phase pressures so equilibration retains the capillary offset.
        const Scalar hydrocarbonPressure = pressFunc->value(depth);
        const bool inWaterZone = (depth > reg.zwoc) && reg.waterPressure.has_value();
        const Scalar press = inWaterZone ? reg.waterPressure->value(depth)
                                         : hydrocarbonPressure;

        fs.setTemperature(Details::evalDepthTable(reg.tempVdTable, depth));
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                fs.setPressure(phaseIdx, hydrocarbonPressure);
                fs.setSaturation(phaseIdx, 0.0);
            }
        }

        // Below the water-oil contact the pore space holds water alone; above
        // it the hydrocarbon leaves room for the connate water only.
        Scalar sWat = 0.0;
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            // The saturation function's own endpoints, per cell: below the
            // contact the water fills what it can, above it only the connate
            // water remains.
            sWat = (depth > reg.zwoc) ? waterLimit(maxWater_, cell, Scalar{1})
                                      : waterLimit(connateWater_, cell, Scalar{0});
            fs.setSaturation(FluidSystem::waterPhaseIdx, sWat);
            if (reg.waterPressure.has_value()) {
                fs.setPressure(FluidSystem::waterPhaseIdx, reg.waterPressure->value(depth));
            }
        }

        // Set a nominal single-phase saturation for the hydrocarbon using the
        // phase represented by the pressure integration. The downstream flash
        // recomputes the phase split from composition, pressure, and temperature.
        fs.setSaturation(inGasZone ? FluidSystem::gasPhaseIdx : reg.nominalPhaseIdx,
                         Scalar{1} - sWat);

        for (int c = 0; c < numComponents; ++c) {
            fs.setMoleFraction(c, z[c]);
        }
        return press;
    }

    CompositionalConfig::EOSType eosType_;
    /// Per-cell scaled water saturation endpoints. The saturation functions are
    /// selected by SATNUM and scaled per cell, so neither can be read off the
    /// equilibration region.
    std::vector<Scalar> connateWater_;
    std::vector<Scalar> maxWater_;
    std::vector<FluidState> fluidStates_;
    std::vector<Scalar> referencePressures_;
};

} // namespace Opm::EQUIL::Comp

#endif // OPM_INIT_STATE_EQUIL_COMP_HPP
