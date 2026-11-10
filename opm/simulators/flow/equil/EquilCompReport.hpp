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
 * \brief Composition versus depth report (RPTSOL SOLVD) of the
 *        compositional equilibration.
 */
#ifndef OPM_EQUIL_COMP_REPORT_HPP
#define OPM_EQUIL_COMP_REPORT_HPP

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <fmt/format.h>

#include <array>
#include <cstddef>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace Opm::EQUIL::Comp {

/// One depth of the composition versus depth report, in SI units.
template <class Scalar, std::size_t NumComponents>
struct CompositionVersusDepthRow
{
    Scalar depth{};
    Scalar poil{};
    Scalar pwat{};
    Scalar pgas{};
    Scalar deno{};
    Scalar denw{};
    Scalar deng{};
    Scalar soil{};
    Scalar swat{};
    Scalar sgas{};
    /// Bubble point of the oil composition, zero when none is reported.
    Scalar psat{};
    /// Saturation pressure that COMPVD states, zero without one.
    Scalar statedPsat{};
    std::array<Scalar, NumComponents> oil{};
    std::array<Scalar, NumComponents> gas{};
    std::array<Scalar, NumComponents> total{};
};

/// The composition versus depth report of one equilibration region as text:
/// a table of pressures, densities, saturations and saturation pressures,
/// followed by the oil, gas and total compositions, one row per depth.
///
/// \param[in] region         One-based equilibration region.
/// \param[in] satTable       One-based saturation table named in the heading.
/// \param[in] rows           Report rows of the region.
/// \param[in] componentNames Name of each component.
/// \param[in] units          Unit system of the report.
template <class Scalar, std::size_t NumComponents>
std::string
formatCompositionVersusDepth(const std::size_t region,
                             const int satTable,
                             const std::vector<CompositionVersusDepthRow<Scalar, NumComponents>>& rows,
                             const std::array<std::string_view, NumComponents>& componentNames,
                             const UnitSystem& units)
{
    using Measure = UnitSystem::measure;
    using Row = CompositionVersusDepthRow<Scalar, NumComponents>;
    std::string report;
    const auto out = std::back_inserter(report);

    fmt::format_to(out, "\nComposition versus depth in equilibration region {} "
                   "(saturation table {})\n\n", region, satTable);
    fmt::format_to(out, "{:>10} {:>10} {:>10} {:>10} {:>9} {:>9} {:>9} "
                   "{:>7} {:>7} {:>7} {:>9} {:>9}\n",
                   "Depth", "Poil", "Pwat", "Pgas", "Deno", "Denw", "Deng",
                   "Soil", "Swat", "Sgas", "Psat", "Ps(obs)");
    const auto length = units.name(Measure::length);
    const auto pressure = units.name(Measure::pressure);
    const auto density = units.name(Measure::density);
    fmt::format_to(out, "{:>10} {:>10} {:>10} {:>10} {:>9} {:>9} {:>9} "
                   "{:>7} {:>7} {:>7} {:>9} {:>9}\n",
                   length, pressure, pressure, pressure, density, density, density,
                   "", "", "", pressure, pressure);
    fmt::format_to(out, "{:->10} {:->10} {:->10} {:->10} {:->9} {:->9} {:->9} "
                   "{:->7} {:->7} {:->7} {:->9} {:->9}\n",
                   "", "", "", "", "", "", "", "", "", "", "", "");
    for (const auto& row : rows) {
        fmt::format_to(out, "{:10.2f} {:10.3f} {:10.3f} {:10.3f} {:9.3f} {:9.3f} {:9.3f} "
                       "{:7.4f} {:7.4f} {:7.4f} {:9.2f} {:9.2f}\n",
                       units.from_si(Measure::length, row.depth),
                       units.from_si(Measure::pressure, row.poil),
                       units.from_si(Measure::pressure, row.pwat),
                       units.from_si(Measure::pressure, row.pgas),
                       units.from_si(Measure::density, row.deno),
                       units.from_si(Measure::density, row.denw),
                       units.from_si(Measure::density, row.deng),
                       row.soil, row.swat, row.sgas,
                       units.from_si(Measure::pressure, row.psat),
                       units.from_si(Measure::pressure, row.statedPsat));
    }

    for (const auto& [phase, fractions] :
             {std::pair{"Oil", &Row::oil},
              std::pair{"Gas", &Row::gas},
              std::pair{"Total", &Row::total}}) {
        fmt::format_to(out, "\n{} composition\n\n{:>10}", phase, "Depth");
        for (const auto name : componentNames) {
            fmt::format_to(out, " {:>8.8}", name);
        }
        fmt::format_to(out, "\n{:>10}\n{:->10}", length, "");
        for (std::size_t c = 0; c < NumComponents; ++c) {
            fmt::format_to(out, " {:->8}", "");
        }
        fmt::format_to(out, "\n");
        for (const auto& row : rows) {
            fmt::format_to(out, "{:10.2f}", units.from_si(Measure::length, row.depth));
            for (const auto fraction : row.*fractions) {
                fmt::format_to(out, " {:8.5f}", fraction);
            }
            fmt::format_to(out, "\n");
        }
    }
    return report;
}

} // namespace Opm::EQUIL::Comp

#endif // OPM_EQUIL_COMP_REPORT_HPP
