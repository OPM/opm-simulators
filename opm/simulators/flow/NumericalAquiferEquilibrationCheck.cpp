/*
  Copyright 2026 Equinor ASA.

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

#include <config.h>
#include <opm/simulators/flow/NumericalAquiferEquilibrationCheck.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/AquiferConfig.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquiferRegionSplit.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquifers.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

namespace {

std::vector<Opm::NumericalAquiferEquilibrationMismatch>
gatherMismatches(const std::vector<Opm::NumericalAquiferEquilibrationMismatch>& local,
                 const Opm::Parallel::Communication& comm)
{
    constexpr int numItems = 4;

    std::vector<double> sendBuffer;
    sendBuffer.reserve(numItems * local.size());
    for (const auto& m : local) {
        sendBuffer.push_back(static_cast<double>(m.aquiferId));
        sendBuffer.push_back(m.aquiferRegion);
        sendBuffer.push_back(m.cellRegion);
        sendBuffer.push_back(m.pressureDifference);
    }

    const auto recvBuffer = Opm::allGatherv(sendBuffer, comm).first;

    std::vector<Opm::NumericalAquiferEquilibrationMismatch> all;
    all.reserve(recvBuffer.size() / numItems);
    for (std::size_t i = 0; i + numItems <= recvBuffer.size(); i += numItems) {
        all.push_back({ static_cast<std::size_t>(recvBuffer[i]),
                        static_cast<int>(recvBuffer[i + 1]),
                        static_cast<int>(recvBuffer[i + 2]),
                        recvBuffer[i + 3] });
    }
    return all;
}

std::string describeMismatches(const std::vector<Opm::NumericalAquiferEquilibrationMismatch>& mismatches,
                               const Opm::UnitSystem& units)
{
    using M = Opm::UnitSystem::measure;

    struct Summary { std::size_t count{0}; double maxDifference{0.0}; };
    std::map<std::tuple<std::size_t, int, int>, Summary> summary;
    for (const auto& m : mismatches) {
        auto& s = summary[{m.aquiferId, m.aquiferRegion, m.cellRegion}];
        ++s.count;
        s.maxDifference = std::max(s.maxDifference, std::abs(m.pressureDifference));
    }

    std::string out;
    for (const auto& [key, s] : summary) {
        const auto& [id, aquiferRegion, cellRegion] = key;
        out += fmt::format("  Aquifer {} (cells in EQLNUM region {}): {} connection(s) to cells in "
                           "EQLNUM region {}, driving an initial flow with pressure differences up to {:.4g} {}\n",
                           id, aquiferRegion, s.count, cellRegion,
                           units.from_si(M::pressure_drop, s.maxDifference),
                           units.name(M::pressure_drop));
    }
    return out;
}

} // Anonymous namespace

namespace Opm {

std::size_t
checkNumericalAquiferEquilibration(const std::vector<NumericalAquiferEquilibrationMismatch>& localMismatches,
                                   const EclipseState& eclState,
                                   const Parallel::Communication& comm,
                                   const bool allowMismatch)
{
    const auto mismatches = gatherMismatches(localMismatches, comm);
    if (mismatches.empty()) {
        return 0;
    }

    if (comm.rank() == 0) {
        std::vector<std::size_t> ids;
        for (const auto& m : mismatches) {
            ids.push_back(m.aquiferId);
        }
        std::ranges::sort(ids);
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

        const auto& units = eclState.getDeckUnitSystem();

        auto message = fmt::format(
            "Numerical aquifer(s) {} connect cells which are initialised from different "
            "equilibration regions (EQLNUM), and which are therefore not in hydrostatic "
            "equilibrium with each other:\n{}"
            "Fluid then flows through the aquifer(s) between these regions from the start "
            "of the simulation, giving an initial transient which may also cause "
            "convergence problems.\n\n",
            fmt::join(ids, ", "), describeMismatches(mismatches, units));

        const auto split = NumericalAquiferRegionSplit {
            eclState.aquifer().numericalAquifers(), ids,
            eclState.getInputGrid(), eclState.globalFieldProps(), eclState.getLgrs()
        };
        message += split.deckChanges(units);

        message += fmt::format(
            "\nTo run the case unchanged, accepting the initial transient, use "
            "--allow-inconsistent-aquifer-equilibration=true.{}",
            allowMismatch ? " This option is set, so the simulation continues." : "");

        if (allowMismatch) {
            OpmLog::warning(message);
        }
        else {
            OpmLog::error(message);
        }
    }

    if (!allowMismatch) {
        // Make sure the report is written before any process stops.
        comm.barrier();
        throw std::domain_error {
            "Numerical aquifers are not in hydrostatic equilibrium "
            "with the cells they are connected to"
        };
    }

    return mismatches.size();
}

} // namespace Opm
