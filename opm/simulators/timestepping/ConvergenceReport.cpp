/*
  Copyright 2024 SINTEF AS.

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

#include <opm/simulators/timestepping/ConvergenceReport.hpp>

#include <cassert>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <fmt/format.h>


namespace Opm
{
    std::string to_string(const ConvergenceReport::ReservoirFailure::Type t)
    {
        using Type = ConvergenceReport::ReservoirFailure::Type;

        const auto type_strings = std::unordered_map<Type, std::string> {
            { Type::Invalid    , "Invalid" },
            { Type::MassBalance, "MB"      },
            { Type::Cnv        , "CNV"     },
        };

        auto strPos = type_strings.find(t);
        assert ((strPos != type_strings.end()) &&
                "Unsupported convergence metric type");

        return strPos->second;
    }


    std::string to_string(const ConvergenceReport::Severity s)
    {
        using S = ConvergenceReport::Severity;
        switch (s) {
        case S::None:
            return "None";
        case S::Normal:
            return "Normal";
        case S::TooLarge:
            return "TooLarge";
        case S::NotANumber:
            return "NotANumber";
        }
        throw std::logic_error("Unknown ConvergenceReport::Severity");
    }


    std::string to_string(const ConvergenceReport::WellFailure::Type t)
    {
        using T = ConvergenceReport::WellFailure::Type;
        switch (t) {
        case T::Invalid:
            return "Invalid";
        case T::MassBalance:
            return "MassBalance";
        case T::Pressure:
            return "Pressure";
        case T::ControlBHP:
            return "ControlBHP";
        case T::ControlTHP:
            return "ControlTHP";
        case T::ControlRate:
            return "ControlRate";
        case T::Unsolvable:
            return "Unsolvable";
        case T::WrongFlowDirection:
            return "WrongFlowDirection";
        }
        throw std::logic_error("Unknown ConvergenceReport::WellFailure::Type");
    }


    std::string to_string(const ConvergenceReport::WellFailure& wf)
    {
        std::string mberror = (wf.type() == ConvergenceReport::WellFailure::Type::MassBalance)
            ? fmt::format(" Severity={} Phase={}", to_string(wf.severity()), wf.phase()) : "";
        return fmt::format("{{ {} {}{} }}", wf.wellName(), to_string(wf.type()), mberror);
    }


} // namespace Opm
