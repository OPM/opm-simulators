/*
  Copyright 2022 Equinor ASA.

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

#ifndef CONVERGENCE_OUTPUT_CONFIGURATION_HPP
#define CONVERGENCE_OUTPUT_CONFIGURATION_HPP

#include <cstddef>
#include <string_view>

namespace Opm {

/// Parse comma separated option strings into a runtime configuration object
/// for whether to output additional convergence information and, if so,
/// what information to output.
///
/// Supported option string values are
///
///   * "none"       -- Dont want any additional convergence output.
///
///   * "steps"      -- Want additional convergence output pertaining to the
///                     converged solution at the end of each timestep.
///
///   * "iterations" -- Want additional convergence output pertaining to each
///                     non-linar ieration in each timestep.
///
/// Option value "none" overrides all other options.  In other words, if the
/// user requests "none", then there will be no additional convergence
/// output, even if there are other options in the option string.
class ConvergenceOutputConfiguration
{
public:
    /// Option values.
    ///
    /// None overrides all other options.  In other words, if the user
    /// requests None, then there will be no additional convergence output,
    /// even if there are other options in the option string.
    enum class Option : unsigned char {
        None = 0,
        Steps = 1 << 1,
        Iterations = 1 << 2,
    };

    /// Constructor
    ///
    /// Parses comma separated option string into runtime configuration
    /// option.
    ///
    /// \param[in] options Comma separated option string.
    ///
    /// \param[in] optionName Name of command line option whose value is \p
    ///    options.  Used as diagnostic information only, and only if
    ///    specified.
    explicit ConvergenceOutputConfiguration(std::string_view options,
                                            std::string_view optionName = "");

    /// Whether or not user wants any additional convergence output at all.
    bool any() const
    {
        return this->flag_ != std::byte{0};
    }

    /// Whether or not user wants specific convergence output.
    ///
    /// \param[in] opt Specific convergence output type.
    bool want(const Option opt) const
    {
        return std::to_integer<int>(this->flag_ & static_cast<std::byte>(opt)) != 0;
    }

private:
    /// Option flags.  Treated as a small bitset.
    std::byte flag_{0};
};

} // namespace Opm

#endif // CONVERGENCE_OUTPUT_CONFIGURATION_HPP
