/*
  Copyright 2025 Equinor ASA

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

#include <algorithm>
#include <stdexcept>
#include <string>

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/linalgparameters.hh>
#include <opm/simulators/linalg/linalgproperties.hh>


namespace Opm::Parameters
{

/**
 * For the CLI options.
 */
struct LinearSolverAccelerator { static constexpr auto value = "cpu"; };

/**
 * \brief Enum class representing the type of linear solver accelerator.
 */
enum class LinearSolverAcceleratorType { GPU = 0, CPU = 1 };

/**
 * \brief Converts a LinearSolverAcceleratorType to a string representation.
 *
 * \param type The LinearSolverAcceleratorType to convert.
 * \return A string representation of the type.
 */
inline std::string
toString(LinearSolverAcceleratorType type)
{
    switch (type) {
    case LinearSolverAcceleratorType::GPU:
        return "gpu";
    case LinearSolverAcceleratorType::CPU:
        return "cpu";
    default:
        OPM_THROW(std::runtime_error, "Unknown LinearSolverAcceleratorType");
    }
}

/**
 * \brief Converts a string representation to a LinearSolverAcceleratorType.
 *
 * \param str The string to convert.
 * \return The corresponding LinearSolverAcceleratorType.
 * \throws std::runtime_error if the string does not match any known type.
 */
inline LinearSolverAcceleratorType
linearSolverAcceleratorTypeFromString(const std::string& str)
{
    std::string lowerStr = str;
    std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(), [](unsigned char c) { return std::tolower(c); });
    if (lowerStr == "gpu") {
        return LinearSolverAcceleratorType::GPU;
    } else if (lowerStr == "cpu") {
        return LinearSolverAcceleratorType::CPU;
    } else {
        OPM_THROW(std::runtime_error, "Unknown LinearSolverAcceleratorType: " + str);
    }
}

/**
 * \brief Converts the CLI option for linear solver accelerator to a LinearSolverAcceleratorType.
 */
inline LinearSolverAcceleratorType
linearSolverAcceleratorTypeFromCLI()
{
    return linearSolverAcceleratorTypeFromString(Parameters::Get<LinearSolverAccelerator>());
}

} // namespace Opm::Parameters
