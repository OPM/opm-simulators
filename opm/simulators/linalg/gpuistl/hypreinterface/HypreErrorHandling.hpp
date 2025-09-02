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

#ifndef OPM_HYPRE_ERROR_HANDLING_HPP
#define OPM_HYPRE_ERROR_HANDLING_HPP

#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_utilities.h>
#include <fmt/core.h>
#include <stdexcept>
#include <string>

namespace Opm::gpuistl
{

/**
 * @brief Exception class for Hypre errors
 */
class HypreError : public std::runtime_error
{
public:
    explicit HypreError(const std::string& msg)
        : std::runtime_error(msg)
    {
    }
};

/**
 * @brief Get a descriptive error message for a Hypre error code
 *
 * @param err The Hypre error code
 * @param expression The Hypre expression that caused the error
 * @param file The source file where the error occurred
 * @param function The function where the error occurred
 * @param line The line number where the error occurred
 * @return std::string The formatted error message
 */
inline std::string
getHypreErrorMessage(
    HYPRE_Int err, const std::string& expression, const std::string& file, const std::string& function, int line)
{
    return fmt::format("Hypre error in expression: {}\nError code: {}\nLocation: {}:{} in function {}",
                       expression,
                       err,
                       file,
                       line,
                       function);
}

/**
 * @brief Safe call wrapper for Hypre functions
 *
 * Checks the return code from Hypre functions and throws a HypreError if an error occurred.
 *
 * @param rc The Hypre return code to check
 * @param expression The expression being evaluated (for error reporting)
 * @param file The source file (typically __FILE__)
 * @param function The function name (typically __func__)
 * @param line The line number (typically __LINE__)
 * @throws HypreError if the return code indicates an error
 */
inline void
hypreSafeCall(
    HYPRE_Int rc, const std::string& expression, const std::string& file, const std::string& function, int line)
{
    if (rc != 0) {
        throw HypreError(getHypreErrorMessage(rc, expression, file, function, line));
    }
}

/**
 * @brief Macro to wrap Hypre function calls with error checking
 *
 * Example usage:
 * @code
 * OPM_HYPRE_SAFE_CALL(HYPRE_BoomerAMGCreate(&solver));
 * @endcode
 */
#define OPM_HYPRE_SAFE_CALL(expr) ::Opm::gpuistl::hypreSafeCall((expr), #expr, __FILE__, __func__, __LINE__)

/**
 * @brief Short form macro for Hypre function calls (for backward compatibility)
 */
#ifndef HYPRE_SAFE_CALL
#define HYPRE_SAFE_CALL(expr) OPM_HYPRE_SAFE_CALL(expr)
#endif

} // namespace Opm::gpuistl

#endif // OPM_HYPRE_ERROR_HANDLING_HPP
