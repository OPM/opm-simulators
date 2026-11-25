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

#ifndef OPM_NON_RETRYABLE_ERROR_HPP
#define OPM_NON_RETRYABLE_ERROR_HPP

#include <stdexcept>
#include <string>

namespace Opm {

/// A failure that a shorter time step cannot cure (e.g. exhausted reserved
/// capacity); adaptive time stepping rethrows it instead of chopping.
class NonRetryableError : public std::runtime_error
{
public:
    explicit NonRetryableError(const std::string& message)
        : std::runtime_error(message)
    {}
};

} // namespace Opm

#endif // OPM_NON_RETRYABLE_ERROR_HPP
