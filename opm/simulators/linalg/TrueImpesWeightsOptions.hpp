/*
  Copyright 2026 Equinor ASA

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

#ifndef OPM_TRUE_IMPES_WEIGHTS_OPTIONS_HEADER_INCLUDED
#define OPM_TRUE_IMPES_WEIGHTS_OPTIONS_HEADER_INCLUDED

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <stdexcept>
#include <string>

namespace Opm::Amg
{

/// Options for the numerical true-IMPES CPR weights, read from the
/// preconditioner sub-tree of the linear solver configuration:
///
///   "weight_scaling": "pore_volume" (default) or "max_abs"
struct TrueImpesWeightsOptions
{
    /// How the weights of each cell are scaled.
    enum class Scaling {
        /// Largest weight of each cell has unit magnitude.
        MaxAbs,
        /// The weighted storage of the mass conservation equations equals the
        /// porosity, i.e. the pressure equation is in reservoir volume units in
        /// every cell, as with the analytic weights.  Black-oil models only;
        /// other models, and cells where the weighted storage is not positive,
        /// use MaxAbs.
        PoreVolume,
    };

    Scaling scaling = Scaling::PoreVolume;

    static std::string toString(const Scaling scaling)
    {
        return scaling == Scaling::MaxAbs ? "max_abs" : "pore_volume";
    }

    /// Store the options in a preconditioner sub-tree under the given prefix.
    void put(PropertyTree& prm, const std::string& prefix) const
    {
        prm.put(prefix + "weight_scaling", toString(scaling));
    }

    /// Read the options from a preconditioner sub-tree, using the defaults
    /// for missing entries.
    static TrueImpesWeightsOptions fromPropertyTree(const PropertyTree& prm)
    {
        TrueImpesWeightsOptions options;
        const auto scaling = prm.get<std::string>("weight_scaling", toString(options.scaling));
        if (scaling == "max_abs") {
            options.scaling = Scaling::MaxAbs;
        }
        else if (scaling == "pore_volume") {
            options.scaling = Scaling::PoreVolume;
        }
        else {
            OPM_THROW(std::invalid_argument,
                      "Unknown preconditioner.weight_scaling '" + scaling
                      + "'. Valid values are 'pore_volume' and 'max_abs'.");
        }
        return options;
    }
};

} // namespace Opm::Amg

#endif // OPM_TRUE_IMPES_WEIGHTS_OPTIONS_HEADER_INCLUDED
