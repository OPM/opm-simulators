/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_SETUP_PARTITIONING_PARAMS_HPP
#define OPM_SETUP_PARTITIONING_PARAMS_HPP

#include <map>
#include <optional>
#include <string>

namespace Opm {

/// Form collection of Zoltan partitioning parameters from named configuration
///
/// \param[in] conf Named Zoltan configuration. Must either be the name of a
/// JSON configuration file with the filename extension ".json", or one of
/// the known configuration names
///
///   -* graph Generates configuration parameters for the "GRAPH"
///            load-balancing method, using the "PHG" graph package.
///
///   -* hypergraph Generates configuration parameters for the "HYPERGRAPH"
///            load-balancing method.
///
///   -* scotch Generates configuration parameters for the "GRAPH"
///            load-balancing method, using the "Scotch" graph package.
///
/// \param[in] edgeSizeThreshold Low-level Zoltan partitioning control
/// parameter for when to omit a hyperedge in a hypergraph.  Fraction in the
/// range [0,1] representing a threshold above which to omit discard
/// hyperedges.  Used for conf="graph" and conf="hypergraph".  Nullopt to
/// use the built-in default value.
///
/// \return Collection of Zoltan partitioning parameters.
std::map<std::string, std::string>
setupZoltanParams(const std::string&           conf,
                  const std::optional<double>& edgeSizeThreshold = {});

/// Form collection of METIS partitioning parameters from named configuration
///
/// \param[in] conf Named METIS configuration. Must either be the name of a
/// JSON configuration file with the filename extension ".json", or the
/// known configuration name "default" which uses the built-in default
/// partitioning parameters.
///
/// \return Collection of METIS partitioning parameters.
std::map<std::string, std::string>
setupMetisParams(const std::string& conf);

} // namespace Opm

#endif // OPM_SETUP_PARTITIONING_PARAMS_HPP
