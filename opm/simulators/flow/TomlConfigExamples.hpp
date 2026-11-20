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

#ifndef OPM_TOML_CONFIG_EXAMPLES_HPP
#define OPM_TOML_CONFIG_EXAMPLES_HPP

#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace Opm::Parameters {

struct WriteTomlExample { static constexpr auto value = ""; };

} // namespace Opm::Parameters

namespace Opm {

class PropertyTree;

/// Register the --write-toml-example run-time parameter.
void registerTomlExampleParameters();

/// Names of the example configuration files that may be requested through
/// --write-toml-example, e.g., "linear-solver" or "zoltan-params".
std::vector<std::string> tomlExampleNames();

/// Write documented example TOML configuration files.
///
/// The examples are generated from the current values of the run-time
/// parameters, e.g., --linear-solver, so the parameter system must be fully
/// set up before calling this function.
///
/// \param[in] spec Comma separated list of example names (see
/// tomlExampleNames()), or "all".
///
/// \param[in] directory Directory into which to write the files.  Each
/// example named NAME is written to the file NAME.toml.
///
/// \return Names of the files written.
///
/// \throws std::invalid_argument if \p spec holds an unknown name, if an
/// example cannot be generated from the current parameters, or if one of
/// the output files exists already.  No files are written in this case.
std::vector<std::filesystem::path>
writeTomlExamples(const std::string& spec,
                  const std::filesystem::path& directory = ".");

/// Documented TOML representation of a linear solver configuration.
///
/// \param[in] prm Linear solver configuration, e.g., as created by
/// setupPropertyTree().
///
/// \param[in] optionName Name of the command line option for which the
/// configuration is intended, e.g., "--linear-solver".
///
/// \param[in] origin Description of where the configuration came from,
/// e.g., the name of the preset.
std::string linearSolverTomlExample(const PropertyTree& prm,
                                    const std::string& optionName,
                                    const std::string& origin);

/// Documented TOML representation of Zoltan partitioning parameters.
///
/// \param[in] params Parameters, e.g., as created by setupZoltanParams().
///
/// \param[in] origin Description of where the parameters came from.
///
/// \param[in] imbalanceTol Value of the IMBALANCE_TOL parameter that Flow
/// passes to Zoltan unless overridden.
std::string zoltanParamsTomlExample(const std::map<std::string, std::string>& params,
                                    const std::string& origin,
                                    double imbalanceTol);

/// Documented TOML representation of METIS partitioning parameters.
///
/// \param[in] params Parameters, e.g., as created by setupMetisParams().
///
/// \param[in] origin Description of where the parameters came from.
std::string metisParamsTomlExample(const std::map<std::string, std::string>& params,
                                   const std::string& origin);

/// Documented TOML example of a Hybrid Newton configuration.
std::string hybridNewtonTomlExample();

} // namespace Opm

#endif // OPM_TOML_CONFIG_EXAMPLES_HPP
