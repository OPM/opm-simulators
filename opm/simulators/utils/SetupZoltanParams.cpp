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

#include <config.h>
#include <opm/simulators/utils/SetupZoltanParams.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <boost/property_tree/json_parser.hpp>
#include <boost/version.hpp>

#include <filesystem>

namespace Opm
{


std::map<std::string,std::string> setupZoltanParams(const std::string& conf)
{
    std::map<std::string,std::string> result;

    if (conf == "scotch") {
        result.emplace("LB_METHOD", "GRAPH");
        result.emplace("GRAPH_PACKAGE", "Scotch");
    } else if (conf == "hypergraph") {
        result.emplace("LB_METHOD", "HYPERGRAPH");
    } else if (conf == "graph") {
        result.emplace("LB_METHOD", "GRAPH");
        result.emplace("GRAPH_PACKAGE", "PHG");
    } else if (conf.size() > 5 && conf.substr(conf.size() - 5, 5) == ".json") {
#if BOOST_VERSION / 100 % 1000 > 48
        if ( !std::filesystem::exists(conf) ) {
            OPM_THROW(std::invalid_argument,
                      "JSON file " + conf + " does not exist.");
        }
        boost::property_tree::ptree tree;
        try {
            boost::property_tree::read_json(conf, tree);
        } catch (boost::property_tree::json_parser::json_parser_error& err) {
            OpmLog::error(err.what());
        }
        for (const auto& node : tree) {
            auto value = node.second.get_value_optional<std::string>();
            if (value)
              result.insert_or_assign(node.first, *value);
        }
#else
        OPM_THROW(std::invalid_argument,
                  "--zoltan-params=file.json not supported with "
                  "boost version. Needs version > 1.48.");
#endif
    } else {
        // No valid configuration option found.
        OPM_THROW(std::invalid_argument,
                  conf + " is not a valid setting for --zoltan-params."
                  " Please use graph, hypergraph or scotch");
    }

    return result;
}

} // namespace Opm
