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

#include <opm/simulators/utils/SetupPartitioningParams.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <boost/version.hpp>

#if BOOST_VERSION / 100 % 1000 > 48
#include <boost/property_tree/json_parser.hpp>
#endif

#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>

#include <fmt/format.h>

namespace {

#if BOOST_VERSION / 100 % 1000 > 48

    std::map<std::string, std::string>
    jsonConfiguration(const std::string&    conf,
                      std::string_view   /* paramName */)
    {
        if (! std::filesystem::exists(conf)) {
            OPM_THROW(std::invalid_argument,
                      fmt::format("JSON file {} does not exist", conf));
        }

        boost::property_tree::ptree tree;
        try {
            boost::property_tree::read_json(conf, tree);
        }
        catch (const boost::property_tree::json_parser::json_parser_error& err) {
            Opm::OpmLog::error(err.what());
        }

        auto result = std::map<std::string, std::string>{};

        for (const auto& node : tree) {
            auto value = node.second.get_value_optional<std::string>();
            if (value) {
                result.insert_or_assign(node.first, *value);
            }
        }

        return result;
    }

#else // ! Boost.PTree has JSON support.

    [[noreturn]] std::map<std::string, std::string>
    jsonConfiguration(const std::string& conf,
                      std::string_view   paramName)
    {
        OPM_THROW(std::invalid_argument,
                  fmt::format("{}=file.json ({}) is not supported in "
                              "current Boost version (1.{}). "
                              "Need Boost version 1.49 or later.",
                              paramName, conf, (BOOST_VERSION / 100) % 1000));
    }

#endif // Boost.PTree has JSON support.

    bool isJsonConfiguration(const std::string& conf)
    {
        const auto ext = std::string_view { ".json" };

        return conf.rfind(ext) == conf.size() - ext.size();
    }

} // Anonymous namespace

// ---------------------------------------------------------------------------

namespace {

    std::map<std::string, std::string>
    zoltanGraphParameters(std::string_view graphPackage)
    {
        return {
            { "LB_METHOD", "GRAPH" },
            { "GRAPH_PACKAGE", graphPackage.data() },
        };
    }

    auto zoltanGraphParameters()
    {
        return zoltanGraphParameters("PHG");
    }

    auto zoltanScotchParameters()
    {
        return zoltanGraphParameters("Scotch");
    }

    std::map<std::string, std::string>
    zoltanHyperGraphParameters()
    {
        return {
            { "LB_METHOD", "HYPERGRAPH" },
        };
    }

} // Anonymous namespace

// ===========================================================================

std::map<std::string, std::string>
Opm::setupZoltanParams(const std::string& conf)
{
    if (conf == "graph") {
        return zoltanGraphParameters();
    }

    else if (conf == "hypergraph") {
        return zoltanHyperGraphParameters();
    }

    else if (conf == "scotch") {
        return zoltanScotchParameters();
    }

    else if (isJsonConfiguration(conf)) {
        return jsonConfiguration(conf, "--zoltan-params");
    }

    else {
        // No valid configuration option found.
        OPM_THROW(std::invalid_argument,
                  fmt::format("{} is not a valid setting for --zoltan-params. "
                              "Please use \"graph\", \"hypergraph\", \"scotch\", "
                              "or a JSON file containing the Zoltan parameters.",
                              conf));
    }
}

std::map<std::string, std::string>
Opm::setupMetisParams(const std::string& conf)
{
    if (conf == "default") {
        return {};
    }

    else if (isJsonConfiguration(conf)) {
        return jsonConfiguration(conf, "--metis-params");
    }

    else {
        // No valid configuration option found.
        OPM_THROW(std::invalid_argument,
                  fmt::format("{} is not a valid setting for --metis-params. "
                              "Please use \"default\" or a JSON file containing "
                              "the METIS parameters.", conf));
    }
}
