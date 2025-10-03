// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE Research AS

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

#ifndef HYBRID_NEWTON_CONFIG_HPP
#define HYBRID_NEWTON_CONFIG_HPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>

#include <opm/simulators/linalg/PropertyTree.hpp>

namespace Opm {

/*!
 * \brief Represents scaling information for a feature.
 *
 * Supports standard (mean/std) and min-max scaling.
 */
struct Scaler {
    enum class Type { None, Standard, MinMax } type = Type::None;
    double mean = 0.0;
    double std  = 1.0;
    double min  = 0.0;
    double max  = 1.0;

    double scale(double raw_value) const {
        switch (type) {
            case Type::Standard:
                return (std == 0.0) ? mean : (raw_value - mean) / std;
            case Type::MinMax: {
                double denom = max - min;
                return (denom == 0.0) ? min : (raw_value - min) / denom;
            }
            case Type::None:
            default:
                return raw_value;
        }
    }

    double unscale(double scaled_value) const {
        switch (type) {
            case Type::Standard: return scaled_value * std + mean;
            case Type::MinMax:   return scaled_value * (max - min) + min;
            case Type::None:
            default:             return scaled_value;
        }
    }
};

/*!
 * \brief Represents a transformation applied to a feature.
 *
 * Supports log, log10, log1p, and no transform. Provides forward
 * and inverse methods.
 */
struct Transform {
    enum class Type { None, Log, Log10, Log1p } type = Type::None;

    Transform() = default;
    explicit Transform(const std::string& name) {
        if      (name == "log10") type = Type::Log10;
        else if (name == "log")   type = Type::Log;
        else if (name == "log1p") type = Type::Log1p;
        else                       type = Type::None;
    }

    double apply(double raw_value) const {
        switch (type) {
            case Type::Log10: return std::log10(raw_value);
            case Type::Log:   return std::log(raw_value);
            case Type::Log1p: return std::log1p(raw_value);
            case Type::None:
            default:          return raw_value;
        }
    }

    double applyInverse(double transformed_value) const {
        switch (type) {
            case Type::Log10: return std::pow(10.0, transformed_value);
            case Type::Log:   return std::exp(transformed_value);
            case Type::Log1p: return std::expm1(transformed_value);
            case Type::None:
            default:          return transformed_value;
        }
    }
};

/*!
 * \brief Metadata for a single feature (input or output).
 *
 * Includes transformation, scaling, delta flag, and actual feature name.
 */
struct FeatureSpec {
    Transform transform;
    Scaler scaler;
    bool is_delta = false;
    std::string actual_name;

    FeatureSpec() = default;
};

/*!
 * \brief Configuration for a Hybrid Newton ML model.
 *
 * Encapsulates model path, cell indices, apply times, input/output
 * features, and validation. Can be constructed from a PropertyTree.
 */
class HybridNewtonConfig {
public:
    std::string model_path;
    std::string cell_indices_file;
    std::vector<int> cell_indices;
    std::size_t n_cells = 0;
    std::vector<double> apply_times;
    std::vector<std::pair<std::string, FeatureSpec>> input_features;
    std::vector<std::pair<std::string, FeatureSpec>> output_features;

    // Default constructor
    HybridNewtonConfig() = default;

    /*!
    * \brief Construct configuration from a PropertyTree.
    *
    * Loads model path, cell indices, apply times, and features from
    * the provided PropertyTree. Throws on missing or invalid entries.
    */
    explicit HybridNewtonConfig(const PropertyTree& model_config)
    {
        model_path = model_config.get<std::string>("model_path", "");
        if (model_path.empty())
            throw std::runtime_error("Missing 'model_path' in HybridNewton config");

        cell_indices_file = model_config.get<std::string>("cell_indices_file", "");
        if (cell_indices_file.empty())
            throw std::runtime_error("Missing 'cell_indices_file' in HybridNewton config");

        // Load cell indices immediately
        cell_indices = loadCellIndicesFromFile(cell_indices_file);
        n_cells = cell_indices.size();

        // Load apply times
        auto applyTimesOpt = model_config.get_child_optional("apply_times");
        if (!applyTimesOpt)
            throw std::runtime_error("Missing 'apply_times' in HybridNewton config");

        for (const auto& key : applyTimesOpt->get_child_keys())
            apply_times.push_back(applyTimesOpt->get_child(key).get<double>(""));

        if (apply_times.empty())
            throw std::runtime_error("'apply_times' must contain at least one value");

        // Parse features
        parseFeatures(model_config, "features.inputs", input_features);
        parseFeatures(model_config, "features.outputs", output_features);
    }

    bool hasInputFeature(const std::string& name) const {
        return std::any_of(input_features.begin(), input_features.end(),
                           [&](const auto& p) { return p.first == name; });
    }

    bool hasOutputFeature(const std::string& name) const {
        return std::any_of(output_features.begin(), output_features.end(),
                           [&](const auto& p) { return p.first == name; });
    }

    /*!
    * \brief Validate feature compatibility with simulator settings.
    *
    * Ensures RS/RV features are only used if composition support is enabled.
    * Throws std::runtime_error otherwise.
    */
    void validateConfig(bool compositionSwitchEnabled) const
    {
        bool hasRsFeature = hasInputFeature("RS") ||
                            hasOutputFeature("RS") ||
                            hasOutputFeature("DELTA_RS");

        bool hasRvFeature = hasInputFeature("RV") ||
                            hasOutputFeature("RV") ||
                            hasOutputFeature("DELTA_RV");

        if ((hasRsFeature || hasRvFeature) && !compositionSwitchEnabled) {
            OPM_THROW(std::runtime_error,
                "HybridNewton: RS or RV features detected but composition support is disabled. "
                "CompositionSwitch must be enabled for RS/RV features."
            );
        }
    }
    
private:
    /*!
    * \brief Load cell indices from a plain text file.
    *
    * Each line must contain one integer index. Lines starting with '#' or
    * empty lines are ignored. Throws if the file is missing, invalid, or empty.
    */
    std::vector<int> loadCellIndicesFromFile(const std::string& filename) const {
        std::vector<int> indices;
        std::ifstream cellFile(filename);
        if (!cellFile.is_open())
            throw std::runtime_error("Cannot open cell indices file: " + filename);

        std::string line;
        int lineNumber = 0;
        while (std::getline(cellFile, line)) {
            ++lineNumber;
            if (line.empty() || line[0] == '#') continue;
            try { indices.push_back(std::stoi(line)); }
            catch (...) {
                throw std::runtime_error("Invalid cell index at line " + std::to_string(lineNumber) +
                                         " in file " + filename + ": " + line);
            }
        }

        if (indices.empty())
            throw std::runtime_error("No valid cell indices found in file: " + filename);

        return indices;
    }

    /*!
    * \brief Parse feature specifications from a PropertyTree.
    *
    * Reads transform, scaling parameters, and delta flag.
    * Stores the results in the given `features` vector.
    *
    * \param pt       PropertyTree containing the model configuration.
    * \param path     Path to the feature subtree ("features.inputs" or "features.outputs").
    * \param features Destination vector of (name, FeatureSpec) pairs.
    */
    void parseFeatures(const PropertyTree& pt, const std::string& path,
                       std::vector<std::pair<std::string, FeatureSpec>>& features) {
        auto subtreeOpt = pt.get_child_optional(path);
        if (!subtreeOpt) return;

        for (const auto& name : subtreeOpt->get_child_keys()) {
            const PropertyTree& ft = subtreeOpt->get_child(name);
            FeatureSpec spec;
            spec.transform = Transform(ft.get<std::string>("feature_engineering", "none"));

            if (auto sOpt = ft.get_child_optional("scaling_params")) {
                const PropertyTree& s = *sOpt;
                if (s.get_child_optional("mean") && s.get_child_optional("std")) {
                    spec.scaler.type = Scaler::Type::Standard;
                    spec.scaler.mean = s.get<double>("mean", 0.0);
                    spec.scaler.std  = s.get<double>("std", 1.0);
                } else if (s.get_child_optional("min") && s.get_child_optional("max")) {
                    spec.scaler.type = Scaler::Type::MinMax;
                    spec.scaler.min  = s.get<double>("min", 0.0);
                    spec.scaler.max  = s.get<double>("max", 1.0);
                } else {
                    spec.scaler.type = Scaler::Type::None;
                }
            } else {
                spec.scaler.type = Scaler::Type::None;
            }

            spec.is_delta = name.compare(0, 6, "DELTA_") == 0;
            spec.actual_name = spec.is_delta ? name.substr(6) : name;

            features.emplace_back(name, std::move(spec));
        }
    }
};

} // namespace Opm

#endif
