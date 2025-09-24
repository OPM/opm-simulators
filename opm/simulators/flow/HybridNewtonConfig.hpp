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

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include <opm/simulators/linalg/PropertyTree.hpp>

namespace Opm {

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

struct FeatureSpec {
    Transform transform;
    Scaler scaler;
    bool is_delta = false;
    std::string actual_name;

    FeatureSpec() = default;
};

class HybridNewtonConfig {
public:
    std::string model_path;
    std::string cell_indices_file;
    std::vector<double> apply_times;
    std::vector<std::pair<std::string, FeatureSpec>> input_features;
    std::vector<std::pair<std::string, FeatureSpec>> output_features;

    HybridNewtonConfig() = default;

    bool hasInputFeature(const std::string& name) const {
        return std::any_of(input_features.begin(), input_features.end(),
                           [&](const auto& p) { return p.first == name; });
    }

    bool hasOutputFeature(const std::string& name) const {
        return std::any_of(output_features.begin(), output_features.end(),
                           [&](const auto& p) { return p.first == name; });
    }

    std::vector<int> loadCellIndicesFromFile() const {
        std::vector<int> indices;
        std::ifstream cellFile(cell_indices_file);
        if (!cellFile.is_open()) {
            throw std::runtime_error("Cannot open cell indices file: " + cell_indices_file);
        }

        std::string line;
        int lineNumber = 0;
        while (std::getline(cellFile, line)) {
            ++lineNumber;
            if (line.empty() || line[0] == '#') continue;
            try { indices.push_back(std::stoi(line)); }
            catch (...) {
                throw std::runtime_error("Invalid cell index at line " + std::to_string(lineNumber) +
                                         " in file " + cell_indices_file + ": " + line);
            }
        }

        if (indices.empty()) {
            throw std::runtime_error("No valid cell indices found in file: " + cell_indices_file);
        }
        return indices;
    }

private:
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
            }

            spec.is_delta = (name.size() >= 6 && std::memcmp(name.data(), "DELTA_", 6) == 0);
            spec.actual_name = spec.is_delta ? name.substr(6) : name; 

            features.emplace_back(name, std::move(spec));
        }
    }

public:
    /*!
     * \brief Parse a JSON config file containing one or multiple HybridNewton configs.
     * \param path Path to the JSON configuration file.
     * \return A vector of HybridNewtonConfig objects, one per model.
     */
    static std::vector<HybridNewtonConfig> parseHybridNewtonConfigs(const std::string& path) {
        std::vector<HybridNewtonConfig> configs;
        PropertyTree pt(path);

        for (const auto& model_key : pt.get_child_keys()) {
            const PropertyTree mt = pt.get_child(model_key);
            HybridNewtonConfig cfg;

            cfg.model_path = mt.get<std::string>("model_path", "");
            if (cfg.model_path.empty()) {
                throw std::runtime_error("Missing 'model_path' for model " + model_key + " in file: " + path);
            }

            cfg.cell_indices_file = mt.get<std::string>("cell_indices_file", "");
            if (cfg.cell_indices_file.empty()) {
                throw std::runtime_error("Missing 'cell_indices_file' for model " + model_key + " in file: " + path);
            }

            auto applyTimesOpt = mt.get_child_optional("apply_times");
            if (!applyTimesOpt) {
                throw std::runtime_error("Missing 'apply_times' for model " + model_key + " in file: " + path);
            }
            for (const auto& key : applyTimesOpt->get_child_keys()) {
                cfg.apply_times.push_back(applyTimesOpt->get_child(key).get<double>(""));
            }

            if (cfg.apply_times.empty()) {
                throw std::runtime_error("'apply_times' must contain at least one value for model " + model_key);
            }

            cfg.parseFeatures(mt, "features.inputs", cfg.input_features);
            cfg.parseFeatures(mt, "features.outputs", cfg.output_features);

            configs.push_back(std::move(cfg));
        }

        if (configs.empty()) {
            throw std::runtime_error("No models found in HybridNewton config file: " + path);
        }

        return configs;
    }
};

} // namespace Opm

#endif
