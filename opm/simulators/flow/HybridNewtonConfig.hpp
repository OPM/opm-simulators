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
#include <unordered_map>
#include <optional>


namespace Opm {

struct Scaler {
    enum class Type { None, Standard, MinMax } type = Type::None;

    // Standard scaler params
    double mean = 0.0;
    double std = 1.0;

    // MinMax scaler params
    double min = 0.0;
    double max = 1.0;

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
            case Type::Standard:
                return scaled_value * std + mean;
            case Type::MinMax:
                return scaled_value * (max - min) + min;
            case Type::None:
            default:
                return scaled_value;
        }
    }
};

struct Transform {
    enum class Type { None, Log, Log10, Log1p } type = Type::None;

    Transform() = default;
    explicit Transform(const std::string& name) {
        if (name == "log10") type = Type::Log10;
        else if (name == "log") type = Type::Log;
        else if (name == "log1p") type = Type::Log1p;
        else type = Type::None;
    }

    double apply(double raw_value) const {
        switch (type) {
            case Type::Log10: return std::log10(raw_value);
            case Type::Log: return std::log(raw_value);
            case Type::Log1p: return std::log1p(raw_value);
            case Type::None:
            default: return raw_value;
        }
    }

    double applyInverse(double transformed_value) const {
        switch (type) {
            case Type::Log10: return std::pow(10.0, transformed_value);
            case Type::Log: return std::exp(transformed_value);
            case Type::Log1p: return std::expm1(transformed_value);
            case Type::None:
            default: return transformed_value;
        }
    }
};


struct FeatureSpec {
    Transform transform;   // replace string transform with Transform struct
    Scaler scaler;

    FeatureSpec() = default;
    FeatureSpec(const std::string& transform_name, const Scaler& scaler_)
      : transform(transform_name), scaler(scaler_) {}
};

struct HybridNewtonConfig {
    std::string model_path;
    std::string cell_indices_file;
    std::vector<double> apply_times;

    std::vector<std::pair<std::string, FeatureSpec>> input_features;
    std::vector<std::pair<std::string, FeatureSpec>> output_features;

    bool hasInputFeature(const std::string& name) const {
        return std::any_of(input_features.begin(), input_features.end(),
                           [&](const auto& p) { return p.first == name; });
    }

    bool hasOutputFeature(const std::string& name) const {
        return std::any_of(output_features.begin(), output_features.end(),
                           [&](const auto& p) { return p.first == name; });
    }
};

} // namespace Opm


#endif