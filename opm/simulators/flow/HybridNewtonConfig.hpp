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

#include <cmath>
#include <string>
#include <vector>

namespace Opm {

class PropertyTree;

/*!
 * \brief Represents scaling information for a feature.
 *
 * Supports standard (mean/std) and min-max scaling.
 */
struct Scaler
{
    enum class Type { None, Standard, MinMax } type = Type::None;
    double mean = 0.0;
    double std  = 1.0;
    double min  = 0.0;
    double max  = 1.0;

    double scale(double raw_value) const
    {
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

    double unscale(double scaled_value) const
    {
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
struct Transform
{
    enum class Type { None, Log, Log10, Log1p } type = Type::None;

    Transform() = default;
    explicit Transform(const std::string& name)
    {
        if      (name == "log10") type = Type::Log10;
        else if (name == "log")   type = Type::Log;
        else if (name == "log1p") type = Type::Log1p;
        else                       type = Type::None;
    }

    double apply(double raw_value) const
    {
        switch (type) {
            case Type::Log10: return std::log10(raw_value);
            case Type::Log:   return std::log(raw_value);
            case Type::Log1p: return std::log1p(raw_value);
            case Type::None:
            default:          return raw_value;
        }
    }

    double applyInverse(double transformed_value) const
    {
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
struct FeatureSpec
{
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
class HybridNewtonConfig
{
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
    explicit HybridNewtonConfig(const PropertyTree& model_config);

    bool hasInputFeature(const std::string& name) const;

    bool hasOutputFeature(const std::string& name) const;

    /*!
    * \brief Validate feature compatibility with simulator settings.
    *
    * Ensures RS/RV features are only used if composition support is enabled.
    * Throws std::runtime_error otherwise.
    */
    void validateConfig(bool compositionSwitchEnabled) const;

private:
    /*!
    * \brief Load cell indices from a plain text file.
    *
    * Each line must contain one integer index. Lines starting with '#' or
    * empty lines are ignored. Throws if the file is missing, invalid, or empty.
    */
    std::vector<int> loadCellIndicesFromFile(const std::string& filename) const;

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
                       std::vector<std::pair<std::string, FeatureSpec>>& features);
};

} // namespace Opm

#endif
