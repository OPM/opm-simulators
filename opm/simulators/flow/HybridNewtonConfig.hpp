/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Statoil ASA.

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



  struct FeatureSpec {
      std::string unit;
      std::string transform = "none"; // "log", "log10", "log1p", etc.
      Scaler scaler;
  };

  struct HybridNewtonConfig {
      std::string model_path;
      std::vector<int> cell_indices;
      std::vector<double> apply_times;

      std::unordered_map<std::string, FeatureSpec> input_features;
      std::unordered_map<std::string, FeatureSpec> output_features;
  };

} // namespace Opm


#endif