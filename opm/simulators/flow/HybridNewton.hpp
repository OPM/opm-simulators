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

#ifndef HYBRID_NEWTON_HPP
#define HYBRID_NEWTON_HPP


#include <opm/simulators/flow/FlowBaseProblemProperties.hpp>
#include <opm/simulators/flow/HybridNewtonConfig.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/ml/ml_model.hpp>

#include <string>         
#include <iostream>
#include <fstream>
#include <boost/property_tree/json_parser.hpp>
#include <vector>

namespace Opm {

/*!
 * \brief Hybrid Newton solver extension for the black-oil model.
 *
 * This class integrates machine learning–based corrections into the
 * standard Newton solver of the black-oil simulator. It uses one or
 * more trained Hybrid Newton models to adjust the initial guess of the
 * nonlinear solver at specific timesteps.
 *
 * Each model is described by a HybridNewtonConfig, which specifies:
 *   - the path to the model file,
 *   - the cells to which the model applies,
 *   - the time points when it should be applied,
 *   - the input and output features with their transformations and scalings.
 *
 * At runtime, the class:
 *   1. Loads the relevant configuration(s),
 *   2. Constructs an input tensor [n_cells x n_input_features],
 *   3. Runs the ML model to produce an output tensor [n_cells x n_output_features],
 *   4. Updates the nonlinear solver’s initial guess accordingly.
 */
template <typename TypeTag>
class BlackOilHybridNewton
{
protected:
    using Simulator  = GetPropType<TypeTag, Properties::Simulator>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { numPhases = FluidSystem::numPhases };

    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    static constexpr bool compositionSwitchEnabled = Indices::compositionSwitchIdx >= 0;

public:
    explicit BlackOilHybridNewton(Simulator& simulator)
        : simulator_(simulator)
        , configsLoaded_(false)
    {}
    void tryApplyHybridNewton()
    {   
        // Check if flag activated 
        if (!Parameters::Get<Parameters::UseHyNe>())
            return;
        
        validateFluidSystem();

        // Load configs if not yet loaded
        if (!configsLoaded_) {
            configs_ = parseHybridNewtonConfigs(Parameters::Get<Parameters::HyNeConfigFile>());  // note: plural
            configsLoaded_ = true;
            
            // Validate all configs
            validateAllConfigs();
        }

        double current_time = simulator_.time();
        // Find and apply all models that should run at this time
        auto applicable_models = getApplicableModels(current_time);
        
        if (!applicable_models.empty()) {
            for (const auto& config : applicable_models) {
                runHybridNewton(config);
            }
        }
    }

private:
    /*!
    * \brief Apply the Hybrid Newton method using a given model configuration.
    *
    * This function loads the cell indices on demand from the configuration's
    * `cell_indices_file`, constructs the input tensor, evaluates the model to
    * produce the output tensor, and applies the result to update the initial
    * guess of the nonlinear solver.
    *
    * \param config The HybridNewtonConfig object containing model path,
    *        feature specifications, and cell indices file reference.
    *
    * \throws std::runtime_error if cell indices cannot be loaded or if
    *         feature extraction or application fails.
    */
    void runHybridNewton(const HybridNewtonConfig& config)
    {
        // Load cell indices from file at runtime
        std::vector<int> cell_indices = loadCellIndicesFromFile(config.cell_indices_file);
        const size_t n_cells = cell_indices.size();

        // Pass cells to each stage
        auto input  = constructInputTensor(config, cell_indices, n_cells);
        auto output = constructOutputTensor(input, config, n_cells);
        updateInitialGuess(output, config, cell_indices, n_cells);
    }



protected:

    void validateFluidSystem()
    {
        const auto& eclState = simulator_.vanguard().eclState();
        const auto& phases = eclState.runspec().phases();
        
        bool hasWater = phases.active(Phase::WATER);
        bool hasGas = phases.active(Phase::GAS);
        bool hasOil = phases.active(Phase::OIL);
        bool hasSolvent = phases.active(Phase::SOLVENT);

        // Check for three-phase black oil system
        if (!(hasWater && hasOil && hasGas && !hasSolvent)) {
            OPM_THROW(std::runtime_error, 
                "HybridNewton: Unsupported fluid system. Only three-phase black oil is supported.");
        }
    }

    std::vector<HybridNewtonConfig> getApplicableModels(double current_time)
    {
        std::vector<HybridNewtonConfig> applicable;
        
        for (const auto& config : configs_) {
            if (shouldApplyHybridNewton(current_time, config)) {
                applicable.push_back(config);
            }
        }
        
        return applicable;
    }

    bool shouldApplyHybridNewton(double current_time, const HybridNewtonConfig& config)
    {
        if (config.apply_times.empty()) {
            return false;
        }

        if (config.apply_times.size() == 1) {
            // Apply exactly at one point in time (with optional tolerance)
            constexpr double tolerance = 1e-6;
            return std::abs(current_time - config.apply_times[0]) < tolerance;
        }

        if (config.apply_times.size() == 2) {
            double start_time = config.apply_times[0];
            double end_time = config.apply_times[1];
            return (current_time >= start_time) && (current_time <= end_time);
        }

        // If apply_times has an unexpected number of entries
        return false;
    }

    /*!
    * \brief Parse the Hybrid Newton configuration file (JSON format).
    *
    * Reads a JSON configuration file containing one or more model
    * specifications. Each specification must include the model path,
    * the path to a file containing cell indices (`cell_indices_file`),
    * the application times, and input/output feature descriptions.
    *
    * \param path Path to the JSON configuration file.
    * \return A vector of HybridNewtonConfig objects populated from the file.
    *
    * \throws std::runtime_error if the file cannot be opened, parsed,
    *         or if no valid model configurations are found.
    */
    std::vector<HybridNewtonConfig> parseHybridNewtonConfigs(const std::string& path)
    {
        std::vector<HybridNewtonConfig> configs;
        boost::property_tree::ptree pt;
        std::ifstream file(path);
        
        if (!file.is_open()) {
            OPM_THROW(std::runtime_error, "Cannot open config file: " + path);
        }
        
        try {
            boost::property_tree::read_json(file, pt);
        } catch (const boost::property_tree::json_parser_error& e) {
            OPM_THROW(std::runtime_error, "Failed to parse JSON: " + std::string(e.what()));
        }
        
        auto parseFeatures = [](const boost::property_tree::ptree& subtree) {
            std::vector<std::pair<std::string, FeatureSpec>> result;
            for (const auto& [key, val] : subtree) {
                FeatureSpec spec;
                spec.transform = Transform(val.get<std::string>("feature_engineering", "none"));
                
                if (auto s = val.get_child_optional("scaling_params")) {
                    if (s->get_child_optional("mean") && s->get_child_optional("std")) {
                        spec.scaler.type = Scaler::Type::Standard;
                        spec.scaler.mean = s->get<double>("mean", 0.0);
                        spec.scaler.std = s->get<double>("std", 1.0);
                    } else if (s->get_child_optional("min") && s->get_child_optional("max")) {
                        spec.scaler.type = Scaler::Type::MinMax;
                        spec.scaler.min = s->get<double>("min", 0.0);
                        spec.scaler.max = s->get<double>("max", 1.0);
                    } else {
                        spec.scaler.type = Scaler::Type::None;
                    }
                } else {
                    spec.scaler.type = Scaler::Type::None;
                }
                
                result.emplace_back(key, std::move(spec));
            }
            return result;
        };
        
        for (const auto& [index, modelTree] : pt) {
            HybridNewtonConfig config;
            
            config.model_path = modelTree.get<std::string>("model_path", "");
            
            // Require cell_indices to be a file path (string)
            config.cell_indices_file = modelTree.get<std::string>("cell_indices_file", "");
            if (config.cell_indices_file.empty()) {
                OPM_THROW(std::runtime_error, 
                    "Model must have 'cell_indices' as a file path (string)");
            }
            
            auto applyTimesTree = modelTree.get_child("apply_times", {});
            for (const auto& item : applyTimesTree) {
                config.apply_times.push_back(item.second.get_value<double>());
            }
            
            config.input_features = parseFeatures(modelTree.get_child("features.inputs"));
            config.output_features = parseFeatures(modelTree.get_child("features.outputs"));
            
            configs.push_back(std::move(config));
        }
        
        if (configs.empty()) {
            OPM_THROW(std::runtime_error, "No models found in config file: " + path);
        }
        
        return configs;
    }

    
    void validateAllConfigs()
    {
        for (const auto& config : configs_) {
            // Check if RS or RV features appear in this config
            bool hasRsFeature = 
                config.hasInputFeature("RS") ||
                config.hasOutputFeature("RS") ||
                config.hasOutputFeature("DELTA_RS");
            
            bool hasRvFeature =
                config.hasInputFeature("RV") ||
                config.hasOutputFeature("RV") ||
                config.hasOutputFeature("DELTA_RV");

            // If RS or RV appear, ensure compositionSwitchEnabled is true
            if ((hasRsFeature || hasRvFeature) && !compositionSwitchEnabled) {
                OPM_THROW(std::runtime_error,
                    "HybridNewton: RS or RV features detected but compositionSwitchEnabled is false. "
                    "Composition support required for RS/RV."
                );
            }
        }
    }

    /*!
    * \brief Construct the input feature tensor for the Hybrid Newton model.
    *
    * For each cell listed in \p cell_indices, this function extracts
    * all input features defined in the configuration, applies the
    * necessary transformations and scaling, and assembles them into
    * a row-major 2D tensor.
    *
    * \param config The HybridNewtonConfig containing feature definitions.
    * \param cell_indices The list of cell indices to process.
    * \param n_cells The number of cells (typically cell_indices.size()).
    * \return A tensor of shape [n_cells x n_input_features], where rows
    *         correspond to cells and columns correspond to input features.
    */
    Opm::ML::Tensor<Evaluation> 
    constructInputTensor(const HybridNewtonConfig& config, const std::vector<int>& cell_indices, size_t n_cells)
    {
        const auto& features = config.input_features;

        // Count scalar and per-cell features
        int num_scalar_features = 0;
        int num_per_cell_features = 0;
        for (const auto& [name, spec] : features) {
            if (name == "TIMESTEP") {
                ++num_scalar_features;
            } else {
                ++num_per_cell_features;
            }
        }

        const auto& unitSyst = simulator_.vanguard().schedule().getUnits();
        // Calculate total input size (feature-major)
        size_t input_size = num_per_cell_features * n_cells + num_scalar_features;
        Opm::ML::Tensor<Evaluation> input(input_size);

        // Track offset for each feature in flat tensor
        size_t offset = 0;

        for (const auto& [name, spec] : features) {
            if (name == "TIMESTEP") {
                // Scalar feature: assign single value at offset
                double raw_value = simulator_.timeStepSize();
                raw_value = spec.transform.apply(raw_value);
                double scaled_value = spec.scaler.scale(raw_value);
                input(offset) = scaled_value;

                offset += 1; // advance offset by 1 for scalar
            } else {
                // Per-cell feature: assign values for each cell
                for (size_t cell_idx = 0; cell_idx < n_cells; ++cell_idx) {
                    const auto& intQuants = simulator_.model().intensiveQuantities(cell_indices[cell_idx], 0);
                    const auto& fs = intQuants.fluidState();

                    double raw_value = 0.0;

                    if (name == "PRESSURE") {
                        raw_value = getValue(fs.pressure(oilPhaseIdx));
                        raw_value = unitSyst.from_si(UnitSystem::measure::pressure, raw_value); 
                    } else if (name == "SWAT") {
                        raw_value = getValue(fs.saturation(waterPhaseIdx));
                    } else if (name == "SGAS") {
                        raw_value = getValue(fs.saturation(gasPhaseIdx));
                    } else if (name == "SOIL") {
                        raw_value = getValue(fs.saturation(oilPhaseIdx));
                    } else if (name == "RS") {
                        raw_value = getValue(fs.Rs());
                        raw_value = unitSyst.from_si(UnitSystem::measure::gas_oil_ratio, raw_value); 
                    } else if (name == "RV") {
                        raw_value = getValue(fs.Rv());
                        raw_value = unitSyst.from_si(UnitSystem::measure::oil_gas_ratio, raw_value); 
                    } else if (name == "PERMX") {
                        const auto& eclState = simulator_.vanguard().eclState();
                        const auto& fp = eclState.fieldProps();
                        auto permX = fp.get_double("PERMX");
                        raw_value = permX[cell_indices[cell_idx]];
                        raw_value = unitSyst.from_si(UnitSystem::measure::permeability, raw_value);
                    } else {
                        OPM_THROW(std::runtime_error, "Unknown feature: " + name);
                    }
                    
                    // Transform and Scaling
                    raw_value = spec.transform.apply(raw_value);
                    double scaled_value = spec.scaler.scale(raw_value);

                    input(offset + cell_idx) = scaled_value;
                }
                offset += n_cells; // advance offset by number of cells
            }
        }

        return input;
    }

    /*!
    * \brief Run the Hybrid Newton model to produce output predictions.
    *
    * Uses the input tensor (prepared by constructInputTensor) together
    * with the model specified in \p config to compute the corresponding
    * output tensor. Scaling and transformations defined in the configuration
    * are applied automatically.
    *
    * \param input The input tensor of shape [n_cells x n_input_features].
    * \param config The HybridNewtonConfig specifying model path and output
    *        feature definitions.
    * \return A tensor of shape [n_cells x n_output_features], where rows
    *         correspond to cells and columns correspond to output features.
    *
    * \throws std::runtime_error if model inference fails or if the
    *         output tensor does not match the expected feature layout.
    */
    Opm::ML::Tensor<Evaluation>
    constructOutputTensor(const Opm::ML::Tensor<Evaluation>& input, const HybridNewtonConfig& config, size_t n_cells)
    {
        const auto& features = config.output_features;
        const int n_features = features.size();
        
        // TODO check is model path exists 
        Opm::ML::NNModel<Evaluation> model;
        model.loadModel(config.model_path);

        Opm::ML::Tensor<Evaluation> output(1, n_cells * n_features);
        model.apply(input, output);

        return output; 
    }


    /*!
    * \brief Update the nonlinear solver's initial guess using model predictions.
    *
    * For each cell listed in \p cell_indices, this function applies the
    * corresponding row from the \p output tensor to update the initial
    * guess of the simulation state. The mapping between tensor columns
    * and state variables is defined by the configuration's output features.
    *
    * \param output The model output tensor of shape
    *        [n_cells x n_output_features], where rows correspond to cells
    *        and columns correspond to output features.
    * \param config The HybridNewtonConfig specifying model path and output
    *        feature definitions.
    * \param cell_indices The list of cells whose initial guesses are updated.
    * \param n_cells The number of cells (should match cell_indices.size()).
    *
    * \throws std::runtime_error if an expected output feature is missing
    *         or if state update fails for any cell.
    */
    void updateInitialGuess(Opm::ML::Tensor<Evaluation>& output, const HybridNewtonConfig& config, const std::vector<int>& cell_indices, size_t n_cells)
    {
        const auto& features = config.output_features;

        bool has_SWAT = false;
        bool has_SOIL = false;
        bool has_SGAS = false;
        bool has_PRESSURE = false;

        // Determine which features are present
        for (const auto& [name, spec] : features) {
            bool is_delta = name.compare(0, 6, "DELTA_") == 0;
            std::string actual_name = is_delta ? name.substr(6) : name;
            if (actual_name == "SWAT")          has_SWAT = true;
            else if (actual_name == "SOIL")     has_SOIL = true;
            else if (actual_name == "SGAS")     has_SGAS = true;
            else if (actual_name == "PRESSURE") has_PRESSURE = true;
        }

        const auto& unitSyst = simulator_.vanguard().schedule().getUnits();

        for (size_t i = 0; i < n_cells; ++i) {
            const int cell_idx = cell_indices[i];
            const auto& intQuants = simulator_.model().intensiveQuantities(cell_idx, /*timeIdx*/0);
            auto fs = intQuants.fluidState();

            int feature_idx = 0;

            // Temp variables per cell
            double sw_val = -1.0;
            double so_val = -1.0;
            double sg_val = -1.0;
            double po_val = -1.0;

            for (const auto& [name, spec] : features) {
                bool is_delta = name.compare(0, 6, "DELTA_") == 0;
                std::string actual_name = is_delta ? name.substr(6) : name;

                auto scaled_value = getValue(output(feature_idx * n_cells + i));
                
                // Inverse scaling
                double raw_value = spec.scaler.unscale(scaled_value);

                // Inverse transform
                raw_value = spec.transform.applyInverse(raw_value);

                if (is_delta) {
                    if (actual_name == "PRESSURE") {
                        raw_value = raw_value  + unitSyst.from_si(UnitSystem::measure::pressure, getValue(fs.pressure(oilPhaseIdx)));
                    } else if (actual_name == "SWAT") {
                        raw_value += getValue(fs.saturation(waterPhaseIdx));
                    } else if (actual_name == "SOIL") {
                        raw_value += getValue(fs.saturation(oilPhaseIdx));
                    } else if (actual_name == "SGAS") {
                        raw_value += getValue(fs.saturation(gasPhaseIdx));
                    } else if (actual_name == "RS") {
                        raw_value = raw_value + unitSyst.from_si(UnitSystem::measure::gas_oil_ratio, getValue(fs.Rs()));
                    } else if (actual_name == "RV") {
                        raw_value = raw_value + unitSyst.from_si(UnitSystem::measure::oil_gas_ratio, getValue(fs.Rv()));
                         
                    } else {
                        OPM_THROW(std::runtime_error, "Unknown delta feature: " + name);
                    }
                }

                if (actual_name == "PRESSURE") {
                    po_val = unitSyst.to_si(UnitSystem::measure::pressure, raw_value);
                    std::cout << "Cell idx: " << cell_idx << " Feature: "<< name << " Value: " << po_val <<std::endl;
                } else if (actual_name == "SWAT") {
                    sw_val = raw_value;
                } else if (actual_name == "SOIL") {
                    so_val = raw_value;
                } else if (actual_name == "SGAS") {
                    sg_val = raw_value;
                } else if (actual_name == "RS") {
                    if constexpr (compositionSwitchEnabled) {
                        raw_value = unitSyst.to_si(UnitSystem::measure::gas_oil_ratio, raw_value); 
                        fs.setRs(raw_value);
                    }
                } else if (actual_name == "RV") {
                    if constexpr (compositionSwitchEnabled) {
                        raw_value = unitSyst.to_si(UnitSystem::measure::oil_gas_ratio, raw_value);
                        fs.setRv(raw_value);
                    }
                } else {
                    OPM_THROW(std::runtime_error, "Unknown output feature: " + actual_name);
                }

                ++feature_idx;
            }

            int sat_count = static_cast<int>(has_SWAT) + static_cast<int>(has_SOIL) + static_cast<int>(has_SGAS);

            if (sat_count >= 2) {
                double sw = sw_val;
                double so = so_val;
                double sg = sg_val;

                if (!has_SWAT) {
                    sw = 1.0 - so - sg;
                } else if (!has_SOIL) {
                    so = 1.0 - sw - sg;
                } else if (!has_SGAS) {
                    sg = 1.0 - sw - so;
                }

                sw = std::max(0.0, std::min(sw, 1.0));
                so = std::max(0.0, std::min(so, 1.0));
                sg = std::max(0.0, std::min(sg, 1.0));

                double sum = sw + so + sg;
                if (sum <= 1e-12) {
                    OPM_THROW(std::runtime_error, "Saturation sum is zero in cell " + std::to_string(cell_idx));
                }

                fs.setSaturation(waterPhaseIdx, sw / sum);
                fs.setSaturation(oilPhaseIdx,   so / sum);
                fs.setSaturation(gasPhaseIdx,   sg / sum);
            }

            if (has_PRESSURE) {
                std::array<Evaluation, numPhases> pC;
                const auto& materialParams = simulator_.problem().materialLawParams(cell_idx);
                MaterialLaw::capillaryPressures(pC, materialParams, fs);

                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) continue;
                    if (phaseIdx == oilPhaseIdx) {
                        fs.setPressure(phaseIdx, po_val);
                    } else {
                        fs.setPressure(phaseIdx, po_val - pC[phaseIdx]);
                    }
                }
            }

            auto& primaryVars = simulator_.model().solution(/*timeIdx*/0)[cell_idx];
            primaryVars.assignNaive(fs);
        }

        simulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx*/0);
    }

    /*!
    * \brief Load a list of cell indices from a plain text file.
    *
    * Reads one integer index per line, ignoring blank lines and
    * lines beginning with `#` (treated as comments). Used to specify
    * which cells are subject to a given Hybrid Newton model.
    *
    * \param filepath Path to the text file containing cell indices.
    * \return A vector of integer cell indices.
    *
    * \throws std::runtime_error if the file cannot be opened,
    *         or if no valid indices are found.
    */
    std::vector<int> loadCellIndicesFromFile(const std::string& filepath) const
    {
        std::vector<int> indices;
        std::ifstream cellFile(filepath);
        
        if (!cellFile.is_open()) {
            OPM_THROW(std::runtime_error, "Cannot open cell indices file: " + filepath);
        }
        
        std::string line;
        int lineNumber = 0;
        while (std::getline(cellFile, line)) {
            ++lineNumber;
            if (line.empty() || line[0] == '#') {
                continue;
            }
            try {
                int cellIndex = std::stoi(line);
                indices.push_back(cellIndex);
            } catch (...) {
                OPM_THROW(std::runtime_error, 
                    "Invalid cell index at line " + std::to_string(lineNumber) +
                    " in file " + filepath + ": " + line);
            }
        }
        
        if (indices.empty()) {
            OPM_THROW(std::runtime_error, "No valid cell indices found in file: " + filepath);
        }
        return indices;
}

    protected:
        Simulator& simulator_;
        std::vector<HybridNewtonConfig> configs_;
        bool configsLoaded_;
    };

} // namespace Opm

#endif // HYBRID_NEWTON_CLASS_HPP
