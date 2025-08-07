#ifndef HYBRID_NEWTON_HPP
#define HYBRID_NEWTON_HPP


#include <opm/simulators/flow/FlowBaseProblemProperties.hpp>
#include <opm/simulators/flow/HybridNewtonConfig.hpp>

#include <opm/ml/ml_model.hpp>

#include <unordered_map>   
#include <string>         
#include <iostream>
#include <fstream>
#include <boost/property_tree/json_parser.hpp>

namespace Opm {

template <typename TypeTag>
class HybridNewton
{
protected:
    using Simulator  = GetPropType<TypeTag, Properties::Simulator>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;

    enum { numPhases = FluidSystem::numPhases };

    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

public:
    explicit HybridNewton(Simulator& simulator)
        : simulator_(simulator)
        , configLoaded_(false)
    {}
    void tryApplyHybridNewton()
    {
        if (!Parameters::Get<Parameters::UseHyNe>())
            return;
        
        if (!configLoaded_) {
            config_ = parseHybridNewtonConfig(Parameters::Get<Parameters::HyNeConfigFile>());
            configLoaded_ = true;
        }

        double current_time = simulator_.time();
        if (shouldApplyHybridNewton(current_time, config_)) {
            std::cout << "Using Hybrid Newton....." << std::endl;
            runHybridNewton(config_);
        }
    }

private: 
    void runHybridNewton(const HybridNewtonConfig& config)
    {
        std::cout << "Running HybridNewton with model path: " << config.model_path << std::endl;

        auto input = constructInputTensor();
        auto output = constructOutputTensor(input);

        updateInitialGuess(output);

        std::cout << "HybridNewton algorithm executed." << std::endl;
    }


protected:
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

        HybridNewtonConfig parseHybridNewtonConfig(const std::string& path)
    {
        HybridNewtonConfig config;
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

        config.model_path = pt.get<std::string>("model_path", "");

        for (const auto& item : pt.get_child("cell_indices", {})) {
            config.cell_indices.push_back(item.second.get_value<int>());
        }

        for (const auto& item : pt.get_child("apply_times", {})) {
            config.apply_times.push_back(item.second.get_value<double>());
        }

        auto parseFeatures = [](const boost::property_tree::ptree& subtree) {
            std::unordered_map<std::string, FeatureSpec> result;
            for (const auto& [key, val] : subtree) {
                FeatureSpec spec;
                spec.unit = val.get<std::string>("unit", "");
                spec.transform = val.get<std::string>("transform", "none");
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
                result[key] = spec;
            }
            return result;
        };

        config.input_features = parseFeatures(pt.get_child("features.inputs"));
        config.output_features = parseFeatures(pt.get_child("features.outputs"));

        return config;
    }

    Opm::ML::Tensor<Evaluation>
    constructInputTensor()
    {   
        const auto& features = config_.input_features;
        const auto& cell_indices = config_.cell_indices;
        const size_t n_cells = cell_indices.size();
        const int n_features = features.size();

        auto& eclState = simulator_.vanguard().eclState();
        auto& fp = eclState.fieldProps();
        auto permX = fp.get_double("PERMX"); 
        
        Opm::ML::Tensor<Evaluation> input(n_cells * n_features);

        int feature_idx = 0;
        for (const auto& [name, spec] : features) {
            // Special case: scalar values
            if (name == "TIMESTEP") {
                double raw_value = simulator_.timeStepSize(); 

                // Apply transform
                if (spec.transform == "log10") {
                    raw_value = std::log10(raw_value);
                } else if (spec.transform == "log") {
                    raw_value = std::log(raw_value);
                } else if (spec.transform == "log1p") {
                    raw_value = std::log1p(raw_value);
                }

                // Apply scaling
                double scaled_value = spec.scaler.scale(raw_value);

                for (size_t i = 0; i < n_cells; ++i) {
                    input(feature_idx * n_cells + i) = scaled_value;
                }

                ++feature_idx;
                continue;
            }
            for (size_t i = 0; i < cell_indices.size(); ++i) {
                double raw_value = 0.0;
                const auto& intQuants = simulator_.model().intensiveQuantities(cell_indices[i], /*timeIdx*/ 0);
                const auto& fs = intQuants.fluidState();

                // Get raw value based on feature name
                if (name == "PRESSURE") {
                    raw_value = fs.pressure(oilPhaseIdx).value();
                } else if (name == "SWAT") {
                    raw_value = fs.saturation(waterPhaseIdx).value();
                } else if (name == "SGAS") {
                    raw_value = fs.saturation(gasPhaseIdx).value();
                } else if (name == "SOIL") {
                    raw_value = fs.saturation(oilPhaseIdx).value();
                } else if (name == "RS") {
                    raw_value = fs.Rs().value();
                } else if (name == "RV") {
                    raw_value = fs.Rv().value();
                } else if (name == "PERMX") {
                    raw_value = permX[cell_indices[i]];
                }
                else {
                    OPM_THROW(std::runtime_error, "Unknown feature: " + name);
                }

                // Apply transform
                if (spec.transform == "log10") {
                    raw_value = std::log10(raw_value);
                } else if (spec.transform == "log") {
                    raw_value = std::log(raw_value);
                } else if (spec.transform == "log1p") {
                    raw_value = std::log1p(raw_value);
                }

                // Apply scaling
                const double mean = spec.scaler.mean;
                const double stddev = spec.scaler.std;
                double scaled_value = (stddev == 0.0) ? mean : (raw_value - mean) / stddev;

                // Store in tensor
                input(feature_idx * n_cells + i) = scaled_value;
            }
            ++feature_idx;
        }

        return input;
    }

    Opm::ML::Tensor<Evaluation>
    constructOutputTensor(const Opm::ML::Tensor<Evaluation>& input)
    {
        const auto& features = config_.output_features;
        const auto& cell_indices = config_.cell_indices;
        const int n_cells = cell_indices.size();
        const int n_features = features.size();
        
        Opm::ML::NNModel<Evaluation> model;
        model.loadModel(config_.model_path);

        Opm::ML::Tensor<Evaluation> output(1, n_cells * n_features);
        model.apply(input, output);

        return output; 
    }

    void updateInitialGuess(Opm::ML::Tensor<Evaluation>& output)
    {
        const auto& features = config_.output_features;
        const auto& cell_indices = config_.cell_indices;
        const size_t n_cells = cell_indices.size();

        // Temporary buffers for saturation and pressure
        std::vector<double> sw_buffer(n_cells, -1.0);
        std::vector<double> so_buffer(n_cells, -1.0);
        std::vector<double> sg_buffer(n_cells, -1.0);
        std::vector<double> po_buffer(n_cells, -1.0);

        bool has_SWAT = false;
        bool has_SOIL = false;
        bool has_SGAS = false;
        bool has_PRESSURE = false;

        int feature_idx = 0;
        for (const auto& [name, spec] : features) {
            bool is_delta = name.compare(0, 6, "DELTA_") == 0;
            std::string actual_name = is_delta ? name.substr(6) : name;

            if (actual_name == "SWAT")          has_SWAT = true;
            else if (actual_name == "SOIL")     has_SOIL = true;
            else if (actual_name == "SGAS")     has_SGAS = true;
            else if (actual_name == "PRESSURE") has_PRESSURE = true; 
            
            for (size_t i = 0; i < n_cells; ++i) {
                const int cell_idx = cell_indices[i];
                auto scaled_value = output(feature_idx * n_cells + i).value();

                // Inverse scaling
                double raw_value = spec.scaler.unscale(scaled_value);
                
                // Inverse transform
                if (spec.transform == "log10")      raw_value = std::pow(10.0, raw_value);
                else if (spec.transform == "log")   raw_value = std::exp(raw_value);
                else if (spec.transform == "log1p") raw_value = std::expm1(raw_value);
                
                // Convert delta feature to absolute 
                if (is_delta) {
                    auto& intQuants = simulator_.model().intensiveQuantities(cell_idx, /*timeIdx*/0);
                    auto fs = intQuants.fluidState();
                    if (actual_name == "PRESSURE") {
                        raw_value += fs.pressure(oilPhaseIdx).value(); 
                    } else if (actual_name == "SWAT") {
                        raw_value += fs.saturation(waterPhaseIdx).value();
                    } else if (actual_name == "SOIL") {
                        raw_value += fs.saturation(oilPhaseIdx).value();
                    } else if (actual_name == "SGAS") {
                        raw_value += fs.saturation(gasPhaseIdx).value();
                    } else if (actual_name == "RS") {
                        raw_value += fs.Rs().value();
                    } else if (actual_name == "RV") {
                        raw_value += fs.Rv().value();
                    } else {
                        OPM_THROW(std::runtime_error, "Unknown delta feature: " + name);
                    }
                }

                // Store in buffers or assign
                if (actual_name == "PRESSURE") {
                    po_buffer[i] = raw_value;
                } else if (actual_name == "SWAT") {
                    sw_buffer[i] = raw_value;
                } else if (actual_name == "SOIL") {
                    so_buffer[i] = raw_value;
                } else if (actual_name == "SGAS") {
                    sg_buffer[i] = raw_value;
                } else if (actual_name == "RS") {
                    // Access mutable fluid state
                    auto& intQuants = simulator_.model().intensiveQuantities(cell_idx, /*timeIdx*/0);
                    auto fs = intQuants.fluidState();
                    fs.setRs(raw_value);
                } else if (actual_name == "RV") {
                    // Access mutable fluid state
                    auto& intQuants = simulator_.model().intensiveQuantities(cell_idx, /*timeIdx*/0);
                    auto fs = intQuants.fluidState();
                    fs.setRv(raw_value);
                } else {
                    OPM_THROW(std::runtime_error, "Unknown output feature: " + actual_name);
                }
            }
            ++feature_idx;
        }

        // After all features processed, set saturations (ensuring sum = 1)
        int sat_count = static_cast<int>(has_SWAT) + static_cast<int>(has_SOIL) + static_cast<int>(has_SGAS);
        if (sat_count >= 2) {
            // compute third and set saturations
            for (size_t i = 0; i < n_cells; ++i) {
                const int cell_idx = cell_indices[i];
                const auto& intQuants = simulator_.model().intensiveQuantities(cell_idx, /*timeIdx*/0);
                auto fs = intQuants.fluidState();

                double sw = sw_buffer[i];
                double so = so_buffer[i];
                double sg = sg_buffer[i];

                int known = int(has_SWAT) + int(has_SOIL) + int(has_SGAS);

                if (known < 2) {
                    std::cerr << "Warning: Insufficient saturation information for cell " << cell_idx << ". Skipping saturation update." << std::endl;
                    continue;
                }

                // Infer missing saturation
                if (!has_SWAT) {
                    sw = 1.0 - so - sg;
                } else if (!has_SOIL) {
                    so = 1.0 - sw - sg;
                } else if (!has_SGAS) {
                    sg = 1.0 - sw - so;
                }

                // Clip to avoid small negative values or overshoots
                sw = std::max(0.0, std::min(sw, 1.0));
                so = std::max(0.0, std::min(so, 1.0));
                sg = std::max(0.0, std::min(sg, 1.0));

                // Renormalize to sum to 1
                double sum = sw + so + sg;
                if (sum <= 1e-12) {
                    OPM_THROW(std::runtime_error, "Saturation sum is zero in cell " + std::to_string(cell_idx));
                }

                fs.setSaturation(waterPhaseIdx, sw / sum);
                fs.setSaturation(oilPhaseIdx,   so / sum);
                fs.setSaturation(gasPhaseIdx,   sg / sum);
            }
        }

        // Now update pressures with fresh saturations (capillary pressure recalculated)
        if (has_PRESSURE) {
            for (size_t i = 0; i < n_cells; ++i) {
                const int cell_idx = cell_indices[i];
                const auto& intQuants = simulator_.model().intensiveQuantities(cell_idx, /*timeIdx*/0);
                auto fs = intQuants.fluidState();
                const double raw_pressure = po_buffer[i];
                std::array<Evaluation, numPhases> pC;
                const auto& materialParams = simulator_.problem().materialLawParams(cell_idx);
                MaterialLaw::capillaryPressures(pC, materialParams, fs);

                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) continue;
                    if (phaseIdx == oilPhaseIdx) {
                        fs.setPressure(phaseIdx, raw_pressure);
                    } else {
                        fs.setPressure(phaseIdx, raw_pressure - pC[phaseIdx]);
                    }
                }
            }
        }

        // Assign to primary variables
        for (size_t i = 0; i < n_cells; ++i) {
            const int cell_idx = cell_indices[i];
            const auto& intQuants = simulator_.model().intensiveQuantities(cell_idx, /*timeIdx*/0);
            auto fs = intQuants.fluidState();
            auto& primaryVars = simulator_.model().solution(/*timeIdx*/0)[cell_idx];
            primaryVars.assignNaive(fs);
        }

        // Update intensive quantities
        simulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx*/0);
    }

protected:
    Simulator& simulator_;
    HybridNewtonConfig config_;
    bool configLoaded_;
};

} // namespace Opm

#endif // HYBRID_NEWTON_CLASS_HPP
