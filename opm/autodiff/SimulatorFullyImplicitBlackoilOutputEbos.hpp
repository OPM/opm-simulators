/*
  Copyright (c) 2014 SINTEF ICT, Applied Mathematics.
  Copyright (c) 2015 IRIS AS

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
#ifndef OPM_SIMULATORFULLYIMPLICITBLACKOILOUTPUTEBOS_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITBLACKOILOUTPUTEBOS_HEADER_INCLUDED
#include <opm/core/grid.h>
#include <opm/core/simulator/SimulatorTimerInterface.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/DataMap.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/output/eclipse/EclipseReader.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/wells/DynamicListEconLimited.hpp>

#include <opm/output/eclipse/EclipseWriter.hpp>

#include <opm/autodiff/Compat.hpp>
#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/ParallelDebugOutput.hpp>

#include <opm/autodiff/WellStateFullyImplicitBlackoilDense.hpp>
#include <opm/autodiff/ThreadHandle.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>


#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <thread>
#include <memory>

#include <boost/filesystem.hpp>

#ifdef HAVE_OPM_GRID
#include <dune/grid/CpGrid.hpp>
#endif
namespace Opm
{
    class BlackoilState;

    /** \brief Wrapper class for VTK, Matlab, and ECL output. */
    class BlackoilOutputWriterEbos
    {

    public:
        // constructor creating different sub writers
        template <class Grid>
        BlackoilOutputWriterEbos(const Grid& grid,
                                 const parameter::ParameterGroup& param,
                                 const EclipseState& eclipseState,
                                 std::unique_ptr<EclipseWriter>&& eclWriter,
                                 const Opm::PhaseUsage &phaseUsage,
                                 const double* permeability );

        /*!
         * \brief Write a blackoil reservoir state to disk for later inspection with
         *        visualization tools like ResInsight. This function will extract the
         *        requested output cell properties specified by the RPTRST keyword
         *        and write these to file.
         */
        template<class Model>
        void writeTimeStep(const SimulatorTimerInterface& timer,
                           const SimulationDataContainer& reservoirState,
                           const Opm::WellState& wellState,
                           const Model& physicalModel,
                           bool substep = false);


        /*!
         * \brief Write a blackoil reservoir state to disk for later inspection with
         *        visualization tools like ResInsight. This function will write all
         *        CellData in simProps to the file as well.
         */
        void writeTimeStepWithCellProperties(
                           const SimulatorTimerInterface& timer,
                           const SimulationDataContainer& reservoirState,
                           const Opm::WellState& wellState,
                           const data::Solution& sol,
                           bool substep = false);

        /*!
         * \brief Write a blackoil reservoir state to disk for later inspection with
         *        visualization tools like ResInsight. This function will not write
         *        any cell properties (e.g., those requested by RPTRST keyword)
         */
        void writeTimeStepWithoutCellProperties(
                           const SimulatorTimerInterface& timer,
                           const SimulationDataContainer& reservoirState,
                           const Opm::WellState& wellState,
                           bool substep = false);

        /*!
         * \brief Write a blackoil reservoir state to disk for later inspection withS
         *        visualization tools like ResInsight. This is the function which does
         *        the actual write to file.
         */
        void writeTimeStepSerial(const SimulatorTimerInterface& timer,
                                 const SimulationDataContainer& reservoirState,
                                 const Opm::WellState& wellState,
                                 const data::Solution& simProps,
                                 bool substep);

        /** \brief return output directory */
        const std::string& outputDirectory() const { return outputDir_; }

        /** \brief return true if output is enabled */
        bool output () const { return output_; }

        void restore(SimulatorTimerInterface& timer,
                     BlackoilState& state,
                     WellStateFullyImplicitBlackoilDense& wellState,
                     const std::string& filename,
                     const int desiredReportStep);


        template <class Grid>
        void initFromRestartFile(const PhaseUsage& phaseusage,
                                 const double* permeability,
                                 const Grid& grid,
                                 SimulationDataContainer& simulatorstate,
                                 WellStateFullyImplicitBlackoilDense& wellstate);

        bool isRestart() const;

    protected:
        const bool output_;
        std::unique_ptr< ParallelDebugOutputInterface > parallelOutput_;

        // Parameters for output.
        const std::string outputDir_;
        const int output_interval_;

        int lastBackupReportStep_;

        std::ofstream backupfile_;
        Opm::PhaseUsage phaseUsage_;
        std::unique_ptr<EclipseWriter> eclWriter_;
        const EclipseState& eclipseState_;

        std::unique_ptr< ThreadHandle > asyncOutput_;
    };


    //////////////////////////////////////////////////////////////
    //
    //  Implementation
    //
    //////////////////////////////////////////////////////////////
    template <class Grid>
    inline
    BlackoilOutputWriterEbos::
    BlackoilOutputWriterEbos(const Grid& grid,
                             const parameter::ParameterGroup& param,
                             const Opm::EclipseState& eclipseState,
                             std::unique_ptr<EclipseWriter>&& eclWriter,
                             const Opm::PhaseUsage &phaseUsage,
                             const double* permeability )
      : output_( param.getDefault("output", true) ),
        parallelOutput_( output_ ? new ParallelDebugOutput< Grid >( grid, eclipseState, phaseUsage.num_phases, permeability ) : 0 ),
        outputDir_( output_ ? param.getDefault("output_dir", std::string("output")) : "." ),
        output_interval_( output_ ? param.getDefault("output_interval", 1): 0 ),
        lastBackupReportStep_( -1 ),
        phaseUsage_( phaseUsage ),
        eclipseState_(eclipseState),
        asyncOutput_()
    {
        // For output.
        if (output_ && parallelOutput_->isIORank() ) {
            eclWriter_ = std::move(eclWriter);

            // Ensure that output dir exists
            boost::filesystem::path fpath(outputDir_);
            try {
                create_directories(fpath);
            }
            catch (...) {
                OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
            }

            // create output thread if enabled and rank is I/O rank
            // async output is enabled by default if pthread are enabled
#if HAVE_PTHREAD
            const bool asyncOutputDefault = false;
#else
            const bool asyncOutputDefault = false;
#endif
            if( param.getDefault("async_output", asyncOutputDefault ) )
            {
#if HAVE_PTHREAD
                asyncOutput_.reset( new ThreadHandle() );
#else
                OPM_THROW(std::runtime_error,"Pthreads were not found, cannot enable async_output");
#endif
            }

            std::string backupfilename = param.getDefault("backupfile", std::string("") );
            if( ! backupfilename.empty() )
            {
                backupfile_.open( backupfilename.c_str() );
            }
        }
    }


    template <class Grid>
    inline void
    BlackoilOutputWriterEbos::
    initFromRestartFile( const PhaseUsage& phaseusage,
                         const double* permeability,
                         const Grid& grid,
                         SimulationDataContainer& simulatorstate,
                         WellStateFullyImplicitBlackoilDense& wellstate)
    {
        // gives a dummy dynamic_list_econ_limited
        DynamicListEconLimited dummy_list_econ_limited;
        WellsManager wellsmanager(eclipseState_,
                                  eclipseState_.getInitConfig().getRestartStep(),
                                  Opm::UgGridHelpers::numCells(grid),
                                  Opm::UgGridHelpers::globalCell(grid),
                                  Opm::UgGridHelpers::cartDims(grid),
                                  Opm::UgGridHelpers::dimensions(grid),
                                  Opm::UgGridHelpers::cell2Faces(grid),
                                  Opm::UgGridHelpers::beginFaceCentroids(grid),
                                  permeability,
                                  dummy_list_econ_limited
                                  // We need to pass the optionaly arguments
                                  // as we get the following error otherwise
                                  // with c++ (Debian 4.9.2-10) 4.9.2 and -std=c++11
                                  // converting to ‘const std::unordered_set<std::basic_string<char> >’ from initializer list would use explicit constructo
                                  , false,
                                  std::vector<double>(),
                                  std::unordered_set<std::string>());

        const Wells* wells = wellsmanager.c_wells();
        wellstate.resize(wells, simulatorstate); //Resize for restart step
        auto restarted = Opm::init_from_restart_file(
                                eclipseState_,
                                Opm::UgGridHelpers::numCells(grid) );

        solutionToSim( restarted.first, phaseusage, simulatorstate );
        wellsToState( restarted.second, phaseusage, wellstate );
    }





    namespace detail {
        template<class Model>
    void getOutputDataEbos(data::Solution& output,
                           const Opm::PhaseUsage& phaseUsage,
                           const Model& model,
                           const RestartConfig& restartConfig,
                           const int reportStepNum,
                           const bool log)
    {
            typedef typename Model::FluidSystem FluidSystem;

            //Get the value of each of the keys
            std::map<std::string, int> outKeywords = restartConfig.getRestartKeywords(reportStepNum);
            for (auto& keyValue : outKeywords) {
                keyValue.second = restartConfig.getKeyword(keyValue.first, reportStepNum);
            }

            const auto& ebosModel = model.ebosSimulator().model();

            //Get shorthands for water, oil, gas
            const int aqua_active = phaseUsage.phase_used[Opm::PhaseUsage::Aqua];
            const int liquid_active = phaseUsage.phase_used[Opm::PhaseUsage::Liquid];
            const int vapour_active = phaseUsage.phase_used[Opm::PhaseUsage::Vapour];

            // extract everything which can possibly be written to disk
            int numCells = ebosModel.numGridDof();

            std::vector<double> pressureOil(numCells);
            std::vector<double> temperature(numCells);

            std::vector<double> satWater(numCells);
            std::vector<double> satGas(numCells);

            std::vector<double> bWater(numCells);
            std::vector<double> bOil(numCells);
            std::vector<double> bGas(numCells);

            std::vector<double> rhoWater(numCells);
            std::vector<double> rhoOil(numCells);
            std::vector<double> rhoGas(numCells);

            std::vector<double> muWater(numCells);
            std::vector<double> muOil(numCells);
            std::vector<double> muGas(numCells);

            std::vector<double> krWater(numCells);
            std::vector<double> krOil(numCells);
            std::vector<double> krGas(numCells);

            std::vector<double> Rs(numCells);
            std::vector<double> Rv(numCells);
            std::vector<double> RsSat(numCells);
            std::vector<double> RvSat(numCells);

            for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
                const auto& intQuants = *ebosModel.cachedIntensiveQuantities(cellIdx, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                pressureOil[cellIdx] = fs.pressure(FluidSystem::oilPhaseIdx).value();

                temperature[cellIdx] = fs.temperature(FluidSystem::oilPhaseIdx).value();

                if (aqua_active) {
                    satWater[cellIdx] = fs.saturation(FluidSystem::waterPhaseIdx).value();
                    bWater[cellIdx] = fs.invB(FluidSystem::waterPhaseIdx).value();
                    rhoWater[cellIdx] = fs.density(FluidSystem::waterPhaseIdx).value();
                    muWater[cellIdx] = fs.viscosity(FluidSystem::waterPhaseIdx).value();
                    krWater[cellIdx] = intQuants.relativePermeability(FluidSystem::waterPhaseIdx).value();
                }
                if (vapour_active) {
                    satGas[cellIdx] = fs.saturation(FluidSystem::gasPhaseIdx).value();
                    bGas[cellIdx] = fs.invB(FluidSystem::gasPhaseIdx).value();
                    rhoGas[cellIdx] = fs.density(FluidSystem::gasPhaseIdx).value();
                    muGas[cellIdx] = fs.viscosity(FluidSystem::gasPhaseIdx).value();
                    krGas[cellIdx] = intQuants.relativePermeability(FluidSystem::gasPhaseIdx).value();
                    Rs[cellIdx] = fs.Rs().value();
                    Rv[cellIdx] = fs.Rv().value();
                    RsSat[cellIdx] = FluidSystem::saturatedDissolutionFactor(fs,
                                                                             FluidSystem::oilPhaseIdx,
                                                                             intQuants.pvtRegionIndex(),
                                                                             /*maxOilSaturation=*/1.0).value();
                    RvSat[cellIdx] = FluidSystem::saturatedDissolutionFactor(fs,
                                                                             FluidSystem::gasPhaseIdx,
                                                                             intQuants.pvtRegionIndex(),
                                                                             /*maxOilSaturation=*/1.0).value();
                }

                // oil is always active
                bOil[cellIdx] = fs.invB(FluidSystem::oilPhaseIdx).value();              
                rhoOil[cellIdx] = fs.density(FluidSystem::oilPhaseIdx).value();               
                muOil[cellIdx] = fs.viscosity(FluidSystem::oilPhaseIdx).value();
                krOil[cellIdx] = intQuants.relativePermeability(FluidSystem::oilPhaseIdx).value();
            }

            /**
             * Oil Pressures
             */
            outKeywords["PRESSURE"] = 0;
            output.insert("PRESSURE",
                       UnitSystem::measure::pressure,
                       std::move(pressureOil),
                       data::TargetType::RESTART_SOLUTION);

            /**
             * Temperatures
             */
            outKeywords["TEMP"] = 0;
            output.insert("TEMP",
                       UnitSystem::measure::temperature,
                       std::move(temperature),
                       data::TargetType::RESTART_SOLUTION);

            /**
             * Water and gas saturation.
             */
            outKeywords["SWAT"] = 0;
            outKeywords["SGAS"] = 0;
            output.insert("SWAT",
                       UnitSystem::measure::identity,
                       std::move(satWater),
                       data::TargetType::RESTART_SOLUTION);
            output.insert("SGAS",
                       UnitSystem::measure::identity,
                       std::move(satGas),
                       data::TargetType::RESTART_SOLUTION);

            /**
             * the dissolution factors
             */
            outKeywords["RS"] = 0;
            outKeywords["RV"] = 0;
            output.insert("RS",
                       UnitSystem::measure::gas_oil_ratio,
                       std::move(Rs),
                       data::TargetType::RESTART_SOLUTION);
            output.insert("RV",
                       UnitSystem::measure::oil_gas_ratio,
                       std::move(Rv),
                       data::TargetType::RESTART_SOLUTION);

            /**
             * Formation volume factors for water, oil, gas
             */
            if (aqua_active && outKeywords["BW"] > 0) {
                outKeywords["BW"] = 0;
                output.insert("BW",
                              Opm::UnitSystem::measure::water_inverse_formation_volume_factor,
                              std::move(bWater),
                              data::TargetType::RESTART_AUXILLARY);

            }
            if (liquid_active && outKeywords["BO"]  > 0) {
                outKeywords["BO"] = 0;
                output.insert("BO",
                           Opm::UnitSystem::measure::oil_inverse_formation_volume_factor,
                           std::move(bOil),
                           data::TargetType::RESTART_AUXILLARY);
            }
            if (vapour_active && outKeywords["BG"] > 0) {
                outKeywords["BG"] = 0;
                output.insert("BG",
                           Opm::UnitSystem::measure::gas_inverse_formation_volume_factor,
                           std::move(bGas),
                           data::TargetType::RESTART_AUXILLARY);
            }

            /**
             * Densities for water, oil gas
             */
            if (outKeywords["DEN"] > 0) {
                outKeywords["DEN"] = 0;
                if (aqua_active) {
                    output.insert("WAT_DEN",
                               Opm::UnitSystem::measure::density,
                               std::move(rhoWater),
                               data::TargetType::RESTART_AUXILLARY);
                }
                if (liquid_active) {
                    output.insert("OIL_DEN",
                               Opm::UnitSystem::measure::density,
                               std::move(rhoOil),
                               data::TargetType::RESTART_AUXILLARY);
                }
                if (vapour_active) {
                output.insert("GAS_DEN",
                           Opm::UnitSystem::measure::density,
                           std::move(rhoGas),
                           data::TargetType::RESTART_AUXILLARY);
                }
            }

            /**
             * Viscosities for water, oil gas
             */
            if (outKeywords["VISC"] > 0) {
                outKeywords["VISC"] = 0;
                if (aqua_active) {
                output.insert("WAT_VISC",
                           Opm::UnitSystem::measure::viscosity,
                           std::move(muWater),
                           data::TargetType::RESTART_AUXILLARY);
                }
                if (liquid_active) {
                output.insert("OIL_VISC",
                           Opm::UnitSystem::measure::viscosity,
                           std::move(muOil),
                           data::TargetType::RESTART_AUXILLARY);
                }
                if (vapour_active) {
                output.insert("GAS_VISC",
                           Opm::UnitSystem::measure::viscosity,
                           std::move(muGas),
                           data::TargetType::RESTART_AUXILLARY);
                }
            }

            /**
             * Relative permeabilities for water, oil, gas
             */
            if (aqua_active && outKeywords["KRW"] > 0) {
                outKeywords["KRW"] = 0;
                output.insert("WATKR",
                           Opm::UnitSystem::measure::permeability,
                           std::move(krWater),
                           data::TargetType::RESTART_AUXILLARY);
            }
            if (liquid_active && outKeywords["KRO"] > 0) {
                outKeywords["KRO"] = 0;
                output.insert("OILKR",
                           Opm::UnitSystem::measure::permeability,
                           std::move(krOil),
                           data::TargetType::RESTART_AUXILLARY);
            }
            if (vapour_active && outKeywords["KRG"] > 0) {
                outKeywords["KRG"] = 0;
                output.insert("GASKR",
                           Opm::UnitSystem::measure::permeability,
                           std::move(krGas),
                           data::TargetType::RESTART_AUXILLARY);
            }

            /**
             * Vaporized and dissolved gas/oil ratio
             */
            if (vapour_active && liquid_active && outKeywords["RSSAT"] > 0) {
                outKeywords["RSSAT"] = 0;
                output.insert("RSSAT",
                           Opm::UnitSystem::measure::gas_oil_ratio,
                           std::move(RsSat),
                           data::TargetType::RESTART_AUXILLARY);
            }
            if (vapour_active && liquid_active && outKeywords["RVSAT"] > 0) {
                outKeywords["RVSAT"] = 0;
                output.insert("RVSAT",
                           Opm::UnitSystem::measure::oil_gas_ratio,
                           std::move(RvSat),
                           data::TargetType::RESTART_AUXILLARY);
            }


            /**
             * Bubble point and dew point pressures
             */
            if (log && vapour_active &&
                    liquid_active && outKeywords["PBPD"] > 0) {
                Opm::OpmLog::warning("Bubble/dew point pressure output unsupported",
                                     "Writing bubble points and dew points (PBPD) to file is unsupported, "
                                     "as the simulator does not use these internally.");
            }

            //Warn for any unhandled keyword
            if (log) {
                for (auto& keyValue : outKeywords) {
                    if (keyValue.second > 0) {
                        std::string logstring = "Keyword '";
                        logstring.append(keyValue.first);
                        logstring.append("' is unhandled for output to file.");
                        Opm::OpmLog::warning("Unhandled output keyword", logstring);
                    }
                }
            }

        }

        /**
         * Checks if the summaryConfig has a keyword with the standardized field, region, or block prefixes.
         */
        inline bool hasFRBKeyword(const SummaryConfig& summaryConfig, const std::string keyword) {
            std::string field_kw = "F" + keyword;
            std::string region_kw = "R" + keyword;
            std::string block_kw = "B" + keyword;
            return summaryConfig.hasKeyword(field_kw)
                    || summaryConfig.hasKeyword(region_kw)
                    || summaryConfig.hasKeyword(block_kw);
        }


        /**
         * Returns the data as asked for in the summaryConfig
         */
        template<class Model>
        void getSummaryData(data::Solution& output,
                            const Opm::PhaseUsage& phaseUsage,
                            const Model& physicalModel,
                            const SummaryConfig& summaryConfig) {

            const typename Model::FIPData& fip = physicalModel.getFIPData();

            //Get shorthands for water, oil, gas
            const int aqua_active = phaseUsage.phase_used[Opm::PhaseUsage::Aqua];
            const int liquid_active = phaseUsage.phase_used[Opm::PhaseUsage::Liquid];
            const int vapour_active = phaseUsage.phase_used[Opm::PhaseUsage::Vapour];

            /**
             * Now process all of the summary config files
             */
            // Water in place
            if (aqua_active && hasFRBKeyword(summaryConfig, "WIP")) {
                output.insert("WIP",
                              Opm::UnitSystem::measure::volume,
                              fip.fip[Model::FIPData::FIP_AQUA],
                              data::TargetType::SUMMARY );
            }
            if (liquid_active) {
                const std::vector<double>& oipl = fip.fip[Model::FIPData::FIP_LIQUID];
                const int size = oipl.size();

                const std::vector<double>& oipg = vapour_active ? fip.fip[Model::FIPData::FIP_VAPORIZED_OIL] : std::vector<double>(size,0.0);
                std::vector<double> oip = oipl;
                if (vapour_active) {
                    oip.insert(oip.end(), oipg.begin(), oipg.end());
                }

                //Oil in place (liquid phase only)
                if (hasFRBKeyword(summaryConfig, "OIPL")) {
                    output.insert("OIPL",
                                  Opm::UnitSystem::measure::volume,
                                  oipl,
                                  data::TargetType::SUMMARY );
                }

                //Oil in place (gas phase only)
                if (hasFRBKeyword(summaryConfig, "OIPG")) {
                    output.insert("OIPG",
                                  Opm::UnitSystem::measure::volume,
                                  oipg,
                                  data::TargetType::SUMMARY );
                }
                // Oil in place (in liquid and gas phases)
                if (hasFRBKeyword(summaryConfig, "OIP")) {
                    output.insert("OIP",
                                  Opm::UnitSystem::measure::volume,
                                  oip,
                                  data::TargetType::SUMMARY );
                }



            }
            if (vapour_active) {
                const std::vector<double>& gipg = fip.fip[Model::FIPData::FIP_VAPOUR];
                const int size = gipg.size();

                const std::vector<double>& gipl= liquid_active ? fip.fip[Model::FIPData::FIP_DISSOLVED_GAS] : std::vector<double>(size,0.0);
                std::vector<double> gip = gipg;
                if (liquid_active) {
                    gip.insert(gip.end(), gipl.begin(), gipl.end());
                }

                // Gas in place (gas phase only)
                if (hasFRBKeyword(summaryConfig, "GIPG")) {
                    output.insert("GIPG",
                                  Opm::UnitSystem::measure::volume,
                                  gipg,
                                  data::TargetType::SUMMARY );
                }
                // Gas in place (liquid phase only)
                if (hasFRBKeyword(summaryConfig, "GIPL")) {
                    output.insert("GIPL",
                                  Opm::UnitSystem::measure::volume,
                                  gipl,
                                  data::TargetType::SUMMARY );
                }
                // Gas in place (in both liquid and gas phases)
                if (hasFRBKeyword(summaryConfig, "GIP")) {
                    output.insert("GIP",
                                  Opm::UnitSystem::measure::volume,
                                  gip,
                                  data::TargetType::SUMMARY );
                }
            }
            // Cell pore volume in reservoir conditions
            if (hasFRBKeyword(summaryConfig, "RPV")) {
                output.insert("RPV",
                              Opm::UnitSystem::measure::volume,
                              fip.fip[Model::FIPData::FIP_PV],
                              data::TargetType::SUMMARY );
            }
            // Pressure averaged value (hydrocarbon pore volume weighted)
            if (summaryConfig.hasKeyword("FPRH") || summaryConfig.hasKeyword("RPRH")) {
                output.insert("PRH",
                              Opm::UnitSystem::measure::pressure,
                              fip.fip[Model::FIPData::FIP_WEIGHTED_PRESSURE],
                              data::TargetType::SUMMARY );
            }
        }

    }




    template<class Model>
    inline void
    BlackoilOutputWriterEbos::
    writeTimeStep(const SimulatorTimerInterface& timer,
                  const SimulationDataContainer& localState,
                  const WellState& localWellState,
                  const Model& physicalModel,
                  bool substep)
    {
        data::Solution cellData{};
        const RestartConfig& restartConfig = eclipseState_.getRestartConfig();
        const SummaryConfig& summaryConfig = eclipseState_.getSummaryConfig();
        const int reportStepNum = timer.reportStepNum();
        bool logMessages = output_ && parallelOutput_->isIORank();

        if( output_ && !parallelOutput_->isParallel() )
        {

            detail::getOutputDataEbos( cellData, phaseUsage_, physicalModel,
                                    restartConfig, reportStepNum, logMessages );
            detail::getSummaryData( cellData, phaseUsage_, physicalModel, summaryConfig );
        }
        else
        {
            if ( logMessages )
            {
                std::map<std::string, int> rstKeywords = restartConfig.getRestartKeywords(reportStepNum);
                std::vector<const char*> keywords =
                    { "WIP", "OIPL", "OIPG", "OIP", "GIPG", "GIPL", "GIP",
                      "RPV", "FRPH", "RPRH"};

                std::ostringstream str;
                str << "Output of restart/summary config not supported in parallel. Requested keywords were ";
                std::size_t no_kw = 0;

                auto func = [&] (const char* kw)
                    {
                        if ( detail::hasFRBKeyword(summaryConfig, kw) )
                        {
                            str << kw << " ";
                            ++ no_kw;
                        }
                    };

                std::for_each(keywords.begin(), keywords.end(), func);

                for (auto& keyValue : rstKeywords)
                {
                        str << keyValue.first << " ";
                        ++ no_kw;
                }

                if ( no_kw )
                {
                    Opm::OpmLog::warning("Unhandled ouput request", str.str());
                }
            }
        }

        writeTimeStepWithCellProperties(timer, localState, localWellState, cellData, substep);
    }
}
#endif
