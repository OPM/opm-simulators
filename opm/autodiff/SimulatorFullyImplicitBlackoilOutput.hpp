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
#ifndef OPM_SIMULATORFULLYIMPLICITBLACKOILOUTPUT_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITBLACKOILOUTPUT_HEADER_INCLUDED
#include <opm/core/grid.h>
#include <opm/core/simulator/SimulatorTimerInterface.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/autodiff/Compat.hpp>
#include <opm/core/utility/DataMap.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/output/eclipse/EclipseReader.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/wells/DynamicListEconLimited.hpp>

#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/eclipse/EclipseWriter.hpp>

#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/ParallelDebugOutput.hpp>

#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/ThreadHandle.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>


#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <thread>

#include <boost/filesystem.hpp>

#ifdef HAVE_OPM_GRID
#include <dune/grid/CpGrid.hpp>
#endif
namespace Opm
{

    class SimulationDataContainer;
    class BlackoilState;

    void outputStateVtk(const UnstructuredGrid& grid,
                        const Opm::SimulationDataContainer& state,
                        const int step,
                        const std::string& output_dir);


    void outputStateMatlab(const UnstructuredGrid& grid,
                           const Opm::SimulationDataContainer& state,
                           const int step,
                           const std::string& output_dir);

    void outputWellStateMatlab(const Opm::WellState& well_state,
                               const int step,
                               const std::string& output_dir);
#ifdef HAVE_OPM_GRID
    void outputStateVtk(const Dune::CpGrid& grid,
                        const Opm::SimulationDataContainer& state,
                        const int step,
                        const std::string& output_dir);
#endif

    template<class Grid>
    void outputStateMatlab(const Grid& grid,
                           const Opm::SimulationDataContainer& state,
                           const int step,
                           const std::string& output_dir)
    {
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        for (const auto& pair : state.cellData())
        {
            const std::string& name = pair.first;
            std::string key;
            if( name == "SURFACEVOL" ) {
                key = "surfvolume";
            }
            else if( name == "RV" ) {
                key = "rv";
            }
            else if( name == "GASOILRATIO" ) {
                key = "rs";
            }
            else { // otherwise skip entry
                continue;
            }
            // set data to datmap
            dm[ key ] = &pair.second;
        }

        std::vector<double> cell_velocity;
        Opm::estimateCellVelocity(AutoDiffGrid::numCells(grid),
                                  AutoDiffGrid::numFaces(grid),
                                  AutoDiffGrid::beginFaceCentroids(grid),
                                  UgGridHelpers::faceCells(grid),
                                  AutoDiffGrid::beginCellCentroids(grid),
                                  AutoDiffGrid::beginCellVolumes(grid),
                                  AutoDiffGrid::dimensions(grid),
                                  state.faceflux(), cell_velocity);
        dm["velocity"] = &cell_velocity;

        // Write data (not grid) in Matlab format
        for (Opm::DataMap::const_iterator it = dm.begin(); it != dm.end(); ++it) {
            std::ostringstream fname;
            fname << output_dir << "/" << it->first;
            boost::filesystem::path fpath = fname.str();
            try {
                create_directories(fpath);
            }
            catch (...) {
                OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
            }
            fname << "/" << std::setw(3) << std::setfill('0') << step << ".txt";
            std::ofstream file(fname.str().c_str());
            if (!file) {
                OPM_THROW(std::runtime_error, "Failed to open " << fname.str());
            }
            file.precision(15);
            const std::vector<double>& d = *(it->second);
            std::copy(d.begin(), d.end(), std::ostream_iterator<double>(file, "\n"));
        }
    }

    class BlackoilSubWriter {
        public:
            BlackoilSubWriter( const std::string& outputDir )
                : outputDir_( outputDir )
        {}

        virtual void writeTimeStep(const SimulatorTimerInterface& timer,
                           const SimulationDataContainer& state,
                           const WellState&,
                           bool /*substep*/ = false) = 0;
        protected:
            const std::string outputDir_;
    };

    template< class Grid >
    class BlackoilVTKWriter : public BlackoilSubWriter {
        public:
            BlackoilVTKWriter( const Grid& grid,
                               const std::string& outputDir )
                : BlackoilSubWriter( outputDir )
                , grid_( grid )
        {}

            void writeTimeStep(const SimulatorTimerInterface& timer,
                    const SimulationDataContainer& state,
                    const WellState&,
                    bool /*substep*/ = false) override
            {
                outputStateVtk(grid_, state, timer.currentStepNum(), outputDir_);
            }

        protected:
            const Grid& grid_;
    };

    template< typename Grid >
    class BlackoilMatlabWriter : public BlackoilSubWriter
    {
        public:
            BlackoilMatlabWriter( const Grid& grid,
                             const std::string& outputDir )
                : BlackoilSubWriter( outputDir )
                , grid_( grid )
        {}

        void writeTimeStep(const SimulatorTimerInterface& timer,
                           const SimulationDataContainer& reservoirState,
                           const WellState& wellState,
                           bool /*substep*/ = false) override
        {
            outputStateMatlab(grid_, reservoirState, timer.currentStepNum(), outputDir_);
            outputWellStateMatlab(wellState, timer.currentStepNum(), outputDir_);
        }

        protected:
            const Grid& grid_;
    };

    /** \brief Wrapper class for VTK, Matlab, and ECL output. */
    class BlackoilOutputWriter
    {

    public:
        // constructor creating different sub writers
        template <class Grid>
        BlackoilOutputWriter(const Grid& grid,
                             const parameter::ParameterGroup& param,
                             const Opm::EclipseState& eclipseState,
                             std::unique_ptr<EclipseWriter>&& eclWriter,
                             const Opm::PhaseUsage &phaseUsage,
                             const double* permeability );

        /** \copydoc Opm::OutputWriter::writeInit */
        void writeInit(const data::Solution& simProps, const NNC& nnc);

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
                           const data::Solution& solution,
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
         * \brief Write a blackoil reservoir state to disk for later inspection with
         *        visualization tools like ResInsight. This is the function which does
         *        the actual write to file.
         */
        void writeTimeStepSerial(const SimulatorTimerInterface& timer,
                                 const SimulationDataContainer& reservoirState,
                                 const Opm::WellState& wellState,
                                 const data::Solution& cellData,
                                 bool substep);

        /** \brief return output directory */
        const std::string& outputDirectory() const { return outputDir_; }

        /** \brief return true if output is enabled */
        bool output () const { return output_; }

        /** \brief Whether this process does write to disk */
        bool isIORank () const
        {
            return parallelOutput_->isIORank();
        }

        void restore(SimulatorTimerInterface& timer,
                     BlackoilState& state,
                     WellStateFullyImplicitBlackoil& wellState,
                     const std::string& filename,
                     const int desiredReportStep);


        template <class Grid>
        void initFromRestartFile(const PhaseUsage& phaseusage,
                                 const double* permeability,
                                 const Grid& grid,
                                 SimulationDataContainer& simulatorstate,
                                 WellStateFullyImplicitBlackoil& wellstate);

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
        std::unique_ptr< BlackoilSubWriter > vtkWriter_;
        std::unique_ptr< BlackoilSubWriter > matlabWriter_;
        std::unique_ptr< EclipseWriter > eclWriter_;
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
    BlackoilOutputWriter::
    BlackoilOutputWriter(const Grid& grid,
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
        if ( output_ )
        {
            if ( param.getDefault("output_vtk",false) )
            {
                vtkWriter_
                    .reset(new BlackoilVTKWriter< Grid >( grid, outputDir_ ));
            }

            if( parallelOutput_->isIORank() ) {

                if ( param.getDefault("output_matlab", false ) )
                {
                    matlabWriter_
                        .reset(new BlackoilMatlabWriter< Grid >( grid, outputDir_ ));
                }

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
    }


    template <class Grid>
    inline void
    BlackoilOutputWriter::
    initFromRestartFile( const PhaseUsage& phaseusage,
                         const double* permeability,
                         const Grid& grid,
                         SimulationDataContainer& simulatorstate,
                         WellStateFullyImplicitBlackoil& wellstate)
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

        /**
         * Converts an ADB::V into a standard vector by copy
         */
        inline std::vector<double> adbVToDoubleVector(const Opm::AutoDiffBlock<double>::V& adb_v) {
            std::vector<double> vec(adb_v.data(), adb_v.data() + adb_v.size());
            return vec;
        }


        /**
         * Converts an ADB into a standard vector by copy
         */
        inline std::vector<double> adbToDoubleVector(const Opm::AutoDiffBlock<double>& adb) {
            return adbVToDoubleVector(adb.value());
        }


        /**
         * Returns the data requested in the restartConfig
         */
        template<class Model>
        void getRestartData(data::Solution& output,
                            const Opm::PhaseUsage& phaseUsage,
                            const Model& physicalModel,
                            const RestartConfig& restartConfig,
                            const int reportStepNum,
                            const bool log) {

            typedef Opm::AutoDiffBlock<double> ADB;

            const typename Model::SimulatorData& sd = physicalModel.getSimulatorData();
            //Get the value of each of the keys for the restart keywords
            std::map<std::string, int> rstKeywords = restartConfig.getRestartKeywords(reportStepNum);
            for (auto& keyValue : rstKeywords) {
                keyValue.second = restartConfig.getKeyword(keyValue.first, reportStepNum);
            }

            //Get shorthands for water, oil, gas
            const int aqua_active = phaseUsage.phase_used[Opm::PhaseUsage::Aqua];
            const int liquid_active = phaseUsage.phase_used[Opm::PhaseUsage::Liquid];
            const int vapour_active = phaseUsage.phase_used[Opm::PhaseUsage::Vapour];

            const int aqua_idx = phaseUsage.phase_pos[Opm::PhaseUsage::Aqua];
            const int liquid_idx = phaseUsage.phase_pos[Opm::PhaseUsage::Liquid];
            const int vapour_idx = phaseUsage.phase_pos[Opm::PhaseUsage::Vapour];


            /**
             * Formation volume factors for water, oil, gas
             */
            if (aqua_active && rstKeywords["BW"] > 0) {
                rstKeywords["BW"] = 0;
                output.insert("1OVERBW",
                              Opm::UnitSystem::measure::water_inverse_formation_volume_factor,
                              adbToDoubleVector(sd.rq[aqua_idx].b),
                              data::TargetType::RESTART_AUXILLARY);
            }
            if (liquid_active && rstKeywords["BO"]  > 0) {
                rstKeywords["BO"] = 0;
                output.insert("1OVERBO",
                              Opm::UnitSystem::measure::oil_inverse_formation_volume_factor,
                              adbToDoubleVector(sd.rq[liquid_idx].b),
                              data::TargetType::RESTART_AUXILLARY);
            }
            if (vapour_active && rstKeywords["BG"] > 0) {
                rstKeywords["BG"] = 0;
                output.insert("1OVERBG",
                              Opm::UnitSystem::measure::gas_inverse_formation_volume_factor,
                              adbToDoubleVector(sd.rq[vapour_idx].b),
                              data::TargetType::RESTART_AUXILLARY);
            }

            /**
             * Densities for water, oil gas
             */
            if (rstKeywords["DEN"] > 0) {
                rstKeywords["DEN"] = 0;
                if (aqua_active) {
                    output.insert("WAT_DEN",
                                  Opm::UnitSystem::measure::density,
                                  adbToDoubleVector(sd.rq[aqua_idx].rho),
                                  data::TargetType::RESTART_AUXILLARY);
                }
                if (liquid_active) {
                    output.insert("OIL_DEN",
                                  Opm::UnitSystem::measure::density,
                                  adbToDoubleVector(sd.rq[liquid_idx].rho),
                                  data::TargetType::RESTART_AUXILLARY);
                }
                if (vapour_active) {
                    output.insert("GAS_DEN",
                                  Opm::UnitSystem::measure::density,
                                  adbToDoubleVector(sd.rq[vapour_idx].rho),
                                  data::TargetType::RESTART_AUXILLARY);
                }
            }

            /**
             * Viscosities for water, oil gas
             */
            if (rstKeywords["VISC"] > 0) {
                rstKeywords["VISC"] = 0;
                if (aqua_active) {
                    output.insert("WAT_VISC",
                                  Opm::UnitSystem::measure::viscosity,
                                  adbToDoubleVector(sd.rq[aqua_idx].mu),
                                  data::TargetType::RESTART_AUXILLARY);
                }
                if (liquid_active) {
                    output.insert("OIL_VISC",
                                  Opm::UnitSystem::measure::viscosity,
                                  adbToDoubleVector(sd.rq[liquid_idx].mu),
                                  data::TargetType::RESTART_AUXILLARY);
                }
                if (vapour_active) {
                    output.insert("GAS_VISC",
                                  Opm::UnitSystem::measure::viscosity,
                                  adbToDoubleVector(sd.rq[vapour_idx].mu),
                                  data::TargetType::RESTART_AUXILLARY);
                }
            }

            /**
             * Relative permeabilities for water, oil, gas
             */
            if (aqua_active && rstKeywords["KRW"] > 0) {
                if (sd.rq[aqua_idx].kr.size() > 0) {
                    rstKeywords["KRW"] = 0;
                    output.insert("WATKR",
                                  Opm::UnitSystem::measure::permeability,
                                  adbToDoubleVector(sd.rq[aqua_idx].kr),
                                  data::TargetType::RESTART_AUXILLARY);
                }
                else {
                    if ( log )
                    {
                        Opm::OpmLog::warning("Empty:WATKR",
                                             "Not emitting empty Water Rel-Perm");
                    }
                }
            }
            if (liquid_active && rstKeywords["KRO"] > 0) {
                if (sd.rq[liquid_idx].kr.size() > 0) {
                    rstKeywords["KRO"] = 0;
                    output.insert("OILKR",
                                  Opm::UnitSystem::measure::permeability,
                                  adbToDoubleVector(sd.rq[liquid_idx].kr),
                                  data::TargetType::RESTART_AUXILLARY);
                }
                else {
                    if ( log )
                    {
                        Opm::OpmLog::warning("Empty:OILKR",
                                             "Not emitting empty Oil Rel-Perm");
                    }
                }
            }
            if (vapour_active && rstKeywords["KRG"] > 0) {
                if (sd.rq[vapour_idx].kr.size() > 0) {
                    rstKeywords["KRG"] = 0;
                    output.insert("GASKR",
                                  Opm::UnitSystem::measure::permeability,
                                  adbToDoubleVector(sd.rq[vapour_idx].kr),
                                  data::TargetType::RESTART_AUXILLARY);
                }
                else {
                    if ( log )
                    {
                        Opm::OpmLog::warning("Empty:GASKR",
                                             "Not emitting empty Gas Rel-Perm");
                    }
                }
            }

            /**
             * Vaporized and dissolved gas/oil ratio
             */
            if (vapour_active && liquid_active && rstKeywords["RSSAT"] > 0) {
                rstKeywords["RSSAT"] = 0;
                output.insert("RSSAT",
                              Opm::UnitSystem::measure::gas_oil_ratio,
                              adbToDoubleVector(sd.rsSat),
                              data::TargetType::RESTART_AUXILLARY);
            }
            if (vapour_active && liquid_active && rstKeywords["RVSAT"] > 0) {
                rstKeywords["RVSAT"] = 0;
                output.insert("RVSAT",
                              Opm::UnitSystem::measure::oil_gas_ratio,
                              adbToDoubleVector(sd.rvSat),
                              data::TargetType::RESTART_AUXILLARY);
            }


            /**
             * Bubble point and dew point pressures
             */
            if (log && vapour_active && 
                liquid_active && rstKeywords["PBPD"] > 0) {
                rstKeywords["PBPD"] = 0;
                Opm::OpmLog::warning("Bubble/dew point pressure output unsupported",
                        "Writing bubble points and dew points (PBPD) to file is unsupported, "
                        "as the simulator does not use these internally.");
            }

            //Warn for any unhandled keyword
            if (log) {
                for (auto& keyValue : rstKeywords) {
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

            typedef Opm::AutoDiffBlock<double> ADB;

            const typename Model::SimulatorData& sd = physicalModel.getSimulatorData();

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
                              adbVToDoubleVector(sd.fip[Model::SimulatorData::FIP_AQUA]),
                              data::TargetType::SUMMARY );
            }
            if (liquid_active) {
                const ADB::V& oipl = sd.fip[Model::SimulatorData::FIP_LIQUID];
                const ADB::V& oipg = vapour_active ? sd.fip[Model::SimulatorData::FIP_VAPORIZED_OIL] : ADB::V();
                const ADB::V& oip = vapour_active ? oipl + oipg : oipl;

                //Oil in place (liquid phase only)
                if (hasFRBKeyword(summaryConfig, "OIPL")) {
                    output.insert("OIPL",
                                  Opm::UnitSystem::measure::volume,
                                  adbVToDoubleVector(oipl),
                                  data::TargetType::SUMMARY );
                }
                //Oil in place (gas phase only)
                if (hasFRBKeyword(summaryConfig, "OIPG")) {
                    output.insert("OIPG",
                                  Opm::UnitSystem::measure::volume,
                                  adbVToDoubleVector(oipg),
                                  data::TargetType::SUMMARY );
                }
                // Oil in place (in liquid and gas phases)
                if (hasFRBKeyword(summaryConfig, "OIP")) {
                    output.insert("OIP",
                                  Opm::UnitSystem::measure::volume,
                                  adbVToDoubleVector(oip),
                                  data::TargetType::SUMMARY );
                }
            }
            if (vapour_active) {
                const ADB::V& gipg = sd.fip[Model::SimulatorData::FIP_VAPOUR];
                const ADB::V& gipl = liquid_active ? sd.fip[Model::SimulatorData::FIP_DISSOLVED_GAS] : ADB::V();
                const ADB::V& gip = liquid_active ? gipg + gipl : gipg;

                // Gas in place (gas phase only)
                if (hasFRBKeyword(summaryConfig, "GIPG")) {
                    output.insert("GIPG",
                                  Opm::UnitSystem::measure::volume,
                                  adbVToDoubleVector(gipg),
                                  data::TargetType::SUMMARY );
                }
                // Gas in place (liquid phase only)
                if (hasFRBKeyword(summaryConfig, "GIPL")) {
                    output.insert("GIPL",
                                  Opm::UnitSystem::measure::volume,
                                  adbVToDoubleVector(gipl),
                                  data::TargetType::SUMMARY );
                }
                // Gas in place (in both liquid and gas phases)
                if (hasFRBKeyword(summaryConfig, "GIP")) {
                    output.insert("GIP",
                                  Opm::UnitSystem::measure::volume,
                                  adbVToDoubleVector(gip),
                                  data::TargetType::SUMMARY );
                }
            }
            // Cell pore volume in reservoir conditions
            if (hasFRBKeyword(summaryConfig, "RPV")) {
                output.insert("RPV",
                              Opm::UnitSystem::measure::volume,
                              adbVToDoubleVector(sd.fip[Model::SimulatorData::FIP_PV]),
                              data::TargetType::SUMMARY );
            }
            // Pressure averaged value (hydrocarbon pore volume weighted)
            if (summaryConfig.hasKeyword("FPRH") || summaryConfig.hasKeyword("RPRH")) {
                output.insert("PRH",
                              Opm::UnitSystem::measure::pressure,
                              adbVToDoubleVector(sd.fip[Model::SimulatorData::FIP_WEIGHTED_PRESSURE]),
                              data::TargetType::SUMMARY );
            }
        }

    }




    template<class Model>
    inline void
    BlackoilOutputWriter::
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

            detail::getRestartData( cellData, phaseUsage_, physicalModel,
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
