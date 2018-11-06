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
#include <opm/grid/UnstructuredGrid.h>
#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/autodiff/Compat.hpp>
#include <opm/core/utility/DataMap.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/core/wells/DynamicListEconLimited.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>

#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>

#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/ParallelDebugOutput.hpp>

#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/ThreadHandle.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/simulators/ensureDirectoryExists.hpp>

#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <thread>
#include <map>

#include <boost/filesystem.hpp>

#ifdef HAVE_OPM_GRID
#include <opm/grid/CpGrid.hpp>
#endif
namespace Opm
{

    class SimulationDataContainer;
    class BlackoilState;

    void outputStateVtk(const UnstructuredGrid& grid,
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
            ensureDirectoryExists(fname.str());
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
                           const WellStateFullyImplicitBlackoil&,
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
                    const WellStateFullyImplicitBlackoil&,
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
                           const WellStateFullyImplicitBlackoil& wellState,
                           bool /*substep*/ = false) override
        {
            outputStateMatlab(grid_, reservoirState, timer.currentStepNum(), outputDir_);
            outputWellStateMatlab(wellState, timer.currentStepNum(), outputDir_);
        }

        protected:
            const Grid& grid_;
    };


    /// Extra data to read/write for OPM restarting
    struct ExtraData
    {
        double suggested_step = -1.0;
    };


    /** \brief Wrapper class for VTK, Matlab, and ECL output. */
    class BlackoilOutputWriter
    {

    public:
        // constructor creating different sub writers
        template <class Grid>
        BlackoilOutputWriter(const Grid& grid,
                             const ParameterGroup& param,
                             const Opm::EclipseState& eclipseState,
                             const Opm::Schedule& schedule,
                             const Opm::SummaryConfig& summaryConfig,
                             std::unique_ptr<EclipseIO>&& eclIO,
                             const Opm::PhaseUsage &phaseUsage);

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
                           const Opm::WellStateFullyImplicitBlackoil& wellState,
                           const Model& physicalModel,
                           const bool substep = false,
                           const double nextstep = -1.0,
                           const SimulatorReport& simulatorReport = SimulatorReport());


        /*!
         * \brief Write a blackoil reservoir state to disk for later inspection with
         *        visualization tools like ResInsight. This function will write all
         *        CellData in simProps to the file as well as the extraRestartData.
         */
        void writeTimeStepWithCellProperties(
                           const SimulatorTimerInterface& timer,
                           const SimulationDataContainer& reservoirState,
                           const data::Solution& cellData,
                           const Opm::WellStateFullyImplicitBlackoil& wellState,
                           const std::map<std::string, double>& miscSummaryData,
                           const RestartValue::ExtraVector& extraRestartData,
                           bool substep = false);

        /*!
         * \brief Write a blackoil reservoir state to disk for later inspection with
         *        visualization tools like ResInsight. This function will not write
         *        any cell properties (e.g., those requested by RPTRST keyword)
         */
        void writeTimeStepWithoutCellProperties(
                           const SimulatorTimerInterface& timer,
                           const SimulationDataContainer& reservoirState,
                           const Opm::WellStateFullyImplicitBlackoil& wellState,
                           const std::map<std::string, double>& miscSummaryData,
                           const RestartValue::ExtraVector& extraRestartData,
                           bool substep = false);

        /*!
         * \brief Write a blackoil reservoir state to disk for later inspection with
         *        visualization tools like ResInsight. This is the function which does
         *        the actual write to file.
         */
        void writeTimeStepSerial(const SimulatorTimerInterface& timer,
                                 const SimulationDataContainer& reservoirState,
                                 const Opm::WellStateFullyImplicitBlackoil& wellState,
                                 const data::Solution& simProps,
                                 const std::map<std::string, double>& miscSummaryData,
                                 const RestartValue::ExtraVector& extraRestartData,
                                 bool substep );

        /** \brief return output directory */
        const std::string& outputDirectory() const { return outputDir_; }

        /** \brief return true if output is enabled */
        bool output () const { return output_; }

        /** \brief Whether this process does write to disk */
        bool isIORank () const
        {
            return parallelOutput_->isIORank();
        }

        template <class Grid, class WellState>
        void initFromRestartFile(const PhaseUsage& phaseUsage,
                                 const Grid& grid,
                                 SimulationDataContainer& simulatorstate,
                                 WellState& wellstate,
                                 ExtraData& extra);

        bool isRestart() const;

        bool requireFIPNUM() const;

    protected:
        const bool output_;
        std::unique_ptr< ParallelDebugOutputInterface > parallelOutput_;

        // Parameters for output.
        const std::string outputDir_;
        const bool restart_double_si_;

        Opm::PhaseUsage phaseUsage_;
        std::unique_ptr< BlackoilSubWriter > vtkWriter_;
        std::unique_ptr< BlackoilSubWriter > matlabWriter_;
        std::unique_ptr< EclipseIO > eclIO_;
        const EclipseState& eclipseState_;
        const Schedule& schedule_;
        const SummaryConfig& summaryConfig_;

        std::unique_ptr< ThreadHandle > asyncOutput_;
        const int* globalCellIdxMap_;
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
                         const ParameterGroup& param,
                         const Opm::EclipseState& eclipseState,
                         const Opm::Schedule& schedule,
                         const Opm::SummaryConfig& summaryConfig,
                         std::unique_ptr<EclipseIO>&& eclIO,
                         const Opm::PhaseUsage &phaseUsage)
        : output_( [ &param ] () -> bool {
                // If output parameter is true or all, then we do output
                const std::string outputString = param.getDefault("output", std::string("all"));
                return ( outputString == "all" ||  outputString == "true" );
            }()
            ),
        parallelOutput_( output_ ? new ParallelDebugOutput< Grid >( grid, eclipseState, schedule, phaseUsage.num_phases, phaseUsage ) : 0 ),
        outputDir_( eclipseState.getIOConfig().getOutputDir() ),
        restart_double_si_( output_ ? param.getDefault("restart_double_si", false) : false ),
        phaseUsage_( phaseUsage ),
        eclipseState_(eclipseState),
        schedule_(schedule),
        summaryConfig_(summaryConfig),
        asyncOutput_(),
        globalCellIdxMap_(Opm::UgGridHelpers::globalCell(grid))
    {
        // For output.
        if ( output_ )
        {
            if ( param.getDefault("output_vtk",false) )
            {
                vtkWriter_
                    .reset(new BlackoilVTKWriter< Grid >( grid, outputDir_ ));
            }

            auto output_matlab = param.getDefault("output_matlab", false );

            if ( parallelOutput_->isParallel() && output_matlab )
            {
                Opm::OpmLog::warning("Parallel Output Config",
                                     "Velocity output for matlab is broken in parallel.");
            }

            if( parallelOutput_->isIORank() ) {

                if ( output_matlab )
                {
                    matlabWriter_
                        .reset(new BlackoilMatlabWriter< Grid >( grid, outputDir_ ));
                }

                eclIO_ = std::move(eclIO);

                // Ensure that output dir exists
                ensureDirectoryExists(outputDir_);
            }

            // create output thread if enabled and rank is I/O rank
            // async output is enabled by default if pthread are enabled
#if HAVE_PTHREAD
            const bool asyncOutputDefault = true;
#else
            const bool asyncOutputDefault = false;
#endif
            if( param.getDefault("async_output", asyncOutputDefault ) )
            {
                const bool isIORank = parallelOutput_ ? parallelOutput_->isIORank() : true;
#if HAVE_PTHREAD
                asyncOutput_.reset( new ThreadHandle( isIORank ) );
#else
                OPM_THROW(std::runtime_error,"Pthreads were not found, cannot enable async_output");
#endif
            }
        }
    }


    template <class Grid, class WellState>
    inline void
    BlackoilOutputWriter::
    initFromRestartFile( const PhaseUsage& phaseUsage,
                         const Grid& grid,
                         SimulationDataContainer& simulatorstate,
                         WellState& wellstate,
                         ExtraData& extra )
    {
        std::vector<RestartKey> solution_keys = {{"PRESSURE" , UnitSystem::measure::pressure},
                                                 {"SWAT" , UnitSystem::measure::identity},
                                                 {"SGAS" , UnitSystem::measure::identity},
                                                 {"TEMP" , UnitSystem::measure::temperature},
                                                 {"RS" , UnitSystem::measure::gas_oil_ratio},
                                                 {"RV" , UnitSystem::measure::oil_gas_ratio},
                                                 {"SOMAX", UnitSystem::measure::identity, false},
                                                 {"PCSWM_OW", UnitSystem::measure::identity, false},
                                                 {"KRNSW_OW", UnitSystem::measure::identity, false},
                                                 {"PCSWM_GO", UnitSystem::measure::identity, false},
                                                 {"KRNSW_GO", UnitSystem::measure::identity, false}};

        std::vector<RestartKey> extra_keys {
          {"OPMEXTRA" , UnitSystem::measure::identity, false}
        };

        if (restart_double_si_) {
            // Avoid any unit conversions, treat restart input as SI units.
            for (auto& elem : solution_keys) {
                elem.dim = UnitSystem::measure::identity;
            }
        }

        // gives a dummy dynamic_list_econ_limited
        DynamicListEconLimited dummy_list_econ_limited;
        WellsManager wellsmanager(eclipseState_,
                                  schedule_,
                                  // The restart step value is used to identify wells present at the given time step.
                                  // Wells that are added at the same time step as RESTART is initiated will not be
                                  // present in a restart file. Use the previous time step to retrieve wells
                                  // that have information written to the restart file.
                                  std::max(eclipseState_.getInitConfig().getRestartStep() - 1, 0),
                                  Opm::UgGridHelpers::numCells(grid),
                                  Opm::UgGridHelpers::globalCell(grid),
                                  Opm::UgGridHelpers::cartDims(grid),
                                  Opm::UgGridHelpers::dimensions(grid),
                                  Opm::UgGridHelpers::cell2Faces(grid),
                                  Opm::UgGridHelpers::beginFaceCentroids(grid),
                                  dummy_list_econ_limited
                                  // We need to pass the optionaly arguments
                                  // as we get the following error otherwise
                                  // with c++ (Debian 4.9.2-10) 4.9.2 and -std=c++11
                                  // converting to ‘const std::unordered_set<std::basic_string<char> >’ from initializer list would use explicit constructo
                                  , false,
                                  std::unordered_set<std::string>());

        const Wells* wells = wellsmanager.c_wells();
        size_t numCells = Opm::UgGridHelpers::numCells(grid);
        wellstate.resize(wells, numCells); //Resize for restart step
        auto restart_values = eclIO_->loadRestart(solution_keys, extra_keys);

        solutionToSim( restart_values, phaseUsage, simulatorstate );
        wellsToState( restart_values.wells, phaseUsage, wellstate );

        if (restart_values.hasExtra("OPMEXTRA")) {
            const std::vector<double>& opmextra = restart_values.getExtra("OPMEXTRA");
            assert(opmextra.size() == 1);
            extra.suggested_step = opmextra[0];
        } else {
            OpmLog::warning("Restart data is missing OPMEXTRA field, restart run may deviate from original run.");
            extra.suggested_step = -1.0;
        }
    }





    namespace detail {


        template <class V>
        void addToSimData( SimulationDataContainer& simData,
                           const std::string& name,
                           const V& vec )
        {
            if (vec.size() == 0) {
                return;
            }

            typedef std::vector< double > OutputVectorType;

            // get data map
            auto& dataMap = simData.cellData();

            // insert name,vector into data map
            dataMap.insert( std::make_pair( name, OutputVectorType( vec.data(), vec.data() + vec.size() ) ) );
        }

        template <class Scalar>
        void addToSimData( SimulationDataContainer& simData,
                           const std::string& name,
                           const AutoDiffBlock<Scalar>& adb )
        {
            // forward value of ADB to output
            addToSimData( simData, name, adb.value() );
        }


        // this method basically converts all Eigen vectors to std::vectors
        // stored in a SimulationDataContainer
        template <class SimulatorData>
        SimulationDataContainer
        convertToSimulationDataContainer( const SimulatorData& sd,
                                          const SimulationDataContainer& localState,
                                          const Opm::PhaseUsage& phaseUsage )
        {
            // copy local state and then add missing data
            SimulationDataContainer simData( localState );

            //Get shorthands for water, oil, gas
            const int aqua_active   = phaseUsage.phase_used[Opm::PhaseUsage::Aqua];
            const int liquid_active = phaseUsage.phase_used[Opm::PhaseUsage::Liquid];
            const int vapour_active = phaseUsage.phase_used[Opm::PhaseUsage::Vapour];

            const int aqua_idx   = phaseUsage.phase_pos[Opm::PhaseUsage::Aqua];
            const int liquid_idx = phaseUsage.phase_pos[Opm::PhaseUsage::Liquid];
            const int vapour_idx = phaseUsage.phase_pos[Opm::PhaseUsage::Vapour];

            // WATER
            if( aqua_active ) {
                addToSimData( simData, "1OVERBW",  sd.rq[aqua_idx].b   );
                addToSimData( simData, "WAT_DEN",  sd.rq[aqua_idx].rho );
                addToSimData( simData, "WAT_VISC", sd.rq[aqua_idx].mu  );
                addToSimData( simData, "WATKR",    sd.rq[aqua_idx].kr  );
            }

            // OIL
            if( liquid_active ) {
                addToSimData( simData, "1OVERBO",  sd.rq[liquid_idx].b   );
                addToSimData( simData, "OIL_DEN",  sd.rq[liquid_idx].rho );
                addToSimData( simData, "OIL_VISC", sd.rq[liquid_idx].mu  );
                addToSimData( simData, "OILKR",    sd.rq[liquid_idx].kr  );
            }

            // GAS
            if( vapour_active ) {
                addToSimData( simData, "1OVERBG",  sd.rq[vapour_idx].b   );
                addToSimData( simData, "GAS_DEN",  sd.rq[vapour_idx].rho );
                addToSimData( simData, "GAS_VISC", sd.rq[vapour_idx].mu  );
                addToSimData( simData, "GASKR",    sd.rq[vapour_idx].kr  );
            }

            // RS and RV
            addToSimData( simData, "RSSAT", sd.rsSat );
            addToSimData( simData, "RVSAT", sd.rvSat );

            addToSimData( simData, "SOMAX", sd.soMax );
            addToSimData( simData, "PBUB", sd.Pb );
            addToSimData( simData, "PDEW", sd.Pd );
            addToSimData( simData, "PCSWMDC_OW", sd.pcswmdc_ow );
            addToSimData( simData, "KRNSWMDC_OW", sd.krnswdc_ow );
            addToSimData( simData, "PCSWMDC_GO", sd.pcswmdc_go );
            addToSimData( simData, "KRNSWMDC_GO", sd.krnswdc_go );

            return simData;
        }

        // in case the data is already in a SimulationDataContainer no
        // conversion is needed
        inline
        SimulationDataContainer&&
        convertToSimulationDataContainer( SimulationDataContainer&& sd,
                                          const SimulationDataContainer& ,
                                          const Opm::PhaseUsage& )
        {
            return std::move( sd );
        }

        /**
         * Returns the data requested in the restartConfig
         * NOTE: Since this function steals data from the SimulationDataContainer (std::move),
         * the variable sd becomes "invalid" after calling this function.
         */
        template<class Model>
        void getRestartData(data::Solution& output,
                            SimulationDataContainer&& sd,
                            const Opm::PhaseUsage& /* phaseUsage */,
                            const Model& /* physicalModel */,
                            const RestartConfig& restartConfig,
                            const int reportStepNum,
                            const bool log)
        {
            //Get the value of each of the keys for the restart keywords
            std::map<std::string, int> rstKeywords = restartConfig.getRestartKeywords(reportStepNum);
            for (auto& keyValue : rstKeywords) {
                keyValue.second = restartConfig.getKeyword(keyValue.first, reportStepNum);
            }

            const bool aqua_active   = sd.hasCellData("1OVERBW");
            const bool liquid_active = sd.hasCellData("1OVERBO");
            const bool vapour_active = sd.hasCellData("1OVERBG");

            assert( aqua_active == (sd.hasCellData("WAT_DEN")  &&
                                    sd.hasCellData("WAT_VISC") &&
                                    sd.hasCellData("WATKR")
                                   )
                  );
            assert( liquid_active == (sd.hasCellData("OIL_DEN")  &&
                                      sd.hasCellData("OIL_VISC") &&
                                      sd.hasCellData("OILKR")
                                   )
                  );
            assert( vapour_active == (sd.hasCellData("GAS_DEN")  &&
                                      sd.hasCellData("GAS_VISC") &&
                                      sd.hasCellData("GASKR")
                                   )
                  );

            /**
             * Formation volume factors for water, oil, gas
             */
            if (aqua_active && rstKeywords["BW"] > 0) {
                rstKeywords["BW"] = 0;
                output.insert("1OVERBW",
                              Opm::UnitSystem::measure::water_inverse_formation_volume_factor,
                              std::move( sd.getCellData("1OVERBW") ),
                              data::TargetType::RESTART_AUXILIARY);
            }
            if (liquid_active && rstKeywords["BO"]  > 0) {
                rstKeywords["BO"] = 0;
                output.insert("1OVERBO",
                              Opm::UnitSystem::measure::oil_inverse_formation_volume_factor,
                              std::move( sd.getCellData("1OVERBO") ),
                              data::TargetType::RESTART_AUXILIARY);
            }
            if (vapour_active && rstKeywords["BG"] > 0) {
                rstKeywords["BG"] = 0;
                output.insert("1OVERBG",
                              Opm::UnitSystem::measure::gas_inverse_formation_volume_factor,
                              std::move( sd.getCellData("1OVERBG") ),
                              data::TargetType::RESTART_AUXILIARY);
            }

            /**
             * Densities for water, oil gas
             */
            if (rstKeywords["DEN"] > 0) {
                rstKeywords["DEN"] = 0;
                if (aqua_active) {
                    output.insert("WAT_DEN",
                                  Opm::UnitSystem::measure::density,
                                  std::move( sd.getCellData("WAT_DEN") ),
                                  data::TargetType::RESTART_AUXILIARY);
                }
                if (liquid_active) {
                    output.insert("OIL_DEN",
                                  Opm::UnitSystem::measure::density,
                                  std::move( sd.getCellData("OIL_DEN") ),
                                  data::TargetType::RESTART_AUXILIARY);
                }
                if (vapour_active) {
                    output.insert("GAS_DEN",
                                  Opm::UnitSystem::measure::density,
                                  std::move( sd.getCellData("GAS_DEN") ),
                                  data::TargetType::RESTART_AUXILIARY);
                }
            }

            /**
             * Viscosities for water, oil gas
             */
            {
                const bool has_vwat = (rstKeywords["VISC"] > 0) || (rstKeywords["VWAT"] > 0);
                const bool has_voil = (rstKeywords["VISC"] > 0) || (rstKeywords["VOIL"] > 0);
                const bool has_vgas = (rstKeywords["VISC"] > 0) || (rstKeywords["VGAS"] > 0);
                rstKeywords["VISC"] = 0;
                if (aqua_active && has_vwat) {
                    output.insert("WAT_VISC",
                                  Opm::UnitSystem::measure::viscosity,
                                  std::move( sd.getCellData("WAT_VISC") ),
                                  data::TargetType::RESTART_AUXILIARY);
                    rstKeywords["VWAT"] = 0;
                }
                if (liquid_active && has_voil) {
                    output.insert("OIL_VISC",
                                  Opm::UnitSystem::measure::viscosity,
                                  std::move( sd.getCellData("OIL_VISC") ),
                                  data::TargetType::RESTART_AUXILIARY);
                    rstKeywords["VOIL"] = 0;
                }
                if (vapour_active && has_vgas) {
                    output.insert("GAS_VISC",
                                  Opm::UnitSystem::measure::viscosity,
                                  std::move( sd.getCellData("GAS_VISC") ),
                                  data::TargetType::RESTART_AUXILIARY);
                    rstKeywords["VGAS"] = 0;
                }
            }

            /**
             * Relative permeabilities for water, oil, gas
             */
            if (aqua_active && rstKeywords["KRW"] > 0) {
                auto& krWater = sd.getCellData("WATKR");
                if (krWater.size() > 0) {
                    rstKeywords["KRW"] = 0;
                    output.insert("WATKR", // WAT_KR ???
                                  Opm::UnitSystem::measure::identity,
                                  std::move( krWater ),
                                  data::TargetType::RESTART_AUXILIARY);
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
                auto& krOil = sd.getCellData("OILKR");
                if (krOil.size() > 0) {
                    rstKeywords["KRO"] = 0;
                    output.insert("OILKR",
                                  Opm::UnitSystem::measure::identity,
                                  std::move( krOil ),
                                  data::TargetType::RESTART_AUXILIARY);
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
                auto& krGas = sd.getCellData("GASKR");
                if (krGas.size() > 0) {
                    rstKeywords["KRG"] = 0;
                    output.insert("GASKR",
                                  Opm::UnitSystem::measure::identity,
                                  std::move( krGas ),
                                  data::TargetType::RESTART_AUXILIARY);
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
                              std::move( sd.getCellData("RSSAT") ),
                              data::TargetType::RESTART_AUXILIARY);
            }
            if (vapour_active && liquid_active && rstKeywords["RVSAT"] > 0) {
                rstKeywords["RVSAT"] = 0;
                output.insert("RVSAT",
                              Opm::UnitSystem::measure::oil_gas_ratio,
                              std::move( sd.getCellData("RVSAT") ),
                              data::TargetType::RESTART_AUXILIARY);
            }


            /**
             * Bubble point and dew point pressures
             */
            if (vapour_active && liquid_active && rstKeywords["PBPD"] > 0) {
                rstKeywords["PBPD"] = 0;
                if (sd.hasCellData("PBUB")) {
                    output.insert("PBUB",
                                  Opm::UnitSystem::measure::pressure,
                                  std::move( sd.getCellData("PBUB") ),
                                  data::TargetType::RESTART_AUXILIARY);
                }
                else if (log) {
                    Opm::OpmLog::warning("Bubble point pressure unavailable", "Output of bubble point pressure requested but not available in this simulator. Ignoring.");
                }

                if (sd.hasCellData("PDEW")) {
                    output.insert("PDEW",
                                  Opm::UnitSystem::measure::pressure,
                                  std::move( sd.getCellData("PDEW") ),
                                  data::TargetType::RESTART_AUXILIARY);
                }
                else if (log) {
                    Opm::OpmLog::warning("Dew point pressure unavailable", "Output of dew point pressure requested but not available in this simulator. Ignoring.");
                }

            }

            if (sd.hasCellData("SOMAX")) {
                output.insert("SOMAX",
                        Opm::UnitSystem::measure::identity,
                        std::move( sd.getCellData("SOMAX") ),
                        data::TargetType::RESTART_AUXILIARY);
            }

            if (sd.hasCellData("PCSWMDC_OW")) {
                output.insert("PCSWM_OW", //FIXME: Eight-long variable name
                        Opm::UnitSystem::measure::identity,
                        std::move( sd.getCellData("PCSWMDC_OW") ),
                        data::TargetType::RESTART_AUXILIARY);
            }
            if (sd.hasCellData("KRNSWMDC_OW")) {
                output.insert("KRNSW_OW",
                        Opm::UnitSystem::measure::identity,
                        std::move( sd.getCellData("KRNSWMDC_OW") ),
                        data::TargetType::RESTART_AUXILIARY);
            }

            if (sd.hasCellData("PCSWMDC_GO")) {
                output.insert("PCSWM_GO", //FIXME: Eight-long variable name
                        Opm::UnitSystem::measure::identity,
                        std::move( sd.getCellData("PCSWMDC_GO") ),
                        data::TargetType::RESTART_AUXILIARY);
            }
            if (sd.hasCellData("KRNSWMDC_GO")) {
                output.insert("KRNSW_GO",
                        Opm::UnitSystem::measure::identity,
                        std::move( sd.getCellData("KRNSWMDC_GO") ),
                        data::TargetType::RESTART_AUXILIARY);
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

            typedef typename Model::FIPDataType FIPDataType;
            typedef typename FIPDataType::VectorType VectorType;

            FIPDataType fd = physicalModel.getFIPData();

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
                              std::move( fd.fip[ FIPDataType::FIP_AQUA ] ),
                              data::TargetType::SUMMARY );
            }
            if (liquid_active) {
                const VectorType& oipl = fd.fip[FIPDataType::FIP_LIQUID];
                VectorType  oip ( oipl );
                const size_t size = oip.size();

                const VectorType& oipg = vapour_active ? fd.fip[FIPDataType::FIP_VAPORIZED_OIL] : VectorType(size, 0.0);
                if( vapour_active )
                {
                    // oip = oipl + oipg
                    for( size_t i=0; i<size; ++ i ) {
                        oip[ i ] += oipg[ i ];
                    }
                }

                //Oil in place (liquid phase only)
                if (hasFRBKeyword(summaryConfig, "OIPL")) {
                    output.insert("OIPL",
                                  Opm::UnitSystem::measure::volume,
                                  std::move( oipl ),
                                  data::TargetType::SUMMARY );
                }
                //Oil in place (gas phase only)
                if (hasFRBKeyword(summaryConfig, "OIPG")) {
                    output.insert("OIPG",
                                  Opm::UnitSystem::measure::volume,
                                  std::move( oipg ),
                                  data::TargetType::SUMMARY );
                }
                // Oil in place (in liquid and gas phases)
                if (hasFRBKeyword(summaryConfig, "OIP") || hasFRBKeyword(summaryConfig, "OE")) {
                    output.insert("OIP",
                                  Opm::UnitSystem::measure::volume,
                                  std::move( oip ),
                                  data::TargetType::SUMMARY );
                }
            }
            if (vapour_active) {
                const VectorType& gipg = fd.fip[ FIPDataType::FIP_VAPOUR];
                VectorType  gip( gipg );
                const size_t size = gip.size();

                const VectorType& gipl = liquid_active ? fd.fip[ FIPDataType::FIP_DISSOLVED_GAS ] : VectorType(size,0.0);
                if( liquid_active )
                {
                    // gip = gipg + gipl
                    for( size_t i=0; i<size; ++ i ) {
                        gip[ i ] += gipl[ i ];
                    }
                }

                // Gas in place (gas phase only)
                if (hasFRBKeyword(summaryConfig, "GIPG")) {
                    output.insert("GIPG",
                                  Opm::UnitSystem::measure::volume,
                                  std::move( gipg ),
                                  data::TargetType::SUMMARY );
                }

                // Gas in place (liquid phase only)
                if (hasFRBKeyword(summaryConfig, "GIPL")) {
                    output.insert("GIPL",
                                  Opm::UnitSystem::measure::volume,
                                  std::move( gipl ),
                                  data::TargetType::SUMMARY );
                }
                // Gas in place (in both liquid and gas phases)
                if (hasFRBKeyword(summaryConfig, "GIP")) {
                    output.insert("GIP",
                                  Opm::UnitSystem::measure::volume,
                                  std::move( gip ),
                                  data::TargetType::SUMMARY );
                }
            }
            // Cell pore volume in reservoir conditions
            if (hasFRBKeyword(summaryConfig, "RPV")) {
                output.insert("RPV",
                              Opm::UnitSystem::measure::volume,
                              std::move( fd.fip[FIPDataType::FIP_PV]),
                              data::TargetType::SUMMARY );
            }
            // Pressure averaged value (hydrocarbon pore volume weighted)
            if (summaryConfig.hasKeyword("FPRH") || summaryConfig.hasKeyword("RPRH")) {
                output.insert("PRH",
                              Opm::UnitSystem::measure::pressure,
                              std::move(fd.fip[FIPDataType::FIP_WEIGHTED_PRESSURE]),
                              data::TargetType::SUMMARY );
            }
        }

    }




    template<class Model>
    inline void
    BlackoilOutputWriter::
    writeTimeStep(const SimulatorTimerInterface& timer,
                  const SimulationDataContainer& localState,
                  const WellStateFullyImplicitBlackoil& localWellState,
                  const Model& physicalModel,
                  const bool substep,
                  const double nextstep,
                  const SimulatorReport& simulatorReport)
    {
        data::Solution localCellData{};
        const RestartConfig& restartConfig = eclipseState_.getRestartConfig();
        const int reportStepNum = timer.reportStepNum();
        bool logMessages = output_ && parallelOutput_->isIORank();
        RestartValue::ExtraVector extraRestartData;
        std::map<std::string, double> miscSummaryData;

        if( output_ )
        {
            {
                // get all data that need to be included in output from the model
                // for flow_legacy and polymer this is a struct holding the data
                // while for flow_ebos a SimulationDataContainer is returned
                // this is addressed in the above specialized methods
                SimulationDataContainer sd =
                    detail::convertToSimulationDataContainer( physicalModel.getSimulatorData(localState), localState, phaseUsage_ );

                localCellData = simToSolution( sd, restart_double_si_, phaseUsage_); // Get "normal" data (SWAT, PRESSURE, ...);

                detail::getRestartData( localCellData, std::move(sd), phaseUsage_, physicalModel,
                                        restartConfig, reportStepNum, logMessages );
                // sd will be invalid after getRestartData has been called
            }
            detail::getSummaryData( localCellData, phaseUsage_, physicalModel, summaryConfig_ );
            assert(!localCellData.empty());

            // Add suggested next timestep to extra data.
            extraRestartData.push_back({{"OPMEXTRA", UnitSystem::measure::identity}, std::vector<double>(1, nextstep)});

            // Add TCPU if simulatorReport is not defaulted.
            const double totalSolverTime = simulatorReport.solver_time;
            if (totalSolverTime != 0.0) {
                miscSummaryData["TCPU"] = totalSolverTime;
            }
        }
        writeTimeStepWithCellProperties(timer, localState, localCellData, physicalModel.wellModel().wellState(localWellState), miscSummaryData, extraRestartData, substep);
    }
}
#endif
