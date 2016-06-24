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
#include <opm/core/utility/Compat.hpp>
#include <opm/core/utility/DataMap.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/output/eclipse/EclipseReader.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/output/Cells.hpp>
#include <opm/output/OutputWriter.hpp>
#include <opm/output/eclipse/EclipseWriter.hpp>

#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/ParallelDebugOutput.hpp>

#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/ThreadHandle.hpp>

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
                             Opm::EclipseStateConstPtr eclipseState,
                             const NNC&,
                             const Opm::PhaseUsage &phaseUsage,
                             const double* permeability );

        /** \copydoc Opm::OutputWriter::writeInit */
        void writeInit(const std::vector<data::CellData>& simProps, const NNC& nnc);

        /** \copydoc Opm::OutputWriter::writeTimeStep */
        void writeTimeStep(const SimulatorTimerInterface& timer,
                           const SimulationDataContainer& reservoirState,
                           const Opm::WellState& wellState,
                           bool substep = false);

        /** \copydoc Opm::OutputWriter::writeTimeStep */
        void writeTimeStepSerial(const SimulatorTimerInterface& timer,
                                 const SimulationDataContainer& reservoirState,
                                 const Opm::WellState& wellState,
                                 bool substep);

        /** \brief return output directory */
        const std::string& outputDirectory() const { return outputDir_; }

        /** \brief return true if output is enabled */
        bool output () const { return output_; }

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
        EclipseStateConstPtr eclipseState_;

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
                         Opm::EclipseStateConstPtr eclipseState,
                         const NNC& nnc,
                         const Opm::PhaseUsage &phaseUsage,
                         const double* permeability )
      : output_( param.getDefault("output", true) ),
        parallelOutput_( output_ ? new ParallelDebugOutput< Grid >( grid, eclipseState, phaseUsage.num_phases, permeability ) : 0 ),
        outputDir_( output_ ? param.getDefault("output_dir", std::string("output")) : "." ),
        output_interval_( output_ ? param.getDefault("output_interval", 1): 0 ),
        lastBackupReportStep_( -1 ),
        phaseUsage_( phaseUsage ),
        vtkWriter_( output_ && param.getDefault("output_vtk",false) ?
                     new BlackoilVTKWriter< Grid >( grid, outputDir_ ) : 0 ),
        matlabWriter_( output_ && parallelOutput_->isIORank() &&
                       param.getDefault("output_matlab", false) ?
                     new BlackoilMatlabWriter< Grid >( grid, outputDir_ ) : 0 ),
        eclWriter_( output_ && parallelOutput_->isIORank() &&
                    param.getDefault("output_ecl", true) ?
                    new EclipseWriter(eclipseState,
                                      parallelOutput_->numCells(),
                                      parallelOutput_->globalCell(),
                                      nnc )
                   : 0 ),
        eclipseState_(eclipseState),
        asyncOutput_()
    {
        // For output.
        if (output_ && parallelOutput_->isIORank() ) {
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
    BlackoilOutputWriter::
    initFromRestartFile( const PhaseUsage& phaseusage,
                         const double* permeability,
                         const Grid& grid,
                         SimulationDataContainer& simulatorstate,
                         WellStateFullyImplicitBlackoil& wellstate)
    {
        WellsManager wellsmanager(eclipseState_,
                                  eclipseState_->getInitConfig()->getRestartStep(),
                                  Opm::UgGridHelpers::numCells(grid),
                                  Opm::UgGridHelpers::globalCell(grid),
                                  Opm::UgGridHelpers::cartDims(grid),
                                  Opm::UgGridHelpers::dimensions(grid),
                                  Opm::UgGridHelpers::cell2Faces(grid),
                                  Opm::UgGridHelpers::beginFaceCentroids(grid),
                                  permeability);

        const Wells* wells = wellsmanager.c_wells();
        wellstate.resize(wells, simulatorstate); //Resize for restart step
        auto restarted = Opm::init_from_restart_file(
                                *eclipseState_,
                                Opm::UgGridHelpers::numCells(grid) );

        solutionToSim( restarted.first, phaseusage, simulatorstate );
        wellsToState( restarted.second, wellstate );
    }


}
#endif
