/*
  Copyright (c) 2014 SINTEF ICT, Applied Mathematics.
  Copyright (c) 2015-2016 IRIS AS

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
#include "config.h"

#include "SimulatorFullyImplicitBlackoilOutput.hpp"

#include <opm/common/data/SimulationDataContainer.hpp>

#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/utility/DataMap.hpp>
#include <opm/autodiff/Compat.hpp>
#include <opm/simulators/vtk/writeVtkData.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>

#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/BackupRestore.hpp>

#include <sstream>
#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>

//For OutputWriterHelper
#include <map>
#include <opm/parser/eclipse/Units/UnitSystem.hpp>


#ifdef HAVE_OPM_GRID
#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <dune/common/version.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>
#endif
namespace Opm
{



    void outputStateVtk(const UnstructuredGrid& grid,
                        const SimulationDataContainer& state,
                        const int step,
                        const std::string& output_dir)
    {
        // Write data in VTK format.
        std::ostringstream vtkfilename;
        vtkfilename << output_dir << "/vtk_files";
        boost::filesystem::path fpath(vtkfilename.str());
        try {
            create_directories(fpath);
        }
        catch (...) {
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        vtkfilename << "/output-" << std::setw(3) << std::setfill('0') << step << ".vtu";
        std::ofstream vtkfile(vtkfilename.str().c_str());
        if (!vtkfile) {
            OPM_THROW(std::runtime_error, "Failed to open " << vtkfilename.str());
        }
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        std::vector<double> cell_velocity;
        Opm::estimateCellVelocity(AutoDiffGrid::numCells(grid),
                                  AutoDiffGrid::numFaces(grid),
                                  AutoDiffGrid::beginFaceCentroids(grid),
                                  AutoDiffGrid::faceCells(grid),
                                  AutoDiffGrid::beginCellCentroids(grid),
                                  AutoDiffGrid::beginCellVolumes(grid),
                                  AutoDiffGrid::dimensions(grid),
                                  state.faceflux(), cell_velocity);
        dm["velocity"] = &cell_velocity;
        Opm::writeVtkData(grid, dm, vtkfile);
    }

    void outputWellStateMatlab(const Opm::WellState& well_state,
                               const int step,
                               const std::string& output_dir)
    {
        Opm::DataMap dm;
        dm["bhp"] = &well_state.bhp();
        dm["wellrates"] = &well_state.wellRates();

        // Write data (not grid) in Matlab format
        for (Opm::DataMap::const_iterator it = dm.begin(); it != dm.end(); ++it) {
            std::ostringstream fname;
            fname << output_dir << "/" << it->first;
            boost::filesystem::path fpath = fname.str();
            try {
                create_directories(fpath);
            }
            catch (...) {
                OPM_THROW(std::runtime_error,"Creating directories failed: " << fpath);
            }
            fname << "/" << std::setw(3) << std::setfill('0') << step << ".txt";
            std::ofstream file(fname.str().c_str());
            if (!file) {
                OPM_THROW(std::runtime_error,"Failed to open " << fname.str());
            }
            file.precision(15);
            const std::vector<double>& d = *(it->second);
            std::copy(d.begin(), d.end(), std::ostream_iterator<double>(file, "\n"));
        }
    }

#if 0
    void outputWaterCut(const Opm::Watercut& watercut,
                        const std::string& output_dir)
    {
        // Write water cut curve.
        std::string fname = output_dir  + "/watercut.txt";
        std::ofstream os(fname.c_str());
        if (!os) {
            OPM_THROW(std::runtime_error, "Failed to open " << fname);
        }
        watercut.write(os);
    }

    void outputWellReport(const Opm::WellReport& wellreport,
                          const std::string& output_dir)
    {
        // Write well report.
        std::string fname = output_dir  + "/wellreport.txt";
        std::ofstream os(fname.c_str());
        if (!os) {
            OPM_THROW(std::runtime_error, "Failed to open " << fname);
        }
        wellreport.write(os);
    }
#endif

#ifdef HAVE_OPM_GRID
    void outputStateVtk(const Dune::CpGrid& grid,
                        const Opm::SimulationDataContainer& state,
                        const int step,
                        const std::string& output_dir)
    {
        // Write data in VTK format.
        std::ostringstream vtkfilename;
        std::ostringstream vtkpath;
        vtkpath << output_dir << "/vtk_files";
        vtkpath << "/output-" << std::setw(3) << std::setfill('0') << step;
        boost::filesystem::path fpath(vtkpath.str());
        try {
            create_directories(fpath);
        }
        catch (...) {
            OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
        }
        vtkfilename << "output-" << std::setw(3) << std::setfill('0') << step;
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
        Dune::VTKWriter<Dune::CpGrid::LeafGridView> writer(grid.leafGridView(), Dune::VTK::nonconforming);
#else
        Dune::VTKWriter<Dune::CpGrid::LeafGridView> writer(grid.leafView(), Dune::VTK::nonconforming);
#endif
        writer.addCellData(state.saturation(), "saturation", state.numPhases());
        writer.addCellData(state.pressure(), "pressure", 1);

        std::vector<double> cell_velocity;
        Opm::estimateCellVelocity(AutoDiffGrid::numCells(grid),
                                  AutoDiffGrid::numFaces(grid),
                                  AutoDiffGrid::beginFaceCentroids(grid),
                                  AutoDiffGrid::faceCells(grid),
                                  AutoDiffGrid::beginCellCentroids(grid),
                                  AutoDiffGrid::beginCellVolumes(grid),
                                  AutoDiffGrid::dimensions(grid),
                                  state.faceflux(), cell_velocity);
        writer.addCellData(cell_velocity, "velocity", Dune::CpGrid::dimension);
        writer.pwrite(vtkfilename.str(), vtkpath.str(), std::string("."), Dune::VTK::ascii);
    }
#endif




    namespace detail {

        struct WriterCall : public ThreadHandle :: ObjectInterface
        {
            BlackoilOutputWriter& writer_;
            std::unique_ptr< SimulatorTimerInterface > timer_;
            const SimulationDataContainer state_;
            const WellStateFullyImplicitBlackoil wellState_;
            data::Solution simProps_;
            const bool substep_;

            explicit WriterCall( BlackoilOutputWriter& writer,
                                 const SimulatorTimerInterface& timer,
                                 const SimulationDataContainer& state,
                                 const WellStateFullyImplicitBlackoil& wellState,
                                 const data::Solution& simProps,
                                 bool substep )
                : writer_( writer ),
                  timer_( timer.clone() ),
                  state_( state ),
                  wellState_( wellState ),
                  simProps_( simProps ),
                  substep_( substep )
            {
            }

            // callback to writer's serial writeTimeStep method
            void run ()
            {
                // write data
                writer_.writeTimeStepSerial( *timer_, state_, wellState_, simProps_, substep_ );
            }
        };
    }




    void
    BlackoilOutputWriter::
    writeTimeStepWithoutCellProperties(
                  const SimulatorTimerInterface& timer,
                  const SimulationDataContainer& localState,
                  const WellStateFullyImplicitBlackoil& localWellState,
                  bool substep)
    {
        data::Solution localCellData{};
        if( output_ )
        {
            localCellData = simToSolution(localState, phaseUsage_); // Get "normal" data (SWAT, PRESSURE, ...);
        }
        writeTimeStepWithCellProperties(timer, localState, localCellData ,
                                        localWellState, substep);
    }





    void
    BlackoilOutputWriter::
    writeTimeStepWithCellProperties(
                  const SimulatorTimerInterface& timer,
                  const SimulationDataContainer& localState,
                  const data::Solution& localCellData,
                  const WellStateFullyImplicitBlackoil& localWellState,
                  bool substep)
    {
        // VTK output (is parallel if grid is parallel)
        if( vtkWriter_ ) {
            vtkWriter_->writeTimeStep( timer, localState, localWellState, false );
        }

        bool isIORank = output_ ;
        if( parallelOutput_ && parallelOutput_->isParallel() )
        {
            // If this is not the initial write and no substep, then the well
            // state used in the computation is actually the one of the last
            // step. We need that well state for the gathering. Otherwise
            // It an exception with a message like "global state does not
            // contain well ..." might be thrown.
            int wellStateStepNumber = ( ! substep && timer.reportStepNum() > 0) ?
                (timer.reportStepNum() - 1) : timer.reportStepNum();
            // collect all solutions to I/O rank
            isIORank = parallelOutput_->collectToIORank( localState, localWellState,
                                                         localCellData,
                                                         wellStateStepNumber );
        }

        const data::Solution& cellData = ( parallelOutput_ && parallelOutput_->isParallel() ) ? parallelOutput_->globalCellData() : localCellData;
        const SimulationDataContainer& state = (parallelOutput_ && parallelOutput_->isParallel() ) ? parallelOutput_->globalReservoirState() : localState;
        const WellStateFullyImplicitBlackoil& wellState  = (parallelOutput_ && parallelOutput_->isParallel() ) ? parallelOutput_->globalWellState() : localWellState;

        // serial output is only done on I/O rank
        if( isIORank )
        {
            if( asyncOutput_ ) {
                // dispatch the write call to the extra thread
                asyncOutput_->dispatch( detail::WriterCall( *this, timer, state, wellState, cellData, substep ) );
            }
            else {
                // just write the data to disk
                writeTimeStepSerial( timer, state, wellState, cellData, substep );
            }
        }
    }



    void
    BlackoilOutputWriter::
    writeTimeStepSerial(const SimulatorTimerInterface& timer,
                        const SimulationDataContainer& state,
                        const WellStateFullyImplicitBlackoil& wellState,
                        const data::Solution& simProps,
                        bool substep)
    {
        // Matlab output
        if( matlabWriter_ ) {
            matlabWriter_->writeTimeStep( timer, state, wellState, substep );
        }

        // ECL output
        if ( eclWriter_ )
        {
            const auto& initConfig = eclipseState_.getInitConfig();
            if (initConfig.restartRequested() && ((initConfig.getRestartStep()) == (timer.currentStepNum()))) {
                std::cout << "Skipping restart write in start of step " << timer.currentStepNum() << std::endl;
            } else {
                // ... insert "extra" data (KR, VISC, ...)
                eclWriter_->writeTimeStep(timer.reportStepNum(),
                                          substep,
                                          timer.simulationTimeElapsed(),
                                          simProps,
                                          wellState.report(phaseUsage_));
            }
        }

        // write backup file
        if( backupfile_.is_open() )
        {
            int reportStep      = timer.reportStepNum();
            int currentTimeStep = timer.currentStepNum();
            if( (reportStep == currentTimeStep || // true for SimulatorTimer
                 currentTimeStep == 0 || // true for AdaptiveSimulatorTimer at reportStep
                 timer.done() ) // true for AdaptiveSimulatorTimer at reportStep
               && lastBackupReportStep_ != reportStep ) // only backup report step once
            {
                // store report step
                lastBackupReportStep_ = reportStep;
                // write resport step number
                backupfile_.write( (const char *) &reportStep, sizeof(int) );

                try {
                    backupfile_ << state;

                    const WellStateFullyImplicitBlackoil& boWellState = static_cast< const WellStateFullyImplicitBlackoil& > (wellState);
                    backupfile_ << boWellState;
                }
                catch ( const std::bad_cast& e )
                {
                }

                backupfile_ << std::flush;
            }
        } // end backup
    }

    void
    BlackoilOutputWriter::
    restore(SimulatorTimerInterface& timer,
            BlackoilState& state,
            WellStateFullyImplicitBlackoil& wellState,
            const std::string& filename,
            const int desiredResportStep )
    {
        std::ifstream restorefile( filename.c_str() );
        if( restorefile )
        {
            std::cout << "============================================================================"<<std::endl;
            std::cout << "Restoring from ";
            if( desiredResportStep < 0 ) {
                std::cout << "last";
            }
            else {
                std::cout << desiredResportStep;
            }
            std::cout << " report step! filename = " << filename << std::endl << std::endl;

            int reportStep;
            restorefile.read( (char *) &reportStep, sizeof(int) );

            const int readReportStep = (desiredResportStep < 0) ?
                std::numeric_limits<int>::max() : desiredResportStep;

            while( reportStep <= readReportStep && ! timer.done() && restorefile )
            {
                restorefile >> state;
                restorefile >> wellState;

                // No per cell data is written for restore steps, but will be
                // for subsequent steps, when we have started simulating
                writeTimeStepWithoutCellProperties( timer, state, wellState );

                // some output
                std::cout << "Restored step " << timer.reportStepNum() << " at day "
                          <<  unit::convert::to(timer.simulationTimeElapsed(),unit::day) << std::endl;

                if( readReportStep == reportStep ) {
                    break;
                }

                // if the stream is not valid anymore we just use the last state read
                if( ! restorefile ) {
                    std::cerr << "Reached EOF, using last state read!" << std::endl;
                    break;
                }

                // try to read next report step
                restorefile.read( (char *) &reportStep, sizeof(int) );

                // if read failed, exit loop
                if( ! restorefile ) {
                    break;
                }

                // next step
                timer.advance();

                if( timer.reportStepNum() != reportStep ) {
                    break;
                }
            }
        }
        else
        {
            std::cerr << "Warning: Couldn't open restore file '" << filename << "'" << std::endl;
        }
    }


    bool BlackoilOutputWriter::isRestart() const {
        const auto& initconfig = eclipseState_.getInitConfig();
        return initconfig.restartRequested();
    }
}
