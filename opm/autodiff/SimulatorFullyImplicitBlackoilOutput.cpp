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
#include "config.h"

#include "SimulatorFullyImplicitBlackoilOutput.hpp"

#include <opm/core/utility/DataMap.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>

#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/BackupRestore.hpp>

#include <sstream>
#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>

#ifdef HAVE_DUNE_CORNERPOINT
#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <dune/common/version.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>
#endif
namespace Opm
{


    void outputStateVtk(const UnstructuredGrid& grid,
                        const SimulatorState& state,
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


    void outputStateMatlab(const UnstructuredGrid& grid,
                           const Opm::SimulatorState& state,
                           const int step,
                           const std::string& output_dir)
    {
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        for( unsigned int i=0; i<state.cellDataNames().size(); ++ i )
        {
            const std::string& name = state.cellDataNames()[ i ];
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
            dm[ key ] = &state.cellData()[ i ];
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

#ifdef HAVE_DUNE_CORNERPOINT
    void outputStateVtk(const Dune::CpGrid& grid,
                        const Opm::SimulatorState& state,
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

    void
    BlackoilOutputWriter::
    writeInit(const SimulatorTimerInterface& timer)
    {
        if( eclWriter_ ) {
            eclWriter_->writeInit(timer);
        }
    }

    void
    BlackoilOutputWriter::
    writeTimeStep(const SimulatorTimerInterface& timer,
                  const SimulatorState& localState,
                  const WellState& localWellState,
                  bool substep)
    {
        // VTK output (is parallel if grid is parallel)
        if( vtkWriter_ ) {
            vtkWriter_->writeTimeStep( timer, localState, localWellState, false );
        }

        bool isIORank = true ;
        if( parallelOutput_ && parallelOutput_->isParallel() )
        {
            // collect all solutions to I/O rank
            isIORank = parallelOutput_->collectToIORank( localState, localWellState );
        }

        const SimulatorState& state = (parallelOutput_ && parallelOutput_->isParallel() ) ? parallelOutput_->globalReservoirState() : localState;
        const WellState& wellState  = (parallelOutput_ && parallelOutput_->isParallel() ) ? parallelOutput_->globalWellState() : localWellState;

        // output is only done on I/O rank
        if( isIORank )
        {
            // Matlab output
            if( matlabWriter_ ) {
                matlabWriter_->writeTimeStep( timer, state, wellState, substep );
            }
            // ECL output
            if ( eclWriter_ ) {
                eclWriter_->writeTimeStep(timer, state, wellState, substep );
            }

            // write backup file
            if( backupfile_ )
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
                        const BlackoilState& boState = dynamic_cast< const BlackoilState& > (state);
                        backupfile_ << boState;

                        const WellStateFullyImplicitBlackoil& boWellState = static_cast< const WellStateFullyImplicitBlackoil& > (wellState);
                        backupfile_ << boWellState;
                    }
                    catch ( const std::bad_cast& e )
                    {

                    }

                    /*
                    const WellStateFullyImplicitBlackoil* boWellState =
                        dynamic_cast< const WellStateFullyImplicitBlackoil* > (&wellState);
                    if( boWellState ) {
                        backupfile_ << (*boWellState);
                    }
                    else
                        OPM_THROW(std::logic_error,"cast to WellStateFullyImplicitBlackoil failed");
                    */
                    backupfile_ << std::flush;
                }
            } // end backup
        } // end isIORank
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

                writeTimeStep( timer, state, wellState );
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
}
