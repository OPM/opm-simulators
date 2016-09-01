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
#include <opm/output/Cells.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/utility/DataMap.hpp>
#include <opm/core/utility/Compat.hpp>
#include <opm/output/vtk/writeVtkData.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>

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


    void outputStateMatlab(const UnstructuredGrid& grid,
                           const Opm::SimulationDataContainer& state,
                           const int step,
                           const std::string& output_dir)
    {
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        for (const auto& pair : state.cellData()) {
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





    void
    BlackoilOutputWriter::
    writeTimeStepSerial(const SimulatorTimerInterface& timer,
                        const SimulationDataContainer& state,
                        const WellState& wellState,
                        const std::vector<data::CellData>& simProps,
                        bool substep)
    {
        // Matlab output
        if( matlabWriter_ ) {
            matlabWriter_->writeTimeStep( timer, state, wellState, substep );
        }

        // ECL output
        if ( eclWriter_ )
        {
            const auto& initConfig = eclipseState_->getInitConfig();
            if (initConfig.restartRequested() && ((initConfig.getRestartStep()) == (timer.currentStepNum()))) {
                std::cout << "Skipping restart write in start of step " << timer.currentStepNum() << std::endl;
            } else {
                /*
                  The simProps vector can be passed to the writeTimestep routine
                  to add more properties to the restart file. Examples of the
                  elements for the simProps vector can be the relative
                  permeabilites KRO, KRG and KRW and the fluxes.

                  Which properties are requested are configured with the RPTRST
                  keyword, which is internalized in the RestartConfig class in
                  EclipseState.
                */

                /*
                  Assuming we already have correctly initialized
                  std::vector<double> instances kro,krw and krg with the oil,
                  water and gas relative permeabilities. Then we can write those
                  to the restart file with:

                     std::vector<data::CellData> simProps;

                     simProps.emplace_back( {"KRO" , UnitSystem::measure::identity , kro} );
                     simProps.emplace_back( {"KRG" , UnitSystem::measure::identity , krg} );
                     simProps.emplace_back( {"KRW" , UnitSystem::measure::identity , krw} );

                */

                eclWriter_->writeTimeStep(timer.reportStepNum(),
                                          substep,
                                          timer.simulationTimeElapsed(),
                                          simToSolution( state, phaseUsage_ ),
                                          wellState.report(),
                                          simProps);
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

                // FIXME: We this should optimally have the proper per cell data to dump
                // Right now it will not dump any per cell data until we start simulating
                std::vector<data::CellData> noData;
                writeTimeStep( timer, state, wellState, noData );

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
        const auto& initconfig = eclipseState_->getInitConfig();
        return initconfig.restartRequested();
    }
}
