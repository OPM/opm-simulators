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
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/miscUtilities.hpp>

#include <opm/autodiff/GridHelpers.hpp>

#include <sstream>
#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>

#ifdef HAVE_DUNE_CORNERPOINT
#include <opm/core/utility/platform_dependent/disable_warnings.h>
#include <dune/common/version.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/core/utility/platform_dependent/reenable_warnings.h>
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
                           const Opm::BlackoilState& state,
                           const int step,
                           const std::string& output_dir)
    {
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        dm["surfvolume"] = &state.surfacevol();
        dm["rs"] = &state.gasoilratio();
        dm["rv"] = &state.rv();
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
                  const SimulatorState& state,
                  const WellState& wellState)
    {
        // VTK output
        if( vtkWriter_ ) {
            vtkWriter_->writeTimeStep( timer, state, wellState );
        }
        // Matlab output
        if( matlabWriter_ ) {
            matlabWriter_->writeTimeStep( timer, state, wellState );
        }
        // ECL output
        if ( eclWriter_ ) {
            eclWriter_->writeTimeStep(timer, state, wellState);
        }
    }
}
