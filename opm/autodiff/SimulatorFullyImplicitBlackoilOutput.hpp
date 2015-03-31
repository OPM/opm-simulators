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
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/DataMap.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/io/OutputWriter.hpp>
#include <opm/core/io/eclipse/EclipseWriter.hpp>

#include <opm/autodiff/GridHelpers.hpp>

#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>

#ifdef HAVE_DUNE_CORNERPOINT
#include <dune/grid/CpGrid.hpp>
#endif
namespace Opm
{

    void outputStateVtk(const UnstructuredGrid& grid,
                        const Opm::SimulatorState& state,
                        const int step,
                        const std::string& output_dir);


    void outputStateMatlab(const UnstructuredGrid& grid,
                           const Opm::BlackoilState& state,
                           const int step,
                           const std::string& output_dir);

    void outputWellStateMatlab(const Opm::WellState& well_state,
                               const int step,
                               const std::string& output_dir);
#ifdef HAVE_DUNE_CORNERPOINT
    void outputStateVtk(const Dune::CpGrid& grid,
                        const Opm::SimulatorState& state,
                        const int step,
                        const std::string& output_dir);
#endif

    template<class Grid>
    void outputStateMatlab(const Grid& grid,
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

    /** Wrapper class satisfying the OutputWriter interface and writing VTK output. */
    template <class Grid>
    class BlackoilVTKWriter : public OutputWriter
    {
    public:
        // constructor taking grid and directory
        BlackoilVTKWriter(const Grid& grid,
                          const std::string& outputDir)
          : grid_( grid ),
            outputDir_( outputDir )
        {}

        /** \copydoc Opm::OutputWriter::writeInit */
        void writeInit(const SimulatorTimerInterface& /* timer */)
        {}

        /** \copydoc Opm::OutputWriter::writeTimeStep */
        void writeTimeStep(const SimulatorTimerInterface& timer,
                           const SimulatorState& state,
                           const WellState& )
        {
            outputStateVtk(grid_, state, timer.currentStepNum(), outputDir_);
        }

    protected:
        const Grid& grid_;
        const std::string outputDir_;
    };

    /** Wrapper class satisfying the OutputWriter interface and writing Matlab output. */
    template <class Grid>
    class BlackoilMatlabWriter : public OutputWriter
    {
    public:
        // constructor taking grid and directory
        BlackoilMatlabWriter(const Grid& grid,
                             const std::string& outputDir)
          : grid_( grid ),
            outputDir_( outputDir )
        {}

        /** \copydoc Opm::OutputWriter::writeInit */
        void writeInit(const SimulatorTimerInterface& /* timer */)
        {}

        /** \copydoc Opm::OutputWriter::writeTimeStep */
        void writeTimeStep(const SimulatorTimerInterface& timer,
                           const SimulatorState& reservoirState,
                           const WellState& wellState)
        {
            const BlackoilState* state =
                dynamic_cast< const BlackoilState* > (&reservoirState);
            if( state ) {
                outputStateMatlab(grid_, *state, timer.currentStepNum(), outputDir_);
            }
            else {
                OPM_THROW(std::logic_error,"BlackoilMatlabWriter only works for BlackoilState");
            }
            outputWellStateMatlab(wellState, timer.currentStepNum(), outputDir_);
        }
    protected:
        const Grid& grid_;
        const std::string outputDir_;
    };

    /** \brief Wrapper class for VTK, Matlab, and ECL output. */
    class BlackoilOutputWriter : public OutputWriter
    {
    public:
        // constructor creating different sub writers
        template <class Grid>
        BlackoilOutputWriter(const Grid& grid,
                             const parameter::ParameterGroup& param,
                             Opm::EclipseStateConstPtr eclipseState,
                             const Opm::PhaseUsage &phaseUsage);

        /** \copydoc Opm::OutputWriter::writeInit */
        void writeInit(const SimulatorTimerInterface &timer);

        /** \copydoc Opm::OutputWriter::writeTimeStep */
        void writeTimeStep(const SimulatorTimerInterface& timer,
                           const SimulatorState& reservoirState,
                           const WellState& wellState);

        /** \brief return output directory */
        const std::string& outputDirectory() const { return outputDir_; }

        /** \brief return true if output is enabled */
        bool output () const { return output_; }

        void restore(SimulatorTimerInterface& timer,
                     BlackoilState& state,
                     WellStateFullyImplicitBlackoil& wellState,
                     const std::string& filename,
                     const int desiredReportStep);

    protected:
        // Parameters for output.
        const bool output_;
        const std::string outputDir_;
        const int output_interval_;

        int lastBackupReportStep_;

        std::ofstream backupfile_;
        std::unique_ptr< OutputWriter  > vtkWriter_;
        std::unique_ptr< OutputWriter  > matlabWriter_;
        std::unique_ptr< EclipseWriter > eclWriter_;
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
                         const Opm::PhaseUsage &phaseUsage )
      : output_( param.getDefault("output", true) ),
        outputDir_( output_ ? param.getDefault("output_dir", std::string("output")) : "." ),
        output_interval_( output_ ? param.getDefault("output_interval", 1): 0 ),
        lastBackupReportStep_( -1 ),
        vtkWriter_( output_ && param.getDefault("output_vtk",false) ?
                     new BlackoilVTKWriter< Grid >( grid, outputDir_ ) : 0 ),
        matlabWriter_( output_ && param.getDefault("output_matlab", false) ?
                     new BlackoilMatlabWriter< Grid >( grid, outputDir_ ) : 0 ),
        eclWriter_( output_ && param.getDefault("output_ecl", true) ?
                    new EclipseWriter(param, eclipseState, phaseUsage,
                                      Opm::UgGridHelpers::numCells( grid ),
                                      Opm::UgGridHelpers::globalCell( grid ) )
                   : 0 )
    {
        // For output.
        if (output_) {
            // Ensure that output dir exists
            boost::filesystem::path fpath(outputDir_);
            try {
                create_directories(fpath);
            }
            catch (...) {
                OPM_THROW(std::runtime_error, "Creating directories failed: " << fpath);
            }

            std::string backupfilename = param.getDefault("backupfile", std::string("") );
            if( ! backupfilename.empty() )
            {
                backupfile_.open( backupfilename.c_str() );
            }
        }
    }
}
#endif
