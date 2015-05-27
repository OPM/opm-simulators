/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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

#ifndef OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_HEADER_INCLUDED

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/SimulatorBase.hpp>

#include <opm/polymer/fullyimplicit/FullyImplicitCompressiblePolymerSolver.hpp>
#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/io/eclipse/EclipseWriter.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>

#include <opm/core/wells/WellsManager.hpp>

#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/grid/ColumnExtract.hpp>
#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/polymer/PolymerInflow.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/WellProductionProperties.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <boost/filesystem.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <numeric>
#include <fstream>
#include <iostream>

namespace Opm
{
    /// Class collecting all necessary components for a two-phase simulation.
    template <class GridT>
    class SimulatorFullyImplicitCompressiblePolymer
        : public SimulatorBase<GridT, SimulatorFullyImplicitCompressiblePolymer<GridT> >
    {
        typedef SimulatorFullyImplicitCompressiblePolymer<GridT> ThisType;
        typedef SimulatorBase<GridT, ThisType> BaseType;

    public:
        /// Initialise from parameters and objects to observe.
        SimulatorFullyImplicitCompressiblePolymer(const parameter::ParameterGroup& param,
                        		                  const UnstructuredGrid& grid,
                                                  const DerivedGeology& geo,
                                   				  BlackoilPropsAdInterface& props,
                                       			  const PolymerPropsAd&    polymer_props,
                                       			  const RockCompressibility* rock_comp_props,
                                                  std::shared_ptr<EclipseState> eclipse_state,
                                                  BlackoilOutputWriter& output_writer,
                                                  Opm::DeckConstPtr& deck,
                                       			  NewtonIterationBlackoilInterface& linsolver,
                                       			  const double* gravity);

        /// Run the simulation.
        /// This will run succesive timesteps until timer.done() is true. It will
        /// modify the reservoir and well states.
        /// \param[in,out] timer       governs the requested reporting timesteps
        /// \param[in,out] state       state of reservoir: pressure, fluxes
        /// \return                    simulation report, with timing data
        SimulatorReport run(SimulatorTimer& timer,
                            PolymerBlackoilState& state);

private:
        Opm::DeckConstPtr deck_;
        const PolymerPropsAd& polymer_props_;

#warning "To remove"
        bool output_;
        int output_interval_;
        bool output_vtk_;
        std::string output_dir_;
        bool check_well_controls_;
        int max_well_control_iterations_;
    };

} // namespace Opm

#include "SimulatorFullyImplicitCompressiblePolymer_impl.hpp"

#endif // OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_HEADER_INCLUDED
