/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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


#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cassert>
#include <opm/core/grid.h>
#include <opm/core/GridManager.hpp>
#include <opm/core/utility/writeVtkData.hpp>
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
#include <opm/core/pressure/IncompTpfa.hpp>
#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/fluid/IncompPropertiesBasic.hpp>

#include <opm/core/transport/reorder/TransportModelTwophase.hpp>

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>

#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/wells/WellCollection.hpp>

/// \page tutorial4 Well controls

/// \page tutorial4
/// \details
/// Main function
/// \code
int main ()
{
    /// \endcode
    /// \page tutorial4
    /// \details
    /// We define the grid. A cartesian grid with 1200 cells.
    /// \code
    int dim = 3;
    int nx = 20;
    int ny = 20;
    int nz = 1;
    double dx = 10.;
    double dy = 10.;
    double dz = 10.;
    using namespace Opm;
    GridManager grid_manager(nx, ny, nz, dx, dy, dz);
    const UnstructuredGrid& grid = *grid_manager.c_grid();
    int num_cells = grid.number_of_cells;
    /// \endcode

    /// \page tutorial4
    /// \details
    /// We define the properties of the fluid.\n 
    /// Number of phases.
    /// \code
    int num_phases = 2;
    using namespace unit;
    using namespace prefix;
    /// \endcode
    /// \page tutorial4
    /// \details density vector (one component per phase).
    /// \code
    std::vector<double> rho(2, 1000.);
    /// \endcode
    /// \page tutorial4
    /// \details viscosity vector (one component per phase).
    /// \code
    std::vector<double> mu(2, 1.*centi*Poise);
    /// \endcode
    /// \page tutorial4
    /// \details porosity and permeability of the rock.
    /// \code
    double porosity = 0.5;
    double k = 10*milli*darcy;
    /// \endcode

    /// \page tutorial4
    /// \details We define the relative permeability function. We use a basic fluid
    /// description and set this function to be linear. For more realistic fluid, the 
    /// saturation function is given by the data.
    /// \code
    SaturationPropsBasic::RelPermFunc rel_perm_func = SaturationPropsBasic::Linear;
    /// \endcode

    /// \page tutorial4
    /// \details We construct a basic fluid with the properties we have defined above.
    /// Each property is constant and hold for all cells.
    /// \code
    IncompPropertiesBasic props(num_phases, rel_perm_func, rho, mu,
                                     porosity, k, dim, num_cells);
    /// \endcode

    
    /// \page tutorial4
    /// \details Gravity parameters. Here, we set zero gravity.
    /// \code
    const double *grav = 0;
    std::vector<double> omega; 
    /// \endcode

    /// \page tutorial4
    /// \details We set up the source term. Positive numbers indicate that the cell is a source,
    /// while negative numbers indicate a sink.
    /// \code
    std::vector<double> src(num_cells, 0.0);
    src[0] = 1.;
    src[num_cells-1] = -1.;
    /// \endcode

    /// \page tutorial4
    /// \details We compute the pore volume
    /// \code
    std::vector<double> porevol;
    Opm::computePorevolume(grid, props.porosity(), porevol);
    /// \endcode

    /// \page tutorial4
    /// \details Set up the transport solver. This is a reordering implicit Euler transport solver.
    /// \code
    const double tolerance = 1e-9;
    const int max_iterations = 30;
    Opm::TransportModelTwophase transport_solver(grid, props, tolerance, max_iterations);
    /// \endcode

    /// \page tutorial4
    /// \details Time integration parameters
    /// \code
    double dt = 0.1*day;
    int num_time_steps = 20;
    /// \endcode


    /// \page tutorial4
    /// \details We define a vector which contains all cell indexes. We use this 
    /// vector to set up parameters on the whole domains.
    std::vector<int> allcells(num_cells);
    for (int cell = 0; cell < num_cells; ++cell) {
        allcells[cell] = cell;
    }


    /// \page tutorial4
    /// \details We set up the boundary conditions. Letting bcs empty is equivalent 
    /// to no flow boundary conditions.
    /// \code
    FlowBCManager bcs;
    /// \endcode

    /// \page tutorial4
    /// \details 
    /// We set up a two-phase state object, and
    /// initialise water saturation to minimum everywhere.
    /// \code
    TwophaseState state;
    state.init(grid, 2);
    state.setFirstSat(allcells, props, TwophaseState::MinSat);
    /// \endcode
    
    /// \page tutorial4
    /// \details This string will contain the name of a VTK output vector.
    /// \code
    std::ostringstream vtkfilename;
    /// \endcode


    /// \page tutorial4
    /// To create wells we need an instance of the PhaseUsage-object
    /// \code
    PhaseUsage phase_usage;
    phase_usage.num_phases = num_phases;
    phase_usage.phase_used[BlackoilPhases::Aqua] = 1;
    phase_usage.phase_used[BlackoilPhases::Liquid] = 1;
    phase_usage.phase_used[BlackoilPhases::Vapour] = 0;
    
    phase_usage.phase_pos[BlackoilPhases::Aqua] = 0;
    phase_usage.phase_pos[BlackoilPhases::Liquid] = 1;
    /// \endcode

    /// \page tutorial4
    /// \details This will contain our well-specific information
    /// \code
    WellCollection well_collection;
    /// \endcode
    
    /// \page tutorial4
    /// \details Create the production specification for our top well group.
    ///          We set a target limit for total reservoir rate, and set the controlling
    ///          mode of the group to be controlled by the reservoir rate.
    /// \code 
    ProductionSpecification well_group_prod_spec;
    well_group_prod_spec.reservoir_flow_max_rate_ = 0.1;
    well_group_prod_spec.control_mode_ = ProductionSpecification::RESV;
    /// \endcode
    
    /// \page tutorial4
    /// \details Create our well group. We hand it an empty injection specification, 
    ///          as we don't want to control its injection target. We use the shared_ptr-type because that's 
    ///          what the interface expects. The first argument is the (unique) name of the group.
    /// \code 
    std::tr1::shared_ptr<WellsGroupInterface> well_group(new WellsGroup("group", well_group_prod_spec, InjectionSpecification(), 
                                                                        phase_usage));
    /// \endcode
    
    /// \page tutorial4
    /// \details We add our well_group to the well_collection
    /// \code
    well_collection.addChild(well_group);
    /// \endcode
    
    
    /// \page tutorial4
    /// \details Create the production specification and Well objects (total 4 wells). We set all our wells to be group controlled. 
    ///          We pass in the string argument \C "group" to set the parent group.
    /// \code
    const int num_wells = 4;
    for (int i = 0; i < num_wells; ++i) {
        std::stringstream well_name;
        well_name << "well" << i;
        ProductionSpecification production_specification;
        production_specification.control_mode_ = ProductionSpecification::GRUP;
        std::tr1::shared_ptr<WellsGroupInterface> well_leaf_node(new WellNode(well_name.str(), production_specification, InjectionSpecification(), 
                                                                              phase_usage));
        well_collection.addChild(well_leaf_node, "group");
        
    }
    /// \endcode
    
    /// \page tutorial4
    /// \details Now we create the C struct to hold our wells (this is to interface with the solver code). For now we
    /// \code
    Wells* wells = create_wells(num_phases, num_wells, num_wells /*number of perforations. We'll only have one perforation per well*/);
    /// \endcode
    
    /// 
    
    /// \page tutorial4
    /// \details We need to add each well to the C API. 
    /// To do this we need to specify the relevant cells the well will be located in (\C well_cells). 
    /// \code
    for (int i = 0; i < num_wells; ++i) {
        const int well_cells = i*nx;
        const double well_index = 1;
        std::stringstream well_name;
        well_name << "well" << i;
        add_well(PRODUCER, 0, 1, NULL, &well_cells, &well_index,
                 well_name.str().c_str(), wells);
    }
    /// \endcode
    
    /// \page tutorial4
    /// \details We need to make the well collection aware of our wells object
    /// \code
    well_collection.setWellsPointer(wells);
    /// \endcode
    /// \page tutorial4
    /// We're not using well controls, just group controls, so we need to apply them.
    /// \code   
    well_collection.applyGroupControls();
    ///\endcode
    
    /// \page tutorial4
    /// \details We set up necessary information for the wells
    /// \code
    WellState well_state;
    well_state.init(wells, state);
    std::vector<double> well_resflowrates_phase;
    std::vector<double> well_surflowrates_phase;
    std::vector<double> fractional_flows;
    /// \endcode

    /// \page tutorial4
    /// \details We set up the pressure solver.
    /// \code
    LinearSolverUmfpack linsolver;
    IncompTpfa psolver(grid, props, linsolver,
                       grav, wells, src, bcs.c_bcs());
    /// \endcode

    
    /// \page tutorial4
    /// \details Loop over the time steps.
    /// \code
    for (int i = 0; i < num_time_steps; ++i) {
        /// \endcode

        /// \page tutorial4
        /// \details We're solving the pressure until the well conditions are met 
        ///          or until we reach the maximum number of iterations.
        /// \code
        const int max_well_iterations = 10;
        int well_iter = 0;
        bool well_conditions_met = false;
        while (!well_conditions_met) {
            /// \endcode

            /// \page tutorial4
            /// \details Solve the pressure equation
            /// \code
            psolver.solve(dt, state, well_state);

            /// \endcode
            /// \page tutorial4
            /// \details We compute the new well rates. Notice that we approximate (wrongly) surfflowsrates := resflowsrate 
            Opm::computeFractionalFlow(props, allcells, state.saturation(), fractional_flows);
            Opm::computePhaseFlowRatesPerWell(*wells, well_state.perfRates(), fractional_flows, well_resflowrates_phase);
            Opm::computePhaseFlowRatesPerWell(*wells, well_state.perfRates(), fractional_flows, well_surflowrates_phase);
            /// \endcode

            /// \page tutorial4
            /// \details We check if the well conditions are met.
            well_conditions_met = well_collection.conditionsMet(well_state.bhp(), well_resflowrates_phase, well_surflowrates_phase);
            ++well_iter;
            if (!well_conditions_met && well_iter == max_well_iterations) {
                THROW("Conditions not met within " << max_well_iterations<< " iterations.");
            }
        }
        /// \endcode

        /// \page tutorial4
        /// \details  Transport solver
        /// \TODO We must call computeTransportSource() here, since we have wells.
        /// \code
        transport_solver.solve(&state.faceflux()[0], &porevol[0], &src[0], dt, state.saturation());
        /// \endcode

        /// \page tutorial4
        /// \details Write the output to file.
        /// \code
        vtkfilename.str(""); 
        vtkfilename << "tutorial4-" << std::setw(3) << std::setfill('0') << i << ".vtu";
        std::ofstream vtkfile(vtkfilename.str().c_str());
        Opm::DataMap dm;
        dm["saturation"] = &state.saturation();
        dm["pressure"] = &state.pressure();
        Opm::writeVtkData(grid, dm, vtkfile);
    }
}
/// \endcode



/// \page tutorial4
/// \section results4 Results.
/// <TABLE>
/// <TR>
/// <TD> \image html tutorial4-000.png </TD>
/// <TD> \image html tutorial4-005.png </TD>
/// </TR>
/// <TR>
/// <TD> \image html tutorial4-010.png </TD>
/// <TD> \image html tutorial4-015.png </TD>
/// </TR>
/// <TR>
/// <TD> \image html tutorial4-019.png </TD>
/// <TD> </TD>
/// </TR>
/// </TABLE>


/// \page tutorial4
/// \details
/// \section completecode4 Complete source code:
/// \include tutorial4.cpp 
/// \include generate_doc_figures.py 
