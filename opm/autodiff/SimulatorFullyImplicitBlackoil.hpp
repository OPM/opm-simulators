/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_SIMULATORFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITBLACKOIL_HEADER_INCLUDED

#include <memory>
#include <vector>

struct UnstructuredGrid;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm
{
    namespace parameter { class ParameterGroup; }
    class BlackoilPropsAdInterface;
    class RockCompressibility;
    class DerivedGeology;
    class NewtonIterationBlackoilInterface;
    class SimulatorTimer;
    class BlackoilState;
    class WellStateFullyImplicitBlackoil;
    class EclipseState;
    class EclipseWriter;
    struct SimulatorReport;

    /// Class collecting all necessary components for a two-phase simulation.
    template<class T>
    class SimulatorFullyImplicitBlackoil
    {
    public:
        /// \brief The type of the grid that we use.
        typedef T Grid;
        /// Initialise from parameters and objects to observe.
        /// \param[in] param       parameters, this class accepts the following:
        ///     parameter (default)            effect
        ///     -----------------------------------------------------------
        ///     output (true)                  write output to files?
        ///     output_dir ("output")          output directoty
        ///     output_interval (1)            output every nth step
        ///     nl_pressure_residual_tolerance (0.0) pressure solver residual tolerance (in Pascal)
        ///     nl_pressure_change_tolerance (1.0)   pressure solver change tolerance (in Pascal)
        ///     nl_pressure_maxiter (10)       max nonlinear iterations in pressure
        ///     nl_maxiter (30)                max nonlinear iterations in transport
        ///     nl_tolerance (1e-9)            transport solver absolute residual tolerance
        ///     num_transport_substeps (1)     number of transport steps per pressure step
        ///     use_segregation_split (false)  solve for gravity segregation (if false,
        ///                                    segregation is ignored).
        ///
        /// \param[in] grid          grid data structure
        /// \param[in] geo           derived geological properties
        /// \param[in] props         fluid and rock properties
        /// \param[in] rock_comp_props     if non-null, rock compressibility properties
        /// \param[in] linsolver     linear solver
        /// \param[in] gravity       if non-null, gravity vector
        /// \param[in] disgas        true for dissolved gas option
        /// \param[in] vapoil        true for vaporized oil option
        /// \param[in] eclipse_state
        /// \param[in] output_writer
        /// \param[in] threshold_pressures_by_face   if nonempty, threshold pressures that inhibit flow
        SimulatorFullyImplicitBlackoil(const parameter::ParameterGroup& param,
                                       const Grid& grid,
                                       const DerivedGeology& geo,
                                       BlackoilPropsAdInterface& props,
                                       const RockCompressibility* rock_comp_props,
                                       NewtonIterationBlackoilInterface& linsolver,
                                       const double* gravity,
                                       const bool disgas,
                                       const bool vapoil,
                                       std::shared_ptr<EclipseState> eclipse_state,
                                       EclipseWriter& output_writer,
                                       const std::vector<double>& threshold_pressures_by_face);

        /// Run the simulation.
        /// This will run succesive timesteps until timer.done() is true. It will
        /// modify the reservoir and well states.
        /// \param[in,out] timer       governs the requested reporting timesteps
        /// \param[in,out] state       state of reservoir: pressure, fluxes
        /// \param[in,out] well_state  state of wells: bhp, perforation rates
        /// \return                    simulation report, with timing data
        SimulatorReport run(SimulatorTimer& timer,
                            BlackoilState& state);

    private:
        class Impl;
        // Using shared_ptr instead of scoped_ptr since scoped_ptr requires complete type for Impl.
        std::shared_ptr<Impl> pimpl_;
    };

} // namespace Opm

#include "SimulatorFullyImplicitBlackoil_impl.hpp"
#endif // OPM_SIMULATORFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
