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

#ifndef OPM_SIMULATORINCOMPTWOPHASE_HEADER_INCLUDED
#define OPM_SIMULATORINCOMPTWOPHASE_HEADER_INCLUDED

#include <boost/shared_ptr.hpp>
#include <vector>

struct UnstructuredGrid;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm
{
    namespace parameter { class ParameterGroup; }
    class IncompPropertiesInterface;
    class RockCompressibility;
    class LinearSolverInterface;
    class SimulatorTimer;
    class TwophaseState;
    class WellState;
    struct SimulatorReport;

    /// Class collecting all necessary components for a two-phase simulation.
    class SimulatorIncompTwophase
    {
    public:
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
        /// \param[in] grid        grid data structure
        /// \param[in] props       fluid and rock properties
        /// \param[in] rock_comp   if non-null, rock compressibility properties
        /// \param[in] wells       if non-null, wells data structure
        /// \param[in] src         source terms
        /// \param[in] bcs         boundary conditions, treat as all noflow if null
        /// \param[in] linsolver   linear solver
        /// \param[in] gravity     if non-null, gravity vector
       SimulatorIncompTwophase(const parameter::ParameterGroup& param,
                               const UnstructuredGrid& grid,
                               const IncompPropertiesInterface& props,
                               const RockCompressibility* rock_comp,
                               const Wells* wells,
                               const std::vector<double>& src,
                               const FlowBoundaryConditions* bcs,
                               LinearSolverInterface& linsolver,
                               const double* gravity);

        /// Run the simulation.
        /// This will run succesive timesteps until timer.done() is true. It will
        /// modify the reservoir and well states.
        /// \param[in,out] timer       governs the requested reporting timesteps
        /// \param[in,out] state       state of reservoir: pressure, fluxes
        /// \param[in,out] well_state  state of wells: bhp, perforation rates
        /// \return                    simulation report, with timing data
        SimulatorReport run(SimulatorTimer& timer,
                            TwophaseState& state,
                            WellState& well_state);

    private:
        class Impl;
        // Using shared_ptr instead of scoped_ptr since scoped_ptr requires complete type for Impl.
        boost::shared_ptr<Impl> pimpl_;
    };

} // namespace Opm

#endif // OPM_SIMULATORTWOPHASE_HEADER_INCLUDED
