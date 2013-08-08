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

#include <memory>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <vector>

struct UnstructuredGrid;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm
{
    namespace parameter { class ParameterGroup; }
    class IncompPropertiesInterface;
    class RockCompressibility;
    class WellsManager;
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
        /// \param[in] grid          grid data structure
        /// \param[in] props         fluid and rock properties
        /// \param[in] rock_comp_props     if non-null, rock compressibility properties
        /// \param[in] well_manager  well manager, may manage no (null) wells
        /// \param[in] src           source terms
        /// \param[in] bcs           boundary conditions, treat as all noflow if null
        /// \param[in] linsolver     linear solver
        /// \param[in] gravity       if non-null, gravity vector
       SimulatorIncompTwophase(const parameter::ParameterGroup& param,
                               const UnstructuredGrid& grid,
                               const IncompPropertiesInterface& props,
                               const RockCompressibility* rock_comp_props,
                               WellsManager& wells_manager,
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

        /// Register a callback for notifications about completed timestep.
        /// The specified method of the object is called when the simulator
        /// has completed a timestep and the state objects are (possibly)
        /// changed.
        /// If you want to know the current timestep, the callback must
        /// also monitor the timer object which was passed to run().
        /// \tparam T        Type of the callback object.
        /// \tparam callback Address of a member function of T which will
        ///                  be invoked when a timestep completes.
        /// \param[in] t     Object which will receive notifications. The
        ///                  lifetime of the object must span over the
        ///                  duration of the simulation, and it is your
        ///                  responsibility to ensure that.
        /// \example
        /// \code{.cpp}
        /// struct Foo {
        ///   void bar () { cout << "Called!" << endl; }
        /// };
        ///
        /// SimulatorIncompTwophase sim (...);
        /// Foo f;
        /// sim.connect_timestep <Foo, &Foo::bar> (f);
        /// sim.run (...);
        /// \endcode
        template <typename T, void (T::*callback)()>
        void connect_timestep (T& t);

    private:
        struct Impl;
        // Using shared_ptr instead of scoped_ptr since scoped_ptr requires complete type for Impl.
        std::shared_ptr<Impl> pimpl_;

        // implementation which is not templated, and can be in library
        void connect_timestep_impl (boost::function0<void> hook);
    };

} // namespace Opm

#include "SimulatorIncompTwophase_impl.hpp"

#endif // OPM_SIMULATORINCOMPTWOPHASE_HEADER_INCLUDED
