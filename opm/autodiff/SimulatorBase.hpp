/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Andreas Lauser

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

#ifndef OPM_SIMULATORBASE_HEADER_INCLUDED
#define OPM_SIMULATORBASE_HEADER_INCLUDED

#include <opm/autodiff/SimulatorFullyImplicitBlackoilOutput.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/BlackoilModel.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/RateConverter.hpp>

#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
#include <opm/core/simulator/AdaptiveSimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/output/vtk/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>

#include <opm/core/props/rock/RockCompressibility.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/AdaptiveTimeStepping.hpp>
#include <opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/WellProductionProperties.hpp>


#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <cstddef>
#include <cassert>
#include <functional>
#include <memory>
#include <numeric>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm
{

    template <class Simulator>
    struct SimulatorTraits;

    /// Class collecting all necessary components for a two-phase simulation.
    template <class Implementation>
    class SimulatorBase
    {
        typedef SimulatorTraits<Implementation> Traits;

    public:
        typedef typename Traits::ReservoirState ReservoirState;
        typedef typename Traits::WellState WellState;
        typedef typename Traits::OutputWriter OutputWriter;
        typedef typename Traits::Grid Grid;
        typedef typename Traits::Solver Solver;
        typedef typename Traits::WellModel WellModel;

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
        /// \param[in] eclipse_state the object which represents an internalized ECL deck
        /// \param[in] output_writer
        /// \param[in] threshold_pressures_by_face   if nonempty, threshold pressures that inhibit flow
        SimulatorBase(const parameter::ParameterGroup& param,
                      const Grid& grid,
                      DerivedGeology& geo,
                      BlackoilPropsAdInterface& props,
                      const RockCompressibility* rock_comp_props,
                      NewtonIterationBlackoilInterface& linsolver,
                      const double* gravity,
                      const bool disgas,
                      const bool vapoil,
                      std::shared_ptr<EclipseState> eclipse_state,
                      OutputWriter& output_writer,
                      const std::vector<double>& threshold_pressures_by_face);

        /// Run the simulation.
        /// This will run succesive timesteps until timer.done() is true. It will
        /// modify the reservoir and well states.
        /// \param[in,out] timer       governs the requested reporting timesteps
        /// \param[in,out] state       state of reservoir: pressure, fluxes
        /// \return                    simulation report, with timing data
        SimulatorReport run(SimulatorTimer& timer,
                            ReservoirState& state);

    protected:
        Implementation& asImpl() { return *static_cast<Implementation*>(this); }
        const Implementation& asImpl() const { return *static_cast<const Implementation*>(this); }

        void handleAdditionalWellInflow(SimulatorTimer& timer,
                                        WellsManager& wells_manager,
                                        WellState& well_state,
                                        const Wells* wells);

        std::unique_ptr<Solver> createSolver(const WellModel& well_model);

        void
        computeRESV(const std::size_t               step,
                    const Wells*                    wells,
                    const BlackoilState&            x,
                    WellState& xw);

        void
        FIPUnitConvert(const UnitSystem& units,
                       std::vector<V>& fip);
        
        V
        FIPTotals(const std::vector<V>& fip, const ReservoirState& state);

        void
        outputFluidInPlace(const V& oip, const V& cip, const UnitSystem& units, const int reg);

        void computeWellPotentials(const Wells*                    wells,
                                   const WellState& xw,
                                   std::vector<double>& well_potentials);

        void updateListEconLimited(const std::unique_ptr<Solver>& solver,
                                   ScheduleConstPtr schedule,
                                   const int current_step,
                                   const Wells* wells,
                                   const WellState& well_state,
                                   DynamicListEconLimited& list_econ_limited) const;


        // Data.
        typedef RateConverter::
        SurfaceToReservoirVoidage< BlackoilPropsAdInterface,
                                   std::vector<int> > RateConverterType;
        typedef typename Traits::Model Model;
        typedef typename Model::ModelParameters ModelParameters;
        typedef typename Solver::SolverParameters SolverParameters;

        const parameter::ParameterGroup param_;
        ModelParameters model_param_;
        SolverParameters solver_param_;

        // Observed objects.
        const Grid& grid_;
        BlackoilPropsAdInterface& props_;
        const RockCompressibility* rock_comp_props_;
        const double* gravity_;
        // Solvers
        DerivedGeology& geo_;
        NewtonIterationBlackoilInterface& solver_;
        // Misc. data
        std::vector<int> allcells_;
        const bool has_disgas_;
        const bool has_vapoil_;
        bool       terminal_output_;
        // eclipse_state
        std::shared_ptr<EclipseState> eclipse_state_;
        // output_writer
        OutputWriter& output_writer_;
        RateConverterType rateConverter_;
        // Threshold pressures.
        std::vector<double> threshold_pressures_by_face_;
        // Whether this a parallel simulation or not
        bool is_parallel_run_;
    };

} // namespace Opm

#include "SimulatorBase_impl.hpp"
#endif // OPM_SIMULATORBASE_HEADER_INCLUDED
