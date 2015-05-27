/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
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

#ifndef OPM_SIMULATORFULLYIMPLICITBLACKOILPOLYMER_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITBLACKOILPOLYMER_HEADER_INCLUDED

#include <opm/autodiff/SimulatorBase.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoilOutput.hpp>
#include <opm/polymer/fullyimplicit/BlackoilPolymerModel.hpp>
#include <opm/polymer/fullyimplicit/WellStateFullyImplicitBlackoilPolymer.hpp>
#include <opm/polymer/fullyimplicit/PolymerBlackoilOutputWriter.hpp>
#include <opm/polymer/PolymerBlackoilState.hpp>
#include <opm/polymer/PolymerInflow.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/RateConverter.hpp>
#include <opm/autodiff/NewtonSolver.hpp>

#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
//#include <opm/core/simulator/AdaptiveSimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/io/vtk/writeVtkData.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>

#include <opm/core/props/rock/RockCompressibility.hpp>

//#include <opm/core/simulator/AdaptiveTimeStepping.hpp>
#include <opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/WellProductionProperties.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

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
    template<class GridT>
    class SimulatorFullyImplicitBlackoilPolymer;

    template<class GridT>
    struct SimulatorTraits<SimulatorFullyImplicitBlackoilPolymer<GridT> >
    {
        typedef WellStateFullyImplicitBlackoilPolymer WellState;
        typedef PolymerBlackoilState ReservoirState;
        typedef PolymerBlackoilOutputWriter OutputWriter;
    };

    /// Class collecting all necessary components for a two-phase simulation.
    template<class GridT>
    class SimulatorFullyImplicitBlackoilPolymer
        : public SimulatorBase<GridT, SimulatorFullyImplicitBlackoilPolymer<GridT> >
    {
        typedef SimulatorFullyImplicitBlackoilPolymer<GridT> ThisType;
        typedef SimulatorBase<GridT, ThisType> BaseType;

    public:
        /// \brief The type of the grid that we use.
        typedef GridT Grid;

        SimulatorFullyImplicitBlackoilPolymer(const parameter::ParameterGroup& param,
                                              const Grid& grid,
                                              const DerivedGeology& geo,
                                              BlackoilPropsAdInterface& props,
                                              const PolymerPropsAd& polymer_props,
                                              const RockCompressibility* rock_comp_props,
                                              NewtonIterationBlackoilInterface& linsolver,
                                              const double* gravity,
                                              const bool disgas,
                                              const bool vapoil,
                                              const bool polymer,
                                              std::shared_ptr<EclipseState> eclipse_state,
                                              PolymerBlackoilOutputWriter& output_writer,
                                              Opm::DeckConstPtr& deck,
                                              const std::vector<double>& threshold_pressures_by_face);

#if 0
        SimulatorReport run(SimulatorTimer& timer,
                            PolymerBlackoilState& state);
#endif

        void computeRESV(const std::size_t step,
                         const Wells* wells,
                         const PolymerBlackoilState& x,
                         WellState& xw)
        {
            BaseType::computeRESV(step, wells, x.blackoilState(), xw);
        }

        void handleAdditionalWellInflow(SimulatorTimer& timer,
                                        WellsManager& wells_manager,
                                        typename BaseType::WellState& well_state,
                                        const Wells* wells)
        {
            // compute polymer inflow
            std::unique_ptr<PolymerInflowInterface> polymer_inflow_ptr;
            if (deck_->hasKeyword("WPOLYMER")) {
                if (wells_manager.c_wells() == 0) {
                    OPM_THROW(std::runtime_error, "Cannot control polymer injection via WPOLYMER without wells.");
                }
                polymer_inflow_ptr.reset(new PolymerInflowFromDeck(deck_, BaseType::eclipse_state_, *wells, Opm::UgGridHelpers::numCells(BaseType::grid_), timer.currentStepNum()));
            } else {
                polymer_inflow_ptr.reset(new PolymerInflowBasic(0.0*Opm::unit::day,
                                                                1.0*Opm::unit::day,
                                                                0.0));
            }
            std::vector<double> polymer_inflow_c(Opm::UgGridHelpers::numCells(BaseType::grid_));
            polymer_inflow_ptr->getInflowValues(timer.simulationTimeElapsed(),
                                                timer.simulationTimeElapsed() + timer.currentStepLength(),
                                                polymer_inflow_c);
            well_state.polymerInflow() = polymer_inflow_c;
        }

    private:
        const PolymerPropsAd& polymer_props_;
        bool has_polymer_;
        DeckConstPtr deck_;
    };

} // namespace Opm

#include "SimulatorFullyImplicitBlackoilPolymer_impl.hpp"
#endif // OPM_SIMULATORFULLYIMPLICITBLACKOILPOLYMER_HEADER_INCLUDED
