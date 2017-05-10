/*
  Copyright 2015 IRIS AS.

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

#ifndef OPM_SIMULATORFULLYIMPLICITBLACKOILSOLVENT_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITBLACKOILSOLVENT_HEADER_INCLUDED

#include <opm/autodiff/SimulatorBase.hpp>
#include <opm/autodiff/SimulatorFullyImplicitBlackoilOutput.hpp>
#include <opm/autodiff/BlackoilSolventModel.hpp>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/SolventPropsAdFromDeck.hpp>
#include <opm/autodiff/RateConverter.hpp>
#include <opm/autodiff/NonlinearSolver.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoilSolvent.hpp>

#include <opm/core/grid.h>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
//#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/core/utility/miscUtilitiesBlackoil.hpp>

#include <opm/core/props/rock/RockCompressibility.hpp>

//#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
//#include <opm/core/transport/reorder/TransportSolverCompressibleTwophaseReorder.hpp>

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
    template <class GridT>
    class SimulatorFullyImplicitBlackoilSolvent;

    class StandardWellsSolvent;

    template<class GridT>
    struct SimulatorTraits<SimulatorFullyImplicitBlackoilSolvent<GridT> >
    {
        typedef WellStateFullyImplicitBlackoilSolvent WellState;
        typedef BlackoilState ReservoirState;
        typedef BlackoilOutputWriter OutputWriter;
        typedef GridT Grid;
        typedef BlackoilSolventModel<Grid> Model;
        typedef NonlinearSolver<Model> Solver;
        typedef StandardWellsSolvent WellModel;

    };

    /// Class collecting all necessary components for a blackoil simulation with polymer
    /// injection.
    template <class GridT>
    class SimulatorFullyImplicitBlackoilSolvent
        : public SimulatorBase<SimulatorFullyImplicitBlackoilSolvent<GridT> >
    {
        typedef SimulatorFullyImplicitBlackoilSolvent<GridT> ThisType;
        typedef SimulatorBase<ThisType> BaseType;

        typedef SimulatorTraits<ThisType> Traits;
        typedef typename Traits::Solver Solver;
        typedef typename Traits::WellModel WellModel;

    public:
        SimulatorFullyImplicitBlackoilSolvent(const ParameterGroup& param,
                                              const GridT& grid,
                                              DerivedGeology& geo,
                                              BlackoilPropsAdFromDeck& props,
                                              const SolventPropsAdFromDeck& solvent_props,
                                              const RockCompressibility* rock_comp_props,
                                              NewtonIterationBlackoilInterface& linsolver,
                                              const double* gravity,
                                              const bool disgas,
                                              const bool vapoil,
                                              std::shared_ptr<EclipseState> eclipse_state,
                                              BlackoilOutputWriter& output_writer,
                                              std::shared_ptr< Deck > deck,
                                              const std::vector<double>& threshold_pressures_by_face,
                                              const bool solvent);

        std::unique_ptr<Solver> createSolver(const WellModel& well_model);

        void handleAdditionalWellInflow(SimulatorTimer& timer,
                                        WellsManager& wells_manager,
                                        typename BaseType::WellState& well_state,
                                        const Wells* wells);

    private:
        bool has_solvent_;
        std::shared_ptr< Deck > deck_;
        SolventPropsAdFromDeck solvent_props_;
        bool is_miscible_;

    };

} // namespace Opm

#include "SimulatorFullyImplicitBlackoilSolvent_impl.hpp"

#endif // OPM_SIMULATORFULLYIMPLICITBLACKOILSOLVENT_HEADER_INCLUDED
