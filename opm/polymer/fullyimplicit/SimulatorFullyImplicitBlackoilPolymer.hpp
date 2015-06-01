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
    template <class GridT>
    class SimulatorFullyImplicitBlackoilPolymer;

    template<class GridT>
    struct SimulatorTraits<SimulatorFullyImplicitBlackoilPolymer<GridT> >
    {
        typedef WellStateFullyImplicitBlackoilPolymer WellState;
        typedef PolymerBlackoilState ReservoirState;
        typedef BlackoilOutputWriter OutputWriter;
        typedef GridT Grid;
        typedef BlackoilPolymerModel<Grid> Model;
        typedef NewtonSolver<Model> Solver;
    };

    /// Class collecting all necessary components for a blackoil simulation with polymer
    /// injection.
    template <class GridT>
    class SimulatorFullyImplicitBlackoilPolymer
        : public SimulatorBase<SimulatorFullyImplicitBlackoilPolymer<GridT> >
    {
        typedef SimulatorFullyImplicitBlackoilPolymer<GridT> ThisType;
        typedef SimulatorBase<ThisType> BaseType;

        typedef SimulatorTraits<ThisType> Traits;
        typedef typename Traits::Solver Solver;

    public:
        SimulatorFullyImplicitBlackoilPolymer(const parameter::ParameterGroup& param,
                                              const GridT& grid,
                                              const DerivedGeology& geo,
                                              BlackoilPropsAdInterface& props,
                                              const PolymerPropsAd& polymer_props,
                                              const RockCompressibility* rock_comp_props,
                                              NewtonIterationBlackoilInterface& linsolver,
                                              const double* gravity,
                                              const bool disgas,
                                              const bool vapoil,
                                              const bool polymer,
                                              const bool plyshlog,
                                              std::shared_ptr<EclipseState> eclipse_state,
                                              BlackoilOutputWriter& output_writer,
                                              Opm::DeckConstPtr& deck,
                                              const std::vector<double>& threshold_pressures_by_face);

        std::unique_ptr<Solver> createSolver(const Wells* wells);

        void handleAdditionalWellInflow(SimulatorTimer& timer,
                                        WellsManager& wells_manager,
                                        typename BaseType::WellState& well_state,
                                        const Wells* wells);

    private:
        const PolymerPropsAd& polymer_props_;
        bool has_polymer_;
        /// flag for PLYSHLOG keywords
        bool has_plyshlog_;
        DeckConstPtr deck_;
    };

} // namespace Opm

#include "SimulatorFullyImplicitBlackoilPolymer_impl.hpp"

#endif // OPM_SIMULATORFULLYIMPLICITBLACKOILPOLYMER_HEADER_INCLUDED
