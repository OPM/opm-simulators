/*
  Copyright 2013, 2015, 2016 SINTEF ICT, Applied Mathematics.
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

#ifndef OPM_SIMULATORSEQUENTIALBLACKOIL_HEADER_INCLUDED
#define OPM_SIMULATORSEQUENTIALBLACKOIL_HEADER_INCLUDED

#include <opm/autodiff/SimulatorBase.hpp>
#include <opm/autodiff/NonlinearSolver.hpp>
#include <opm/autodiff/BlackoilSequentialModel.hpp>
#include <opm/autodiff/StandardWells.hpp>

namespace Opm {

template <class GridT>
class SimulatorSequentialBlackoil;
class StandardWells;

template <class GridT>
struct SimulatorTraits<SimulatorSequentialBlackoil<GridT> >
{
    typedef WellStateFullyImplicitBlackoil WellState;
    typedef BlackoilState ReservoirState;
    typedef BlackoilOutputWriter OutputWriter;
    typedef GridT Grid;
    typedef BlackoilSequentialModel<Grid, StandardWells> Model;
    typedef NonlinearSolver<Model> Solver;
    typedef StandardWells WellModel;
};

/// a simulator for the blackoil model
template <class GridT>
class SimulatorSequentialBlackoil
    : public SimulatorBase<SimulatorSequentialBlackoil<GridT> >
{
    typedef SimulatorBase<SimulatorSequentialBlackoil<GridT> > Base;
public:
    // forward the constructor to the base class
    SimulatorSequentialBlackoil(const ParameterGroup& param,
                                const typename Base::Grid& grid,
                                DerivedGeology& geo,
                                BlackoilPropsAdFromDeck& props,
                                const RockCompressibility* rock_comp_props,
                                NewtonIterationBlackoilInterface& linsolver,
                                const double* gravity,
                                const bool disgas,
                                const bool vapoil,
                                std::shared_ptr<EclipseState> eclipse_state,
                                BlackoilOutputWriter& output_writer,
                                const std::vector<double>& threshold_pressures_by_face)
    : Base(param, grid, geo, props, rock_comp_props, linsolver, gravity, disgas, vapoil,
           eclipse_state, output_writer, threshold_pressures_by_face,
           // names of deactivated wells in parallel run
           std::unordered_set<std::string>())
    {}
};

} // namespace Opm

#endif // OPM_SIMULATORSEQUENTIALBLACKOIL_HEADER_INCLUDED
