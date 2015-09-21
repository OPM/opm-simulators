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

#ifndef OPM_SIMULATORFULLYIMPLICITBLACKOILMULTISEGMENT_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITBLACKOILMULTISEGMENT_HEADER_INCLUDED

#include "SimulatorBase.hpp"

#include "NewtonSolver.hpp"

#include <opm/autodiff/BlackoilMultiSegmentModel.hpp>
#include <opm/autodiff/WellStateMultiSegment.hpp>

namespace Opm {

template <class GridT>
class SimulatorFullyImplicitBlackoilMultiSegment;

template <class GridT>
struct SimulatorTraits<SimulatorFullyImplicitBlackoilMultiSegment<GridT> >
{
    typedef WellStateMultiSegment WellState;
    typedef BlackoilState ReservoirState;
    typedef BlackoilOutputWriter OutputWriter;
    typedef GridT Grid;
    typedef BlackoilMultiSegmentModel<Grid> Model;
    typedef NewtonSolver<Model> Solver;
};

/// a simulator for the blackoil model
template <class GridT>
class SimulatorFullyImplicitBlackoilMultiSegment
    : public SimulatorBase<SimulatorFullyImplicitBlackoilMultiSegment<GridT> >
{
    typedef SimulatorBase<SimulatorFullyImplicitBlackoilMultiSegment<GridT> > Base;
    typedef SimulatorFullyImplicitBlackoilMultiSegment<GridT> ThisType;
    typedef SimulatorTraits<ThisType> Traits;
    typedef typename Traits::ReservoirState ReservoirState;
    typedef typename Traits::WellState WellState;
public:
    // forward the constructor to the base class
    SimulatorFullyImplicitBlackoilMultiSegment(const parameter::ParameterGroup& param,
                                   const GridT& grid,
                                   const DerivedGeology& geo,
                                   BlackoilPropsAdInterface& props,
                                   const RockCompressibility* rock_comp_props,
                                   NewtonIterationBlackoilInterface& linsolver,
                                   const double* gravity,
                                   const bool disgas,
                                   const bool vapoil,
                                   std::shared_ptr<EclipseState> eclipse_state,
                                   BlackoilOutputWriter& output_writer,
                                   const std::vector<double>& threshold_pressures_by_face)
    : Base(param, grid, geo, props, rock_comp_props, linsolver, gravity, disgas, vapoil,
           eclipse_state, output_writer, threshold_pressures_by_face)
    {}

    SimulatorReport run(SimulatorTimer& timer,
                        ReservoirState& state);

protected:
    using Base::output_writer_;
    using Base::param_;
    using Base::solver_;
    using Base::terminal_output_;
    using Base::eclipse_state_;
    using Base::grid_;
    using Base::props_;
    using Base::is_parallel_run_;
    using Base::allcells_;
};

} // namespace Opm

#include "SimulatorFullyImplicitBlackoilMultiSegment_impl.hpp"

#endif // OPM_SIMULATORFULLYIMPLICITBLACKOILMULTISEGMENT_HEADER_INCLUDED
