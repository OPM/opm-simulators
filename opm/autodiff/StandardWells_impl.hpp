/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

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


#include <opm/autodiff/StandardWells.hpp>




namespace Opm
{


    StandardWells::
    WellOps::WellOps(const Wells* wells)
      : w2p(),
        p2w(),
        well_cells()
    {
        if( wells )
        {
            w2p = Eigen::SparseMatrix<double>(wells->well_connpos[ wells->number_of_wells ], wells->number_of_wells);
            p2w = Eigen::SparseMatrix<double>(wells->number_of_wells, wells->well_connpos[ wells->number_of_wells ]);

            const int        nw   = wells->number_of_wells;
            const int* const wpos = wells->well_connpos;

            typedef Eigen::Triplet<double> Tri;

            std::vector<Tri> scatter, gather;
            scatter.reserve(wpos[nw]);
            gather .reserve(wpos[nw]);

            for (int w = 0, i = 0; w < nw; ++w) {
                for (; i < wpos[ w + 1 ]; ++i) {
                    scatter.push_back(Tri(i, w, 1.0));
                    gather .push_back(Tri(w, i, 1.0));
                }
            }

            w2p.setFromTriplets(scatter.begin(), scatter.end());
            p2w.setFromTriplets(gather .begin(), gather .end());

            well_cells.assign(wells->well_cells, wells->well_cells + wells->well_connpos[wells->number_of_wells]);
        }
    }





    StandardWells::StandardWells(const Wells* wells_arg)
      : wells_(wells_arg)
      , wops_(wells_arg)
      , well_perforation_densities_(Vector())
      , well_perforation_pressure_diffs_(Vector())
    {
    }





    const Wells& StandardWells::wells() const
    {
        assert(wells_ != 0);
        return *(wells_);
    }





    bool StandardWells::wellsActive() const
    {
        return wells_active_;
    }





    void StandardWells::setWellsActive(const bool wells_active)
    {
        wells_active_ = wells_active;
    }





    bool StandardWells::localWellsActive() const
    {
        return wells_ ? (wells_->number_of_wells > 0 ) : false;
    }





    const StandardWells::WellOps&
    StandardWells::wellOps() const
    {
        return wops_;
    }





    Vector& StandardWells::wellPerforationDensities()
    {
        return well_perforation_densities_;
    }





    const Vector&
    StandardWells::wellPerforationDensities() const
    {
        return well_perforation_densities_;
    }





    Vector&
    StandardWells::wellPerforationPressureDiffs()
    {
        return well_perforation_pressure_diffs_;
    }





    const Vector&
    StandardWells::wellPerforationPressureDiffs() const
    {
        return well_perforation_pressure_diffs_;
    }




    template<class SolutionState, class WellState>
    void
    StandardWells::
    computePropertiesForWellConnectionPressures(const SolutionState& state,
                                                const WellState& xw,
                                                const BlackoilPropsAdInterface& fluid,
                                                const std::vector<bool>& active,
                                                const std::vector<PhasePresence>& pc,
                                                std::vector<double>& b_perf,
                                                std::vector<double>& rsmax_perf,
                                                std::vector<double>& rvmax_perf,
                                                std::vector<double>& surf_dens_perf)
    {
        const int nperf = wells().well_connpos[wells().number_of_wells];
        const int nw = wells().number_of_wells;

        // Compute the average pressure in each well block
        const Vector perf_press = Eigen::Map<const Vector>(xw.perfPress().data(), nperf);
        Vector avg_press = perf_press*0;
        for (int w = 0; w < nw; ++w) {
            for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                const double p_above = perf == wells().well_connpos[w] ? state.bhp.value()[w] : perf_press[perf - 1];
                const double p_avg = (perf_press[perf] + p_above)/2;
                avg_press[perf] = p_avg;
            }
        }

        const std::vector<int>& well_cells = wellOps().well_cells;

        // Use cell values for the temperature as the wells don't knows its temperature yet.
        const ADB perf_temp = subset(state.temperature, well_cells);

        // Compute b, rsmax, rvmax values for perforations.
        // Evaluate the properties using average well block pressures
        // and cell values for rs, rv, phase condition and temperature.
        const ADB avg_press_ad = ADB::constant(avg_press);
        std::vector<PhasePresence> perf_cond(nperf);
        // const std::vector<PhasePresence>& pc = phaseCondition();
        for (int perf = 0; perf < nperf; ++perf) {
            perf_cond[perf] = pc[well_cells[perf]];
        }
        const PhaseUsage& pu = fluid.phaseUsage();
        DataBlock b(nperf, pu.num_phases);
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            const Vector bw = fluid.bWat(avg_press_ad, perf_temp, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Aqua]) = bw;
        }
        assert(active[Oil]);
        const Vector perf_so =  subset(state.saturation[pu.phase_pos[Oil]].value(), well_cells);
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            const ADB perf_rs = subset(state.rs, well_cells);
            const Vector bo = fluid.bOil(avg_press_ad, perf_temp, perf_rs, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Liquid]) = bo;
            // const V rssat = fluidRsSat(avg_press, perf_so, well_cells);
            const Vector rssat = fluid.rsSat(ADB::constant(avg_press), ADB::constant(perf_so), well_cells).value();
            rsmax_perf.assign(rssat.data(), rssat.data() + nperf);
        } else {
            rsmax_perf.assign(nperf, 0.0);
        }
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            const ADB perf_rv = subset(state.rv, well_cells);
            const Vector bg = fluid.bGas(avg_press_ad, perf_temp, perf_rv, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg;
            // const V rvsat = fluidRvSat(avg_press, perf_so, well_cells);
            const Vector rvsat = fluid.rvSat(ADB::constant(avg_press), ADB::constant(perf_so), well_cells).value();
            rvmax_perf.assign(rvsat.data(), rvsat.data() + nperf);
        } else {
            rvmax_perf.assign(nperf, 0.0);
        }
        // b is row major, so can just copy data.
        b_perf.assign(b.data(), b.data() + nperf * pu.num_phases);

        // Surface density.
        // The compute density segment wants the surface densities as
        // an np * number of wells cells array
        Vector rho = superset(fluid.surfaceDensity(0 , well_cells), Span(nperf, pu.num_phases, 0), nperf*pu.num_phases);
        for (int phase = 1; phase < pu.num_phases; ++phase) {
            rho += superset(fluid.surfaceDensity(phase , well_cells), Span(nperf, pu.num_phases, phase), nperf*pu.num_phases);
        }
        surf_dens_perf.assign(rho.data(), rho.data() + nperf * pu.num_phases);

    }


}
