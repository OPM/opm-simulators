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


    template <class SolutionState, class WellState>
    StandardWells<SolutionState, WellState>::
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





    template <class SolutionState, class WellState>
    StandardWells<SolutionState, WellState>::
    StandardWells(const Wells* wells_arg)
      : wells_(wells_arg)
      , wops_(wells_arg)
      , well_perforation_densities_(Vector())
      , well_perforation_pressure_diffs_(Vector())
    {
    }





    template <class SolutionState, class WellState>
    const Wells&
    StandardWells<SolutionState, WellState>::
    wells() const
    {
        assert(wells_ != 0);
        return *(wells_);
    }





    template <class SolutionState, class WellState>
    bool
    StandardWells<SolutionState, WellState>::
    wellsActive() const
    {
        return wells_active_;
    }





    template <class SolutionState, class WellState>
    void
    StandardWells<SolutionState, WellState>::
    setWellsActive(const bool wells_active)
    {
        wells_active_ = wells_active;
    }





    template <class SolutionState, class WellState>
    bool
    StandardWells<SolutionState, WellState>::
    localWellsActive() const
    {
        return wells_ ? (wells_->number_of_wells > 0 ) : false;
    }





    template <class SolutionState, class WellState>
    const typename StandardWells<SolutionState, WellState>::WellOps&
    StandardWells<SolutionState, WellState>::
    wellOps() const
    {
        return wops_;
    }





    template <class SolutionState, class WellState>
    Vector&
    StandardWells<SolutionState, WellState>::
    wellPerforationDensities()
    {
        return well_perforation_densities_;
    }





    template <class SolutionState, class WellState>
    const Vector&
    StandardWells<SolutionState, WellState>::
    wellPerforationDensities() const
    {
        return well_perforation_densities_;
    }





    template <class SolutionState, class WellState>
    Vector&
    StandardWells<SolutionState, WellState>::
    wellPerforationPressureDiffs()
    {
        return well_perforation_pressure_diffs_;
    }





    template <class SolutionState, class WellState>
    const Vector&
    StandardWells<SolutionState, WellState>::
    wellPerforationPressureDiffs() const
    {
        return well_perforation_pressure_diffs_;
    }

}
