/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2013 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include "ImpesTPFAAD.hpp"

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>

#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <opm/core/simulator/initState.hpp>

#include <algorithm>

namespace {
    template <class Geology, class Vector>
    class DerivedGeology {
    public:
        typedef Vector V;

        DerivedGeology(const UnstructuredGrid& grid,
                       const Geology&          geo)
            : pvol_ (grid.number_of_cells)
            , trans_(grid.number_of_faces)
        {
            // Pore volume
            const typename Vector::Index nc = grid.number_of_cells;
            std::transform(grid.cell_volumes, grid.cell_volumes + nc,
                           geo.porosity(), pvol_.data(),
                           std::multiplies<double>());

            // Transmissibility
            Vector htrans(grid.cell_facepos[nc]);
            UnstructuredGrid* ug = const_cast<UnstructuredGrid*>(& grid);
            tpfa_htrans_compute(ug, geo.permeability(), htrans.data());
            tpfa_trans_compute (ug, htrans.data()     , trans_.data());
        }

        const Vector& poreVolume()       const { return pvol_ ; }
        const Vector& transmissibility() const { return trans_; }

    private:
        Vector pvol_ ;
        Vector trans_;
    };
}

int
main(int argc, char* argv[])
{
    const Opm::parameter::ParameterGroup param(argc, argv, false);
    const Opm::GridManager               gm(3, 3);

    const UnstructuredGrid*              g  = gm.c_grid();
    const int                            nc = g->number_of_cells;
    const Opm::BlackoilPropertiesBasic   props(param, 2, nc);

    typedef AutoDiff::ForwardBlock<double>      ADB;
    typedef Opm::BlackoilPropertiesInterface    Geology;
    typedef DerivedGeology<Geology, ADB::V>     GeoProps;
    typedef Opm::BlackoilPropertiesInterface    BOFluid;
    typedef Opm::ImpesTPFAAD<BOFluid, GeoProps> PSolver;

    GeoProps geo(*g, props);
    PSolver  ps (*g, props, geo);

    Opm::BlackoilState state;
    initStateBasic(*g, props, param, 0.0, state);
    initBlackoilSurfvol(*g, props, state);

    ps.solve(1.0, state);

    return 0;
}
