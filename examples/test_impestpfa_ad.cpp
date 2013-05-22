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

#include <config.h>

#define HACK_INCOMPRESSIBLE_GRAVITY 0

#include <opm/autodiff/ImpesTPFAAD.hpp>
#include <opm/autodiff/BlackoilPropsAd.hpp>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>

#include <opm/core/props/BlackoilPropertiesBasic.hpp>

#include <opm/core/pressure/tpfa/trans_tpfa.h>

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>

#include <opm/core/simulator/initState.hpp>

#include <opm/core/wells.h>
// #include <opm/core/WellsManager.hpp>

#include <algorithm>

namespace {
    template <class Geology, class Vector>
    class DerivedGeology {
    public:
        typedef Vector V;

        DerivedGeology(const UnstructuredGrid& grid,
                       const Geology&          geo ,
                       const double*           grav = 0)
            : pvol_ (grid.number_of_cells)
            , trans_(grid.number_of_faces)
            , gpot_ (grid.cell_facepos[ grid.number_of_cells ])
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

            if (grav != 0) {
                const typename Vector::Index nd = grid.dimensions;

                for (typename Vector::Index c = 0; c < nc; ++c) {
                    const double* const cc = & grid.cell_centroids[c*nd + 0];

                    const int* const p = grid.cell_facepos;
                    for (int i = p[c]; i < p[c + 1]; ++i) {
                        const int f = grid.cell_faces[i];

                        const double* const fc = & grid.face_centroids[f*nd + 0];

                        for (typename Vector::Index d = 0; d < nd; ++d) {
                            gpot_[i] += grav[d] * (fc[d] - cc[d]);
                        }
                    }
                }
            }
        }

        const Vector& poreVolume()       const { return pvol_ ; }
        const Vector& transmissibility() const { return trans_; }
        const Vector& gravityPotential() const { return gpot_ ; }

    private:
        Vector pvol_ ;
        Vector trans_;
        Vector gpot_ ;
    };
}





int
main(int argc, char* argv[])
{
    const Opm::parameter::ParameterGroup param(argc, argv, false);
    const Opm::GridManager               gm(5, 5);

    const UnstructuredGrid*              g  = gm.c_grid();
    const int                            nc = g->number_of_cells;
    const Opm::BlackoilPropertiesBasic   oldprops(param, 2, nc);
    const Opm::BlackoilPropsAd           props(oldprops);

    typedef AutoDiff::ForwardBlock<double>      ADB;
    typedef Opm::BlackoilPropertiesInterface    Geology;
    typedef DerivedGeology<Geology, ADB::V>     GeoProps;
    typedef Opm::BlackoilPropsAd    BOFluid;
    typedef Opm::ImpesTPFAAD<GeoProps> PSolver;

    Wells* wells = create_wells(2, 2, 5);
    const double inj_frac[] = { 1.0, 0.0 };
    const double prod_frac[] = { 0.0, 0.0 };
    const int num_inj = 3;
    const int inj_cells[num_inj] = { 0, 1, 2 };
    const int num_prod = 2;
    const int prod_cells[num_prod] = { 20, 21 };
    const double WI[3] = { 1e-12 };
    bool ok = add_well(INJECTOR, 0.0, num_inj, inj_frac, inj_cells, WI, "Inj", wells);
    ok = ok && add_well(PRODUCER, 0.0, num_prod, prod_frac, prod_cells, WI, "Prod", wells);
    ok = ok && append_well_controls(BHP, 500.0*Opm::unit::barsa, 0, 0, wells);
    ok = ok && append_well_controls(BHP, 200.0*Opm::unit::barsa, 0, 1, wells);
    double oildistr[2] = { 0.0, 1.0 };
    // ok = ok && append_well_controls(SURFACE_RATE, 8.64297e-05, oildistr, 1, wells);
    if (!ok) {
        THROW("Something went wrong with well init.");
    }
    set_current_control(0, 0, wells);
    set_current_control(1, 0, wells);

    double grav[] = { /*1.0*/ 0.0, 0.0 };
    GeoProps geo(*g, oldprops, grav);
    Opm::LinearSolverFactory linsolver(param);
    PSolver  ps (*g, props, geo, *wells, linsolver);

    Opm::BlackoilState state;
    initStateBasic(*g, oldprops, param, 0.0, state);
    initBlackoilSurfvol(*g, oldprops, state);
    Opm::WellState well_state;
    well_state.init(wells, state);

    ps.solve(1.0, state, well_state);

    std::cout << "Cell pressure:" << std::endl;
    std::copy(state.pressure().begin(), state.pressure().end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "Face flux:" << std::endl;
    std::copy(state.faceflux().begin(), state.faceflux().end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "Well bhp pressure:" << std::endl;
    std::copy(well_state.bhp().begin(), well_state.bhp().end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    return 0;
}
