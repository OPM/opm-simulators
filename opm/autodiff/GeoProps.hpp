/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_GEOPROPS_HEADER_INCLUDED
#define OPM_GEOPROPS_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <Eigen/Eigen>

namespace Opm
{

    class DerivedGeology
    {
    public:
        typedef Eigen::ArrayXd Vector;

        /// Construct contained derived geological properties
        /// from grid and property information.
        template <class Props>
        DerivedGeology(const UnstructuredGrid& grid,
                       const Props&            props ,
                       const double*           grav = 0)
            : pvol_ (grid.number_of_cells)
            , trans_(grid.number_of_faces)
            , gpot_ (Vector::Zero(grid.cell_facepos[ grid.number_of_cells ], 1))
            , z_(grid.number_of_cells)
        {
            // Pore volume
            const typename Vector::Index nc = grid.number_of_cells;
            std::transform(grid.cell_volumes, grid.cell_volumes + nc,
                           props.porosity(), pvol_.data(),
                           std::multiplies<double>());

            // Transmissibility
            Vector htrans(grid.cell_facepos[nc]);
            UnstructuredGrid* ug = const_cast<UnstructuredGrid*>(& grid);
            tpfa_htrans_compute(ug, props.permeability(), htrans.data());
            tpfa_trans_compute (ug, htrans.data()     , trans_.data());

            // Compute z coordinates
            for (int c = 0; c<nc; ++c){
                z_[c] = grid.cell_centroids[c*3 + 2];
            }


            // Gravity potential
            std::fill(gravity_, gravity_ + 3, 0.0);
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
                std::copy(grav, grav + nd, gravity_);
            }
        }

        const Vector& poreVolume()       const { return pvol_ ; }
        const Vector& transmissibility() const { return trans_; }
        const Vector& gravityPotential() const { return gpot_ ; }
        const Vector& z()                const  {return z_;}
        const double* gravity()          const { return gravity_; }

    private:
        Vector pvol_ ;
        Vector trans_;
        Vector gpot_ ;
        Vector z_;
        double gravity_[3]; // Size 3 even if grid is 2-dim.
    };
}




#endif // OPM_GEOPROPS_HEADER_INCLUDED
