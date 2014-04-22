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
#include <opm/autodiff/GridHelpers.hpp>
//#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/core/pressure/tpfa/TransTpfa.hpp>
#include <Eigen/Eigen>

namespace Opm
{

    /// Class containing static geological properties that are
    /// derived from grid and petrophysical properties:
    ///   - pore volume
    ///   - transmissibilities
    ///   - gravity potentials
    class DerivedGeology
    {
    public:
        typedef Eigen::ArrayXd Vector;

        /// Construct contained derived geological properties
        /// from grid and property information.
        template <class Props, class Grid>
        DerivedGeology(const Grid&             grid,
                       const Props&            props ,
                       const double*           grav = 0)
            : pvol_ (Opm::AutoDiffGrid::numCells(grid))
            , trans_(Opm::AutoDiffGrid::numFaces(grid))
            , gpot_ (Vector::Zero(Opm::AutoDiffGrid::cell2Faces(grid).noEntries(), 1))
            , z_(Opm::AutoDiffGrid::numCells(grid))
        {
            using namespace Opm::AutoDiffGrid;
            // Pore volume
            const typename Vector::Index nc = numCells(grid);
            std::transform(beginCellVolumes(grid), endCellVolumes(grid),
                           props.porosity(), pvol_.data(),
                           std::multiplies<double>());

            // Transmissibility
            Vector htrans(numCellFaces(grid));
            Grid* ug = const_cast<Grid*>(& grid);
            tpfa_htrans_compute(ug, props.permeability(), htrans.data());
            tpfa_trans_compute (ug, htrans.data()     , trans_.data());

            // Compute z coordinates
            for (int c = 0; c<nc; ++c){
                z_[c] = cellCentroid(grid, c)[2];
            }


            // Gravity potential
            std::fill(gravity_, gravity_ + 3, 0.0);
            if (grav != 0) {
                const typename Vector::Index nd = dimensions(grid);
                typedef typename ADCell2FacesTraits<Grid>::Type Cell2Faces;
                Cell2Faces c2f=cell2Faces(grid);

                for (typename Vector::Index c = 0; c < nc; ++c) {
                    const double* const cc = cellCentroid(grid, c);

                    typename Cell2Faces::row_type faces=c2f[c];
                    typedef typename Cell2Faces::row_type::iterator Iter;

                    for (Iter f=faces.begin(), end=faces.end(); f!=end; ++f) {
                        const double* const fc = faceCentroid(grid, *f);

                        for (typename Vector::Index d = 0; d < nd; ++d) {
                            gpot_[f-faces.begin()] += grav[d] * (fc[d] - cc[d]);
                        }
                    }
                }
                std::copy(grav, grav + nd, gravity_);
            }
        }

        const Vector& poreVolume()       const { return pvol_   ;}
        const Vector& transmissibility() const { return trans_  ;}
        const Vector& gravityPotential() const { return gpot_   ;}
        const Vector& z()                const { return z_      ;}
        const double* gravity()          const { return gravity_;}

    private:
        Vector pvol_ ;
        Vector trans_;
        Vector gpot_ ;
        Vector z_;
        double gravity_[3]; // Size 3 even if grid is 2-dim.
    };
}




#endif // OPM_GEOPROPS_HEADER_INCLUDED
