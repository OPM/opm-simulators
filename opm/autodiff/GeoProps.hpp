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

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include "disable_warning_pragmas.h"

#include <Eigen/Eigen>

#include "reenable_warning_pragmas.h"


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
        DerivedGeology(const Grid&              grid,
                       const Props&             props ,
                       Opm::EclipseStateConstPtr eclState,
                       const double*            grav = 0)
            : pvol_ (Opm::AutoDiffGrid::numCells(grid))
            , trans_(Opm::AutoDiffGrid::numFaces(grid))
            , gpot_ (Vector::Zero(Opm::AutoDiffGrid::cell2Faces(grid).noEntries(), 1))
            , z_(Opm::AutoDiffGrid::numCells(grid))
        {
            auto multipliers = eclState->getTransMult();

            int numCells = AutoDiffGrid::numCells(grid);
            int numFaces = AutoDiffGrid::numFaces(grid);
            const int *cartDims = AutoDiffGrid::cartDims(grid);
            int numCartesianCells =
                cartDims[0]
                * cartDims[1]
                * cartDims[2];

            // get the pore volume multipliers from the EclipseState
            std::vector<double> multpv(numCartesianCells, 1.0);
            if (eclState->hasDoubleGridProperty("MULTPV")) {
                multpv = eclState->getDoubleGridProperty("MULTPV")->getData();
            }

            // get the net-to-gross cell thickness from the EclipseState
            std::vector<double> ntg(numCartesianCells, 1.0);
            if (eclState->hasDoubleGridProperty("NTG")) {
                multpv = eclState->getDoubleGridProperty("NTG")->getData();
            }

            // Pore volume
            const auto &poreVolumes = pvol_.data();
            for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
                int cartesianCellIdx = AutoDiffGrid::globalCell(grid)[cellIdx];
                poreVolumes[cellIdx] =
                    props.porosity()[cellIdx]
                    * multpv[cartesianCellIdx]
                    * ntg[cartesianCellIdx]
                    * AutoDiffGrid::cellVolume(grid, cellIdx);
            }

            // Transmissibility
            Vector htrans(AutoDiffGrid::numCellFaces(grid));
            Grid* ug = const_cast<Grid*>(& grid);
            tpfa_htrans_compute(ug, props.permeability(), htrans.data());

            // multiply the half-face transmissibilities with the appropriate NTG values
            for (int cellIdx = 0; cellIdx < numCells; ++cellIdx) {
                for (int cellFaceIdx = grid.cell_facepos[cellIdx];
                     cellFaceIdx < grid.cell_facepos[cellIdx + 1];
                     ++ cellFaceIdx)
                {
                    int cartesianCellIdx = grid.global_cell[cellIdx];

                    // std::cout << "cellFaceIdx: " << cellFaceIdx << ", cellIdx: " << cellIdx << "\n";
                    switch (grid.cell_facetag[cellFaceIdx]) {
                    case 0: // left
                    case 1: // right
                    case 2: // front
                    case 3: // back
                        htrans.data()[cellFaceIdx] *= ntg[cartesianCellIdx];
                        break;
                    }
                }
            }

            // combine the half-face transmissibilites into the final face
            // transmissibilites.
            tpfa_trans_compute (ug, htrans.data()     , trans_.data());

            // multiply the face transmissibilities with the appropriate transmissibility
            // multipliers
            for (int faceIdx = 0; faceIdx < numFaces; faceIdx++) {
                // get the two cells adjacent to the face and ignore all boundary faces
                int insideCompressedCellIdx = grid.face_cells[2*faceIdx + 0];
                int outsideCompressedCellIdx = grid.face_cells[2*faceIdx + 1];
                if (insideCompressedCellIdx < 0 || outsideCompressedCellIdx  < 0)
                    continue;

                // convert the compressed cell indices to indices of the logically
                // cartesian grid
                int insideCartesianCellIdx = AutoDiffGrid::globalCell(grid)[insideCompressedCellIdx];
                int outsideCartesianCellIdx = AutoDiffGrid::globalCell(grid)[outsideCompressedCellIdx];

                // find out the orientation of the face in the Cartesian grid from the
                // inside cell's perspective. Note that -- as an unlikely corner case --
                // the face could exhibit a different Cartesian orientation for the
                // outside cell. We don't consider this case here, though...
                int insideFaceTag = -1;
                for (int cellFaceIdx = grid.cell_facepos[insideCompressedCellIdx];
                     cellFaceIdx < grid.cell_facepos[insideCompressedCellIdx + 1];
                     ++ cellFaceIdx)
                {
                    if (faceIdx == grid.cell_faces[cellFaceIdx]) {
                        insideFaceTag = grid.cell_facetag[cellFaceIdx];
                        break;
                    }
                }
                // we need a face tag!
                assert(insideFaceTag >= 0);

                double multiplier = 1.0;
                switch (insideFaceTag) {
                case 0: // left
                    multiplier *= multipliers->getMultiplier(insideCartesianCellIdx, FaceDir::XMinus);
                    multiplier *= multipliers->getMultiplier(outsideCartesianCellIdx, FaceDir::XPlus);
                    break;

                case 1: // right
                    multiplier *= multipliers->getMultiplier(insideCartesianCellIdx, FaceDir::XPlus);
                    multiplier *= multipliers->getMultiplier(outsideCartesianCellIdx, FaceDir::XMinus);
                    break;

                case 2: // front
                    multiplier *= multipliers->getMultiplier(insideCartesianCellIdx, FaceDir::YMinus);
                    multiplier *= multipliers->getMultiplier(outsideCartesianCellIdx, FaceDir::YPlus);
                    break;
                case 3: // back
                    multiplier *= multipliers->getMultiplier(insideCartesianCellIdx, FaceDir::YPlus);
                    multiplier *= multipliers->getMultiplier(outsideCartesianCellIdx, FaceDir::YMinus);
                    break;

                case 4: // top
                    multiplier *= multipliers->getMultiplier(insideCartesianCellIdx, FaceDir::ZMinus);
                    multiplier *= multipliers->getMultiplier(outsideCartesianCellIdx, FaceDir::ZPlus);
                    break;
                case 5: // bottom
                    multiplier *= multipliers->getMultiplier(insideCartesianCellIdx, FaceDir::ZPlus);
                    multiplier *= multipliers->getMultiplier(outsideCartesianCellIdx, FaceDir::ZMinus);
                    break;

                default:
                    // all cells need to be hexahedrons in logically cartesian
                    // coordinates!
                    assert(false);
                }

                trans_.data()[faceIdx] *= multiplier;
            }

            // Compute z coordinates
            for (int c = 0; c<numCells; ++c){
                z_[c] = AutoDiffGrid::cellCentroid(grid, c)[2];
            }


            // Gravity potential
            std::fill(gravity_, gravity_ + 3, 0.0);
            if (grav != 0) {
                const typename Vector::Index nd = AutoDiffGrid::dimensions(grid);
                typedef typename AutoDiffGrid::ADCell2FacesTraits<Grid>::Type Cell2Faces;
                Cell2Faces c2f=AutoDiffGrid::cell2Faces(grid);

                for (typename Vector::Index c = 0; c < numCells; ++c) {
                    const double* const cc = AutoDiffGrid::cellCentroid(grid, c);

                    typename Cell2Faces::row_type faces=c2f[c];
                    typedef typename Cell2Faces::row_type::iterator Iter;

                    for (Iter f=faces.begin(), end=faces.end(); f!=end; ++f) {
                        const double* const fc = AutoDiffGrid::faceCentroid(grid, *f);

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
