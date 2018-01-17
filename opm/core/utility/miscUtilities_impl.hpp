#include <opm/core/grid/GridHelpers.hpp>
#include <opm/core/wells.h>
#include <opm/core/props/rock/RockCompressibility.hpp>

namespace Opm
{
    /// @brief Estimates a scalar cell velocity from face fluxes.
    /// @param[in]  number_of_cells      The number of cells of the grid
    /// @param[in]  begin_face_centroids Iterator pointing to first face centroid.
    /// @param[in]  face_cells           Mapping from faces to connected cells.
    /// @param[in]  dimensions           The dimensions of the grid.
    /// @param[in]  begin_cell_centroids Iterator pointing to first cell centroid.
    /// @param[in]  face_flux            signed per-face fluxes
    /// @param[out] cell_velocity        the estimated velocities.
    template<class CC, class FC, class FC1, class CV>
    void estimateCellVelocity(int number_of_cells,
                              int number_of_faces,
                              FC begin_face_centroids,
                              FC1 face_cells,
                              CC begin_cell_centroids,
                              CV begin_cell_volumes,
                              int dimension,
                              const std::vector<double>& face_flux,
                              std::vector<double>& cell_velocity)
    {
        cell_velocity.clear();
        cell_velocity.resize(number_of_cells*dimension, 0.0);
        for (int face = 0; face < number_of_faces; ++face) {
            int c[2] = { face_cells(face, 0), face_cells(face, 1) };
            FC fc = UgGridHelpers::increment(begin_face_centroids, face, dimension);
            double flux = face_flux[face];
            for (int i = 0; i < 2; ++i) {
                if (c[i] >= 0) {
                    CC cc = UgGridHelpers::increment(begin_cell_centroids, c[i], dimension);
                    for (int d = 0; d < dimension; ++d) {
                        double v_contrib = UgGridHelpers::getCoordinate(fc, d) - UgGridHelpers::getCoordinate(cc, d);
                        v_contrib *= flux/begin_cell_volumes[c[i]];
                        cell_velocity[c[i]*dimension + d] += (i == 0) ? v_contrib : -v_contrib;
                    }
                }
            }
        }
    }

    template<class T>
    void computePorevolume(int number_of_cells,
                           T begin_cell_volume,
                           const double* porosity,
                           std::vector<double>& porevol)
    {
        porevol.resize(number_of_cells);
        std::transform(porosity, porosity + number_of_cells,
                       begin_cell_volume,
                       porevol.begin(),
                       std::multiplies<double>());
    }

    /// @brief Computes pore volume of all cells in a grid, with rock compressibility effects.
    /// @param[in]  number_of_cells The number of cells of the grid.
    /// @param[in]  porosity  array of grid.number_of_cells porosity values
    /// @param[in]  rock_comp rock compressibility properties
    /// @param[in]  pressure  pressure by cell
    /// @param[out] porevol   the pore volume by cell.
    template<class T>
    void computePorevolume(int number_of_cells,
                           T begin_cell_volumes,
                           const double* porosity,
                           const RockCompressibility& rock_comp,
                           const std::vector<double>& pressure,
                           std::vector<double>& porevol)
    {
        porevol.resize(number_of_cells);
        for (int i = 0; i < number_of_cells; ++i) {
            porevol[i] = porosity[i]*begin_cell_volumes[i]*rock_comp.poroMult(pressure[i]);
        }
    }

    template<class T>
    void computeWDP(const Wells& wells, int number_of_cells, T begin_cell_centroids, const std::vector<double>& saturations,
                    const double* densities, const double gravity, const bool per_grid_cell,
                    std::vector<double>& wdp)
    {
        const int nw = wells.number_of_wells;
        const size_t np = per_grid_cell ?
            saturations.size()/number_of_cells
            : saturations.size()/wells.well_connpos[nw];
        // Simple for now:
        for (int i = 0; i < nw; i++) {
            double depth_ref = wells.depth_ref[i];
            for (int j = wells.well_connpos[i]; j < wells.well_connpos[i + 1]; j++) {
                int cell = wells.well_cells[j];

                // Is this correct wrt. depth_ref?
                double cell_depth = UgGridHelpers
                    ::getCoordinate(UgGridHelpers::increment(begin_cell_centroids, cell, 3), 2);

                double saturation_sum = 0.0;
                for (size_t p = 0; p < np; p++) {
                    if (!per_grid_cell) {
                        saturation_sum += saturations[j * np + p];
                    } else {
                        saturation_sum += saturations[np * cell + p];
                    }
                }
                if (saturation_sum == 0) {
                    saturation_sum = 1.0;
                }
                double density = 0.0;
                for (size_t p = 0; p < np; p++) {
                    if (!per_grid_cell) {
                        density += saturations[j * np + p] * densities[p] / saturation_sum;
                    } else {
                        // Is this a smart way of doing it?
                        density += saturations[np * cell + p] * densities[p] / saturation_sum;
                    }
                }

                // Is the sign correct?
                wdp.push_back(density * (cell_depth - depth_ref) * gravity);
            }
        }
    }
}
