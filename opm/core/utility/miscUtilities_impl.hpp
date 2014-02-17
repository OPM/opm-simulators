#include <opm/core/grid/GridHelpers.hpp>
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
}
