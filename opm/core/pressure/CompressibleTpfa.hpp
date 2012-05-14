/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_COMPRESSIBLETPFA_HEADER_INCLUDED
#define OPM_COMPRESSIBLETPFA_HEADER_INCLUDED


#include <vector>

struct UnstructuredGrid;
struct cfs_tpfa_res_data;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm
{

    class LinearSolverInterface;

    /// Encapsulating a tpfa pressure solver for the compressible-fluid case.
    /// Supports gravity, wells and simple sources as driving forces.
    /// Below we use the shortcuts D for the number of dimensions, N
    /// for the number of cells and F for the number of faces.
    class CompressibleTpfa
    {
    public:
	/// Construct solver.
	/// \param[in] g             A 2d or 3d grid.
	/// \param[in] permeability  Array of permeability tensors, the array
	///                          should have size N*D^2, if D == g.dimensions
	///                          and N == g.number_of_cells.
	/// \param[in] gravity       Gravity vector. If nonzero, the array should
	///                          have D elements.
        /// \param[in] wells         The wells argument. Will be used in solution, 
        ///                          is ignored if NULL
        /// \param[in] num_phases    Must be 2 or 3.
	CompressibleTpfa(const UnstructuredGrid& g,
                         const double* permeability,
                         const double* gravity,
                         const LinearSolverInterface& linsolver,
                         const struct Wells* wells,
                         const int num_phases);

	/// Destructor.
	~CompressibleTpfa();

        void solve();

    private:
        void computeDynamicData();
        void assemble();
        void solveIncrement();

	void computeResults(std::vector<double>& pressure,
                            std::vector<double>& faceflux,
                            std::vector<double>& well_bhp,
                            std::vector<double>& well_rate);

        // ------ Data that will remain unmodified after construction. ------
	const UnstructuredGrid& grid_;
        const LinearSolverInterface& linsolver_;
	std::vector<double> htrans_;
	std::vector<double> trans_ ;
        const Wells* wells_;   // Outside may modify controls (only) between calls to solve().

        // ------ Internal data for the cfs_tpfa_res solver. ------
	struct cfs_tpfa_res_data* h_;

        // ------ Data that will be modified for every solve. ------

        // ------ Data that will be modified for every solver iteration. ------
        // The update to be applied to the pressures (cell and bhp).
        std::vector<double> pressure_increment_;



        // Gravity and capillary contributions (per face).
        std::vector<double> gravcapf_;


    };

} // namespace Opm


#endif // OPM_COMPRESSIBLETPFA_HEADER_INCLUDED
