/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_ANISOTROPICEIKONAL_HEADER_INCLUDED
#define OPM_ANISOTROPICEIKONAL_HEADER_INCLUDED

#include <opm/core/utility/SparseTable.hpp>
#include <vector>

struct UnstructuredGrid;

namespace Opm
{
    /// A solver for the anisotropic eikonal equation:
    ///    \f[ || \nabla u^T M^{-1}(x) \nabla u || = 1 \qquad x \in \Omega \f]
    /// where M(x) is a symmetric positive definite matrix.
    /// The boundary conditions are assumed to be
    ///    \f[ u(x) = 0 \qquad x \in \partial\Omega \f].
    class AnisotropicEikonal2d
    {
    public:
	/// Construct solver.
        /// \param[in] grid      A 2d grid.
        explicit AnisotropicEikonal2d(const UnstructuredGrid& grid);

        /// Solve the eikonal equation.
	/// \param[in]  metric            Array of metric tensors, M, for each cell.
        /// \param[in]  startcells        Array of cells where u = 0 at the centroid.
        /// \param[out] solution          Array of solution to the eikonal equation.
        void solve(const double* metric,
		   const std::vector<int>& startcells,
		   std::vector<double>& solution);
    private:
	const UnstructuredGrid& grid_;
	SparseTable<int> cell_neighbours_;
	typedef std::pair<double, int> ValueAndCell;
	std::vector<ValueAndCell> considered_;
	std::vector<char> is_considered_;

	double computeValue(const int cell) const;

	const ValueAndCell& topConsidered() const;
	void pushConsidered(const ValueAndCell& vc);
	void popConsidered();
    };

} // namespace Opm


#endif // OPM_ANISOTROPICEIKONAL_HEADER_INCLUDED
