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
#include <boost/heap/fibonacci_heap.hpp>

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
        // Grid and topology.
	const UnstructuredGrid& grid_;
	SparseTable<int> cell_neighbours_;

        // Keep track of accepted cells.
	std::vector<char> is_accepted_;
        std::set<int> accepted_front_;

        // Keep track of considered cells.
	typedef std::pair<double, int> ValueAndCell;
        typedef boost::heap::compare<std::greater<ValueAndCell>> Comparator;
        typedef boost::heap::fibonacci_heap<ValueAndCell, Comparator> Heap;
        Heap considered_;
        typedef Heap::handle_type HeapHandle;
        std::map<int, HeapHandle> considered_handles_;
	std::vector<char> is_considered_;

        bool isClose(const int c1, const int c2, const double* metric) const;
	double computeValue(const int cell, const double* metric, const double* solution) const;
	double computeFromLine(const int cell, const int from, const double* metric, const double* solution) const;
	double computeFromTri(const int cell, const int n0, const int n1, const double* metric, const double* solution) const;

	const ValueAndCell& topConsidered() const;
	void pushConsidered(const ValueAndCell& vc);
	void popConsidered();
    };

} // namespace Opm


#endif // OPM_ANISOTROPICEIKONAL_HEADER_INCLUDED
