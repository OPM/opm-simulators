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

#ifndef OPM_GRAVITYCOLUMNSOLVERPOLYMER_HEADER_INCLUDED
#define OPM_GRAVITYCOLUMNSOLVERPOLYMER_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <vector>
#include <map>

namespace Opm
{

    /// Class for doing gravity segregation (only),
    /// on a vertical column of cells.
    template <class Model>
    class GravityColumnSolverPolymer
    {
    public:
	/// Note: the model will be changed since it stores computed
	/// quantities in itself, such as mobilities.
	GravityColumnSolverPolymer(Model& model,
			    const UnstructuredGrid& grid,
			    const double tol,
			    const int maxit);

	/// \param[in] columns         for each column (with logical cartesian indices as key),
	///                            contains the cells on which to solve the segregation
	///                            problem. For each column, its cells must be in a single
	///                            vertical column, and ordered
	///                            (direction doesn't matter).
	void solve(const std::pair<std::vector<int>, std::vector<std::vector<int> > >& columns,
		   const double dt,
		   std::vector<double>& s,
		   std::vector<double>& c,
		   std::vector<double>& cmax);

    private:
	void solveSingleColumn(const std::vector<int>& column_cells,
			       const double dt,
			       std::vector<double>& s,
			       std::vector<double>& c,
			       std::vector<double>& cmax,
			       std::vector<double>& sol_vec
 			       );
	Model& model_;
	const UnstructuredGrid& grid_;
	const double tol_;
	const int maxit_;
};

} // namespace Opm

#include <opm/polymer/GravityColumnSolverPolymer_impl.hpp>

#endif // OPM_GRAVITYCOLUMNSOLVERPOLYMER_HEADER_INCLUDED
