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

#ifndef OPM_TRANSPORTMODELINTERFACE_HEADER_INCLUDED
#define OPM_TRANSPORTMODELINTERFACE_HEADER_INCLUDED

#include <vector>

struct UnstructuredGrid;

namespace Opm
{

    /// Interface for reordering transport models.
    /// A transport model must provide the solveSingleCell() and
    /// solveMultiCell methods, and is expected to implement a solve()
    /// method that will have an interface geared to the model's
    /// needs. (The solve() method is therefore not virtual in this
    /// class.) The reorderAndTransport() method is provided as an aid
    /// to implementing solve() in subclasses, together with the
    /// sequence() and components() methods for accessing the ordering.
    class TransportModelInterface
    {
    public:
	virtual ~TransportModelInterface() {}
    private:
	virtual void solveSingleCell(const int cell) = 0;
	virtual void solveMultiCell(const int num_cells, const int* cells) = 0;
    protected:
	void reorderAndTransport(const UnstructuredGrid& grid, const double* darcyflux);
        const std::vector<int>& sequence() const;
        const std::vector<int>& components() const;
    private:
        std::vector<int> sequence_;
        std::vector<int> components_;
    };


} // namespace Opm

#endif // OPM_TRANSPORTMODELINTERFACE_HEADER_INCLUDED
