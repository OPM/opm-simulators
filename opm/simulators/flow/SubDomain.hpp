/*
  Copyright 2021 Total SE

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

#ifndef OPM_SUBDOMAIN_HEADER_INCLUDED
#define OPM_SUBDOMAIN_HEADER_INCLUDED

#include <opm/grid/common/SubGridPart.hpp>

#include <utility>
#include <vector>

namespace Opm
{
    //! \brief Solver approach for NLDD.
    enum class DomainSolveApproach {
        Jacobi,
        GaussSeidel
    };

    //! \brief Measure to use for domain ordering.
    enum class DomainOrderingMeasure {
        AveragePressure,
        MaxPressure,
        Residual
    };

    /// Representing a part of a grid, in a way suitable for performing
    /// local solves.
    template <class Grid>
    struct SubDomain
    {
        // The index of a subdomain is arbitrary, but can be used by the
        // solvers to keep track of well locations etc.
        int index;
        // The set of cells that make up a SubDomain, stored as cell indices
        // in the local numbering of the current MPI rank.
        std::vector<int> cells;
        // Flag for each cell of the current MPI rank, true if the cell is part
        // of the subdomain. If empty, assumed to be all true. Not required for
        // all nonlinear solver algorithms.
        std::vector<bool> interior;
        // Enables subdomain solves and linearization using the generic linearization
        // approach (i.e. FvBaseLinearizer as opposed to TpfaLinearizer).
        Dune::SubGridPart<Grid> view;
        // Constructor that moves from its argument.
        SubDomain(const int i, std::vector<int>&& c, std::vector<bool>&& in, Dune::SubGridPart<Grid>&& v)
            : index(i), cells(std::move(c)), interior(std::move(in)), view(std::move(v))
        {}
    };

} // namespace Opm


#endif // OPM_SUBDOMAIN_HEADER_INCLUDED
