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

#include <fmt/format.h>

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <string_view>
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

    inline DomainOrderingMeasure domainOrderingMeasureFromString(const std::string_view measure)
    {
        if (measure == "residual") {
            return DomainOrderingMeasure::Residual;
        } else if (measure == "maxpressure") {
            return DomainOrderingMeasure::MaxPressure;
        } else if (measure == "averagepressure") {
            return DomainOrderingMeasure::AveragePressure;
        } else {
            throw std::runtime_error(fmt::format(fmt::runtime("Invalid domain ordering '{}' specified"), measure));
        }
    }

    /// Representing a part of a grid, in a way suitable for performing
    /// local solves.
    ///
    /// A subdomain owns a set of cells (its interior) and may in addition
    /// include overlap cells owned by neighbouring subdomains on the same MPI
    /// rank. Overlap cells are unknowns of the local nonlinear problem, but
    /// convergence, ownership of wells and the result of the local solve are
    /// defined by the interior cells only.
    struct SubDomainIndices
    {
        // The index of a subdomain is arbitrary, but can be used by the
        // solvers to keep track of well locations etc.
        int index;
        // All cells of the local problem (interior and overlap), sorted, stored
        // as cell indices in the local numbering of the current MPI rank.
        std::vector<int> cells;
        // The cells owned by this subdomain, sorted. Equal to cells when there
        // is no overlap.
        std::vector<int> interior_cells;
        // Flag for each cell of the current MPI rank, true if the cell is owned
        // by the subdomain. If empty, assumed to be all true. Not required for
        // all nonlinear solver algorithms.
        std::vector<bool> interior;
        // Flag indicating if this domain should be skipped during solves
        bool skip;
        // Enables subdomain solves and linearization using the generic linearization
        // approach (i.e. FvBaseLinearizer as opposed to TpfaLinearizer).
        SubDomainIndices(const int i, std::vector<int>&& c, std::vector<bool>&& in, bool s)
            : SubDomainIndices(i, std::move(c), {}, std::move(in), s)
        {}
        // The interior cells c and the overlap cells ov may be given in any
        // order, but must not have cells in common. Both are sorted here and
        // merged into cells: extracting the local matrix
        // (Details::extractMatrix) requires sorted cells, and wells are
        // assigned to subdomains by binary search in interior_cells.
        SubDomainIndices(const int i, std::vector<int>&& c, std::vector<int>&& ov,
                         std::vector<bool>&& in, bool s)
            : index(i), interior_cells(std::move(c)), interior(std::move(in)), skip(s)
        {
            // The interior cells normally arrive sorted, so check first.
            if (!std::ranges::is_sorted(interior_cells)) {
                std::ranges::sort(interior_cells);
            }
            if (ov.empty()) {
                cells = interior_cells;
            } else {
                std::ranges::sort(ov);
                cells.reserve(interior_cells.size() + ov.size());
                std::ranges::merge(interior_cells, ov, std::back_inserter(cells));
            }
        }

        bool hasOverlap() const
        {
            return cells.size() != interior_cells.size();
        }
    };

    /// Representing a part of a grid, in a way suitable for performing
    /// local solves.
    template <class Grid>
    struct SubDomain : public SubDomainIndices
    {
        // View of all cells of the local problem. Every entity of the view
        // reports Dune::InteriorEntity, including overlap cells, since all of
        // them are solved for; use SubDomainIndices::interior for ownership.
        Dune::SubGridPart<Grid> view;
        // Constructor that moves from its argument.
        SubDomain(const int i, std::vector<int>&& c, std::vector<bool>&& in, Dune::SubGridPart<Grid>&& v, bool s)
            : SubDomainIndices(i, std::move(c), std::move(in), s)
            , view(std::move(v))
        {}
        SubDomain(const int i, std::vector<int>&& c, std::vector<int>&& ov,
                  std::vector<bool>&& in, Dune::SubGridPart<Grid>&& v, bool s)
            : SubDomainIndices(i, std::move(c), std::move(ov), std::move(in), s)
            , view(std::move(v))
        {}
    };

} // namespace Opm


#endif // OPM_SUBDOMAIN_HEADER_INCLUDED
