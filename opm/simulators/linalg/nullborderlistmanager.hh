// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::Linear::NullBorderListCreator
 */
#ifndef EWOMS_NULL_BORDER_LIST_MANAGER_HH
#define EWOMS_NULL_BORDER_LIST_MANAGER_HH

#include "overlaptypes.hh"

#include <opm/material/common/Unused.hpp>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

#include <algorithm>

namespace Opm {
namespace Linear {
/*!
 * \brief This is a grid manager which does not create any border list.
 *
 * This means that discretizations using this grid manager cannot be
 * used for parallel computations!
 */
template <class GridView, class DofMapper>
class NullBorderListCreator
{
public:
    NullBorderListCreator(const GridView& gridView,
                          const DofMapper& map OPM_UNUSED)
    {
        if (gridView.comm().size() > 1)
            throw std::runtime_error("The used model is not usable for parallel computations");
    }

    // Access to the border list.
    const BorderList& borderList() const
    { return borderList_; }

private:
    BorderList borderList_;
};

} // namespace Linear
} // namespace Opm

#endif
