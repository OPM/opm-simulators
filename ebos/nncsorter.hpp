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

#ifndef EWOMS_EBOS_NNCSORTER_HPP
#define EWOMS_EBOS_NNCSORTER_HPP


#include <opm/parser/eclipse/EclipseState/Grid/NNC.hpp>

#include <vector>

namespace Ewoms
{
/// \brief Scale NNC data wit informtion form EDITNNC and sort it.
/// \param nncData The NNC data as provided by the deck.
/// \param editnncData The EDITNNC data as provided by the deck.
/// \return A lexicographically sorted vector of the scaled NNC data.
///         For each entry entry.cell1<entry.cell2 will hold for convenience.
std::vector<Opm::NNCdata> sortNncAndApplyEditnnc(const std::vector<Opm::NNCdata>& nncData, std::vector<Opm::NNCdata> editnncData,
                                                 bool log = true);
}
#endif
