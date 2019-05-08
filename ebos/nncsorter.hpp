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

#include <opm/common/OpmLog/OpmLog.hpp>

#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>

namespace Ewoms {

/// \brief Scale NNC data wit informtion form EDITNNC and sort it.
/// \param nncData The NNC data as provided by the deck.
/// \param editnncData The EDITNNC data as provided by the deck.
/// \return A lexicographically sorted vector of the scaled NNC data.
///         For each entry entry.cell1<entry.cell2 will hold for convenience.
inline std::vector<Opm::NNCdata> sortNncAndApplyEditnnc(const std::vector<Opm::NNCdata>& nncData,
                                                        std::vector<Opm::NNCdata> editnncData,
                                                        bool doLog = true)
{
    auto nncLess =
        [](const Opm::NNCdata& d1, const Opm::NNCdata& d2) {
            return
                (d1.cell1 < d2.cell1)
                || (d1.cell1 == d2.cell1 && d1.cell2 < d2.cell2);
        };

    auto makeCell1LessCell2 =
        [](const Opm::NNCdata& entry) {
            if (entry.cell2 < entry.cell1)
                return Opm::NNCdata(entry.cell2, entry.cell1, entry.trans);
            else
                return entry;
        };

    // We need to make sure that for each entry cell1<=cell2 holds. Otherwise sorting
    // will not make the search more accurate if the engineer chooses to define NNCs
    // differently.
    std::vector<Opm::NNCdata> nncDataCopy(nncData);
    std::transform(nncDataCopy.begin(), nncDataCopy.end(), nncDataCopy.begin(), makeCell1LessCell2);
    std::transform(editnncData.begin(), editnncData.end(), editnncData.begin(), makeCell1LessCell2);
    std::sort(nncDataCopy.begin(), nncDataCopy.end(), nncLess);
    auto candidate = nncDataCopy.begin();

    for (const auto& edit: editnncData) {
        auto printNncWarning =
            [](int c1, int c2) {
                std::ostringstream sstr;
                sstr << "Cannot edit NNC from " << c1 << " to " << c2
                     << " as it does not exist";
                Opm::OpmLog::warning(sstr.str());
            };

        if (candidate == nncDataCopy.end() && doLog) {
            // no more NNCs left
            printNncWarning(edit.cell1, edit.cell2);
            continue;
        }

        if (candidate->cell1 != edit.cell1 || candidate->cell2 != edit.cell2) {
            candidate = std::lower_bound(nncDataCopy.begin(), nncDataCopy.end(), Opm::NNCdata(edit.cell1, edit.cell2, 0), nncLess);
            if (candidate == nncDataCopy.end() && doLog) {
                // no more NNCs left
                printNncWarning(edit.cell1, edit.cell2);
                continue;
            }
        }

        auto firstCandidate = candidate;
        while (candidate != nncDataCopy.end()
               && candidate->cell1 == edit.cell1
               && candidate->cell2 == edit.cell2)
        {
            candidate->trans *= edit.trans;
            ++candidate;
        }

        // start with first match in next iteration to catch case where next
        // EDITNNC is for same pair.
        candidate = firstCandidate;
    }

    return nncDataCopy;
}

}
#endif
