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
 *
 * \copydoc Opm::BlackOilDynamicPhases
 */
#ifndef OPM_BLACK_OIL_DYNAMIC_PHASES_HH
#define OPM_BLACK_OIL_DYNAMIC_PHASES_HH

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief The set of active phases, known only at run time.
 *
 * Read by the dynamic index classes (e.g. BlackOilDynTwoPhaseIndices). It is a
 * process-wide setting that must be set once, from the deck, before the first
 * primary variables or well structures are created. Flow does this in
 * Main::setupVanguard().
 */
class BlackOilDynamicPhases
{
public:
    static void set(const bool oil, const bool water, const bool gas)
    {
        oil_ = oil;
        water_ = water;
        gas_ = gas;
    }

    static bool oilActive() { return oil_; }
    static bool waterActive() { return water_; }
    static bool gasActive() { return gas_; }

private:
    // Default: oil-water until set.
    inline static bool oil_ = true;
    inline static bool water_ = true;
    inline static bool gas_ = false;
};

} // namespace Opm

#endif
