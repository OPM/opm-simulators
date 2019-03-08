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
 * \copydoc Opm::BlackOilDefaultIndexTraits
 */
#ifndef OPM_BLACK_OIL_DEFAULT_INDEX_TRAITS_HPP
#define OPM_BLACK_OIL_DEFAULT_INDEX_TRAITS_HPP

namespace Opm {

/*!
 * \brief The class which specifies the default phase and component indices for the
 *        black-oil fluid system.
 */
class BlackOilDefaultIndexTraits
{
public:
    //! Index of the water phase
    static const unsigned waterPhaseIdx = 0;
    //! Index of the oil phase
    static const unsigned oilPhaseIdx = 1;
    //! Index of the gas phase
    static const unsigned gasPhaseIdx = 2;

    //! Index of the oil component
    static const unsigned oilCompIdx = 0;
    //! Index of the water component
    static const unsigned waterCompIdx = 1;
    //! Index of the gas component
    static const unsigned gasCompIdx = 2;
};

} // namespace Opm

#endif
