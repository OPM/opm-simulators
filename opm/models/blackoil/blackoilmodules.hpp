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
 * \brief Contains classes extending the black-oil model.
 * \detail This file holds dummy definitions,
 *         actual implementation is in separate headers.
 */
#ifndef OPM_BLACK_OIL_MODULES_HPP
#define OPM_BLACK_OIL_MODULES_HPP

namespace Opm {

#define DECLARE_MODULE(T) \
    template<class TypeTag, bool enable> class T##Module; \
    template<class TypeTag, bool enable> class T##IntensiveQuantities;  \
    template<class TypeTag, bool enable> class T##ExtensiveQuantities; \
    template<class TypeTag> struct T##Params; \
    template<class TypeTag> class T##IntensiveQuantities<TypeTag, false> {}; \
    template<class TypeTag> class T##ExtensiveQuantities<TypeTag, false> {};

DECLARE_MODULE(BlackOilBioeffects)
DECLARE_MODULE(BlackOilBrine)

#undef DECLARE_MODULE

} // namespace Opm

#endif
