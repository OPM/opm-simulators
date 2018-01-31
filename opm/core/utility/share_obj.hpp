/*
  Copyright (c) 2013 Uni Research AS

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

#ifndef OPM_SHARE_OBJ_HPP
#define OPM_SHARE_OBJ_HPP

#include <memory>  // shared_ptr

namespace Opm {

/// Custom deleter that does nothing
inline void no_delete (void const *) { }

/*!
 * Share pointer of a local object.
 *
 * Use this wrapper when an interface needs a shared_ptr, but you
 * want to pass an object that has local storage (and you know
 * that the shared_ptr client doesn't need it outside of the scope).
 *
 * \example
 * \code{.cpp}
 *  Foo obj;
 *  std::shared_ptr <Foo> ptr = share_obj (obj);
 * \endcode
 */
template <typename T> std::shared_ptr <T> share_obj (T& t) {
    return std::shared_ptr <T> (&t, no_delete);
}

} // namespace Opm

#endif /* OPM_SHARE_OBJ_HPP */
