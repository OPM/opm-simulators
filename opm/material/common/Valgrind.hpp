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
 * \brief Some templates to wrap the valgrind client request macros
 */
#ifndef OPM_VALGRIND_HPP
#define OPM_VALGRIND_HPP

#include <opm/material/common/Unused.hpp>

#if HAVE_VALGRIND
#include <valgrind/memcheck.h>
#endif

#if HAVE_VALGRIND
#define OPM_VALGRIND_OPTIM_UNUSED OPM_OPTIM_UNUSED
#else
#define OPM_VALGRIND_OPTIM_UNUSED OPM_UNUSED
#endif

namespace Opm {
namespace Valgrind {
/*!
 * \ingroup Valgrind
 * \brief Returns whether the program is running under Valgrind or not.
 */
inline bool IsRunning()
{
#if !defined NDEBUG && HAVE_VALGRIND
    return RUNNING_ON_VALGRIND;
#else
    return false;
#endif
}

/*!
 * \ingroup Valgrind
 * \brief Make valgrind complain if any of the memory occupied by an object
 *        is undefined.
 *
 * Please note that this does not check whether the destinations of an
 * object's pointers or references are defined. Also, for performance
 * reasons the compiler might insert "padding bytes" between within
 * the objects which leads to false positives.
 *
 * Example:
 *
 * \code
 * int i;
 * Valgrind::CheckDefined(i); // Valgrind complains!
 * \endcode
 *
 * \tparam T The type of the object which ought to be checked
 *
 * \param value the object which valgrind should check
 *
 * \return true iff there are no undefined bytes in the memory
 *         occupied by the object.
 */
template <class T>
inline bool CheckDefined(const T& value OPM_VALGRIND_OPTIM_UNUSED)
{
#if !defined NDEBUG && HAVE_VALGRIND
    auto tmp = VALGRIND_CHECK_MEM_IS_DEFINED(&value, sizeof(T));
    return tmp == 0;
#else
    return true;
#endif
}



/*!
 * \ingroup Valgrind
 * \brief Make valgrind complain if any of the memory occupied by an object
 *        is not addressable.
 *
 * Example:
 *
 * \code
 * int* i = NULL;
 * Valgrind::CheckAddressable(*i); // Valgrind complains!
 * \endcode
 *
 * \tparam T The type of the object which ought to be checked
 *
 * \param value the object which valgrind should check
 *
 * \return true iff there are no unadressable bytes in the memory
 *         occupied by the object.
 */
template <class T>
inline bool CheckAddressable(const T& value OPM_VALGRIND_OPTIM_UNUSED)
{
#if !defined NDEBUG && HAVE_VALGRIND
    auto tmp = VALGRIND_CHECK_MEM_IS_ADDRESSABLE(&value, sizeof(T));
    return tmp == 0;
#else
    return true;
#endif
}

/*!
 * \ingroup Valgrind
 * \brief Make valgrind complain if any of the the memory occupied
 *        by a C-style array objects is undefined.
 *
 * Please note that this does not check whether the destinations of an
 * object's pointers or references are defined. Also, for performance
 * reasons the compiler might insert "padding bytes" between within
 * the objects which leads to false positives.
 *
 * Example:
 *
 * \code
 * int i[2];
 * Valgrind::CheckDefined(i, 2); // Valgrind complains!
 * \endcode
 *
 * \tparam T The type of the object which ought to be checked
 *
 * \param value Pointer to the first object of the array.
 * \param size The size of the array in number of objects
 *
 * \return true iff there are no undefined bytes in the memory
 *         occupied by the array.
 */
template <class T>
inline bool CheckDefined(const T* value OPM_VALGRIND_OPTIM_UNUSED, int size OPM_VALGRIND_OPTIM_UNUSED)
{
#if !defined NDEBUG && HAVE_VALGRIND
    auto tmp = VALGRIND_CHECK_MEM_IS_DEFINED(value, size*sizeof(T));
    return tmp == 0;
#else
    return true;
#endif
}

/*!
 * \ingroup Valgrind
 * \brief Make the memory on which an object resides undefined in
 *        valgrind runs.
 *
 * Example:
 *
 * \code
 * int i = 0;
 * Valgrind::SetUndefined(i);
 * Valgrind::CheckDefined(i); // Valgrind complains!
 * \endcode
 *
 * \tparam T The type of the object which ought to be set to undefined
 *
 * \param value The object which's memory valgrind should be told is undefined
 */
template <class T>
inline void SetUndefined(const T &value OPM_VALGRIND_OPTIM_UNUSED)
{
#if !defined NDEBUG && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_UNDEFINED(&value, sizeof(T));
#endif
}

/*!
 * \ingroup Valgrind
 * \brief Make the memory on which an array of object resides
 *        undefined in valgrind runs.
 *
 * Example:
 *
 * \code
 * int i[3] = {0, 1, 3};
 * Valgrind::SetUndefined(&i[1], 2);
 * Valgrind::CheckDefined(i, 3); // Valgrind complains!
 * \endcode
 *
 * \tparam T The type of the object which ought to be set to undefined
 *
 * \param value Pointer to the first object of the array.
 * \param size The size of the array in number of objects
 */
template <class T>
inline void SetUndefined(const T* value OPM_VALGRIND_OPTIM_UNUSED, int size OPM_VALGRIND_OPTIM_UNUSED)
{
#if !defined NDEBUG && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_UNDEFINED(value, size*sizeof(T));
#endif
}

/*!
 * \ingroup Valgrind
 * \brief Make the memory on which an object resides defined.
 *
 * Example:
 *
 * \code
 * int i;
 * Valgrind::SetDefined(i);
 * Valgrind::CheckDefined(i); // Valgrind does not complain!
 * \endcode
 *
 * \tparam T The type of the object which valgrind should consider as defined
 *
 * \param value The object which's memory valgrind should consider as defined
 */
template <class T>
inline void SetDefined(const T& value OPM_VALGRIND_OPTIM_UNUSED)
{
#if !defined NDEBUG && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_DEFINED(&value, sizeof(T));
#endif
}

/*!
 * \ingroup Valgrind
 * \brief Make the memory on which a C-style array of objects resides
 *        defined.
 *
 * Example:
 *
 * \code
 * int i[3];
 * Valgrind::SetDefined(i, 3);
 * Valgrind::CheckDefined(i, 3); // Valgrind does not complain!
 * \endcode
 *
 * \tparam T The type of the object which valgrind should consider as defined
 *
 * \param value Pointer to the first object of the array.
 * \param n The size of the array in number of objects
 */
template <class T>
inline void SetDefined(const T *value OPM_VALGRIND_OPTIM_UNUSED, int n OPM_VALGRIND_OPTIM_UNUSED)
{
#if !defined NDEBUG && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_DEFINED(value, n*sizeof(T));
#endif
}

/*!
 * \ingroup Valgrind
 * \brief Make valgrind complain if an object's memory is accessed.
 *
 * Example:
 *
 * \code
 * int i = 1;
 * Valgrind::SetNoAccess(i);
 * int j = i; // Valgrind complains!
 * \endcode
 *
 * \tparam T The type of the object which valgrind should complain if accessed
 *
 * \param value The object which's memory valgrind should complain if accessed
 */
template <class T>
inline void SetNoAccess(const T &value OPM_VALGRIND_OPTIM_UNUSED)
{
#if !defined NDEBUG && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_NOACCESS(&value, sizeof(T));
#endif
}

/*!
 * \ingroup Valgrind
 * \brief Make valgrind complain if the memory of a C-style array of
 *        objects is accessed.
 *
 * Example:
 *
 * \code
 * int i[3] = {0, 1, 2};
 * Valgrind::SetNoAccess(i, 2);
 * int j = i[1]; // Valgrind complains!
 * \endcode
 *
 * \param value Pointer to the first object of the array.
 * \param size The size of the array in number of objects
 */
template <class T>
inline void SetNoAccess(const T *value OPM_VALGRIND_OPTIM_UNUSED, int size OPM_VALGRIND_OPTIM_UNUSED)
{
#if !defined NDEBUG && HAVE_VALGRIND
    VALGRIND_MAKE_MEM_NOACCESS(value, size*sizeof(T));
#endif
}

}} // namespace Valgrind, Opm

#endif
