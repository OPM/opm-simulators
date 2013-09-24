// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Some templates to wrap the valgrind macros
 */
#ifndef OPM_VALGRIND_HH
#define OPM_VALGRIND_HH

#include <opm/core/utility/Unused.hpp>

#if ! HAVE_VALGRIND && ! defined(DOXYGEN)
namespace Valgrind
{
bool boolBlubb(bool value) { return value; }
void voidBlubb() { }

#define SetUndefined(t) voidBlubb()
#define SetDefined(t) voidBlubb()
#define CheckDefined(t) boolBlubb(true)
#define SetNoAccess(t) voidBlubb()
#define IsRunning() boolBlubb(false)
}

#else

#include <valgrind/memcheck.h>

namespace Valgrind
{
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
 * Opm::Valgrind::CheckDefined(i); // Valgrind complains!
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
inline bool CheckDefined(const T &value)
{
#if !defined NDEBUG && HAVE_VALGRIND
    unsigned int tmp = VALGRIND_CHECK_MEM_IS_DEFINED(&value, sizeof(T));
    return tmp == 0;
#else
    return true;
#endif
}
/*!
 * \ingroup Valgrind
 *
 *  * \brief Make valgrind complain if any of the the memory occupied
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
 * Opm::Valgrind::CheckDefined(i, 2); // Valgrind complains!
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
inline bool CheckDefined(const T *value, int size)
{
#if !defined NDEBUG && HAVE_VALGRIND
    unsigned int tmp = VALGRIND_CHECK_MEM_IS_DEFINED(value, size*sizeof(T));
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
 * Opm::Valgrind::SetUndefined(i);
 * Opm::Valgrind::CheckDefined(i); // Valgrind complains!
 * \endcode
 *
 * \tparam T The type of the object which ought to be set to undefined
 *
 * \param value The object which's memory valgrind should be told is undefined
 */
template <class T>
inline void SetUndefined(const T &value)
{
#if !defined NDEBUG && HAVE_VALGRIND
    auto OPM_UNUSED result = VALGRIND_MAKE_MEM_UNDEFINED(&value, sizeof(T));
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
 * Opm::Valgrind::SetUndefined(&i[1], 2);
 * Opm::Valgrind::CheckDefined(i, 3); // Valgrind complains!
 * \endcode
 *
 * \tparam T The type of the object which ought to be set to undefined
 *
 * \param value Pointer to the first object of the array.
 * \param size The size of the array in number of objects
 */
template <class T>
inline void SetUndefined(const T *value, int size)
{
#if !defined NDEBUG && HAVE_VALGRIND
    auto OPM_UNUSED result = VALGRIND_MAKE_MEM_UNDEFINED(value, size*sizeof(T));
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
 * Opm::Valgrind::SetDefined(i);
 * Opm::Valgrind::CheckDefined(i); // Valgrind does not complain!
 * \endcode
 *
 * \tparam T The type of the object which valgrind should consider as defined
 *
 * \param value The object which's memory valgrind should consider as defined
 */
template <class T>
inline void SetDefined(const T &value)
{
#if !defined NDEBUG && HAVE_VALGRIND
    auto OPM_UNUSED result = VALGRIND_MAKE_MEM_DEFINED(&value, sizeof(T));
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
 * Opm::Valgrind::SetDefined(i, 3);
 * Opm::Valgrind::CheckDefined(i, 3); // Valgrind does not complain!
 * \endcode
 *
 * \tparam T The type of the object which valgrind should consider as defined
 *
 * \param value Pointer to the first object of the array.
 * \param n The size of the array in number of objects
 */
template <class T>
inline void SetDefined(const T *value, int n)
{
#if !defined NDEBUG && HAVE_VALGRIND
    auto OPM_UNUSED result = VALGRIND_MAKE_MEM_DEFINED(value, n*sizeof(T));
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
 * Opm::Valgrind::SetNoAccess(i);
 * int j = i; // Valgrind complains!
 * \endcode
 *
 * \tparam T The type of the object which valgrind should complain if accessed
 *
 * \param value The object which's memory valgrind should complain if accessed
 */
template <class T>
inline void SetNoAccess(const T &value)
{
#if !defined NDEBUG && HAVE_VALGRIND
    auto OPM_UNUSED result = VALGRIND_MAKE_MEM_NOACCESS(&value, sizeof(T));
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
 * Opm::Valgrind::SetNoAccess(i, 2);
 * int j = i[1]; // Valgrind complains!
 * \endcode
 *
 * \param value Pointer to the first object of the array.
 * \param size The size of the array in number of objects
 *
 * \param value The object which's memory valgrind should complain if accessed
 */
template <class T>
inline void SetNoAccess(const T *value, int size)
{
#if !defined NDEBUG && HAVE_VALGRIND
    auto OPM_UNUSED result = VALGRIND_MAKE_MEM_NOACCESS(value, size*sizeof(T));
#endif
}

}

#endif

#endif
