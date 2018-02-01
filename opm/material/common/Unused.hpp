// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
/*!
 * \file
 * \brief Provides the OPM_UNUSED macro
 *
 * This macro can be used to mark variables as "potentially unused" which suppresses some
 * bogus compiler warnings. If the compiler does not support this, the macro is a no-op.
 */
#ifndef OPM_UNUSED_HH
#define OPM_UNUSED_HH

#ifndef HAS_ATTRIBUTE_UNUSED
#define OPM_UNUSED
#else
#define OPM_UNUSED __attribute__((unused))
#endif

#ifdef NDEBUG
#define OPM_DEBUG_UNUSED
#define OPM_OPTIM_UNUSED OPM_UNUSED
#else
#define OPM_DEBUG_UNUSED OPM_UNUSED
#define OPM_OPTIM_UNUSED
#endif

#endif
