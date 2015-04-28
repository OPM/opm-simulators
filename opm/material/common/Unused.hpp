/*
  Copyright (C) 2010-2013 by Andreas Lauser

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

#endif
