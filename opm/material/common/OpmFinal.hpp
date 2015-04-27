/*
  Copyright (C) 2015 by Andreas Lauser

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
 * \brief This file provides a wrapper around the "final" C++-2011 statement.
 *
 * The "final" C++-2011 statement specifies that a virtual method cannot be overloaded by
 * derived classes anymore. This allows the compiler to de-virtualize calls to such
 * methods and is this an optimization. (it also prevents the programmer from creating a
 * new method instead of overloading an existing one, i.e., it is nice to have from a
 * code quality perspective.) Since not all compilers which must be supported by OPM
 * feature the "final" qualifier, this method provides a wrapper macro around it.
 */
#ifndef OPM_FINAL_HPP
#define OPM_FINAL_HPP

#if HAVE_FINAL
#define OPM_FINAL final
#else
#define OPM_FINAL /* nothing; compiler does not recognize "final" */
#endif

#endif
