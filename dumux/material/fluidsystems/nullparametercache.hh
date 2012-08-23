// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2009-2011 by Bernd Flemisch                               *
 *   Copyright (C) 2010 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Melanie Darcis                                    *
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
 * \ingroup Fluidsystems
 *
 * \brief The a parameter cache which does nothing
 */
#ifndef DUMUX_NULL_PARAMETER_CACHE_HH
#define DUMUX_NULL_PARAMETER_CACHE_HH

#include "parametercachebase.hh"

namespace Dumux
{
/*!
 * \ingroup Fluidsystems
 * \brief The a parameter cache which does nothing
 */
class NullParameterCache : public ParameterCacheBase<NullParameterCache>
{
public:
    NullParameterCache()
    {}
};

} // end namepace

#endif
