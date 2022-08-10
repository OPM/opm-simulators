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
 * \brief This file ensures that ebos can be compiled in the presence of dune-fem
 *
 * It implements a few work-arounds for some incompatibilities with the Dune grid
 * interface of CpGrid. A better way would be to make CpGrid conforming.
 */
#ifndef EWOMS_FEM_CPGRID_COMPAT_HH
#define EWOMS_FEM_CPGRID_COMPAT_HH

#if HAVE_DUNE_FEM
#include <dune/common/version.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/compatibility.hh>
#include <dune/fem/io/streams/streams.hh>

namespace Dune {
namespace cpgrid {
template <int codim>
class Entity;

template <int codim>
class EntityPointer;
}

// specialization of dune-fem compatiblity functions for CpGrid, since CpGrid does not use the interface classes.
namespace Fem {

////////////////////////////////////////////////////////////
//
//  GridEntityAccess for CpGrid entities
//
////////////////////////////////////////////////////////////
template<int codim>
struct GridEntityAccess<Dune::cpgrid::Entity<codim> >
{
    using EntityType = Dune::cpgrid::Entity<codim>;
    using GridEntityType = EntityType;

    static const GridEntityType& gridEntity(const EntityType& entity)
    { return entity; }
};

} // namespace Fem

} // end namespace Dune

#endif // #if HAVE_DUNE_FEM

#endif // EWOMS_FEM_CPGRID_COMPAT_HH
