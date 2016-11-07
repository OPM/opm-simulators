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
#ifndef EWOMS_DUNE_COMPATIBILITY_HH
#define EWOMS_DUNE_COMPATIBILITY_HH

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/compatibility.hh>
#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  namespace cpgrid
  {
    template <int codim>
    class Entity;

    template <int codim>
    class EntityPointer;

  }

  // specialization of dune-fem compatiblity functions for CpGrid, since CpGrid does not use the interface classes.
  namespace Fem
  {

    ////////////////////////////////////////////////////////////
    //
    //  make_entity for CpGrid entities
    //
    ////////////////////////////////////////////////////////////
    template <int codim>
    inline Dune::cpgrid::Entity< codim > make_entity ( const Dune::cpgrid::EntityPointer< codim >& entityPointer )
    {
      return *entityPointer;
    }

    template <int codim>
    inline Dune::cpgrid::Entity<codim> make_entity ( Dune::cpgrid::Entity<codim> entity )
    {
      return std::move( entity );
    }

    ////////////////////////////////////////////////////////////
    //
    //  GridEntityAccess for CpGrid entities
    //
    ////////////////////////////////////////////////////////////
    template< int codim >
    struct GridEntityAccess< Dune::cpgrid::Entity< codim > >
    {

      typedef Dune::cpgrid::Entity< codim >   EntityType;
      typedef EntityType                      GridEntityType;

      static const GridEntityType& gridEntity ( const EntityType& entity )
      {
        return entity;
      }
    };

    ////////////////////////////////////////////////////////////
    //
    //  operator << and operator >> for __float128
    //
    ////////////////////////////////////////////////////////////
#if HAVE_QUAD
    template< class Traits >
    inline OutStreamInterface< Traits > &
      operator<< ( OutStreamInterface< Traits >& out,
                   const __float128 value )
    {
      double val = double( value );
      out.writeDouble( val );
      return out;
    }

    template< class Traits >
    inline InStreamInterface< Traits > &
      operator>> ( InStreamInterface< Traits >& in,
                   __float128& value )
    {
      double val;
      in.readDouble( val );
      value = val;
      return in;
    }
#endif

  } // namespace Fem

} // end namespace Dune

#endif // #if HAVE_DUNE_FEM

#endif // #ifndef EWOMS_DUNE_COMPATIBILITY_HH
