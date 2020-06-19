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
#ifndef EWOMS_COPYRESTRICTPROLONG_HH
#define EWOMS_COPYRESTRICTPROLONG_HH

#include <opm/material/common/Unused.hpp>

#if HAVE_DUNE_FEM
#include <dune/fem/space/common/restrictprolonginterface.hh>
#endif

namespace Opm
{
    template < class Grid, class Container >
    class CopyRestrictProlong;

    template < class Grid, class Container >
    struct CopyRestrictProlongTraits
    {
      using DomainFieldType = typename Grid::ctype;
      using RestProlImp = CopyRestrictProlong< Grid, Container > ;
    };

    template< class Grid, class Container >
    class CopyRestrictProlong
#if HAVE_DUNE_FEM
    : public Dune::Fem::RestrictProlongInterfaceDefault< CopyRestrictProlongTraits< Grid, Container > >
#endif
    {
      using ThisType = CopyRestrictProlong< Grid, Container >;

      Container& container_;
    public:
      using DomainFieldType = typename Grid::ctype;

      explicit CopyRestrictProlong( Container& container )
        : container_( container )
      {}

      /** \brief explicit set volume ratio of son and father
       *
       *  \param[in]  weight  volume of son / volume of father
       *
       *  \note If this ratio is set, it is assume to be constant.
       */
      template <class Field>
      void setFatherChildWeight(const Field& weight OPM_UNUSED) const
      {}

      //! restrict data to father
      template< class Entity >
      void restrictLocal ( const Entity& father, const Entity& son, bool initialize ) const
      {
        container_.resize();
        assert( container_.codimension() == 0 );
        if( initialize )
        {
          // copy values from son to father (only once)
          container_[ father ] = container_[ son ];
        }
      }

      //! restrict data to father
      template< class Entity, class LocalGeometry >
      void restrictLocal(const Entity& father,
                         const Entity& son,
                         const LocalGeometry& geometryInFather OPM_UNUSED,
                         bool initialize) const
      { restrictLocal(father, son, initialize); }

      //! prolong data to children
      template< class Entity >
      void prolongLocal(const Entity& father,
                        const Entity& son,
                        bool initialize OPM_UNUSED) const
      {
        container_.resize();
        assert( container_.codimension() == 0 );
        // copy values from father to son all sons
        container_[ son ] = container_[ father ];
      }

      //! prolong data to children
      template< class Entity, class LocalGeometry >
      void prolongLocal(const Entity& father,
                        const Entity& son,
                        const LocalGeometry& geometryInFather OPM_UNUSED,
                        bool initialize) const
      { prolongLocal(father, son, initialize); }

      /** \brief add discrete function to communicator
       *  \param[in]  comm  Communicator to add the discrete functions to
       */
      template< class Communicator >
      void addToList(Communicator& comm OPM_UNUSED)
      {
        // TODO
      }

      /** \brief add discrete function to load balancer
       *  \param[in]  lb LoadBalancer to add the discrete functions to
       */
      template< class LoadBalancer >
      void addToLoadBalancer(LoadBalancer& lb OPM_UNUSED)
      {
        // TODO
      }

    };


    class EmptyRestrictProlong;

    struct EmptyRestrictProlongTraits
    {
      using DomainFieldType = double               ;
      using RestProlImp = EmptyRestrictProlong ;
    };

    class EmptyRestrictProlong
#if HAVE_DUNE_FEM
    : public Dune::Fem::RestrictProlongInterfaceDefault< EmptyRestrictProlongTraits >
#endif
    {
      using ThisType = EmptyRestrictProlong;

    public:
      /** \brief explicit set volume ratio of son and father
       *
       *  \param[in]  weight  volume of son / volume of father
       *
       *  \note If this ratio is set, it is assume to be constant.
       */
      template <class Field>
      void setFatherChildWeight(const Field& weight OPM_UNUSED) const
      { }

      //! restrict data to father
      template< class Entity >
      void restrictLocal(const Entity& father OPM_UNUSED,
                         const Entity& son OPM_UNUSED,
                         bool initialize OPM_UNUSED) const
      { }

      //! restrict data to father
      template< class Entity, class LocalGeometry >
      void restrictLocal(const Entity& father OPM_UNUSED,
                         const Entity& son OPM_UNUSED,
                         const LocalGeometry& geometryInFather OPM_UNUSED,
                         bool initialize OPM_UNUSED) const
      { }

      //! prolong data to children
      template< class Entity >
      void prolongLocal(const Entity& father OPM_UNUSED,
                        const Entity& son OPM_UNUSED,
                        bool initialize OPM_UNUSED) const
      { }

      //! prolong data to children
      template< class Entity, class LocalGeometry >
      void prolongLocal(const Entity& father OPM_UNUSED,
                        const Entity& son OPM_UNUSED,
                        const LocalGeometry& geometryInFather OPM_UNUSED,
                        bool initialize OPM_UNUSED) const
      { }

      /** \brief add discrete function to communicator
       *  \param[in]  comm  Communicator to add the discrete functions to
       */
      template< class Communicator >
      void addToList(Communicator& comm OPM_UNUSED)
      { }

      /** \brief add discrete function to load balancer
       *  \param[in]  lb LoadBalancer to add the discrete functions to
       */
      template< class LoadBalancer >
      void addToLoadBalancer(LoadBalancer& lb OPM_UNUSED)
      { }
    };

} // namespace Opm

#endif
