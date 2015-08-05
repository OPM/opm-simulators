#ifndef EWOMS_DUNE_COMPATIBILITY_HH
#define EWOMS_DUNE_COMPATIBILITY_HH

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/misc/compatibility.hh>

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

    template< int codim >
    struct GridEntityAccess< Dune::cpgrid::Entity< codim > >
    {

      typedef Dune::cpgrid::Entity< codim >   EntityType;
      typedef EntityType                      GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return entity;
      }
    };

  } // namespace Fem

} // end namespace Dune

#endif // #if HAVE_DUNE_FEM

#endif // #ifndef EWOMS_DUNE_COMPATIBILITY_HH
