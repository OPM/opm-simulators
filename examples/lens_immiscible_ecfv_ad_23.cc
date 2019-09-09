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
 * \brief Two-phase test for the immiscible model which uses the element-centered finite
 *        volume discretization in conjunction with automatic differentiation
 */
#include "config.h"

#include "lens_immiscible_ecfv_ad.hh"

#include <dune/grid/geometrygrid.hh>
#include <dune/grid/io/file/dgfparser/dgfgeogrid.hh>

BEGIN_PROPERTIES

// Use Dune-grid's GeometryGrid< YaspGrid >
SET_PROP(LensProblemEcfvAd, Grid )
{
  template< class ctype, unsigned int dim, unsigned int dimworld >
  class IdentityCoordFct
    : public Dune::AnalyticalCoordFunction
      < ctype, dim, dimworld, IdentityCoordFct< ctype, dim, dimworld > >
  {
    typedef IdentityCoordFct< ctype, dim, dimworld > This;
    typedef Dune::AnalyticalCoordFunction< ctype, dim, dimworld, This > Base;

  public:
    typedef typename Base :: DomainVector DomainVector;
    typedef typename Base :: RangeVector  RangeVector;

    template< typename... Args >
    IdentityCoordFct( Args&... )
    {}

    RangeVector operator()(const DomainVector& x) const
    {
      RangeVector y;
      evaluate( x, y );
      return y;
    }

    void evaluate( const DomainVector &x, RangeVector &y ) const
    {
      y = 0;
      for( unsigned int i = 0; i<dim; ++i )
        y[ i ] = x[ i ];
    }

  };

  typedef Dune::YaspGrid< 2 > MyYaspGrid;

public:
  //typedef MyYaspGrid type;
  typedef Dune::GeometryGrid< MyYaspGrid,
                              IdentityCoordFct< typename MyYaspGrid::ctype,
                                                MyYaspGrid::dimension,
                                                MyYaspGrid::dimensionworld+1> >  type;
};

END_PROPERTIES

#include <ewoms/common/start.hh>

int main(int argc, char **argv)
{
    typedef TTAG(LensProblemEcfvAd) ProblemTypeTag;
    return Opm::start<ProblemTypeTag>(argc, argv);
}
