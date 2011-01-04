/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Copyright (C) 2010 by Benjamin Faigle                                   *
 *   Copyright (C) 2010 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 *
 * \brief Provides defaults for all available components.
 */
#ifndef DUMUX_DEFAULT_COMPONENTS_HH
#define DUMUX_DEFAULT_COMPONENTS_HH

#include <dumux/common/propertysystem.hh>

#include <dumux/material/components/ch4.hh>
#include <dumux/material/components/simpleco2.hh>
#include <dumux/material/components/h2.hh>
#include <dumux/material/components/o2.hh>
#include <dumux/material/components/oil.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/brine.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dune/common/stdstreams.hh>

namespace Dumux
{
namespace Properties
{
//! defines the components which are being used by the fluid system by
//! default and how they are initialized
NEW_PROP_TAG(DefaultComponents);

//! defines the components which are actually being used by the fluid
//! system
NEW_PROP_TAG(Components);

NEW_PROP_TAG(Scalar);

SET_PROP_DEFAULT(DefaultComponents)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef Dumux::H2O<Scalar> H2O_IAPWS;

public:
    typedef Dumux::TabulatedComponent<Scalar, H2O_IAPWS> H2O;
    typedef Dumux::N2<Scalar> N2;
    typedef Dumux::O2<Scalar> O2;
    typedef Dumux::H2<Scalar> H2;
    typedef Dumux::CH4<Scalar> CH4;
    typedef Dumux::SimpleCO2<Scalar> SimpleCO2;
    typedef Dumux::SimpleH2O<Scalar> SimpleH2O;
    typedef Dumux::Brine<Scalar, H2O> Brine;

    static void init()
    {
        int nT = 100;
        int nP = 200;
        Dune::dinfo << "Initializing tables for the H2O fluid properties ("
                    << nT*nP
                    << " entries).\n";
        H2O::init(273.15, 623.15, nT, -10, 20e6, nP);;
    }
};

SET_PROP_DEFAULT(Components)
    : public GET_PROP(TypeTag, PTAG(DefaultComponents))
{};

}; // namespace Properties
}; // namespace Dumux

#endif
