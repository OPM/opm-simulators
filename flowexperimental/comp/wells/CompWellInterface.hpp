/*
  Copyright 2024, SINTEF Digital

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

#ifndef OPM_COMP_WELLINTERFACE_HPP
#define OPM_COMP_WELLINTERFACE_HPP

#include <opm/models/utils/propertysystem.hh>

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <string>

namespace Opm
{

template <typename TypeTag> // TODO: do we need to use TypeTag here?
class CompWellInterface
{
public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    CompWellInterface(const Well& well,
                      const int index_of_well);

protected:

    const Well& well_ecl_;
    int index_of_well_{-1};

    Scalar reference_depth_ {}; // TODO: we might not need it since it in well_ecl_.

    // std::string name_;

};


} // end of namespace Opm

#include "CompWellInterface_impl.hpp"

#endif // OPM_COMP_WELLINTERFACE_HPP
