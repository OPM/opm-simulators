/*
  Copyright 2020 Equinor ASA.

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

#ifndef OPM_PY_MATERIAL_STATE_HEADER_INCLUDED
#define OPM_PY_MATERIAL_STATE_HEADER_INCLUDED

#include <opm/models/utils/propertysystem.hh>

#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Opm::Pybind
{
    template <class TypeTag>
    class PyMaterialState {
        using Simulator = GetPropType<TypeTag, Opm::Properties::Simulator>;
        using Problem = GetPropType<TypeTag, Opm::Properties::Problem>;
        using Model = GetPropType<TypeTag, Opm::Properties::Model>;
        using ElementContext = GetPropType<TypeTag, Opm::Properties::ElementContext>;
        using FluidSystem = GetPropType<TypeTag, Opm::Properties::FluidSystem>;
        using Indices = GetPropType<TypeTag, Opm::Properties::Indices>;
        using GridView = GetPropType<TypeTag, Opm::Properties::GridView>;

    public:
        PyMaterialState(Simulator *ebos_simulator)
            : ebos_simulator_(ebos_simulator) { }

        std::unique_ptr<double []> getCellVolumes( std::size_t *size);
        std::unique_ptr<double []> getPorosity( std::size_t *size);
        void setPorosity(const double *poro, std::size_t size);
    private:
        Simulator *ebos_simulator_;
    };

}
#include "PyMaterialState_impl.hpp"

#endif // OPM_PY_MATERIAL_STATE_HEADER_INCLUDED
