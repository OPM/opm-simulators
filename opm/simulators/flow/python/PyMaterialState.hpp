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

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <cstddef>

namespace Opm::Pybind
{
    template <class TypeTag>
    class PyMaterialState {
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Problem = GetPropType<TypeTag, Properties::Problem>;
        using Model = GetPropType<TypeTag, Properties::Model>;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;
        using GridView = GetPropType<TypeTag, Properties::GridView>;

    public:
        PyMaterialState(Simulator* simulator)
            : simulator_(simulator)
        {}

        std::vector<double> getCellVolumes();
        std::vector<double> getPorosity();
        void setPorosity(const double* poro, std::size_t size);

    private:
        Simulator* simulator_;
    };

}
#include "PyMaterialState_impl.hpp"

#endif // OPM_PY_MATERIAL_STATE_HEADER_INCLUDED
