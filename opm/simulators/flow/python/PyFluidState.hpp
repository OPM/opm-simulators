/*
  Copyright 2023 Equinor ASA.

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

#ifndef OPM_PY_FLUID_STATE_HEADER_INCLUDED
#define OPM_PY_FLUID_STATE_HEADER_INCLUDED

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
    class PyFluidState {
        using Simulator = GetPropType<TypeTag, Opm::Properties::Simulator>;
        using Problem = GetPropType<TypeTag, Opm::Properties::Problem>;
        using Model = GetPropType<TypeTag, Opm::Properties::Model>;
        using ElementContext = GetPropType<TypeTag, Opm::Properties::ElementContext>;
        using FluidSystem = GetPropType<TypeTag, Opm::Properties::FluidSystem>;
        using Indices = GetPropType<TypeTag, Opm::Properties::Indices>;
        using GridView = GetPropType<TypeTag, Opm::Properties::GridView>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

        enum class VariableType {
            // Primary variables: Sw, Sg, po, pg, Rs, Rv
            Sw, Sg, So, pw, pg, po, Rs, Rv, rho_w, rho_g, rho_o, T
        };
    public:
        PyFluidState(Simulator *ebos_simulator);
        std::unique_ptr<double []> getFluidStateVariable(
                                   const std::string &name, std::size_t *size ) const;
        std::unique_ptr<int []> getPrimaryVarMeaning(
                                   const std::string &variable, std::size_t *size) const;
        std::map<std::string, int> getPrimaryVarMeaningMap(const std::string &variable) const;
        std::unique_ptr<double []> getPrimaryVariable(
                                   const std::string &idx_name, std::size_t *size ) const;
        void setPrimaryVariable(const std::string &idx_name, const double *data, std::size_t size);

    private:
        std::size_t getPrimaryVarIndex_(const std::string &idx_name) const;
        int getVariableMeaning_(PrimaryVariables &primary_vars, const std::string &variable) const;
        VariableType getVariableType_(const std::string &name) const;
        template <class FluidState> double getVariableValue_(
            FluidState &fs, VariableType var_type, const std::string &name) const;
        void variableNotFoundError_(const std::string &name) const;

        Simulator *ebos_simulator_;
    };
}
#include "PyFluidState_impl.hpp"

#endif // OPM_PY_FLUID_STATE_HEADER_INCLUDED

