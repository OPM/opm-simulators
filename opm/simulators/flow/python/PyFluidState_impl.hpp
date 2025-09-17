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
#ifndef OPM_PY_FLUID_STATE_IMPL_HEADER_INCLUDED
#define OPM_PY_FLUID_STATE_IMPL_HEADER_INCLUDED

// Improve IDE experience
#ifndef OPM_PY_FLUID_STATE_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/flow/python/PyFluidState.hpp>
#endif

#include <opm/material/common/MathToolbox.hpp>

#include <stdexcept>

#include <fmt/format.h>

namespace Opm::Pybind {

template <class TypeTag>
PyFluidState<TypeTag>::
PyFluidState(Simulator* simulator)
    : simulator_(simulator)
{}

// Public methods alphabetically sorted
// ------------------------------------

template <class TypeTag>
std::vector<int>
PyFluidState<TypeTag>::
getPrimaryVarMeaning(const std::string& variable) const
{
    Model& model = this->simulator_->model();
    auto& sol = model.solution(/*timeIdx*/0);
    const auto size = model.numGridDof();
    std::vector<int> array(size);
    for (unsigned dof_idx = 0; dof_idx < size; ++dof_idx) {
        auto primary_vars = sol[dof_idx];
        array[dof_idx] = getVariableMeaning_(primary_vars, variable);
    }
    return array;
}

template <class TypeTag>
std::map<std::string, int>
PyFluidState<TypeTag>::
getPrimaryVarMeaningMap(const std::string& variable) const
{
    if (variable.compare("pressure") == 0) {
        return {{ "Po", static_cast<int>(PrimaryVariables::PressureMeaning::Po) },
                { "Pw", static_cast<int>(PrimaryVariables::PressureMeaning::Pw) },
                { "Pg", static_cast<int>(PrimaryVariables::PressureMeaning::Pg) }};
    }
    else if (variable.compare("water") == 0) {
        return {{ "Sw", static_cast<int>(PrimaryVariables::WaterMeaning::Sw) },
                { "Rvw", static_cast<int>(PrimaryVariables::WaterMeaning::Rvw) },
                { "Rsw", static_cast<int>(PrimaryVariables::WaterMeaning::Rsw) },
                { "Disabled", static_cast<int>(PrimaryVariables::WaterMeaning::Disabled) }};
    }
    else if (variable.compare("gas") == 0) {
        return {{ "Sg", static_cast<int>(PrimaryVariables::GasMeaning::Sg) },
                { "Rs", static_cast<int>(PrimaryVariables::GasMeaning::Rs) },
                { "Rv", static_cast<int>(PrimaryVariables::GasMeaning::Rv) },
                { "Disabled", static_cast<int>(PrimaryVariables::GasMeaning::Disabled) }};
    }
    else if (variable.compare("brine") == 0) {
        return {{ "Cs", static_cast<int>(PrimaryVariables::BrineMeaning::Cs) },
                { "Sp", static_cast<int>(PrimaryVariables::BrineMeaning::Sp) },
                { "Disabled", static_cast<int>(PrimaryVariables::BrineMeaning::Disabled) }};
    }
    else {
        const std::string msg = fmt::format(
            "Unknown variable meaning '{}': Expected pressure, water, gas, or brine", variable);
        throw std::runtime_error(msg);
    }
}

/* Meaning of the primary variables: Sw, Sg, po, pg, Rs, Rv
 * 1. Sw_po_Sg -> threephase case
 * 2. Sw_po_Rs -> water + oil case
 * 3. Sw_pg_Rv -> water + gas case
 */

/* Variables:
   Sw = Water saturation,
   So = Oil saturation,
   Sg = Gas saturation,
   pw = Water pressure,
   po = Oil pressure,
   pg = Gas pressure,
   Rs = The solution gas oil ratio: The amount of gas dissolved in the oil
   Rv = The oil vaporization factor of the gas phase
   invB = The inverse formation volume factor of a fluid phase
   rho_w = Water density,
   rho_o = Oil density,
   rho_g = Gas density,
   mu_w = Water viscosity,
   mu_o = Oil viscosity,
   mu_g = Gas viscosity,
   kr_w = Water relperm,
   kr_o = Oil relperm,
   kr_g = Gas relperm,
 */
template <class TypeTag>
std::vector<double>
PyFluidState<TypeTag>::
getFluidStateVariable(const std::string& name) const
{
    Model& model = this->simulator_->model();
    const auto size = model.numGridDof();
    std::vector<double> array(size);
    const auto& grid_view = this->simulator_->vanguard().gridView();
    /* NOTE: grid_view.size(0) should give the same value as
     *  model.numGridDof()
     */
    ElementContext elem_ctx(*this->simulator_);
    auto var_type = getVariableType_(name);
    for (const auto& elem : elements(grid_view, Dune::Partitions::interior)) {
        elem_ctx.updatePrimaryStencil(elem);
        elem_ctx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        for (unsigned dof_idx = 0; dof_idx < elem_ctx.numPrimaryDof(/*timeIdx=*/0); ++dof_idx) {
            const auto& int_quants = elem_ctx.intensiveQuantities(dof_idx, /*timeIdx=*/0);
            const auto& fs = int_quants.fluidState();
            unsigned global_dof_idx = elem_ctx.globalSpaceIndex(dof_idx, /*timeIdx=*/0);
            array[global_dof_idx] = getVariableValue_(fs, var_type, name);
        }
    }
    return array;
}

template <class TypeTag>
std::vector<double>
PyFluidState<TypeTag>::
getPrimaryVariable(const std::string& idx_name) const
{
    std::size_t primary_var_idx = getPrimaryVarIndex_(idx_name);
    Model& model = this->simulator_->model();
    auto& sol = model.solution(/*timeIdx*/0);
    const auto size = model.numGridDof();
    std::vector<double> array(size);
    for (unsigned dof_idx = 0; dof_idx < size; ++dof_idx) {
        auto primary_vars = sol[dof_idx];
        array[dof_idx] = primary_vars[primary_var_idx];
    }
    return array;
}

template <class TypeTag>
void
PyFluidState<TypeTag>::
setPrimaryVariable(const std::string& idx_name,
                   const double* data,
                   std::size_t size)
{
    const std::size_t primary_var_idx = getPrimaryVarIndex_(idx_name);
    Model& model = this->simulator_->model();
    auto& sol = model.solution(/*timeIdx*/0);
    const auto model_size = model.numGridDof();
    if (model_size != size) {
        const std::string msg = fmt::format(
            "Cannot set primary variable. Expected array of size {} but got array of size: {}",
            model_size, size);
        throw std::runtime_error(msg);
    }
    for (unsigned dof_idx = 0; dof_idx < size; ++dof_idx) {
        auto& primary_vars = sol[dof_idx];
        primary_vars[primary_var_idx] = data[dof_idx];
    }
}

// Private methods alphabetically sorted
// -------------------------------------

template <class TypeTag>
std::size_t
PyFluidState<TypeTag>::
getPrimaryVarIndex_(const std::string& idx_name) const
{
    if (idx_name.compare("pressure") == 0) {
        return Indices::pressureSwitchIdx;
    }
    else if (idx_name.compare("water_saturation") == 0) {
        return Indices::waterSwitchIdx;
    }
    else if (idx_name.compare("composition") == 0) {
        return Indices::compositionSwitchIdx;
    }
    else {
        const std::string msg = fmt::format("Unknown primary variable index name: {}", idx_name);
        throw std::runtime_error(msg);
    }
}

template <class TypeTag>
int
PyFluidState<TypeTag>::
getVariableMeaning_(PrimaryVariables& primary_vars,
                    const std::string& variable) const
{
    if (variable.compare("pressure") == 0) {
        return static_cast<int>(primary_vars.primaryVarsMeaningPressure());
    }
    else if(variable.compare("water") == 0) {
        return static_cast<int>(primary_vars.primaryVarsMeaningWater());
    }
    else if (variable.compare("gas") == 0) {
        return static_cast<int>(primary_vars.primaryVarsMeaningGas());
    }
    else if (variable.compare("brine") == 0) {
        return static_cast<int>(primary_vars.primaryVarsMeaningBrine());
    }
    else {
        const std::string msg = fmt::format(
            "Unknown variable meaning '{}': Expected pressure, water, gas, or brine", variable);
        throw std::runtime_error(msg);
    }
}

template <class TypeTag>
typename PyFluidState<TypeTag>::VariableType
PyFluidState<TypeTag>::
getVariableType_(const std::string& name) const
{
    static std::map<std::string, VariableType> variable_type_map{
       {"Sw", VariableType::Sw},
       {"Sg", VariableType::Sg},
       {"So", VariableType::So},
       {"pw", VariableType::pw},
       {"pg", VariableType::pg},
       {"po", VariableType::po},
       {"Rs", VariableType::Rs},
       {"Rv", VariableType::Rv},
       {"rho_w", VariableType::rho_w},
       {"rho_g", VariableType::rho_g},
       {"rho_o", VariableType::rho_o},
       {"T", VariableType::T}
    };

    if (variable_type_map.count(name) == 0) {
        variableNotFoundError_(name);
    }
    return variable_type_map.at(name);
}

template <class TypeTag>
template <class FluidState>
double
PyFluidState<TypeTag>::
getVariableValue_(FluidState& fs,
                  VariableType var_type,
                  const std::string& name) const
{
    switch(var_type) {
    case VariableType::pw :
        return getValue(fs.pressure(FluidSystem::waterPhaseIdx));
    case VariableType::pg :
        return getValue(fs.pressure(FluidSystem::gasPhaseIdx));
    case VariableType::po :
        return getValue(fs.pressure(FluidSystem::oilPhaseIdx));
    case VariableType::rho_w :
        return getValue(fs.density(FluidSystem::waterPhaseIdx));
    case VariableType::rho_g :
        return getValue(fs.density(FluidSystem::gasPhaseIdx));
    case VariableType::rho_o :
        return getValue(fs.density(FluidSystem::oilPhaseIdx));
    case VariableType::Rs :
        return getValue(fs.Rs());
    case VariableType::Rv :
        return getValue(fs.Rv());
    case VariableType::Sw :
        return getValue(fs.saturation(FluidSystem::waterPhaseIdx));
    case VariableType::Sg :
        return getValue(fs.saturation(FluidSystem::gasPhaseIdx));
    case VariableType::So :
        return getValue(fs.saturation(FluidSystem::oilPhaseIdx));
    case VariableType::T :
        return getValue(fs.temperature(FluidSystem::waterPhaseIdx));
    default:
        variableNotFoundError_(name);
        return 0.0; // unreachable, shut up compiler
    }
}

template <class TypeTag>
void
PyFluidState<TypeTag>::
variableNotFoundError_(const std::string& name) const
{
    const std::string msg = fmt::format("Access to variable '{}' is not implemented yet!", name);
    throw std::runtime_error(msg);
}

} // namespace Opm::Pybind

#endif // OPM_PY_FLUID_STATE_IMPL_HEADER_INCLUDED
