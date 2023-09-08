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

#include <fmt/format.h>

namespace Opm::Pybind {

template <class TypeTag>
std::unique_ptr<double []>
PyMaterialState<TypeTag>::
getCellVolumes( std::size_t *size)
{
    Model &model = this->ebos_simulator_->model();
    *size = model.numGridDof();
    auto array = std::make_unique<double []>(*size);
    for (unsigned dof_idx = 0; dof_idx < *size; ++dof_idx) {
        array[dof_idx] = model.dofTotalVolume(dof_idx);
    }
    return array;
}

template <class TypeTag>
std::unique_ptr<double []>
PyMaterialState<TypeTag>::
getPorosity( std::size_t *size)
{
    Problem &problem = this->ebos_simulator_->problem();
    Model &model = this->ebos_simulator_->model();
    *size = model.numGridDof();
    auto array = std::make_unique<double []>(*size);
    for (unsigned dof_idx = 0; dof_idx < *size; ++dof_idx) {
        array[dof_idx] = problem.referencePorosity(dof_idx, /*timeIdx*/0);
    }
    return array;
}

template <class TypeTag>
void
PyMaterialState<TypeTag>::
setPorosity(const double *poro, std::size_t size)
{
    Problem &problem = this->ebos_simulator_->problem();
    Model &model = this->ebos_simulator_->model();
    auto model_size = model.numGridDof();
    if (model_size != size) {
        const std::string msg = fmt::format(
            "Cannot set porosity. Expected array of size: {}, got array of size: ",
            model_size, size);
        throw std::runtime_error(msg);
    }
    for (unsigned dof_idx = 0; dof_idx < size; ++dof_idx) {
        problem.setPorosity(poro[dof_idx], dof_idx);
    }
}
} //namespace Opm::Pybind
