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
#ifndef OPM_PY_MATERIAL_STATE_IMPL_HEADER_INCLUDED
#define OPM_PY_MATERIAL_STATE_IMPL_HEADER_INCLUDED

// Improve IDE experience
#ifndef OPM_PY_MATERIAL_STATE_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/flow/python/PyMaterialState.hpp>
#endif

#include <fmt/format.h>

namespace Opm::Pybind {

template <class TypeTag>
std::vector<double>
PyMaterialState<TypeTag>::
getCellVolumes()
{
    Model& model = this->simulator_->model();
    const auto size = model.numGridDof();
    std::vector<double> array(size);
    for (unsigned dof_idx = 0; dof_idx < size; ++dof_idx) {
        array[dof_idx] = model.dofTotalVolume(dof_idx);
    }
    return array;
}

template <class TypeTag>
std::vector<double>
PyMaterialState<TypeTag>::
getPorosity()
{
    Problem& problem = this->simulator_->problem();
    Model& model = this->simulator_->model();
    const auto size = model.numGridDof();
    std::vector<double> array(size);
    for (unsigned dof_idx = 0; dof_idx < size; ++dof_idx) {
        array[dof_idx] = problem.referencePorosity(dof_idx, /*timeIdx*/0);
    }
    return array;
}

template <class TypeTag>
void
PyMaterialState<TypeTag>::
setPorosity(const double* poro, std::size_t size)
{
    Problem& problem = this->simulator_->problem();
    Model& model = this->simulator_->model();
    const auto model_size = model.numGridDof();
    if (model_size != size) {
        const std::string msg = fmt::format(
            "Cannot set porosity. Expected array of size: {}, got array of size: ",
            model_size, size
        );
        throw std::runtime_error(msg);
    }
    for (unsigned dof_idx = 0; dof_idx < size; ++dof_idx) {
        problem.setPorosity(poro[dof_idx], dof_idx);
    }
}

} // namespace Opm::Pybind

#endif // OPM_PY_MATERIAL_STATE_IMPL_HEADER_INCLUDED
