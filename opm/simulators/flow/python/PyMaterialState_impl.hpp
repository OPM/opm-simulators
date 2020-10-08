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

namespace Opm::Pybind {

template <class TypeTag>
std::unique_ptr<double []>
PyMaterialState<TypeTag>::
getCellVolumes( std::size_t *size)
{
    Model &model = ebosSimulator_->model();
    *size = model.numGridDof();
    auto array = std::make_unique<double []>(*size);
    for (unsigned dofIdx = 0; dofIdx < *size; ++dofIdx) {
        array[dofIdx] = model.dofTotalVolume(dofIdx);
    }
    return array;
}

template <class TypeTag>
std::unique_ptr<double []>
PyMaterialState<TypeTag>::
getPorosity( std::size_t *size)
{
    Problem &problem = ebosSimulator_->problem();
    Model &model = ebosSimulator_->model();
    *size = model.numGridDof();
    auto array = std::make_unique<double []>(*size);
    for (unsigned dofIdx = 0; dofIdx < *size; ++dofIdx) {
        array[dofIdx] = problem.referencePorosity(dofIdx, /*timeIdx*/0);
    }
    return array;
}

template <class TypeTag>
void
PyMaterialState<TypeTag>::
setPorosity(const double *poro, std::size_t size)
{
    Problem &problem = ebosSimulator_->problem();
    Model &model = ebosSimulator_->model();
    auto model_size = model.numGridDof();
    if (model_size != size) {
        std::ostringstream message;
        message << "Cannot set porosity. Expected array of size: "
                << model_size << ", got array of size: " << size;
        throw std::runtime_error(message.str());
    }
    for (unsigned dofIdx = 0; dofIdx < size; ++dofIdx) {
        problem.setPorosity(poro[dofIdx], dofIdx);
    }
}
} //namespace Opm::Pybind
