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
/**
 * \file
 *
 * \copydoc Opm::TemperatureModel
 */
#ifndef OPM_GENERIC_TEMPERATURE_MODEL_HPP
#define OPM_GENERIC_TEMPERATURE_MODEL_HPP

#include <dune/istl/bcrsmatrix.hh>

#include <opm/grid/common/CartesianIndexMapper.hpp>

#include <opm/models/blackoil/blackoilmodel.hh>

#include <opm/simulators/linalg/matrixblock.hh>

#include <opm/input/eclipse/EclipseState/Phase.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Opm {

class EclipseState;
class Well;

template<class Grid, class GridView, class DofMapper, class Stencil, class FluidSystem, class Scalar>
class GenericTemperatureModel {
public:
    using EnergyMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<Scalar, 1, 1>>;
    using EnergyVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    static constexpr int dimWorld = Grid::dimensionworld;

    bool doTemp()
    {
        return doTemp_;
    }

    const Scalar temperature(size_t globalIdx) const {
        return temperature_[globalIdx];
    }

protected:
    GenericTemperatureModel(const GridView& gridView,
                           const EclipseState& eclState,
                           const CartesianIndexMapper& cartMapper,
                           const DofMapper& dofMapper);

    /*!
     * \brief Initialize all internal data structures needed by the temperature module
     */
    void doInit(std::size_t numGridDof);

    bool linearSolve_(const EnergyMatrix& M, EnergyVector& x, EnergyVector& b);

    void syncOverlap_();

    const GridView& gridView_;
    const EclipseState& eclState_;
    const CartesianIndexMapper& cartMapper_;
    const DofMapper& dofMapper_;

    EnergyVector energyVector_;
    std::unique_ptr<EnergyMatrix> energyMatrix_;
    std::vector<Scalar> temperature_;
    bool doTemp_{false};

};

} // namespace Opm

#endif // OPM_GENERIC_TEMPERATURE_MODEL_HPP
