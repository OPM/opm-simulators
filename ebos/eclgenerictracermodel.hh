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
 * \copydoc Opm::EclTracerModel
 */
#ifndef EWOMS_ECL_GENERIC_TRACER_MODEL_HH
#define EWOMS_ECL_GENERIC_TRACER_MODEL_HH

#include <dune/istl/bcrsmatrix.hh>

#include <opm/grid/common/CartesianIndexMapper.hpp>

#include <opm/models/blackoil/blackoilmodel.hh>

#include <opm/simulators/linalg/matrixblock.hh>

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

template<class Grid, class GridView, class DofMapper, class Stencil, class Scalar>
class EclGenericTracerModel {
public:
    using TracerMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<Scalar, 1, 1>>;
    using TracerVector = Dune::BlockVector<Dune::FieldVector<Scalar,1>>;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    static constexpr int dimWorld = Grid::dimensionworld;
    /*!
     * \brief Return the number of tracers considered by the tracerModel.
     */
    int numTracers() const;

    /*!
     * \brief Return the tracer name
     */
    const std::string& name(int tracerIdx) const;
    std::string fname(int tracerIdx) const;


    /*!
     * \brief Return the tracer concentration for tracer index and global DofIdx
     */
    Scalar tracerConcentration(int tracerIdx, int globalDofIdx) const;
    void setTracerConcentration(int tracerIdx, int globalDofIdx, Scalar value);

    /*!
    * \brief Return well tracer rates
    */
    const std::map<std::pair<std::string, std::string>, double>&
    getWellTracerRates() const {return wellTracerRate_;}

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(tracerConcentration_);
        serializer(wellTracerRate_);
    }

protected:
    EclGenericTracerModel(const GridView& gridView,
                          const EclipseState& eclState,
                          const CartesianIndexMapper& cartMapper,
                          const DofMapper& dofMapper,
                          const std::function<std::array<double,dimWorld>(int)> centroids);

    /*!
     * \brief Initialize all internal data structures needed by the tracer module
     */
    void doInit(bool rst,
                std::size_t numGridDof,
                std::size_t gasPhaseIdx,
                std::size_t oilPhaseIdx,
                std::size_t waterPhaseIdx);

    bool linearSolve_(const TracerMatrix& M, TracerVector& x, TracerVector& b);

    bool linearSolveBatchwise_(const TracerMatrix& M, std::vector<TracerVector>& x, std::vector<TracerVector>& b);

    double currentConcentration_(const Well& eclWell, const std::string& name) const;

    const GridView& gridView_;
    const EclipseState& eclState_;
    const CartesianIndexMapper& cartMapper_;
    const DofMapper& dofMapper_;

    std::vector<int> tracerPhaseIdx_;
    std::vector<Dune::BlockVector<Dune::FieldVector<Scalar, 1>>> tracerConcentration_;
    std::unique_ptr<TracerMatrix> tracerMatrix_;
    std::vector<Dune::BlockVector<Dune::FieldVector<Scalar, 1>>> storageOfTimeIndex1_;

    // <wellName, tracerIdx> -> wellRate
    std::map<std::pair<std::string, std::string>, double> wellTracerRate_;
    /// \brief Function returning the cell centers
    std::function<std::array<double,dimWorld>(int)> centroids_;
};

} // namespace Opm

#endif
