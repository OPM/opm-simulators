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
#ifndef EWOMS_ECL_GENERIC_TRACER_MODEL_IMPL_HH
#define EWOMS_ECL_GENERIC_TRACER_MODEL_IMPL_HH

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>

#include <ebos/eclgenerictracermodel.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TracerVdTable.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTracerProperties.hpp>

#include <opm/models/discretization/ecfv/ecfvstencil.hh>

#include <opm/simulators/linalg/ilufirstelement.hh>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>

#include <fmt/format.h>

#include <array>
#include <functional>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

namespace Opm {

#if HAVE_MPI
template<class M, class V>
struct TracerSolverSelector
{
    using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;
    using TracerOperator = Dune::OverlappingSchwarzOperator<M, V, V, Comm>;
    using type = Dune::FlexibleSolver<TracerOperator>;
};

template<class Vector, class Grid, class Matrix>
std::tuple<std::unique_ptr<Dune::OverlappingSchwarzOperator<Matrix,Vector,Vector,
                                                            Dune::OwnerOverlapCopyCommunication<int,int>>>,
           std::unique_ptr<typename TracerSolverSelector<Matrix,Vector>::type>>
createParallelFlexibleSolver(const Grid&, const Matrix&, const PropertyTree&)
{
    OPM_THROW(std::logic_error, "Grid not supported for parallel Tracers.");
    return {nullptr, nullptr};
}

template<class Vector, class Matrix>
std::tuple<std::unique_ptr<Dune::OverlappingSchwarzOperator<Matrix,Vector,Vector,
                                                            Dune::OwnerOverlapCopyCommunication<int,int>>>,
           std::unique_ptr<typename TracerSolverSelector<Matrix,Vector>::type>>
createParallelFlexibleSolver(const Dune::CpGrid& grid, const Matrix& M, const PropertyTree& prm)
{
        using TracerOperator = Dune::OverlappingSchwarzOperator<Matrix,Vector,Vector,
                                                                Dune::OwnerOverlapCopyCommunication<int,int>>;
        using TracerSolver = Dune::FlexibleSolver<TracerOperator>;
        const auto& cellComm = grid.cellCommunication();
        auto op = std::make_unique<TracerOperator>(M, cellComm);
        auto dummyWeights = [](){ return Vector();};
        return {std::move(op), std::make_unique<TracerSolver>(*op, cellComm, prm, dummyWeights, 0)};
}
#endif

template<class Grid, class GridView, class DofMapper, class Stencil, class Scalar>
EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
EclGenericTracerModel(const GridView& gridView,
                      const EclipseState& eclState,
                      const CartesianIndexMapper& cartMapper,
                      const DofMapper& dofMapper,
                      const std::function<std::array<double,dimWorld>(int)> centroids)
    : gridView_(gridView)
    , eclState_(eclState)
    , cartMapper_(cartMapper)
    , dofMapper_(dofMapper)
    , centroids_(centroids)
{
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
Scalar EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
tracerConcentration(int tracerIdx, int globalDofIdx) const
{
    if (tracerConcentration_.empty())
        return 0.0;

    return tracerConcentration_[tracerIdx][globalDofIdx];
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
void EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
setTracerConcentration(int tracerIdx, int globalDofIdx, Scalar value)
{
    this->tracerConcentration_[tracerIdx][globalDofIdx] = value;
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
int EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
numTracers() const
{
    return this->eclState_.tracer().size();
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
std::string EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
fname(int tracerIdx) const
{
    return this->eclState_.tracer()[tracerIdx].fname();
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
double EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
currentConcentration_(const Well& eclWell, const std::string& name) const
{
    return eclWell.getTracerProperties().getConcentration(name);
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
const std::string& EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
name(int tracerIdx) const
{
    return this->eclState_.tracer()[tracerIdx].name;
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
void EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
doInit(bool rst, std::size_t numGridDof,
       std::size_t gasPhaseIdx, std::size_t oilPhaseIdx, std::size_t waterPhaseIdx)
{
    const auto& tracers = eclState_.tracer();

    if (tracers.size() == 0)
        return; // tracer treatment is supposed to be disabled

    // retrieve the number of tracers from the deck
    const std::size_t numTracers = tracers.size();
    tracerConcentration_.resize(numTracers);
    storageOfTimeIndex1_.resize(numTracers);

    // the phase where the tracer is
    tracerPhaseIdx_.resize(numTracers);
    for (std::size_t tracerIdx = 0; tracerIdx < numTracers; tracerIdx++) {
        const auto& tracer = tracers[tracerIdx];

        if (tracer.phase == Phase::WATER)
            tracerPhaseIdx_[tracerIdx] = waterPhaseIdx;
        else if (tracer.phase == Phase::OIL)
            tracerPhaseIdx_[tracerIdx] = oilPhaseIdx;
        else if (tracer.phase == Phase::GAS)
            tracerPhaseIdx_[tracerIdx] = gasPhaseIdx;

        tracerConcentration_[tracerIdx].resize(numGridDof);
        storageOfTimeIndex1_[tracerIdx].resize(numGridDof);

        if (rst)
            continue;


        //TBLK keyword
        if (tracer.free_concentration.has_value()){
            const auto& free_concentration = tracer.free_concentration.value();
            int tblkDatasize = free_concentration.size();
            if (tblkDatasize < cartMapper_.cartesianSize()){
                throw std::runtime_error("Wrong size of TBLK for" + tracer.name);
            }
            for (std::size_t globalDofIdx = 0; globalDofIdx < numGridDof; ++globalDofIdx) {
                int cartDofIdx = cartMapper_.cartesianIndex(globalDofIdx);
                tracerConcentration_[tracerIdx][globalDofIdx] = free_concentration[cartDofIdx];
            }
        }
        //TVDPF keyword
        else if (tracer.free_tvdp.has_value()) {
            const auto& free_tvdp = tracer.free_tvdp.value();
            for (std::size_t globalDofIdx = 0; globalDofIdx < numGridDof; ++globalDofIdx) {
                tracerConcentration_[tracerIdx][globalDofIdx] =
                    free_tvdp.evaluate("TRACER_CONCENTRATION",
                                       centroids_(globalDofIdx)[2]);
            }
        } else
            throw std::logic_error(fmt::format("Can not initialize tracer: {}", tracer.name));
    }

    // allocate matrix for storing the Jacobian of the tracer residual
    tracerMatrix_ = std::make_unique<TracerMatrix>(numGridDof, numGridDof, TracerMatrix::random);

    // find the sparsity pattern of the tracer matrix
    using NeighborSet = std::set<unsigned>;
    std::vector<NeighborSet> neighbors(numGridDof);

    Stencil stencil(gridView_, dofMapper_);
    for (const auto& elem : elements(gridView_)) {
        stencil.update(elem);

        for (unsigned primaryDofIdx = 0; primaryDofIdx < stencil.numPrimaryDof(); ++primaryDofIdx) {
            unsigned myIdx = stencil.globalSpaceIndex(primaryDofIdx);

            for (unsigned dofIdx = 0; dofIdx < stencil.numDof(); ++dofIdx) {
                unsigned neighborIdx = stencil.globalSpaceIndex(dofIdx);
                neighbors[myIdx].insert(neighborIdx);
            }
        }
    }

    // allocate space for the rows of the matrix
    for (unsigned dofIdx = 0; dofIdx < numGridDof; ++ dofIdx)
        tracerMatrix_->setrowsize(dofIdx, neighbors[dofIdx].size());
    tracerMatrix_->endrowsizes();

    // fill the rows with indices. each degree of freedom talks to
    // all of its neighbors. (it also talks to itself since
    // degrees of freedom are sometimes quite egocentric.)
    for (unsigned dofIdx = 0; dofIdx < numGridDof; ++ dofIdx) {
        typename NeighborSet::iterator nIt = neighbors[dofIdx].begin();
        typename NeighborSet::iterator nEndIt = neighbors[dofIdx].end();
        for (; nIt != nEndIt; ++nIt)
            tracerMatrix_->addindex(dofIdx, *nIt);
    }
    tracerMatrix_->endindices();
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
bool EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
linearSolve_(const TracerMatrix& M, TracerVector& x, TracerVector& b)
{
    x = 0.0;
    Scalar tolerance = 1e-2;
    int maxIter = 100;

    int verbosity = 0;
    PropertyTree prm;
    prm.put("maxiter", maxIter);
    prm.put("tol", tolerance);
    prm.put("verbosity", verbosity);
    prm.put("solver", std::string("bicgstab"));
    prm.put("preconditioner.type", std::string("ParOverILU0"));

#if HAVE_MPI
    if(gridView_.grid().comm().size() > 1)
    {
        auto [tracerOperator, solver] =
            createParallelFlexibleSolver<TracerVector>(gridView_.grid(), M, prm);
        (void) tracerOperator;

        Dune::InverseOperatorResult result;
        solver->apply(x, b, result);

        // return the result of the solver
        return result.converged;
    }
    else
    {
#endif
        using TracerSolver = Dune::BiCGSTABSolver<TracerVector>;
        using TracerOperator = Dune::MatrixAdapter<TracerMatrix,TracerVector,TracerVector>;
        using TracerScalarProduct = Dune::SeqScalarProduct<TracerVector>;
        using TracerPreconditioner = Dune::SeqILU< TracerMatrix,TracerVector,TracerVector>;

        TracerOperator tracerOperator(M);
        TracerScalarProduct tracerScalarProduct;
        TracerPreconditioner tracerPreconditioner(M, 0, 1); // results in ILU0

        TracerSolver solver (tracerOperator, tracerScalarProduct,
                             tracerPreconditioner, tolerance, maxIter,
                             verbosity);

        Dune::InverseOperatorResult result;
        solver.apply(x, b, result);

        // return the result of the solver
        return result.converged;
#if HAVE_MPI
    }
#endif
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
bool  EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
linearSolveBatchwise_(const TracerMatrix& M, std::vector<TracerVector>& x, std::vector<TracerVector>& b)
{
    Scalar tolerance = 1e-2;
    int maxIter = 100;

    int verbosity = 0;
    PropertyTree prm;
    prm.put("maxiter", maxIter);
    prm.put("tol", tolerance);
    prm.put("verbosity", verbosity);
    prm.put("solver", std::string("bicgstab"));
    prm.put("preconditioner.type", std::string("ParOverILU0"));

#if HAVE_MPI
    if(gridView_.grid().comm().size() > 1)
    {
        auto [tracerOperator, solver] =
            createParallelFlexibleSolver<TracerVector>(gridView_.grid(), M, prm);
        (void) tracerOperator;
        bool converged = true;
        for (std::size_t nrhs = 0; nrhs < b.size(); ++nrhs) {
            x[nrhs] = 0.0;
            Dune::InverseOperatorResult result;
            solver->apply(x[nrhs], b[nrhs], result);
            converged = (converged && result.converged);
        }
        return converged;
    }
    else
    {
#endif
        using TracerSolver = Dune::BiCGSTABSolver<TracerVector>;
        using TracerOperator = Dune::MatrixAdapter<TracerMatrix,TracerVector,TracerVector>;
        using TracerScalarProduct = Dune::SeqScalarProduct<TracerVector>;
        using TracerPreconditioner = Dune::SeqILU< TracerMatrix,TracerVector,TracerVector>;

        TracerOperator tracerOperator(M);
        TracerScalarProduct tracerScalarProduct;
        TracerPreconditioner tracerPreconditioner(M, 0, 1); // results in ILU0

        TracerSolver solver (tracerOperator, tracerScalarProduct,
                             tracerPreconditioner, tolerance, maxIter,
                             verbosity);

        bool converged = true;
        for (std::size_t nrhs = 0; nrhs < b.size(); ++nrhs) {
            x[nrhs] = 0.0;
            Dune::InverseOperatorResult result;
            solver.apply(x[nrhs], b[nrhs], result);
            converged = (converged && result.converged);
        }

        // return the result of the solver
        return converged;
#if HAVE_MPI
    }
#endif
}


} // namespace Opm
#endif
