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
#ifndef OPM_GENERIC_TEMPERATURE_MODEL_IMPL_HPP
#define OPM_GENERIC_TEMPERATURE_MODEL_IMPL_HPP

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
//#include <opm/input/eclipse/EclipseState/Tables/TemperatureVdTable.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
//#include <opm/input/eclipse/Schedule/Well/WellTemperatureProperties.hpp>

#include <opm/models/discretization/ecfv/ecfvstencil.hh>

#include <opm/simulators/flow/GenericTemperatureModel.hpp>
#include <opm/simulators/linalg/ilufirstelement.hh>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>

#include <fmt/format.h>

#include <array>
#include <functional>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

namespace Opm {

#if HAVE_MPI
template<class M, class V>
struct EnergySolverSelector
{
    using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;
    using EnergyOperator = Dune::OverlappingSchwarzOperator<M, V, V, Comm>;
    using type = Dune::FlexibleSolver<EnergyOperator>;
};

template<class Vector, class Grid, class Matrix>
std::tuple<std::unique_ptr<Dune::OverlappingSchwarzOperator<Matrix,Vector,Vector,
                                                            Dune::OwnerOverlapCopyCommunication<int,int>>>,
           std::unique_ptr<typename EnergySolverSelector<Matrix,Vector>::type>>
createParallelFlexibleSolver(const Grid&, const Matrix&, const PropertyTree&)
{
    OPM_THROW(std::logic_error, "Grid not supported for parallel Temperatures.");
    return {nullptr, nullptr};
}

template<class Vector, class Matrix>
std::tuple<std::unique_ptr<Dune::OverlappingSchwarzOperator<Matrix,Vector,Vector,
                                                            Dune::OwnerOverlapCopyCommunication<int,int>>>,
           std::unique_ptr<typename EnergySolverSelector<Matrix,Vector>::type>>
createParallelFlexibleSolver(const Dune::CpGrid& grid, const Matrix& M, const PropertyTree& prm)
{
        using EnergyOperator = Dune::OverlappingSchwarzOperator<Matrix,Vector,Vector,
                                                                Dune::OwnerOverlapCopyCommunication<int,int>>;
        using EnergySolver = Dune::FlexibleSolver<EnergyOperator>;
        const auto& cellComm = grid.cellCommunication();
        auto op = std::make_unique<EnergyOperator>(M, cellComm);
        auto dummyWeights = [](){ return Vector();};
        return {std::move(op), std::make_unique<EnergySolver>(*op, cellComm, prm, dummyWeights, 0)};
}
#endif

template<class Grid, class GridView, class DofMapper, class Stencil, class FluidSystem, class Scalar>
GenericTemperatureModel<Grid,GridView,DofMapper,Stencil,FluidSystem,Scalar>::
GenericTemperatureModel(const GridView& gridView,
                   const EclipseState& eclState,
                   const CartesianIndexMapper& cartMapper,
                   const DofMapper& dofMapper)
    : gridView_(gridView)
    , eclState_(eclState)
    , cartMapper_(cartMapper)
    , dofMapper_(dofMapper)
{
}


template<class Grid, class GridView, class DofMapper, class Stencil, class FluidSystem, class Scalar>
void GenericTemperatureModel<Grid,GridView,DofMapper,Stencil,FluidSystem,Scalar>::
doInit(bool rst, std::size_t numGridDof)
{
    doTemp_ = eclState_.getSimulationConfig().isTemp();

    temperature_.resize(numGridDof);
    energyVector_.resize(numGridDof);
    // allocate matrix for storing the Jacobian of the temperature residual
    energyMatrix_ = std::make_unique<EnergyMatrix>(numGridDof, numGridDof, EnergyMatrix::random);

    // find the sparsity pattern of the temperature matrix
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
    for (unsigned dofIdx = 0; dofIdx < numGridDof; ++ dofIdx) {
        energyMatrix_->setrowsize(dofIdx, neighbors[dofIdx].size());
    }
    energyMatrix_->endrowsizes();

    // fill the rows with indices. each degree of freedom talks to
    // all of its neighbors. (it also talks to itself since
    // degrees of freedom are sometimes quite egocentric.)
    for (unsigned dofIdx = 0; dofIdx < numGridDof; ++ dofIdx) {
        typename NeighborSet::iterator nIt = neighbors[dofIdx].begin();
        typename NeighborSet::iterator nEndIt = neighbors[dofIdx].end();
        for (; nIt != nEndIt; ++nIt) {
            energyMatrix_->addindex(dofIdx, *nIt);
        }
    }
    energyMatrix_->endindices();
}

template<class Grid, class GridView, class DofMapper, class Stencil, class FluidSystem, class Scalar>
bool GenericTemperatureModel<Grid,GridView,DofMapper,Stencil,FluidSystem,Scalar>::
linearSolve_(const EnergyMatrix& M, EnergyVector& x, EnergyVector& b)
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
        auto [energyOperator, solver] =
            createParallelFlexibleSolver<EnergyVector>(gridView_.grid(), M, prm);
        (void) energyOperator;

        Dune::InverseOperatorResult result;
        solver->apply(x, b, result);

        // return the result of the solver
        return result.converged;
    }
    else
    {
#endif
        using EnergySolver = Dune::BiCGSTABSolver<EnergyVector>;
        using EnergyOperator = Dune::MatrixAdapter<EnergyMatrix,EnergyVector,EnergyVector>;
        using EnergyScalarProduct = Dune::SeqScalarProduct<EnergyVector>;
        using EnergyPreconditioner = Dune::SeqILU< EnergyMatrix,EnergyVector,EnergyVector>;

        EnergyOperator energyOperator(M);
        EnergyScalarProduct energyScalarProduct;
        EnergyPreconditioner energyPreconditioner(M, 0, 1); // results in ILU0

        EnergySolver solver (energyOperator, energyScalarProduct,
                                  energyPreconditioner, tolerance, maxIter,
                                  verbosity);

        Dune::InverseOperatorResult result;
        solver.apply(x, b, result);

        // return the result of the solver
        return result.converged;
#if HAVE_MPI
    }
#endif
}

template<class Grid, class GridView, class DofMapper, class Stencil, class FluidSystem, class Scalar>
void GenericTemperatureModel<Grid,GridView,DofMapper,Stencil,FluidSystem,Scalar>::
syncOverlap_()
{
#if HAVE_MPI
    // syncronize the solution on the ghost and overlap elements
    using GhostSyncHandle = GridCommHandleGhostSync<Dune::FieldVector<Scalar, 1>,
                                                    EnergyVector,
                                                    DofMapper,
                                                    /*commCodim=*/0>;

    auto ghostSync = GhostSyncHandle(this->energyVector_,
                                     this->dofMapper_);
    gridView_.communicate(ghostSync,
                          Dune::InteriorBorder_All_Interface,
                          Dune::ForwardCommunication);
#endif
}

} // namespace Opm

#endif // OPM_GENERIC_TEMPERATURE_MODEL_IMPL_HPP
