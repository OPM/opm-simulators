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

#include <config.h>
#include <ebos/eclgenerictracermodel.hh>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/models/discretization/ecfv/ecfvstencil.hh>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Runspec.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TracerVdTable.hpp>

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#include <ebos/femcpgridcompat.hh>
#endif

#include <iostream>
#include <set>
#include <stdexcept>

namespace Opm {

template<class Grid, class GridView, class DofMapper, class Stencil, class Scalar>
EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
EclGenericTracerModel(const GridView& gridView,
                      const EclipseState& eclState,
                      const CartesianIndexMapper& cartMapper,
                      const DofMapper& dofMapper)
    : gridView_(gridView)
    , eclState_(eclState)
    , cartMapper_(cartMapper)
    , dofMapper_(dofMapper)
{
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
const std::string& EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
tracerName(int tracerIdx) const
{
    if (tracerNames_.empty())
        throw std::logic_error("This method should never be called when there are no tracers in the model");

    return tracerNames_[tracerIdx];
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
doInit(bool enabled, size_t numGridDof,
       size_t gasPhaseIdx, size_t oilPhaseIdx, size_t waterPhaseIdx)
{
    const auto& tracers = eclState_.tracer();
    const auto& comm = gridView_.comm();

    if (tracers.size() == 0)
        return; // tracer treatment is supposed to be disabled

    if (!enabled) {
        if (gridView_.comm().rank() == 0) {
            OpmLog::warning("Keyword TRACERS has only experimental support, and is hence ignored.\n"
                            "The experimental tracer model can still be used, but must be set explicitely.\n"
                            "To use tracers, set the command line option: --enable-tracer-model=true"
                            "\n");
        }
        return; // Tracer transport must be enabled by the user
    }

    if (comm.size() > 1) {
        tracerNames_.resize(0);
        if (comm.rank() == 0)
            std::cout << "Warning: The tracer model currently does not work for parallel runs\n"
                      << std::flush;
        return;
    }

    // retrieve the number of tracers from the deck
    const size_t numTracers = tracers.size();
    tracerNames_.resize(numTracers);
    tracerConcentration_.resize(numTracers);
    storageOfTimeIndex1_.resize(numTracers);

    // the phase where the tracer is
    tracerPhaseIdx_.resize(numTracers);
    size_t tracerIdx = 0;
    for (const auto& tracer : tracers) {
        tracerNames_[tracerIdx] = tracer.name;
        if (tracer.phase == Phase::WATER)
            tracerPhaseIdx_[tracerIdx] = waterPhaseIdx;
        else if (tracer.phase == Phase::OIL)
            tracerPhaseIdx_[tracerIdx] = oilPhaseIdx;
        else if (tracer.phase == Phase::GAS)
            tracerPhaseIdx_[tracerIdx] = gasPhaseIdx;

        tracerConcentration_[tracerIdx].resize(numGridDof);
        storageOfTimeIndex1_[tracerIdx].resize(numGridDof);


        //TBLK keyword
        if (!tracer.concentration.empty()){
            int tblkDatasize = tracer.concentration.size();
            if (tblkDatasize < cartMapper_.cartesianSize()){
                throw std::runtime_error("Wrong size of TBLK for" + tracer.name);
            }
            for (size_t globalDofIdx = 0; globalDofIdx < numGridDof; ++globalDofIdx){
                int cartDofIdx = cartMapper_.cartesianIndex(globalDofIdx);
                tracerConcentration_[tracerIdx][globalDofIdx] = tracer.concentration[cartDofIdx];
            }
        }
        //TVDPF keyword
        else {
            const auto& eclGrid = eclState_.getInputGrid();

            for (size_t globalDofIdx = 0; globalDofIdx < numGridDof; ++globalDofIdx){
                int cartDofIdx = cartMapper_.cartesianIndex(globalDofIdx);
                const auto& center = eclGrid.getCellCenter(cartDofIdx);
                tracerConcentration_[tracerIdx][globalDofIdx] = tracer.tvdpf.evaluate("TRACER_CONCENTRATION", center[2]);
            }
        }
        ++tracerIdx;
    }

    // initial tracer concentration
    tracerConcentrationInitial_ = tracerConcentration_;

    // residual of tracers
    tracerResidual_.resize(numGridDof);

    // allocate matrix for storing the Jacobian of the tracer residual
    tracerMatrix_ = new TracerMatrix(numGridDof, numGridDof, TracerMatrix::random);

    // find the sparsity pattern of the tracer matrix
    using NeighborSet = std::set<unsigned>;
    std::vector<NeighborSet> neighbors(numGridDof);

    Stencil stencil(gridView_, dofMapper_);
    auto elemIt = gridView_.template begin<0>();
    const auto elemEndIt = gridView_.template end<0>();
    for (; elemIt != elemEndIt; ++elemIt) {
        const auto& elem = *elemIt;
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

    const int sizeCartGrid = cartMapper_.cartesianSize();
    cartToGlobal_.resize(sizeCartGrid);
    for (unsigned i = 0; i < numGridDof; ++i) {
        int cartIdx = cartMapper_.cartesianIndex(i);
        cartToGlobal_[cartIdx] = i;
    }
}

template<class Grid,class GridView, class DofMapper, class Stencil, class Scalar>
bool EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>::
linearSolve_(const TracerMatrix& M, TracerVector& x, TracerVector& b)
{
#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2,7)
    Dune::FMatrixPrecision<Scalar>::set_singular_limit(1.e-30);
    Dune::FMatrixPrecision<Scalar>::set_absolute_limit(1.e-30);
#endif
    x = 0.0;
    Scalar tolerance = 1e-2;
    int maxIter = 100;

    int verbosity = 0;
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
}

#if HAVE_DUNE_FEM
template class EclGenericTracerModel<Dune::CpGrid,
                                     Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>,
                                     Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>>,
                                     Opm::EcfvStencil<double,Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>,false,false>,
                                     double>;
#else
template class EclGenericTracerModel<Dune::CpGrid,
                                     Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                     Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,Dune::Impl::MCMGFailLayout>,
                                     Opm::EcfvStencil<double,Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,false,false>,
                                     double>;
#endif

template class EclGenericTracerModel<Dune::PolyhedralGrid<3,3,double>,
                                     Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>,Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>,Dune::Impl::MCMGFailLayout>,
                                     Opm::EcfvStencil<double, Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>,false,false>,
                                     double>;

} // namespace Opm
