/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_GRAPHWELL_EQUATIONS_HEADER_INCLUDED
#define OPM_GRAPHWELL_EQUATIONS_HEADER_INCLUDED

#include <flowexperimental/graphwell/GraphWellTopology.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

#if HAVE_SUITESPARSE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif

#include <memory>
#include <stdexcept>
#include <vector>

namespace Opm {

/// \brief Well-internal linear system for the GraphWell model.
///
/// The DOF split (segment mass-conservation DOFs vs connection flux DOFs) maps onto
/// a 2x2 block system, stored as a Dune::MultiTypeBlockMatrix:
///   D = [ Dss (NP x NP)   Dsc (NP x 1) ]    segment rows    (mass equations)
///       [ Dcs (1 x NP)    Dcc (1 x 1)  ]    connection rows (momentum / control eqs)
/// with reservoir couplings B (in) and C (out) only on the segment rows.
///
/// The Schur complement onto the reservoir mirrors MultisegmentWellEquations: the well
/// system is eliminated with a direct solve. Because Dune's UMFPack is specialised for
/// plain BCRSMatrix (not MultiTypeBlockMatrix), D is flattened to a scalar BCRSMatrix
/// for factorisation; the sparsity is built once and only values are recopied.
template<class Scalar, int NumPhases, int NumResEq>
class GraphWellEquations
{
public:
    static constexpr int NP = NumPhases;

    using SegBlock = Dune::FieldMatrix<Scalar, NP, NP>;
    using SCBlock  = Dune::FieldMatrix<Scalar, NP, 1>;
    using CSBlock  = Dune::FieldMatrix<Scalar, 1, NP>;
    using CCBlock  = Dune::FieldMatrix<Scalar, 1, 1>;
    using DSS = Dune::BCRSMatrix<SegBlock>;
    using DSC = Dune::BCRSMatrix<SCBlock>;
    using DCS = Dune::BCRSMatrix<CSBlock>;
    using DCC = Dune::BCRSMatrix<CCBlock>;

    using Row0 = Dune::MultiTypeBlockVector<DSS, DSC>;
    using Row1 = Dune::MultiTypeBlockVector<DCS, DCC>;
    using DMatrix = Dune::MultiTypeBlockMatrix<Row0, Row1>;

    using SegVector  = Dune::BlockVector<Dune::FieldVector<Scalar, NP>>;
    using ConnVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;
    using BVectorWell = Dune::MultiTypeBlockVector<SegVector, ConnVector>;

    // reservoir coupling (segment rows x perforation columns)
    using OffDiagBlock = Dune::FieldMatrix<Scalar, NP, NumResEq>;
    using OffDiagMat = Dune::BCRSMatrix<OffDiagBlock>;
    using BVector = Dune::BlockVector<Dune::FieldVector<Scalar, NumResEq>>;

    // flattened scalar system for the direct solve
    using FlatBlock = Dune::FieldMatrix<Scalar, 1, 1>;
    using FlatMatrix = Dune::BCRSMatrix<FlatBlock>;
    using FlatVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;

    void init(const GraphWellTopology<Scalar>& topo);

    void clear();

    // -- assembly access ------------------------------------------------------
    DSS& Dss() { using namespace Dune::Indices; return D_[_0][_0]; }
    DSC& Dsc() { using namespace Dune::Indices; return D_[_0][_1]; }
    DCS& Dcs() { using namespace Dune::Indices; return D_[_1][_0]; }
    DCC& Dcc() { using namespace Dune::Indices; return D_[_1][_1]; }
    OffDiagMat& B() { return duneB_; }
    OffDiagMat& C() { return duneC_; }
    const OffDiagMat& B() const { return duneB_; }
    const OffDiagMat& C() const { return duneC_; }
    SegVector& resSeg() { using namespace Dune::Indices; return resWell_[_0]; }
    ConnVector& resConn() { using namespace Dune::Indices; return resWell_[_1]; }
    const BVectorWell& residual() const { return resWell_; }

    // -- direct solve & Schur -------------------------------------------------
    void createSolver();
    BVectorWell solve() const;                 //!< inv(D) * residual
    BVectorWell solve(const BVectorWell& rhs) const;

    //! Ax -= C^T inv(D) [B x; 0]   (Schur action on the reservoir)
    void apply(const BVector& x, BVector& Ax) const;
    //! r  -= C^T inv(D) residual   (well contribution to reservoir residual)
    void apply(BVector& r) const;
    //! xw = inv(D) (residual - [B x; 0])
    void recoverSolutionWell(const BVector& x, BVectorWell& xw) const;

    int numSegments() const { return nseg_; }
    int numConnections() const { return nconn_; }
    int numPerforations() const { return nperf_; }

    // -- inspection helpers (testing / debugging) -----------------------------
    //! Dense N x N image of the well Jacobian (N = NP*nseg + nconn).
    std::vector<std::vector<Scalar>> denseJacobian() const;
    //! Flattened residual vector (segment block first, then connections).
    std::vector<Scalar> flatResidual() const;

private:
    int flatSeg(int s, int eq) const { return NP * s + eq; }
    int flatConn(int c) const { return NP * nseg_ + c; }
    int flatSize() const { return NP * nseg_ + nconn_; }

    void fillFlatValues();
    FlatVector flatten(const BVectorWell& v) const;
    BVectorWell unflatten(const FlatVector& v) const;

    DMatrix D_;
    OffDiagMat duneB_;
    OffDiagMat duneC_;
    BVectorWell resWell_;

    FlatMatrix flat_;
#if HAVE_SUITESPARSE_UMFPACK
    mutable std::shared_ptr<Dune::UMFPack<FlatMatrix>> solver_;
#endif

    int nseg_{0};
    int nconn_{0};
    int nperf_{0};
    std::vector<std::vector<int>> seg_perfs_;   // for B/C row->perf column mapping
};

} // namespace Opm

#include <flowexperimental/graphwell/GraphWellEquations_impl.hpp>

#endif // OPM_GRAPHWELL_EQUATIONS_HEADER_INCLUDED
