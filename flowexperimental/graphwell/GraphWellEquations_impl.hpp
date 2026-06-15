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

#ifndef OPM_GRAPHWELL_EQUATIONS_IMPL_HEADER_INCLUDED
#define OPM_GRAPHWELL_EQUATIONS_IMPL_HEADER_INCLUDED

#include <dune/istl/io.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <cmath>

namespace Opm {

template<class Scalar, int NP_, int NumResEq>
void GraphWellEquations<Scalar, NP_, NumResEq>::init(const GraphWellTopology<Scalar>& topo)
{
    using namespace Dune::Indices;
    nseg_ = topo.numSegments();
    nconn_ = topo.numConnections();
    nperf_ = topo.numPerforations();
    constexpr int surf = GraphWellTopology<Scalar>::surface_node;

    // ---- Dss : segment x segment ----
    Dune::MatrixIndexSet pss(nseg_, nseg_);
    for (int s = 0; s < nseg_; ++s)
        pss.add(s, s);
    for (int c = 0; c < nconn_; ++c) {
        const auto& conn = topo.connection(c);
        if (conn.up != surf && conn.down != surf) {
            pss.add(conn.up, conn.down);
            pss.add(conn.down, conn.up);
        }
    }
    pss.exportIdx(Dss());

    // ---- Dsc : segment x connection ----
    Dune::MatrixIndexSet psc(nseg_, nconn_);
    for (int c = 0; c < nconn_; ++c) {
        const auto& conn = topo.connection(c);
        if (conn.up != surf) psc.add(conn.up, c);
        if (conn.down != surf) psc.add(conn.down, c);
    }
    psc.exportIdx(Dsc());

    // ---- Dcs : connection x segment ----
    Dune::MatrixIndexSet pcs(nconn_, nseg_);
    for (int c = 0; c < nconn_; ++c) {
        const auto& conn = topo.connection(c);
        if (conn.up != surf) pcs.add(c, conn.up);
        if (conn.down != surf) pcs.add(c, conn.down);
    }
    pcs.exportIdx(Dcs());

    // ---- Dcc : connection x connection ----
    // Diagonal plus all pairs of connections sharing a segment. The off-diagonal
    // entries carry the reverse-flow pressure-equation and acceleration coupling
    // (a connection's momentum equation can depend on the flux of neighbouring
    // connections through the shared segment).
    Dune::MatrixIndexSet pcc(nconn_, nconn_);
    for (int c = 0; c < nconn_; ++c)
        pcc.add(c, c);
    for (int s = 0; s < nseg_; ++s) {
        const auto& incident = topo.connectionsOfSegment(s);
        for (int ci : incident)
            for (int cj : incident)
                pcc.add(ci, cj);
    }
    pcc.exportIdx(Dcc());

    // ---- B / C : segment x perforation ----
    seg_perfs_.assign(nseg_, {});
    for (int s = 0; s < nseg_; ++s)
        seg_perfs_[s] = topo.perforationsOfSegment(s);
    Dune::MatrixIndexSet pbc(nseg_, nperf_);
    for (int s = 0; s < nseg_; ++s)
        for (int p : seg_perfs_[s])
            pbc.add(s, p);
    pbc.exportIdx(duneB_);
    pbc.exportIdx(duneC_);

    // ---- residual ----
    resWell_[_0].resize(nseg_);
    resWell_[_1].resize(nconn_);

    // ---- flattened scalar pattern (built once) ----
    const int N = flatSize();
    Dune::MatrixIndexSet pflat(N, N);
    auto addBlock = [&](int row0, int nr, int col0, int nc) {
        for (int i = 0; i < nr; ++i)
            for (int j = 0; j < nc; ++j)
                pflat.add(row0 + i, col0 + j);
    };
    for (int s = 0; s < nseg_; ++s) addBlock(flatSeg(s, 0), NP, flatSeg(s, 0), NP); // diag seg
    for (int s = 0; s < nseg_; ++s) {                                              // Dcc shared-segment
        const auto& incident = topo.connectionsOfSegment(s);
        for (int ci : incident)
            for (int cj : incident)
                addBlock(flatConn(ci), 1, flatConn(cj), 1);
    }
    for (int c = 0; c < nconn_; ++c) {
        const auto& conn = topo.connection(c);
        addBlock(flatConn(c), 1, flatConn(c), 1); // Dcc diag
        if (conn.up != surf) {
            addBlock(flatSeg(conn.up, 0), NP, flatConn(c), 1);   // Dsc
            addBlock(flatConn(c), 1, flatSeg(conn.up, 0), NP);   // Dcs
        }
        if (conn.down != surf) {
            addBlock(flatSeg(conn.down, 0), NP, flatConn(c), 1);
            addBlock(flatConn(c), 1, flatSeg(conn.down, 0), NP);
        }
        if (conn.up != surf && conn.down != surf) {
            addBlock(flatSeg(conn.up, 0), NP, flatSeg(conn.down, 0), NP);   // Dss off-diag
            addBlock(flatSeg(conn.down, 0), NP, flatSeg(conn.up, 0), NP);
        }
    }
    pflat.exportIdx(flat_);

    clear();
}

template<class Scalar, int NP_, int NumResEq>
void GraphWellEquations<Scalar, NP_, NumResEq>::clear()
{
    Dss() = 0; Dsc() = 0; Dcs() = 0; Dcc() = 0;
    duneB_ = 0; duneC_ = 0;
    using namespace Dune::Indices;
    resWell_[_0] = 0;
    resWell_[_1] = 0;
}

template<class Scalar, int NP_, int NumResEq>
void GraphWellEquations<Scalar, NP_, NumResEq>::fillFlatValues()
{
    flat_ = 0;
    // Dss
    for (auto row = Dss().begin(); row != Dss().end(); ++row) {
        const int s = row.index();
        for (auto col = row->begin(); col != row->end(); ++col) {
            const int s2 = col.index();
            for (int i = 0; i < NP; ++i)
                for (int j = 0; j < NP; ++j)
                    flat_[flatSeg(s, i)][flatSeg(s2, j)] = (*col)[i][j];
        }
    }
    // Dsc
    for (auto row = Dsc().begin(); row != Dsc().end(); ++row) {
        const int s = row.index();
        for (auto col = row->begin(); col != row->end(); ++col) {
            const int c = col.index();
            for (int i = 0; i < NP; ++i)
                flat_[flatSeg(s, i)][flatConn(c)] = (*col)[i][0];
        }
    }
    // Dcs
    for (auto row = Dcs().begin(); row != Dcs().end(); ++row) {
        const int c = row.index();
        for (auto col = row->begin(); col != row->end(); ++col) {
            const int s = col.index();
            for (int j = 0; j < NP; ++j)
                flat_[flatConn(c)][flatSeg(s, j)] = (*col)[0][j];
        }
    }
    // Dcc
    for (auto row = Dcc().begin(); row != Dcc().end(); ++row) {
        const int c = row.index();
        for (auto col = row->begin(); col != row->end(); ++col)
            flat_[flatConn(c)][flatConn(col.index())] = (*col)[0][0];
    }
}

template<class Scalar, int NP_, int NumResEq>
typename GraphWellEquations<Scalar, NP_, NumResEq>::FlatVector
GraphWellEquations<Scalar, NP_, NumResEq>::flatten(const BVectorWell& v) const
{
    using namespace Dune::Indices;
    FlatVector f(flatSize());
    f = 0;
    for (int s = 0; s < nseg_; ++s)
        for (int eq = 0; eq < NP; ++eq)
            f[flatSeg(s, eq)] = v[_0][s][eq];
    for (int c = 0; c < nconn_; ++c)
        f[flatConn(c)] = v[_1][c][0];
    return f;
}

template<class Scalar, int NP_, int NumResEq>
typename GraphWellEquations<Scalar, NP_, NumResEq>::BVectorWell
GraphWellEquations<Scalar, NP_, NumResEq>::unflatten(const FlatVector& f) const
{
    using namespace Dune::Indices;
    BVectorWell v;
    v[_0].resize(nseg_);
    v[_1].resize(nconn_);
    for (int s = 0; s < nseg_; ++s)
        for (int eq = 0; eq < NP; ++eq)
            v[_0][s][eq] = f[flatSeg(s, eq)];
    for (int c = 0; c < nconn_; ++c)
        v[_1][c][0] = f[flatConn(c)];
    return v;
}

template<class Scalar, int NP_, int NumResEq>
void GraphWellEquations<Scalar, NP_, NumResEq>::createSolver()
{
#if HAVE_SUITESPARSE_UMFPACK
    if constexpr (std::is_same_v<Scalar, double>) {
        fillFlatValues();
        solver_ = std::make_shared<Dune::UMFPack<FlatMatrix>>(flat_, 0);
    } else {
        throw std::runtime_error("GraphWellEquations: UMFPack solve requires double.");
    }
#else
    throw std::runtime_error("GraphWellEquations: built without SuiteSparse/UMFPACK.");
#endif
}

template<class Scalar, int NP_, int NumResEq>
typename GraphWellEquations<Scalar, NP_, NumResEq>::BVectorWell
GraphWellEquations<Scalar, NP_, NumResEq>::solve(const BVectorWell& rhs) const
{
#if HAVE_SUITESPARSE_UMFPACK
    if (!solver_)
        throw std::runtime_error("GraphWellEquations::solve called before createSolver().");
    FlatVector b = flatten(rhs);
    FlatVector x(flatSize());
    x = 0;
    Dune::InverseOperatorResult res;
    solver_->apply(x, b, res);
    // Guard against a singular well system (UMFPack can return inf/nan silently).
    for (const auto& xi : x)
        if (!std::isfinite(xi[0]))
            OPM_THROW_NOLOG(NumericalProblem,
                            "GraphWellEquations: inf/nan from UMFPack (singular well system)");
    return unflatten(x);
#else
    throw std::runtime_error("GraphWellEquations: built without SuiteSparse/UMFPACK.");
#endif
}

template<class Scalar, int NP_, int NumResEq>
typename GraphWellEquations<Scalar, NP_, NumResEq>::BVectorWell
GraphWellEquations<Scalar, NP_, NumResEq>::solve() const
{
    return solve(resWell_);
}

template<class Scalar, int NP_, int NumResEq>
void GraphWellEquations<Scalar, NP_, NumResEq>::apply(const BVector& x, BVector& Ax) const
{
    using namespace Dune::Indices;
    SegVector Bx(nseg_);
    Bx = 0;
    duneB_.mv(x, Bx);
    BVectorWell rhs;
    rhs[_0] = Bx;
    rhs[_1].resize(nconn_);
    rhs[_1] = 0;
    const BVectorWell invDBx = solve(rhs);
    duneC_.mmtv(invDBx[_0], Ax);
}

template<class Scalar, int NP_, int NumResEq>
void GraphWellEquations<Scalar, NP_, NumResEq>::apply(BVector& r) const
{
    using namespace Dune::Indices;
    const BVectorWell invDres = solve(resWell_);
    duneC_.mmtv(invDres[_0], r);
}

template<class Scalar, int NP_, int NumResEq>
std::vector<std::vector<Scalar>>
GraphWellEquations<Scalar, NP_, NumResEq>::denseJacobian() const
{
    const int N = flatSize();
    std::vector<std::vector<Scalar>> M(N, std::vector<Scalar>(N, Scalar{0}));
    auto cthis = const_cast<GraphWellEquations*>(this);
    for (auto row = cthis->Dss().begin(); row != cthis->Dss().end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (int i = 0; i < NP; ++i)
                for (int j = 0; j < NP; ++j)
                    M[flatSeg(row.index(), i)][flatSeg(col.index(), j)] = (*col)[i][j];
    for (auto row = cthis->Dsc().begin(); row != cthis->Dsc().end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (int i = 0; i < NP; ++i)
                M[flatSeg(row.index(), i)][flatConn(col.index())] = (*col)[i][0];
    for (auto row = cthis->Dcs().begin(); row != cthis->Dcs().end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (int j = 0; j < NP; ++j)
                M[flatConn(row.index())][flatSeg(col.index(), j)] = (*col)[0][j];
    for (auto row = cthis->Dcc().begin(); row != cthis->Dcc().end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            M[flatConn(row.index())][flatConn(col.index())] = (*col)[0][0];
    return M;
}

template<class Scalar, int NP_, int NumResEq>
std::vector<Scalar>
GraphWellEquations<Scalar, NP_, NumResEq>::flatResidual() const
{
    using namespace Dune::Indices;
    std::vector<Scalar> r(flatSize(), Scalar{0});
    for (int s = 0; s < nseg_; ++s)
        for (int eq = 0; eq < NP; ++eq)
            r[flatSeg(s, eq)] = resWell_[_0][s][eq];
    for (int c = 0; c < nconn_; ++c)
        r[flatConn(c)] = resWell_[_1][c][0];
    return r;
}

template<class Scalar, int NP_, int NumResEq>
void GraphWellEquations<Scalar, NP_, NumResEq>::
recoverSolutionWell(const BVector& x, BVectorWell& xw) const
{
    using namespace Dune::Indices;
    SegVector Bx(nseg_);
    Bx = 0;
    duneB_.mv(x, Bx);
    BVectorWell rhs = resWell_;
    rhs[_0] -= Bx;
    xw = solve(rhs);
}

} // namespace Opm

#endif // OPM_GRAPHWELL_EQUATIONS_IMPL_HEADER_INCLUDED
