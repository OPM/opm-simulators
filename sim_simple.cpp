/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#include "AutoDiffBlock.hpp"
#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <Eigen/UmfPackSupport>

/*
  Equations for incompressible two-phase flow.

  Using s and p as variables:

  PV (s_i - s0_i) / dt + sum_{j \in U(i)} f(s_j) v_{ij} + sum_{j in D(i) f(s_i) v_{ij} = qw_i

  where

  v_{ij} = totmob_ij T_ij (p_i - p_j)


  Pressure equation:

  sum_{j \in N(i)} totmob_ij T_ij (p_i - p_j) = q_i

*/

/// Contains vectors and sparse matrices that represent subsets or
/// operations on (AD or regular) vectors of data.
struct HelperOps
{
    typedef AutoDiff::ForwardBlock<double>::M M;
    typedef AutoDiff::ForwardBlock<double>::V V;

    /// A list of internal faces.
    Eigen::Array<int, Eigen::Dynamic, 1> internal_faces;

    /// Extract for each face the difference of its adjacent cells'values.
    M ngrad;
    /// Extract for each face the average of its adjacent cells' values.
    M caver;
    /// Extract for each cell the sum of its adjacent faces' (signed) values.
    M div;

    /// Constructs all helper vectors and matrices.
    HelperOps(const UnstructuredGrid& grid)
    {
        const int nc = grid.number_of_cells;
        const int nf = grid.number_of_faces;
        // Define some neighbourhood-derived helper arrays.
        typedef Eigen::Array<int, Eigen::Dynamic, 1> OneColInt;
        typedef Eigen::Array<bool, Eigen::Dynamic, 1> OneColBool;
        typedef Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor> TwoColInt;
        typedef Eigen::Array<bool, Eigen::Dynamic, 2, Eigen::RowMajor> TwoColBool;
        TwoColInt nb = Eigen::Map<TwoColInt>(grid.face_cells, nf, 2);
        // std::cout << "nb = \n" << nb << std::endl;
        TwoColBool nbib = nb >= 0;
        OneColBool ifaces = nbib.rowwise().all();
        const int num_internal = ifaces.cast<int>().sum();
        // std::cout << num_internal << " internal faces." << std::endl;
        TwoColInt nbi(num_internal, 2);
        internal_faces.resize(num_internal);
        int fi = 0;
        for (int f = 0; f < nf; ++f) {
            if (ifaces[f]) {
                internal_faces[fi] = f;
                nbi.row(fi) = nb.row(f);
                ++fi;
            }
        }
        // std::cout << "nbi = \n" << nbi << std::endl;
        // Create matrices.
        ngrad.resize(num_internal, nc);
        caver.resize(num_internal, nc);
        typedef Eigen::Triplet<double> Tri;
        std::vector<Tri> ngrad_tri;
        std::vector<Tri> caver_tri;
        ngrad_tri.reserve(2*num_internal);
        caver_tri.reserve(2*num_internal);
        for (int i = 0; i < num_internal; ++i) {
            ngrad_tri.emplace_back(i, nbi(i,0), 1.0);
            ngrad_tri.emplace_back(i, nbi(i,1), -1.0);
            caver_tri.emplace_back(i, nbi(i,0), 0.5);
            caver_tri.emplace_back(i, nbi(i,1), 0.5);
        }
        ngrad.setFromTriplets(ngrad_tri.begin(), ngrad_tri.end());
        caver.setFromTriplets(caver_tri.begin(), caver_tri.end());
        div = ngrad.transpose();
    }
};



int main()
{
    typedef AutoDiff::ForwardBlock<double> ADB;
    typedef ADB::V V;
    typedef ADB::M M;

    Opm::time::StopWatch clock;
    clock.start();
    Opm::GridManager gm(50, 50, 10);
    const UnstructuredGrid& grid = *gm.c_grid();
    using namespace Opm::unit;
    using namespace Opm::prefix;
    Opm::IncompPropertiesBasic props(2, Opm::SaturationPropsBasic::Quadratic,
                                     { 1000.0, 800.0 },
                                     { 1.0*centi*Poise, 5.0*centi*Poise },
                                     0.2, 100*milli*darcy,
                                     grid.dimensions, grid.number_of_cells);
    std::vector<double> htrans(grid.cell_facepos[grid.number_of_cells]);
    tpfa_htrans_compute((UnstructuredGrid*)&grid, props.permeability(), htrans.data());
    // std::vector<double> trans(grid.number_of_faces);
    V trans_all(grid.number_of_faces);
    tpfa_trans_compute((UnstructuredGrid*)&grid, htrans.data(), trans_all.data());
    const int nc = grid.number_of_cells;
    std::vector<int> allcells(nc);
    for (int i = 0; i < nc; ++i) {
        allcells[i] = i;
    }
    std::cerr << "Opm core " << clock.secsSinceLast() << std::endl;

    // Define neighbourhood-derived operator matrices.
    HelperOps ops(grid);
    const int num_internal = ops.internal_faces.size();
    V transi(num_internal);
    for (int fi = 0; fi < num_internal; ++fi) {
        transi[fi] = trans_all[ops.internal_faces[fi]];
    }
    std::cerr << "Topology matrices " << clock.secsSinceLast() << std::endl;

    typedef AutoDiff::ForwardBlock<double> ADB;
    typedef ADB::V V;

    // q
    V q(nc);
    q.setZero();
    q[0] = 1.0;
    q[nc-1] = -1.0;

    // s - this is explicit now
    typedef Eigen::Array<double, Eigen::Dynamic, 2, Eigen::RowMajor> TwoCol;
    TwoCol s(nc, 2);
    s.leftCols<1>().setZero();
    s.rightCols<1>().setOnes();

    // totmob - explicit as well
    TwoCol kr(nc, 2);
    props.relperm(nc, s.data(), allcells.data(), kr.data(), 0);
    V krw = kr.leftCols<1>();
    V kro = kr.rightCols<1>();
    const double* mu = props.viscosity();
    V totmob = krw/mu[0] + kro/mu[1];
    V totmobf = (ops.caver*totmob.matrix()).array();

    // Mobility-weighted transmissibilities per internal face.
    // Still explicit, and no upwinding!
    V mobtransf = totmobf*transi;

    std::cerr << "Property arrays " << clock.secsSinceLast() << std::endl;

    // Initial pressure.
    V p0(nc,1);
    p0.fill(200*Opm::unit::barsa);

    // First actual AD usage: defining pressure variable.
    std::vector<int> block_pattern = { nc };
    // Could actually write { nc } instead of block_pattern below,
    // but we prefer a named variable since we will repeat it.
    ADB p = ADB::variable(0, p0, block_pattern);
    ADB ngradp = ops.ngrad*p;
    // We want flux = totmob*trans*(p_i - p_j) for the ij-face.
    // We only need to multiply mobtransf and pdiff_face,
    // but currently multiplication with constants is not in,
    // so we define an AD constant to multiply with.
    ADB mobtransf_ad = ADB::constant(mobtransf, block_pattern);
    ADB flux = mobtransf_ad*ngradp;
    ADB residual = ops.div*flux - ADB::constant(q, block_pattern);
    std::cerr << "Construct AD residual " << clock.secsSinceLast() << std::endl;

    // std::cout << div << pdiff_face;
    // std::cout << div*pdiff_face;
    // std::cout << q << std::endl;
    // std::cout << residual << std::endl;

    // It's the residual we want to be zero. We know it's linear in p,
    // so we just need a single linear solve. Since we have formulated
    // ourselves with a residual and jacobian we do this with a single
    // Newton step (hopefully easy to extend later):
    //   p = p0 - J(p0) \ R(p0)
    // Where R(p0) and J(p0) are contained in residual.value() and
    // residual.derived()[0].

    Eigen::UmfPackLU<M> solver;
    M matr = residual.derivative()[0];
    matr.coeffRef(0,0) *= 2.0;
    matr.makeCompressed();
    solver.compute(residual.derivative()[0]);
    // if (solver.info() != Eigen::Succeeded) {
    //     std::cerr << "Decomposition error!\n";
    //     return 1;
    // }
    Eigen::VectorXd x = solver.solve(residual.value().matrix());
    // if (solver.info() != Eigen::Succeeded) {
    //     std::cerr << "Solve failure!\n";
    //     return 1;
    // }
    V p_new = p0 - x.array();
    std::cerr << "Solve " << clock.secsSinceLast() << std::endl;
    std::cout << p_new << std::endl;
}
