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

#include <config.h>

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/GridManager.hpp>
#include <opm/core/props/IncompPropertiesBasic.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/grid/utility/StopWatch.hpp>
#include <opm/grid/transmissibility/trans_tpfa.h>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#if HAVE_SUITESPARSE_UMFPACK_H
#include <Eigen/UmfPackSupport>
#else
#include <Eigen/IterativeLinearSolvers>
#endif
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <iostream>
#include <cstdlib>

/*
  Equations for incompressible two-phase flow.

  Using s and p as variables:

  PV (s_i - s0_i) / dt + sum_{j \in U(i)} f(s_j) v_{ij} + sum_{j in D(i) f(s_i) v_{ij} = qw_i

  where

  v_{ij} = totmob_ij T_ij (p_i - p_j)


  Pressure equation:

  sum_{j \in N(i)} totmob_ij T_ij (p_i - p_j) = q_i

*/


template <class ADB>
std::vector<ADB>
phaseMobility(const Opm::IncompPropertiesInterface& props,
              const std::vector<int>& cells,
              const typename ADB::V& sw)
{
    typedef Eigen::Array<double, Eigen::Dynamic, 2, Eigen::RowMajor> TwoCol;
    typedef Eigen::Array<double, Eigen::Dynamic, 4, Eigen::RowMajor> FourCol;
    typedef Eigen::SparseMatrix<double> S;
    typedef typename ADB::V V;
    typedef typename ADB::M M;
    const int nc = props.numCells();
    TwoCol s(nc, 2);
    s.leftCols<1>() = sw;
    s.rightCols<1>() = 1.0 - s.leftCols<1>();
    TwoCol kr(nc, 2);
    FourCol dkr(nc, 4);
    props.relperm(nc, s.data(), cells.data(), kr.data(), dkr.data());
    V krw = kr.leftCols<1>();
    V kro = kr.rightCols<1>();
    V dkrw = dkr.leftCols<1>();  // Left column is top-left of dkr/ds 2x2 matrix.
    V dkro = -dkr.rightCols<1>(); // Right column is bottom-right of dkr/ds 2x2 matrix.
    S krwjac(nc,nc);
    S krojac(nc,nc);
    auto sizes = Eigen::ArrayXi::Ones(nc);
    krwjac.reserve(sizes);
    krojac.reserve(sizes);
    for (int c = 0; c < nc; ++c) {
        krwjac.insert(c,c) = dkrw(c);
        krojac.insert(c,c) = dkro(c);
    }
    const double* mu = props.viscosity();
    std::vector<M> dmw = { M(krwjac)/mu[0] };
    std::vector<M> dmo = { M(krojac)/mu[1] };

    std::vector<ADB> pmobc = { ADB::function(krw / mu[0], std::move(dmw)) ,
                               ADB::function(kro / mu[1], std::move(dmo)) };
    return pmobc;
}

/// Returns fw(sw).
template <class ADB>
ADB
fluxFunc(const std::vector<ADB>& m)
{
    assert (m.size() == 2);

    ADB f = m[0] / (m[0] + m[1]);

    return f;
}


int main()
try
{
    typedef Opm::AutoDiffBlock<double> ADB;
    typedef ADB::V V;
    typedef Eigen::SparseMatrix<double> S;

    Opm::time::StopWatch clock;
    clock.start();
    const Opm::GridManager gm(3,3);//(50, 50, 10);
    const UnstructuredGrid& grid = *gm.c_grid();
    using namespace Opm::unit;
    using namespace Opm::prefix;
    // const Opm::IncompPropertiesBasic props(2, Opm::SaturationPropsBasic::Linear,
    //                                        { 1000.0, 800.0 },
    //                                        { 1.0*centi*Poise, 5.0*centi*Poise },
    //                                        0.2, 100*milli*darcy,
    //                                        grid.dimensions, grid.number_of_cells);
    // const Opm::IncompPropertiesBasic props(2, Opm::SaturationPropsBasic::Linear,
    //                                        { 1000.0, 1000.0 },
    //                                        { 1.0, 1.0 },
    //                                        1.0, 1.0,
    //                                        grid.dimensions, grid.number_of_cells);
    const Opm::IncompPropertiesBasic props(2, Opm::SaturationPropsBasic::Linear,
                                           { 1000.0, 1000.0 },
                                           { 1.0, 30.0 },
                                           1.0, 1.0,
                                           grid.dimensions, grid.number_of_cells);
    V htrans(grid.cell_facepos[grid.number_of_cells]);
    tpfa_htrans_compute(const_cast<UnstructuredGrid*>(&grid), props.permeability(), htrans.data());
    V trans_all(grid.number_of_faces);
    // tpfa_trans_compute(const_cast<UnstructuredGrid*>(&grid), htrans.data(), trans_all.data());
    const int nc = grid.number_of_cells;
    std::vector<int> allcells(nc);
    for (int i = 0; i < nc; ++i) {
        allcells[i] = i;
    }
    std::cerr << "Opm core " << clock.secsSinceLast() << std::endl;

    // Define neighbourhood-derived operator matrices.
    const Opm::HelperOps ops(grid);
    const int num_internal = ops.internal_faces.size();
    std::cerr << "Topology matrices " << clock.secsSinceLast() << std::endl;

    typedef Opm::AutoDiffBlock<double> ADB;
    typedef ADB::V V;

    // q
    V q(nc);
    q.setZero();
    q[0] = 1.0;
    q[nc-1] = -1.0;

    // s0 - this is explicit now
    typedef Eigen::Array<double, Eigen::Dynamic, 2, Eigen::RowMajor> TwoCol;
    TwoCol s0(nc, 2);
    s0.leftCols<1>().setZero();
    s0.rightCols<1>().setOnes();

    // totmob - explicit as well
    TwoCol kr(nc, 2);
    props.relperm(nc, s0.data(), allcells.data(), kr.data(), 0);
    const V krw = kr.leftCols<1>();
    const V kro = kr.rightCols<1>();
    const double* mu = props.viscosity();
    const V totmob = krw/mu[0] + kro/mu[1];

    // Moved down here because we need total mobility.
    tpfa_eff_trans_compute(const_cast<UnstructuredGrid*>(&grid), totmob.data(),
                           htrans.data(), trans_all.data());
    // Still explicit, and no upwinding!
    V mobtransf(num_internal);
    for (int fi = 0; fi < num_internal; ++fi) {
        mobtransf[fi] = trans_all[ops.internal_faces[fi]];
    }
    std::cerr << "Property arrays " << clock.secsSinceLast() << std::endl;

    // Initial pressure.
    V p0(nc,1);
    p0.fill(200*Opm::unit::barsa);

    // First actual AD usage: defining pressure variable.
    const std::vector<int> bpat = { nc };
    // Could actually write { nc } instead of bpat below,
    // but we prefer a named variable since we will repeat it.
    const ADB p = ADB::variable(0, p0, bpat);
    const ADB ngradp = ops.ngrad*p;
    // We want flux = totmob*trans*(p_i - p_j) for the ij-face.
    const ADB flux = mobtransf*ngradp;
    const ADB residual = ops.div*flux - q;
    std::cerr << "Construct AD residual " << clock.secsSinceLast() << std::endl;

    // It's the residual we want to be zero. We know it's linear in p,
    // so we just need a single linear solve. Since we have formulated
    // ourselves with a residual and jacobian we do this with a single
    // Newton step (hopefully easy to extend later):
    //   p = p0 - J(p0) \ R(p0)
    // Where R(p0) and J(p0) are contained in residual.value() and
    // residual.derived()[0].

#if HAVE_SUITESPARSE_UMFPACK_H
    typedef Eigen::UmfPackLU<S> LinSolver;
#else
    typedef Eigen::BiCGSTAB<S>  LinSolver;
#endif  // HAVE_SUITESPARSE_UMFPACK_H

    LinSolver solver;
    S pmatr;
    residual.derivative()[0].toSparse(pmatr);
    pmatr.coeffRef(0,0) *= 2.0;
    pmatr.makeCompressed();
    solver.compute(pmatr);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Pressure/flow Jacobian decomposition error\n";
        return EXIT_FAILURE;
    }
    // const Eigen::VectorXd dp = solver.solve(residual.value().matrix());
    ADB::V residual_v = residual.value();
    const V dp = solver.solve(residual_v.matrix()).array();
    if (solver.info() != Eigen::Success) {
        std::cerr << "Pressure/flow solve failure\n";
        return EXIT_FAILURE;
    }
    const V p1 = p0 - dp;
    std::cerr << "Solve " << clock.secsSinceLast() << std::endl;
    // std::cout << p1 << std::endl;

    // ------ Transport solve ------

    // Now we'll try to do a transport step as well.
    // Residual formula is
    //   R_w = s_w - s_w^0 + dt/pv * (div v_w)
    // where
    //   v_w = f_w v
    // and f_w is (for now) based on averaged mobilities, not upwind.

    double res_norm = 1e100;
    const V sw0 = s0.leftCols<1>();
    // V sw1 = sw0;
    V sw1 = 0.5*V::Ones(nc,1);
    const V ndp = (ops.ngrad * p1.matrix()).array();
    const V dflux = mobtransf * ndp;
    const Opm::UpwindSelector<double> upwind(grid, ops, dflux);
    const V pv = Eigen::Map<const V>(props.porosity(), nc, 1)
        * Eigen::Map<const V>(grid.cell_volumes, nc, 1);
    const double dt = 0.0005;
    const V dtpv = dt/pv;
    const V qneg = q.min(V::Zero(nc,1));
    const V qpos = q.max(V::Zero(nc,1));

    std::cout.setf(std::ios::scientific);
    std::cout.precision(16);

    int it = 0;
    do {
        const ADB sw = ADB::variable(0, sw1, bpat);
        const std::vector<ADB> pmobc = phaseMobility<ADB>(props, allcells, sw.value());
        const std::vector<ADB> pmobf = upwind.select(pmobc);
        const ADB fw_cell = fluxFunc(pmobc);
        const ADB fw_face = fluxFunc(pmobf);
        const ADB flux1 = fw_face * dflux;
        const ADB qtr_ad = qpos + fw_cell*qneg;
        const ADB transport_residual = sw - sw0 + dtpv*(ops.div*flux1 - qtr_ad);
        res_norm = transport_residual.value().matrix().norm();
        std::cout << "res_norm[" << it << "] = "
                  << res_norm << std::endl;

        S smatr;
        transport_residual.derivative()[0].toSparse(smatr);
        smatr.makeCompressed();
        solver.compute(smatr);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Transport Jacobian decomposition error\n";
            return EXIT_FAILURE;
        }
        ADB::V transport_residual_v = transport_residual.value();
        const V ds = solver.solve(transport_residual_v.matrix()).array();
        if (solver.info() != Eigen::Success) {
            std::cerr << "Transport solve failure\n";
            return EXIT_FAILURE;
        }
        sw1 = sw.value() - ds;
        std::cerr << "Solve for s[" << it << "]: "
                  << clock.secsSinceLast() << '\n';
        sw1 = sw1.min(V::Ones(nc,1)).max(V::Zero(nc,1));

        it += 1;
    } while (res_norm > 1e-7);

    std::cout << "Saturation solution:\n"
              << "function s1 = solution\n"
              << "s1 = [\n" << sw1 << "\n];\n";
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

