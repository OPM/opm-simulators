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

#include <opm/autodiff/TransportSolverTwophaseAd.hpp>
#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/props/IncompPropertiesInterface.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <iostream>



namespace Opm
{

    /// Construct solver.
    /// \param[in] grid       A 2d or 3d grid.
    /// \param[in] props      Rock and fluid properties.
    /// \param[in] linsolver  Linear solver for Newton-Raphson scheme.
    /// \param[in] gravity    Gravity vector (null for no gravity).
    /// \param[in] param      Parameters for the solver.
    TransportSolverTwophaseAd::TransportSolverTwophaseAd(const UnstructuredGrid& grid,
                                                         const IncompPropertiesInterface& props,
                                                         const LinearSolverInterface& linsolver,
                                                         const double* gravity,
                                                         const parameter::ParameterGroup& param)
        : grid_(grid),
          props_(props),
          linsolver_(linsolver),
          ops_(grid),
          gravity_(0.0),
          tol_(param.getDefault("nl_tolerance", 1e-9)),
          maxit_(param.getDefault("nl_maxiter", 30))
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid_);
        allcells_.resize(nc);
        for (int i = 0; i < nc; ++i) {
            allcells_[i] = i;
        }
        if (gravity && gravity[dimensions(grid_) - 1] != 0.0) {
            gravity_ = gravity[dimensions(grid_) - 1];
            for (int dd = 0; dd < dimensions(grid_) - 1; ++dd) {
                if (gravity[dd] != 0.0) {
                    OPM_THROW(std::runtime_error, "TransportSolverTwophaseAd: can only handle gravity aligned with last dimension");
                }
            }
            V htrans(grid.cell_facepos[grid.number_of_cells]);
            tpfa_htrans_compute(const_cast<UnstructuredGrid*>(&grid), props.permeability(), htrans.data());
            V trans(numFaces(grid_));
            tpfa_trans_compute(const_cast<UnstructuredGrid*>(&grid), htrans.data(), trans.data());
            transi_ = subset(trans, ops_.internal_faces);
        }
    }




    // Virtual destructor.
    TransportSolverTwophaseAd::~TransportSolverTwophaseAd()
    {
    }





    namespace
    {

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
            // In dkr, columns col(0..3) are:
            //    dkrw/dsw  dkro/dsw  dkrw/dso  dkrw/dso  <-- partial derivatives, really.
            // If we want the derivatives with respect to some variable x,
            // we must apply the chain rule:
            //    dkrw/dx = dkrw/dsw*dsw/dx + dkrw/dso*dso/dx.
            // If x is sw as in our case we are left with.
            //    dkrw/dsw = col(0) - col(2)
            //    dkro/dsw = col(1) - col(3)
            V dkrw = dkr.leftCols<1>() - dkr.rightCols<2>().leftCols<1>();
            V dkro = dkr.leftCols<2>().rightCols<1>() - dkr.rightCols<1>();
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

    } // anonymous namespace


    /// Solve for saturation at next timestep.
    /// Note that this only performs advection by total velocity, and
    /// no gravity segregation.
    /// \param[in]      porevolume   Array of pore volumes.
    /// \param[in]      source       Transport source term. For interpretation see Opm::computeTransportSource().
    /// \param[in]      dt           Time step.
    /// \param[in, out] state        Reservoir state. Calling solve() will read state.faceflux() and
    ///                              read and write state.saturation().
    void TransportSolverTwophaseAd::solve(const double* porevolume,
                                          const double* source,
                                          const double dt,
                                          TwophaseState& state)
    {
        using namespace Opm::AutoDiffGrid;
        typedef Eigen::Array<double, Eigen::Dynamic, 2, Eigen::RowMajor> TwoCol;
        typedef Eigen::Map<const V> Vec;
        const int nc = numCells(grid_);
        const TwoCol s0 = Eigen::Map<const TwoCol>(state.saturation().data(), nc, 2);
        double res_norm = 1e100;
        const V sw0 = s0.leftCols<1>();
        // sw1 is the object that will be changed every Newton iteration.
        // V sw1 = sw0;
        V sw1 = 0.5*V::Ones(nc,1);
        const V dflux_all = Vec(state.faceflux().data(), numFaces(grid_), 1);
        const int num_internal = ops_.internal_faces.size();
        V dflux = subset(dflux_all, ops_.internal_faces);

        // Upwind selection of mobilities by phase.
        // We have that for a phase P
        //   v_P = lambda_P K (-grad p + rho_P g grad z)
        // and we assume that this has the same direction as
        //   dh_P = -grad p + rho_P g grad z.
        // This may not be true for arbitrary anisotropic situations,
        // but for scalar lambda and using TPFA it holds.
        const V p1 = Vec(state.pressure().data(), nc, 1);
        const V ndp = (ops_.ngrad * p1.matrix()).array();
        const V z = cellCentroidsZToEigen(grid_);
        const V ndz = (ops_.ngrad * z.matrix()).array();
        assert(num_internal == ndp.size());
        const double* density = props_.density();
        const V dhw = ndp - ndz*(gravity_*density[0]);
        const V dho = ndp - ndz*(gravity_*density[1]);
        const UpwindSelector<double> upwind_w(grid_, ops_, dhw);
        const UpwindSelector<double> upwind_o(grid_, ops_, dho);

        // Compute more explicit and constant terms used in the equations.
        const V pv = Vec(porevolume, nc, 1);
        const V dtpv = dt/pv;
        const V q = Vec(source, nc, 1);
        const V qneg = q.min(V::Zero(nc,1));
        const V qpos = q.max(V::Zero(nc,1));
        const double gfactor = gravity_*(density[0] - density[1]);
        const V gravflux = (gravity_ == 0.0) ? V(V::Zero(num_internal, 1))
            : ndz*transi_*gfactor;

        // Block pattern for variables.
        // Primary variables:
        //    sw : one per cell
        std::vector<int> bpat = { nc };

        // Newton-Raphson loop.
        int it = 0;
        do {
            // Assemble linear system.
            const ADB sw = ADB::variable(0, sw1, bpat);
            const std::vector<ADB> pmobc = phaseMobility<ADB>(props_, allcells_, sw.value());
            const ADB fw_cell = fluxFunc(pmobc);
            const std::vector<ADB> pmobf = { upwind_w.select(pmobc[0]),
                                             upwind_o.select(pmobc[1])  };
            const ADB fw_face = fluxFunc(pmobf);
            const ADB flux = fw_face * (dflux - pmobf[1]*gravflux);
            // const ADB fw_face = upwind_w.select(fw_cell);
            // const ADB flux = fw_face * dflux;
            const ADB qtr_ad = qpos + fw_cell*qneg;
            const ADB transport_residual = sw - sw0 + dtpv*(ops_.div*flux - qtr_ad);
            res_norm = transport_residual.value().matrix().norm();
            std::cout << "Residual l2-norm = " << res_norm << std::endl;

            // Solve linear system.
            Eigen::SparseMatrix<double, Eigen::RowMajor> smatr;
            transport_residual.derivative()[0].toSparse(smatr);
            assert(smatr.isCompressed());
            V ds(nc);
            LinearSolverInterface::LinearSolverReport rep
                = linsolver_.solve(nc, smatr.nonZeros(),
                                   smatr.outerIndexPtr(), smatr.innerIndexPtr(), smatr.valuePtr(),
                                   transport_residual.value().data(), ds.data());
            if (!rep.converged) {
                OPM_THROW(LinearSolverProblem, "Linear solver convergence error in TransportSolverTwophaseAd::solve()");
            }

            // Update (possible clamp) sw1.
            sw1 = sw.value() - ds;
            sw1 = sw1.min(V::Ones(nc,1)).max(V::Zero(nc,1));
            it += 1;
        } while (res_norm > tol_ && it < maxit_);

        // Write to output data structure.
        Eigen::Map<TwoCol> sref(state.saturation().data(), nc, 2);
        sref.leftCols<1>() = sw1;
        sref.rightCols<1>() = 1.0 - sw1;
    }


} // namespace Opm
