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

#include "TransportSolverTwophaseAd.hpp"
#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <iostream>



namespace Opm
{

    /// Construct solver.
    /// \param[in] grid       A 2d or 3d grid.
    /// \param[in] props      Rock and fluid properties.
    /// \param[in] linsolver  Linear solver for Newton-Raphson scheme.
    /// \param[in] param      Parameters for the solver.
    TransportSolverTwophaseAd::TransportSolverTwophaseAd(const UnstructuredGrid& grid,
                                                         const IncompPropertiesInterface& props,
                                                         const LinearSolverInterface& linsolver,
                                                         const parameter::ParameterGroup& param)
        : grid_(grid),
          props_(props),
          linsolver_(linsolver),
          ops_(grid),
          tol_(param.getDefault("nl_tolerance", 1e-9)),
          maxit_(param.getDefault("nl_maxiter", 30))
    {
        const int nc = grid_.number_of_cells;
        allcells_.resize(nc);
        for (int i = 0; i < nc; ++i) {
            allcells_[i] = i;
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
            M krwjac(nc,nc);
            M krojac(nc,nc);
            auto sizes = Eigen::ArrayXi::Ones(nc);
            krwjac.reserve(sizes);
            krojac.reserve(sizes);
            for (int c = 0; c < nc; ++c) {
                krwjac.insert(c,c) = dkrw(c);
                krojac.insert(c,c) = dkro(c);
            }
            const double* mu = props.viscosity();
            std::vector<M> dmw = { krwjac/mu[0] };
            std::vector<M> dmo = { krojac/mu[1] };

            std::vector<ADB> pmobc = { ADB::function(krw / mu[0], dmw) ,
                                       ADB::function(kro / mu[1], dmo) };
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
        typedef Eigen::Array<double, Eigen::Dynamic, 2, Eigen::RowMajor> TwoCol;
        typedef Eigen::Map<const V> Vec;
        const int nc = grid_.number_of_cells;
        const TwoCol s0 = Eigen::Map<const TwoCol>(state.saturation().data(), nc, 2);
        double res_norm = 1e100;
        const V sw0 = s0.leftCols<1>();
        // V sw1 = sw0;
        V sw1 = 0.5*V::Ones(nc,1);
        const V p1 = Vec(state.pressure().data(), nc, 1);
        const V ndp = (ops_.ngrad * p1.matrix()).array();
        const V dflux_all = Vec(state.faceflux().data(), grid_.number_of_faces, 1);
        const int num_internal = ops_.internal_faces.size();
        V dflux(num_internal);
        for (int fi = 0; fi < num_internal; ++fi) {
            dflux[fi] = dflux_all[ops_.internal_faces[fi]];
        }
        const UpwindSelectorTotalFlux<double> upwind(grid_, ops_, dflux);
        const V pv = Vec(porevolume, nc, 1);
        const V dtpv = dt/pv;
        const V q = Vec(source, nc, 1);
        const V qneg = q.min(V::Zero(nc,1));
        const V qpos = q.max(V::Zero(nc,1));

        // Block pattern for variables.
        // Primary variables:
        //    sw : one per cell
        std::vector<int> bpat = { nc };

        // Newton-Raphson loop.
        int it = 0;
        do {
            // Assemble linear system/
            const ADB sw = ADB::variable(0, sw1, bpat);
            const std::vector<ADB> pmobc = phaseMobility<ADB>(props_, allcells_, sw.value());
            const std::vector<ADB> pmobf = upwind.select(pmobc);
            const ADB fw_cell = fluxFunc(pmobc);
            const ADB fw_face = fluxFunc(pmobf);
            const ADB flux1 = fw_face * dflux;
            const ADB qtr_ad = qpos + fw_cell*qneg;
            const ADB transport_residual = sw - sw0 + dtpv*(ops_.div*flux1 - qtr_ad);
            res_norm = transport_residual.value().matrix().norm();

            // Solve linear system.
            Eigen::SparseMatrix<double, Eigen::RowMajor> smatr = transport_residual.derivative()[0];
            ASSERT(smatr.isCompressed());
            V ds(nc);
            LinearSolverInterface::LinearSolverReport rep
                = linsolver_.solve(nc, smatr.nonZeros(),
                                   smatr.outerIndexPtr(), smatr.innerIndexPtr(), smatr.valuePtr(),
                                   transport_residual.value().data(), ds.data());
            if (!rep.converged) {
                THROW("Linear solver convergence error in TransportSolverTwophaseAd::solve()");
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
