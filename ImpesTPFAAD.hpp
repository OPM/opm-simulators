/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2013 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_IMPESTPFAAD_HEADER_INCLUDED
#define OPM_IMPESTPFAAD_HEADER_INCLUDED

#include "AutoDiffBlock.hpp"

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>

#include <algorithm>
#include <vector>

#include <boost/shared_ptr.hpp>

struct UnstructuredGrid;

namespace {
    std::vector<int>
    buildAllCells(const int nc)
    {
        std::vector<int> all_cells(nc);

        for (int c = 0; c < nc; ++c) { all_cells[c] = c; }

        return all_cells;
    }
}

namespace Opm {
    class LinearSolverInterface;

    template <typename Scalar, class BOProps>
    class PressureDependentFluidData {
    public:
        typedef AutoDiff::ForwardBlock<Scalar> ADB;

        PressureDependentFluidData(const int      nc,
                                   const BOProps& props)
            : nc_   (nc)
            , np_   (props.numPhases())
            , cells_(buildAllCells(nc))
            , props_(props)
            , A_    (nc_, np_ * np_)
            , dA_   (nc_, np_ * np_)
            , mu_   (nc_, np_      )
            , dmu_  (nc_, np_      )
            , kr_   (nc_, np_      )
            , zero_ (ADB::V::Zero(nc, 1))
            , one_  (ADB::V::Ones(nc, 1))
        {
        }

        void
        computeSatQuant(const BlackoilState& state)
        {
            const std::vector<double>& s = state.saturation();

            double* dkrds = 0;  // Ignore rel-perm derivatives
            props_.relperm(nc_, & s[0], & cells_[0],
                           kr_.data(), dkrds);
        }

        void
        computePressQuant(const BlackoilState& state)
        {
            const std::vector<double>& p = state.pressure();
            const std::vector<double>& z = state.surfacevol();

            props_.matrix   (nc_, & p[0], & z[0], & cells_[0],
                             A_ .data(), dA_ .data());

            props_.viscosity(nc_, & p[0], & z[0], & cells_[0],
                             mu_.data(), dmu_.data());
        }

        ADB
        fvf(const int p) const
        {
            assert (0 <= p  );
            assert (p <  np_);

            typename ADB::V A   = A_ .block(0, p * (np_ + 1), nc_, 1);
            typename ADB::V dA  = dA_.block(0, p * (np_ + 1), nc_, 1);
            typename ADB::M jac = dA.matrix().asDiagonal();

            return one_ / ADB::function(A, jac);
        }

        typename ADB::V
        phaseRelPerm(const int p) const
        {
            typename ADB::V kr = kr_.block(0, p, nc_, 1);

            return kr;
        }

        ADB
        phaseViscosity(const int p) const
        {
            assert (0 <= p  );
            assert (p <  np_);

            typename ADB::V mu  = mu_ .block(0, p, nc_, 1);
            typename ADB::V dmu = dmu_.block(0, p, nc_, 1);
            typename ADB::M jac = dmu.matrix().asDiagonal();

            return ADB::function(mu, jac);
        }

    private:
        const int nc_;
        const int np_;

        const std::vector<int> cells_;
        const BOProps&         props_;

        typedef Eigen::Array<Scalar,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor> DerivedQuant;

        // Pressure dependent quantities (essentially B and \mu)
        DerivedQuant A_  ;
        DerivedQuant dA_ ;
        DerivedQuant mu_ ;
        DerivedQuant dmu_;

        // Saturation dependent quantities (rel-perm only)
        DerivedQuant kr_;

        const typename ADB::V zero_;
        const typename ADB::V one_ ;
    };


    template <class BOProps>
    class ImpesTPFAAD {
    public:
        ImpesTPFAAD(const UnstructuredGrid& grid ,
                    const BOProps&          props)
            : grid_     (grid)
            , pdepfdata_(grid.number_of_cells, props)
        {
        }

    private:
        // Disallow copying and assignment
        ImpesTPFAAD(const ImpesTPFAAD& rhs);
        ImpesTPFAAD& operator=(const ImpesTPFAAD& rhs);

        typedef PressureDependentFluidData<double,BOProps> PDepFData;

        const UnstructuredGrid& grid_;
        PDepFData               pdepfdata_;
    };
} // namespace Opm

#endif  /* OPM_IMPESTPFAAD_HEADER_INCLUDED */
