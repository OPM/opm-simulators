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

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>

struct UnstructuredGrid;
struct Wells;

namespace Opm {

    class DerivedGeology;
    class LinearSolverInterface;
    class BlackoilState;
    class WellState;

    /// Class for solving black-oil impes problems.
    /// Current known limitations:
    ///   - pressure solve only
    ///   - no miscibility
    ///   - no gravity in wells or crossflow
    class ImpesTPFAAD
    {
    public:
        /// Construct impes solver.
        ImpesTPFAAD(const UnstructuredGrid&         grid,
                    const BlackoilPropsAdInterface& fluid,
                    const DerivedGeology&           geo,
                    const Wells&                    wells,
                    const LinearSolverInterface&    linsolver);

        /// Solve forward in time.
        /// Currently this will only modify
        /// state.pressure(), state.faceflux(), well_state.bhp()
        /// and well_state.wellRates().
        void solve(const double   dt,
                   BlackoilState& state,
                   WellState&     well_state);
    private:
        // Disallow copying and assignment
        ImpesTPFAAD(const ImpesTPFAAD& rhs);
        ImpesTPFAAD& operator=(const ImpesTPFAAD& rhs);

        // Types
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;
        typedef Eigen::Array<double,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor> DataBlock;
        enum { Water = BlackoilPropsAdInterface::Water,
               Oil = BlackoilPropsAdInterface::Oil,
               Gas = BlackoilPropsAdInterface::Gas };

        // Data
        const UnstructuredGrid&      grid_;
        const BlackoilPropsAdInterface&               fluid_;
        const DerivedGeology&        geo_ ;
        const Wells&                 wells_;
        const LinearSolverInterface& linsolver_;
        HelperOps                    ops_;
        const M                      grav_;
        ADB                          cell_residual_;
        std::vector<ADB>             well_flow_residual_;
        ADB                          well_residual_;
        ADB                          total_residual_;
        std::vector<V>               kr_;
        std::vector<V>               well_kr_;
        ADB                          qs_;
        V                            well_perf_dp_;

        // Methods for assembling and solving.
        void computeExplicitData(const double         dt,
                                 const BlackoilState& state,
                                 const WellState& well_state);
        void assemble(const double         dt,
                      const BlackoilState& state,
                      const WellState& well_state);
        void solveJacobianSystem(BlackoilState& state,
                                 WellState& well_state) const;
        double residualNorm() const;
        void computeFluxes(BlackoilState& state, WellState& well_state) const;

        // Fluid interface forwarding calls to correct methods of fluid_.
        V fluidMu(const int phase, const V& p, const std::vector<int>& cells) const;
        ADB fluidMu(const int phase, const ADB& p, const std::vector<int>& cells) const;
        V fluidFvf(const int phase, const V& p, const std::vector<int>& cells) const;
        ADB fluidFvf(const int phase, const ADB& p, const std::vector<int>& cells) const;
        V fluidRho(const int phase, const V& p, const std::vector<int>& cells) const;
        ADB fluidRho(const int phase, const ADB& p, const std::vector<int>& cells) const;
        V fluidKr(const int phase) const;
        V fluidKrWell(const int phase) const;
    };


} // namespace Opm

#endif  /* OPM_IMPESTPFAAD_HEADER_INCLUDED */
