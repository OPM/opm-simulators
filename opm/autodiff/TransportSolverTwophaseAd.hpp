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

#ifndef OPM_TRANSPORTSOLVERTWOPHASEAD_HEADER_INCLUDED
#define OPM_TRANSPORTSOLVERTWOPHASEAD_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/transport/TransportSolverTwophaseInterface.hpp>
#include <vector>

struct UnstructuredGrid;

namespace Opm
{

    class IncompPropertiesInterface;
    class LinearSolverInterface;
    namespace parameter { class ParameterGroup; }

    /// Implements an implicit transport solver for incompressible two-phase flow,
    /// using automatic differentiation.
    class TransportSolverTwophaseAd : public TransportSolverTwophaseInterface
    {
    public:
        /// Construct solver.
        /// \param[in] grid       A 2d or 3d grid.
        /// \param[in] props      Rock and fluid properties.
        /// \param[in] linsolver  Linear solver for Newton-Raphson scheme.
        /// \param[in] gravity    Gravity vector (null for no gravity).
        /// \param[in] param      Parameters for the solver.
        TransportSolverTwophaseAd(const UnstructuredGrid& grid,
                                  const IncompPropertiesInterface& props,
                                  const LinearSolverInterface& linsolver,
                                  const double* gravity,
                                  const parameter::ParameterGroup& param);

        // Virtual destructor.
        virtual ~TransportSolverTwophaseAd();

        /// Solve for saturation at next timestep.
        /// Note that this only performs advection by total velocity, and
        /// no gravity segregation.
        /// \param[in]      porevolume   Array of pore volumes.
        /// \param[in]      source       Transport source term. For interpretation see Opm::computeTransportSource().
        /// \param[in]      dt           Time step.
        /// \param[in, out] state        Reservoir state. Calling solve() will read state.faceflux() and
        ///                              read and write state.saturation().
        virtual void solve(const double* porevolume,
                           const double* source,
                           const double dt,
                           TwophaseState& state);

    private:
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;

    private:
        const UnstructuredGrid& grid_;
        const IncompPropertiesInterface& props_;
        const LinearSolverInterface& linsolver_;
        const HelperOps ops_;
        double gravity_;
        double tol_;
        int maxit_;
        std::vector<int> allcells_;
        V transi_;
    };

} // namespace Opm

#endif // OPM_TRANSPORTSOLVERTWOPHASEAD_HEADER_INCLUDED
